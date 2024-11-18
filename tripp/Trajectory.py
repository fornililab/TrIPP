"""
    @release_date  : $release_date
    @version       : $release_version
    @author        : Christos Matsingos, Ka Fu Man

    This file is part of the TrIPP software
    (https://github.com/fornililab/TrIPP).
    Copyright (c) 2024 Christos Matsingos, Ka Fu Man and Arianna Fornili.

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, version 3.

    This program is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program. If not, see <http://www.gnu.org/licenses/>.
"""

import os
import multiprocessing as mp
from tripp._create_mda_universe_ import (
    create_mda_universe,
    create_propka_compatible_universe,
)
from tripp._pka_iterator_ import pka_iterator
from tripp._sort_pka_df_ import sort_pka_df, sort_buriedness_df
from tripp._detect_disulphide_bonds_ import detect_disulphide_bonds 
from datetime import datetime
from tripp._generate_trajectory_log_ import log_header, trajectory_log
from tripp._edit_pdb_ import find_mutation
import logging
import glob


class Trajectory:
    """
    Main class of TrIPP. Calling this class creates an iterable object of
    sliced trajectories which are then used with the run method to run the
    analysis. The arguments taken are a trajectory file (formats supported
    by MDAnalysis), a topology file (usually a PDB file but can be all
    formats supported by MDAnalaysis) and the number of CPU cores to be
    used to run the analysis.

    Parameters
    ----------
    trajectory_file: str
        The path of the file containing the trajectory. The same formats
        permited by MDAnalysis can be used.
    topology_file: str
        The path of the file containing the topology. The same formats
        allowed by MDAnalysis can be used.
    output_directory: str
        The directory you want to save the output files.
    output_prefix: str
        The output prefix for output files.
    cpu_core_number: int, default=-1
        The number of cpu cores used for the calculation.
        If cpu_core_number=-1, all available cores are used.
    hetatm_resid: int or list of int
    
    """
    def __init__(
        self,
        topology_file,
        trajectory_file,
        output_directory,
        output_prefix,
        cpu_core_number=-1,
        hetatm_resid=None,
        custom_terminal_oxygens=None,
        custom_resname_correction=None,
    ):  
        self.output_directory = output_directory
        self.output_prefix = output_prefix
        if not os.path.isdir(output_directory):
            # Make directory if not present
            os.makedirs(output_directory) 
        else:
            # Remove files named .temp* in the directory before proceeding.
            [os.remove(file) for file in glob.glob(f'{output_directory}/.temp*')]
        if os.path.isfile(f'{output_directory}/{output_prefix}.log'):
            os.remove(f'{output_directory}/{output_prefix}.log')
            
        self.logger = logging.getLogger()
        for handler in self.logger.handlers[:]:
            self.logger.removeHandler(handler)
        self.logger.setLevel(logging.INFO)
        handler = logging.FileHandler(f'{output_directory}/{output_prefix}.log',
                                      'a')
        handler.setLevel(logging.INFO)
        formatter = logging.Formatter('%(message)s')
        handler.setFormatter(formatter)
        self.logger.addHandler(handler)
        self.logger.info(log_header())
        self.logger.info('Trajectory class initialised:')
        
        self.trajectory_file = trajectory_file
        self.topology_file = topology_file

        if cpu_core_number == -1:
            self.cpu_core_number = os.cpu_count()

        else:
            self.cpu_core_number = cpu_core_number

        self.universe = create_mda_universe(
            topology_file=self.topology_file,
            trajectory_file=self.trajectory_file
        )

        self.corrected_universe = create_propka_compatible_universe(
            self.universe,
            hetatm_resid,
            custom_terminal_oxygens,
            custom_resname_correction,
        )

        frames_nr = len(self.universe.trajectory)
        slices_nr = self.cpu_core_number
        slice_length = frames_nr // slices_nr
        remainder = frames_nr % slices_nr
        slices = []
        start_frame = 0

        for i in range(slices_nr):
            if i == slices_nr - 1:
                end_frame = start_frame + slice_length + remainder
            else:
                end_frame = start_frame + slice_length
            slices.append([start_frame, end_frame])
            start_frame = end_frame

        self.trajectory_slices = slices

    def run(
        self,
        extract_buriedness_data=True,
        chain="A",
        mutation=None,
        disulphide_bond_detection=True,
        optargs=[],
    ):
        r"""
        Function to perform PROPKA after initialising the Trajectory class

        Parameters
        ----------
        extract_buriedness_data: bool, default=True
            If set to True both data on buriedness and pKa will be extracted.
            If set to False only pKa data will be extracted.
        chain: str, default='A'
            The chain ID to perform PROPKA prediction.
        mutation: int or list of int, default=None
            Peform pseudomutation of residues to alanine.
            The residue number or list of residue numbers (MDAnalysis resid)
            can be provided.
        disulphide_bond_detection: bool, default=True
            If set to True detects all disulphide bonds present in the
            topology file and does not provide pKa values for them.
        optargs: list of str, default=[]
            PROPKA prediction can be run with optional arguments as indicated
            in their documentation. Each flag is string separated by a comma.
            For example, `["-k","--pH=7.2"]` keeps the proton in your topology
            and setting pH 7.2 for ie: stability calculations, respectively.
        """
        start = datetime.today().strftime("%Y-%m-%d %H:%M:%S")
            
        pool = mp.Pool(self.cpu_core_number)
        # Create jobs
        jobs = []
        for trajectory_slice in self.trajectory_slices:
            # Create asynchronous jobs that will be submitted once a
            # processor is ready
            job = pool.apply_async(
                pka_iterator,
                args=(
                    trajectory_slice,
                    self.corrected_universe,
                    self.output_directory,
                    mutation,
                    chain,
                    extract_buriedness_data,
                    optargs,
                ),
            )
            jobs.append(job)
        # Submit jobs
        results = [job.get() for job in jobs]
        pool.close()
        pool.join()
        
        pka_dfs, buriedness_dfs = zip(*results)

        # Detect disulphide bond and remove it from the pKa and buriedness CSV
        # if set to True.
        if disulphide_bond_detection:
            disulphide_cysteines_list = detect_disulphide_bonds(self.topology_file)
        else:
            disulphide_cysteines_list = []
        
        # Combine the temporary pka csv and sort it according to Time [ps] column
        sort_pka_df(
            output_directory=self.output_directory,
            output_prefix=self.output_prefix,
            pka_dfs=pka_dfs,
            disulphide_cysteines_list=disulphide_cysteines_list,
        )
        if extract_buriedness_data:
            # Combine the temporary buriedness csv and sort it according to Time [ps] column
            sort_buriedness_df(
                output_directory=self.output_directory,
                output_prefix=self.output_prefix,
                buriedness_dfs=buriedness_dfs,
                disulphide_cysteines_list=disulphide_cysteines_list,
            )

        end = datetime.today().strftime("%Y-%m-%d %H:%M:%S")
        # Find the corresponding mutation and retreive it to print it in the log.
        if mutation is not None:
            mutation_arg = find_mutation(mutation, self.corrected_universe)
        else:
            mutation_arg = mutation
        # Log all parameter used for the run and also pKa statistic table.
        trajectory_log(
            self.output_directory,
            self.output_prefix,
            extract_buriedness_data,
            chain,
            mutation_arg,
            disulphide_cysteines_list,
            optargs,
            self.cpu_core_number,
            self.trajectory_slices,
            start,
            end,
        )