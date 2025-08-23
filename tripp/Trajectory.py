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
from tripp._sort_pka_df_ import output_df
from datetime import datetime
from tripp._generate_trajectory_log_ import log_header, trajectory_log
import logging
import glob

class Trajectory:
    """
    Main class of TrIPP. Calling this class creates an iterable object of
    sliced trajectories that are then used with the run method to perform the
    analyses.

    Parameters
    ----------
    trajectory_file: str
        The path of the file containing the trajectory. The same formats
        permited by MDAnalysis can be used.
    topology_file: str
        The path of the file containing the topology. The same formats
        allowed by MDAnalysis can be used.
    output_directory: str
        The directory where output files will be saved.
    output_prefix: str
        The output file prefix.
    cpu_core_number: int, default=-1
        The number of cpu cores used for the calculation.
        If cpu_core_number=-1, all available cores are used.
    hetatm_resname: str, list, default=None
        PDB residue name(s) for non-protein molecules that we want to 
        be taken into account in the pKa calculation. Their record type
        will be set to 'HETATM'. See log for info if an error is raised.
    custom_terminal_oxygens: list, default=None
        PROPKA only recognizes C-terminal oxygen atoms named O and OXT. 
        If different, either provide a list of length 2 containing the 
        names of your C-terminal oxygen atoms, or set to None for correction 
        according to our pre-defined dictionary 
        (see _correction_dictionary_.py for the pre-defined names).
    custom_resname_correction: dict, default=None
        dictionary of custom protein residue names not included in
        the hard-coded TrIPP dictionary (tripp._correction_dictionary_.py). 
        Can be given as e.g. {'XXX':'ASP'}, where 'XXX' is the residue
        name in the PDB file and 'ASP' is the corresponding PROPKA name.

    """
    def __init__(
        self,
        topology_file,
        trajectory_file,
        output_directory,
        output_prefix,
        cpu_core_number=-1,
        hetatm_resname=None,
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
        self.logger.setLevel(logging.INFO)
        handler = logging.FileHandler(f'{output_directory}/{output_prefix}.log',
                                      'a')
        handler.setLevel(logging.INFO)
        formatter = logging.Formatter('%(message)s')
        handler.setFormatter(formatter)
        self.logger.addHandler(handler)
        propka_logger = logging.getLogger('propka') # Prevent loggers from PROPKA logging into our root logger.
        propka_logger.propagate = False 
        
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
            trajectory_file=self.trajectory_file,
        )

        self.corrected_universe = create_propka_compatible_universe(
            self.universe,
            hetatm_resname,
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
            if i < remainder:
                end_frame = start_frame + slice_length + 1
            else:
                end_frame = start_frame + slice_length
            slices.append([start_frame, end_frame])
            start_frame = end_frame

        self.trajectory_slices = slices

    def run(
        self,
        extract_buriedness_data=True,
        mutation_selections=None,
        save_disulphide_pka=False,
        optargs=[],
    ):
        """
        Function to run PROPKA after initialising the Trajectory class

        Parameters
        ----------
        extract_buriedness_data: bool, default=True
            If set to True, both buried ratios and pKa values will be extracted.
            If set to False, only pKa values will be extracted.
        mutation_selections: str, default=None
            Peform pseudomutation of residues to alanine.
            Selection is based on MDAnalysis syntax. For multi-chain systems,
            please make sure you include the chainID in the selection. Double
            mutations can also be performed.
            e.g.: chainID A and resid 2 3
        save_disulphide_pka: bool, default=False
            If set to False, pKa and buried ratio values for cysteines 
            forming a disulphide bond (pKa set by PROPKA to 99.99) will 
            not be saved in the CSV files. 
        optargs: list of str, default=[]
            PROPKA predictions can be run with optional arguments
            (see https://propka.readthedocs.io/en/latest/command.html). 
            For example, if optargs is set to `["-k"]`, propka will run with the -k flag 
            (protons from the input file are kept).
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
                    mutation_selections,
                    optargs),
            )
            jobs.append(job)
        # Submit jobs
        results = [job.get() for job in jobs]
        pool.close()
        pool.join()
        
        data, log_contents = zip(*results)
        self.data = data
        log_contents = list(filter(bool,log_contents))
        if log_contents:
            propka_warning_logger = logging.getLogger('proka_warning')
            propka_warning_logger.propagate = False
            propka_warning_handler = logging.FileHandler(f'{self.output_directory}/{self.output_prefix}_propka_warnings.log','w')
            propka_warning_handler.setLevel(logging.WARNING)
            propka_warning_handler.setFormatter(logging.Formatter('%(message)s'))
            propka_warning_logger.addHandler(propka_warning_handler)
            propka_warning_logger.warning(log_contents[0])
            propka_warning_handler.close()
            propka_warning_logger.removeHandler(propka_warning_handler)

        
        # Combine the temporary pka csv and sort it according to Time [ps] column
        disulphide_cys_col = output_df(
            output_directory=self.output_directory,
            output_prefix=self.output_prefix,
            data=data,
            save_disulphide_pka=save_disulphide_pka,
            extract_buriedness_data=extract_buriedness_data
        )

        end = datetime.today().strftime("%Y-%m-%d %H:%M:%S")

        # Log all parameter used for the run and also pKa statistic table.
        trajectory_log(
            self.output_directory,
            self.output_prefix,
            extract_buriedness_data,
            mutation_selections,
            disulphide_cys_col,
            optargs,
            self.cpu_core_number,
            self.trajectory_slices,
            start,
            end,
        )
        for handler in self.logger.handlers[:]:
            handler.close()
            self.logger.removeHandler(handler)