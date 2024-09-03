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
from tripp._create_mda_universe_ import create_mda_universe,create_propka_compatible_universe
from tripp._pka_iterator_ import pka_iterator 
from tripp._sort_pka_df_ import sort_pka_df 
from datetime import datetime
from tripp._generate_trajectory_log_ import trajectory_log

class Trajectory: 

    """
    Main class of TrIPP. Calling this class creates an iterable object of sliced 
    trajectories which are then used with the run method to run the analysis. The 
    arguments taken are a trajectory file (formats supported by MDAnalysis), a topology 
    file (usually a PDB file but can be all formats supported by MDAnalaysis) and the 
    number of CPU cores to be used to run the analysis. 

    The class takes as input: 

    trajectory_file: str 
    The path of the file containing the trajectory. The same formats permited by 
    MDAnalysis can be used. 

    topology_file: str 
    The path of the file containing the topology. The same formats allowed by 
    MDAnalysis can be used. 

    cpu_core_number: int, default=-1 
    The number of cpu cores used for the calculation. If cpu_core_number=-1, all 
    available cores are used. 
    """

    def __init__(self, trajectory_file, topology_file, cpu_core_number=-1): 

        self.trajectory_file = trajectory_file 
        self.topology_file = topology_file 

        if cpu_core_number == -1: 
            self.cpu_core_number = os.cpu_count() 

        else: 
            self.cpu_core_number = cpu_core_number 

        self.universe = create_mda_universe(topology_file=self.topology_file, trajectory_file=self.trajectory_file)
        self.corrected_universe = create_propka_compatible_universe(self.universe)
        frames_nr = len(self.universe.trajectory)
        slices_nr = self.cpu_core_number 
        slice_length = round(frames_nr/slices_nr) 
        slices = [] 
        start_frame = 0 
        end_frame = slice_length 

        for i in range(slices_nr): 

            if i == slices_nr - 1: # Slice_length might not get to the end of trajectory, using this condition make sure the last trunk include the end of the trajectory
                slices.append([start_frame, frames_nr]) 

            else: 
                slices.append([start_frame, end_frame]) 
                start_frame+=slice_length 
                end_frame+=slice_length 

        self.trajectory_slices = slices

    def run(self,
            output_directory,
            output_prefix, 
            extract_buridness_data, 
            chain, 
            mutation, 
            disulphide_bond_detection,
            optargs):
        start = datetime.today().strftime('%Y-%m-%d %H:%M:%S')
        
        if os.path.isdir(output_directory) != True:
            os.makedirs(output_directory)
            
        pool = mp.Pool(self.cpu_core_number)
        # Create jobs
        jobs = []
        for index, _ in enumerate(self.trajectory_slices):
            # Create asynchronous jobs that will be submitted once a processor is ready
            job = pool.apply_async(pka_iterator, args=(self.trajectory_slices, 
                                                       index,
                                                       self.corrected_universe,
                                                       output_directory,
                                                       mutation,
                                                       chain,
                                                       extract_buridness_data,
                                                       optargs))
            jobs.append(job)
        # Submit jobs
        results = [job.get() for job in jobs]
        pool.close()
        pool.join()
        # Post-processing the data
        sort_pka_df(cores=self.cpu_core_number, 
                    topology_file=self.topology_file,
                    output_directory=output_directory,
                    output_prefix=output_prefix, 
                    extract_buridness_data=extract_buridness_data, 
                    mutation=mutation, 
                    disulphide_bond_detection=disulphide_bond_detection) 
        # Generate pKa summary
        # pka_statistics_table(output_directory,output_prefix)
        
        end = datetime.today().strftime('%Y-%m-%d %H:%M:%S')
        # Logging
        trajectory_log(output_directory,
                       output_prefix, 
                       extract_buridness_data,
                       chain, 
                       mutation, 
                       disulphide_bond_detection,
                       optargs,
                       self.cpu_core_number,
                       self.trajectory_slices,
                       start,
                       end)
