"""
    @release_date  : $release_date
    @version       : $release_version
    @author        : Christos Matsingos
    
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

import MDAnalysis as mda
import propka 
from propka import run
import numpy as np 
import os 
import multiprocessing as mp 
from tripp.edit_pdb import edit_pdb 
from tripp.edit_pdb import mutate 
import pandas as pd 

class Trajectory: 

    """
    Main class of TrIPP. Calling this class creates an iterable object of sliced trajectories which are 
    then used with the run method to run the analysis. The arguments taken are a trajectory file (formats supported by MDAnalysis), 
    a topology file (usually a PDB file but can be all formats supported by MDAnalaysis) and the number 
    of CPU cores to be used to run the analyss. 
    """

    def __init__(self, trajectory_file, topology_file, cpu_core_number): 

        self.trajectory_file = trajectory_file 
        self.topology_file = topology_file 
        self.cpu_core_number = cpu_core_number 
        self.universe = mda.Universe(self.topology_file, self.trajectory_file) 
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
        
    def calculate_pka(self, output_file, extract_surface_data=False, chain='A', mutation=None, core=None): 
        
        def extract_data(file, chain, time): 

            pkafile = open(file, 'r') 
            data_pka_list = [] 
            data_surf_list = [] 
            for line in pkafile: 
                line_processed = line.rstrip() 
                line_list = line_processed.strip().split() 
                if len(line_list) > 15 and line_list[2] == chain: 
                    data_pka_list.append([line_list[0]+line_list[1], line_list[3]]) 
                    data_surf_list.append([line_list[0]+line_list[1], line_list[4]]) 
            
            time_array = np.array([['Time [ps]'], [time]]) 
            data_pka_array = np.array(data_pka_list).T 
            data_pka = np.concatenate((time_array, data_pka_array), axis=1).tolist() 
            
            data_surf_array = np.array(data_surf_list).T 
            data_surf = np.concatenate((time_array, data_surf_array), axis=1).tolist() 
            
            return data_pka, data_surf
        
        if type(mutation) == int: 
            out = f'{output_file}_{mutation}' 
            temp_name = f'temp_{mutation}_{core}' 
        elif type(mutation) == list: 
            out = f'{output_file}_{"_".join(map(str, mutation))}' 
            temp_name = f'temp_{"_".join(map(str, mutation))}_{core}' 
        else: 
            out = output_file 
            temp_name = f'temp_{core}' 

        def pka_iterator(core): 
            
            start = self.trajectory_slices[core][0] 
            end = self.trajectory_slices[core][1]

            for index, ts in enumerate(self.universe.trajectory[start:end]):
                #Check if chainID is empty or not, if so default chain A for the whole system.
                if '' in self.universe._topology.chainIDs.values:
                    print('Your topology file contains no chain identity. Will add chain A for your whole system by default')
                    self.universe._topology.chainIDs.values = np.full(len(self.universe._topology.chainIDs.values),'A',dtype=str)
                with mda.Writer(f'{temp_name}.pdb') as w:
                    w.write(self.universe)
                edit_pdb(f'{temp_name}.pdb')
                if mutation != None: 
                    mutate(temp_name, mutation) 
                run.single(f'{temp_name}.pdb')
                time = ts.time
                #Writing pKa csv
                header = ','.join(extract_data(f'{temp_name}.pka', chain=chain, time=time)[0][0]) 
                data = ','.join(extract_data(f'{temp_name}.pka', chain=chain, time=time)[0][1]).replace('*', '')
                if index == 0 and core == 0: 
                    f = open(f'{out}_pka.csv', "w")
                    f.write(header+'\n')
                    f.write(data+'\n')
                    f.close() 
                else: 
                    f = open(f'{out}_pka.csv', "a")
                    f.write(data+'\n')
                    f.close()
                #Writing buridness csv if extract_surface_data set to true.
                if extract_surface_data == True: 
                    header = ','.join(extract_data(f'{temp_name}.pka', chain=chain, time=time)[1][0]) 
                    data = ','.join(extract_data(f'{temp_name}.pka', chain=chain, time=time)[1][1]).replace('*', '') 
                    if index == 0 and core == 0: 
                        f = open(f'{out}_surf.csv', "w")
                        f.write(header+'\n')
                        f.write(data+'\n')
                        f.close() 
                    else: 
                        f = open(f'{out}_surf.csv', "a")
                        f.write(data+'\n')
                        f.close() 
            
            os.remove(f'{temp_name}.pdb') 
            os.remove(f'{temp_name}.pka')
        
        pka_iterator(core) 
    #Sorting pKa data, if extract_surface_data is set to true, buridness data is sorted also.
    def sort_data(self,output_file, extract_surface_data, mutation):
        if type(mutation) == int: 
            out = f'{output_file}_{mutation}' 
        elif type(mutation) == list: 
            out = f'{output_file}_{"_".join(map(str, mutation))}' 
        else: 
            out = output_file 
        df_pka = pd.read_csv(f'{out}_pka.csv')
        df_pka = df_pka.sort_values('Time [ps]')
        df_pka.to_csv(f'{out}_pka.csv', index=False)
        if extract_surface_data == True:
            df_surf = pd.read_csv(f'{out}_surf.csv')
            df_surf = df_surf.sort_values('Time [ps]')
            df_surf.to_csv(f'{out}_surf.csv', index=False)
            
    def loop_function(self, output_file, index, extract_surface_data, chain, mutation): 
        self.calculate_pka(output_file, extract_surface_data=extract_surface_data, chain=chain, mutation=mutation, core=index) 
    
    def run(self, output_file, extract_surface_data, chain, mutation): 
        pool = mp.Pool(self.cpu_core_number)
        # Create jobs
        jobs = []
        for index, item in enumerate(self.trajectory_slices):
            # Create asynchronous jobs that will be submitted once a processor is ready
            job = pool.apply_async(self.loop_function, args=(output_file, index, extract_surface_data, chain, mutation,))
            jobs.append(job)
        # Submit jobs
        results = [job.get() for job in jobs]
        pool.close()
        self.sort_data(output_file, extract_surface_data, mutation) #Sorting the data only once after all calculations are done, rather than at the end of each job.
