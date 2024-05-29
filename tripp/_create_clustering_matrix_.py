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

import MDAnalysis as mda 
import numpy as np 
from tripp._calculate_charge_center_distances_ import calculate_charge_center_distances 
from tripp._determine_charge_center_ import determine_charge_center 
from scipy.stats import zscore 
from tripp._create_mda_universe_ import create_mda_universe 

def create_clustering_matrix(trajectory_file, topology_file, pka_df, residues, include_distances=True): 

    """
    Function that generates the clustering_matrix. The clustering_matrix 
    is a 2D numpy array, where columns correspond to normalized pKa 
    values, extracted from pka_df, and normalized interresidue 
    distances, calculated using the calculate_charge_center_distances 
    function. Normalization is done using the z-score. 
    """

    if type(trajectory_file) == str: 
        trajectory_dict = {'Unnamed Trajectory' : create_mda_universe(topology_file=topology_file, trajectory_file=trajectory_file)} 
    
    elif type(trajectory_file) == dict: 
        trajectory_dict = {} 
        for key in trajectory_file.keys(): 
            traj = trajectory_file[key] 
            trajectory_dict[key] = create_mda_universe(topology_file=topology_file, trajectory_file=traj) 
    
    def extract_pka_distances(trajectory_name, universe, df_traj): 

        pka = [] 
        times = [] 
        frames = [] 
        trajectory_names = [] 

        if include_distances == True: 
            d = [] 
            for ts in universe.trajectory: 
                times.append([ts.time]) 
                frames.append([ts.frame]) 
                trajectory_names.append(trajectory_name) 
                positions = [] 
                pkas = [] 
                for residue_index in residues: 
                    charge_center, residue_identifier = determine_charge_center(universe, residue_index) 
                    positions.append(charge_center) 
                    pka_residue = df_traj.loc[ts.time, residue_identifier] 
                    pkas.append(pka_residue) 
                distances = calculate_charge_center_distances(positions) 
                d.append(distances) 
                pka.append(np.array(pkas)) 
            
            pka = np.array(pka) 
            d = np.array(d) 
            
            clustering_matrix = np.concatenate((zscore(d, axis=None), zscore(pka, axis=None)), axis=1) 
        
        elif include_distances == False: 
            for ts in universe.trajectory: 
                times.append([ts.time]) 
                frames.append([ts.frame]) 
                trajectory_names.append(trajectory_name) 
                pkas = [] 
                for residue_index in residues: 
                    charge_center, residue_identifier = determine_charge_center(universe, residue_index) 
                    pka_residue = df_traj.loc[ts.time, residue_identifier] 
                    pkas.append(pka_residue) 
                pka.append(np.array(pkas)) 
            
            pka = np.array(pka) 
            
            clustering_matrix = zscore(pka, axis=None) 
        
        times = np.array(times) 
        frames = np.array(frames) 
        trajectory_names = np.array(trajectory_names) 

        return clustering_matrix, times, frames, trajectory_names
    
    for index, key in enumerate(trajectory_dict.keys()): 
        df_traj = pka_df[pka_df['Trajectories']==key]
        if index == 0: 
            clustering_matrix, times, frames, trajectory_names = extract_pka_distances(key, trajectory_dict[key], df_traj) 
        
        else: 
            partial_clustering_matrix, partial_times, partial_frames, partial_trajectory_names = extract_pka_distances(key, trajectory_dict[key], df_traj) 
            clustering_matrix = np.concatenate((clustering_matrix, partial_clustering_matrix), axis=0) 
            times = np.concatenate((times, partial_times), axis=0) 
            frames = np.concatenate((frames, partial_frames), axis=0) 
            trajectory_names = np.concatenate((trajectory_names, partial_trajectory_names), axis=0) 

    return clustering_matrix, times, frames, trajectory_names, trajectory_dict 