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

def create_clustering_matrix(universe, pka_df, residues, include_distances=True): 


    """
    Function that generates the clustering_matrix. The clustering_matrix 
    is a 2D numpy array, where columns correspond to normalized pKa 
    values, extracted from pka_df, and normalized interresidue 
    distances, calculated using the calculate_charge_center_distances 
    function. Normalization is done using the z-score. 
    """
    

    def determine_charge_center(universe, residue_index): 
        selection = universe.residues[residue_index-1] 
        residue_type = selection.resname

        #Charge center is determined as in PROPKA3 
        if residue_type in ['ARG', 'ARGN', 'CARG', 'NARG']: 
            atom_coordinates = universe.select_atoms(f'resid {residue_index} and name CZ').positions 
            charge_center = atom_coordinates[0] 
            residue_identifier = f'ARG{residue_index}' 

        elif residue_type in ['ASP', 'ASPH', 'ASPP', 'CASF', 'CASP', 'NASP', 'ASF', 'ASH']: 
            atom_coordinates = universe.select_atoms(f'resid {residue_index} and name OD1 OD2').positions 
            charge_center = np.mean(atom_coordinates, axis=0) 
            residue_identifier = f'ASP{residue_index}' 

        elif residue_type in ['CYS', 'CCYS', 'CCYX', 'CYS1', 'CYS2', 'CYSH', 'NCYS', 'NCYX', 'CYM', 'CYN', 'CYX']: 
            atom_coordinates = universe.select_atoms(f'resid {residue_index} and name SG').positions 
            charge_center = atom_coordinates[0] 
            residue_identifier = f'CYS{residue_index}'
        
        elif residue_type in ['GLU', 'CGLU', 'GLUH', 'GLUP', 'NGLU', 'PGLU', 'GLH']: 
            atom_coordinates = universe.select_atoms(f'resid {residue_index} and name OE1 OE2').positions 
            charge_center = np.mean(atom_coordinates, axis=0) 
            residue_identifier = f'GLU{residue_index}' 

        elif residue_type in ['HIS', 'CHID', 'CHIE', 'CHIP', 'HIS1', 'HIS2', 'HISA', 'HISB', 'HISD', 'HISE', 'HISH', 'NHID', 'NHIE', 'NHIP', 'HID', 'HIE', 'HIP', 'HSD', 'HSE', 'HSP']: 
            atom_coordinates = universe.select_atoms(f'resid {residue_index} and name CG CD2 ND1 CE1 NE2').positions 
            charge_center = np.mean(atom_coordinates, axis=0) 
            residue_identifier = f'HIS{residue_index}' 

        elif residue_type in ['LYS', 'CLYS', 'LYSH', 'NLYS', 'LYN', 'LSN']: 
            atom_coordinates = universe.select_atoms(f'resid {residue_index} and name NZ').positions 
            charge_center = atom_coordinates[0] 
            residue_identifier = f'LYS{residue_index}' 

        elif residue_type in ['TYR', 'CTYR', 'NTYR']: 
            atom_coordinates = universe.select_atoms(f'resid {residue_index} and OH' ).positions 
            charge_center = atom_coordinates[0] 
            residue_identifier = f'TYR{residue_index}' 
        
        else: 
            print(f'Residue type {residue_type} not recognized by TrIPP to be able to change protonation state') 
        
        return charge_center, residue_identifier  
    
    def normalize(arr): 
        arr = np.array(arr) 
        mean = np.mean(arr) 
        std = np.std(arr) 
        normalized_array = (arr-mean)/std 
        return normalized_array 
    
    pka = [] 
    times = [] 
    frames = [] 
    if include_distances == True: 

        d = [] 
        for ts in universe.trajectory: 

            times.append([ts.time]) 
            frames.append([ts.frame]) 
            positions = [] 
            pkas = [] 
            for residue_index in residues: 
                charge_center, residue_identifier = determine_charge_center(universe, residue_index) 
                positions.append(charge_center) 
                pka_residue = pka_df.loc[ts.time, residue_identifier] 
                pkas.append(pka_residue) 
            distances = calculate_charge_center_distances(positions) 
            d.append(distances) 
            pka.append(pkas) 
        
        clustering_matrix = np.concatenate((normalize(d), normalize(pka)), axis=1) 
    
    elif include_distances == False: 
        
        for ts in universe.trajectory: 

            times.append([ts.time]) 
            frames.append([ts.frame]) 
            pkas = [] 
            for residue_index in residues: 
                charge_center, residue_identifier = determine_charge_center(universe, residue_index) 
                pka_residue = pka_df.loc[ts.time, residue_identifier] 
                pkas.append(pka_residue) 
            pka.append(pkas) 
        
        clustering_matrix = normalize(pka) 
    
    times = np.array(times) 
    frames = np.array(frames) 

    return clustering_matrix, times, frames 