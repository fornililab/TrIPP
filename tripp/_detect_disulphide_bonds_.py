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

def detect_disulphide_bonds(universe): 

    residues = universe.residues
    cysteines = [] 
    sulphur_positions = [] 

    for residue in residues: 
        residue_type = residue.resname
        resindex = residue.resindex

        if residue_type in ['CYS', 'CCYS', 'CCYX', 'CYS1', 'CYS2', 'CYSH', 'NCYS', 'NCYX', 'CYM', 'CYN', 'CYX']: 
            atom = universe.select_atoms(f'resindex {resindex} and name SG')
            atom_coordinates = atom.positions
            resid = atom.resids[0]
            chain = atom.chainIDs[0]
            sulphur_atom = atom_coordinates[0] 
            residue_identifier = f'CYS{resid}:{chain}'
            cysteines.append(residue_identifier) 
            sulphur_positions.append(sulphur_atom) 
    
    cysteines = np.array(cysteines) 
    sulphur_positions = np.array(sulphur_positions) 
    distance_matrix = np.zeros(shape=(len(cysteines), len(cysteines))) 

    for i in range(len(cysteines)): 
        for j in range(len(cysteines)): 
            if i == j: 
                distance_matrix[i, j] = np.inf 
            else: 
                distance_matrix[i, j] = np.sqrt(np.sum(np.square(sulphur_positions[i] - sulphur_positions[j]))) 
    
    disulphide_bond_cutoff = 3 
    distance_mask = distance_matrix < disulphide_bond_cutoff

    cysteine_indices = np.nonzero(distance_mask)[0] 

    disulphide_cysteines = cysteines[cysteine_indices] 
    
    return disulphide_cysteines