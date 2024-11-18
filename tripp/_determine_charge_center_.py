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

def determine_charge_center(universe, resid):
    selection = universe.select_atoms(f'resid {resid}')
    residue_type = selection.residues.resnames

    # Charge center is determined as in PROPKA3
    if residue_type in ['ARG', 'ARGN', 'CARG', 'NARG']:
        atom_coordinates = universe.select_atoms(f'resid {resid} and name CZ').positions
        charge_center = atom_coordinates[0]
        residue_identifier = f'ARG{resid}'

    elif residue_type in ['ASP', 'ASPH', 'ASPP', 'CASF', 'CASP', 'NASP', 'ASF', 'ASH']:
        atom_coordinates = universe.select_atoms(f'resid {resid} and name OD1 OD2').positions
        charge_center = np.mean(atom_coordinates, axis=0)
        residue_identifier = f'ASP{resid}'

    elif residue_type in ['CYS', 'CCYS', 'CCYX', 'CYS1', 'CYS2', 'CYSH', 'NCYS', 'NCYX', 'CYM', 'CYN', 'CYX']:
        atom_coordinates = universe.select_atoms(f'resid {resid} and name SG').positions
        charge_center = atom_coordinates[0]
        residue_identifier = f'CYS{resid}'

    elif residue_type in ['GLU', 'CGLU', 'GLUH', 'GLUP', 'NGLU', 'PGLU', 'GLH']:
        atom_coordinates = universe.select_atoms(f'resid {resid} and name OE1 OE2').positions
        charge_center = np.mean(atom_coordinates, axis=0) 
        residue_identifier = f'GLU{resid}'

    elif residue_type in ['HIS', 'CHID', 'CHIE', 'CHIP', 'HIS1', 'HIS2', 'HISA', 'HISB', 'HISD', 'HISE', 'HISH', 'NHID', 'NHIE', 'NHIP', 'HID', 'HIE', 'HIP', 'HSD', 'HSE', 'HSP']:
        atom_coordinates = universe.select_atoms(f'resid {resid} and name CG CD2 ND1 CE1 NE2').positions
        charge_center = np.mean(atom_coordinates, axis=0)
        residue_identifier = f'HIS{resid}'

    elif residue_type in ['LYS', 'CLYS', 'LYSH', 'NLYS', 'LYN', 'LSN']:
        atom_coordinates = universe.select_atoms(f'resid {resid} and name NZ').positions
        charge_center = atom_coordinates[0]
        residue_identifier = f'LYS{resid}'

    elif residue_type in ['TYR', 'CTYR', 'NTYR']:
        atom_coordinates = universe.select_atoms(f'resid {resid} and name OH').positions
        charge_center = atom_coordinates[0]
        residue_identifier = f'TYR{resid}'

    else:
        raise Exception(f'Residue {residue_type}{resid} is not recognized by TrIPP to be able to determine charge centre')

    return charge_center, residue_identifier