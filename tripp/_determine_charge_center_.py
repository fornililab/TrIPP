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

def determine_charge_center(universe, selection):
    """
    Determines the charge center of a residue in a MDAnalysis universe.
    Parameters
    ----------
    universe: MDAnalysis.universe
    
    selection: str
        A string (MDAnalysis selection syntax) to select the residue for
        charge center calculation.
    Returns
    -------
    charge_center: np.ndarray
        The charge center coordinates for the selected residue.
    residue_identifier: str
        A string identifier for the residue, formatted as 'RESID:CHAINID'.
        The following residue types are recognized:
        - ARG, ASP, CYS, GLU, HIS, LYS, TYR
    Raises
    ------
    Exception
        If the residue type is not recognized by TrIPP, an exception is raised.
    """
    ag = universe.select_atoms(selection)
    residue_type = ag.residues.resnames
    resid = ag.residues.resids[0]
    resindex = ag.residues.resindices[0]
    chain = ag.atoms.chainIDs[0]

    # Charge center is determined as in PROPKA3
    if residue_type in ['ARG', 'ARGN', 'CARG', 'NARG']:
        atom_coordinates = universe.select_atoms(f'resindex {resindex} and name CZ').positions
        charge_center = atom_coordinates[0]
        residue_identifier = f'ARG{resid}:{chain}'

    elif residue_type in ['ASP', 'ASPH', 'ASPP', 'CASF', 'CASP', 'NASP', 'ASF', 'ASH']:
        atom_coordinates = universe.select_atoms(f'resindex {resindex} and name OD1 OD2').positions
        charge_center = np.mean(atom_coordinates, axis=0)
        residue_identifier = f'ASP{resid}:{chain}'

    elif residue_type in ['CYS', 'CCYS', 'CCYX', 'CYS1', 'CYS2', 'CYSH', 'NCYS', 'NCYX', 'CYM', 'CYN', 'CYX']:
        atom_coordinates = universe.select_atoms(f'resindex {resindex} and name SG').positions
        charge_center = atom_coordinates[0]
        residue_identifier = f'CYS{resid}:{chain}'

    elif residue_type in ['GLU', 'CGLU', 'GLUH', 'GLUP', 'NGLU', 'PGLU', 'GLH']:
        atom_coordinates = universe.select_atoms(f'resindex {resindex} and name OE1 OE2').positions
        charge_center = np.mean(atom_coordinates, axis=0) 
        residue_identifier = f'GLU{resid}:{chain}'

    elif residue_type in ['HIS', 'CHID', 'CHIE', 'CHIP', 'HIS1', 'HIS2', 'HISA', 'HISB', 'HISD', 'HISE', 'HISH', 'NHID', 'NHIE', 'NHIP', 'HID', 'HIE', 'HIP', 'HSD', 'HSE', 'HSP']:
        atom_coordinates = universe.select_atoms(f'resindex {resindex} and name CG CD2 ND1 CE1 NE2').positions
        charge_center = np.mean(atom_coordinates, axis=0)
        residue_identifier = f'HIS{resid}:{chain}'

    elif residue_type in ['LYS', 'CLYS', 'LYSH', 'NLYS', 'LYN', 'LSN']:
        atom_coordinates = universe.select_atoms(f'resindex {resindex} and name NZ').positions
        charge_center = atom_coordinates[0]
        residue_identifier = f'LYS{resid}:{chain}'

    elif residue_type in ['TYR', 'CTYR', 'NTYR']:
        atom_coordinates = universe.select_atoms(f'resindex {resindex} and name OH').positions
        charge_center = atom_coordinates[0]
        residue_identifier = f'TYR{resid}:{chain}'

    else:
        raise Exception(f'Residue {residue_type}{resid}:{chain} is not recognized by TrIPP: unable to determine charge centre')

    if np.any(np.isnan(charge_center)):
        raise Exception(f'Unable to determine charge center for {residue_type}{resid}:{chain} due to missing atoms in the residue.')

    return charge_center, residue_identifier