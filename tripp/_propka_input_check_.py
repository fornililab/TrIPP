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
import logging
import numpy as np


logger = logging.getLogger(__name__)
    
def check_resname_HETATM(non_protein_ag):
    """
    Check if the resnames in the non-protein atoms universe are compatible with PROPKA.
    Parameters
    ----------
    non_protein_ag : MDAnalysis.AtomGroup
        The AtomGroup containing non-protein atoms, such as ligands or water.
        Note that non-protein selection was defined by MDAnalysis selection algebra.
    Raises
    ------
    NameError
        If there are resnames that are not recognized by PROPKA, and their record type
        is not 'HETATM', an exception will be raised.
    """
    compatible_resnames = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 
                           'GLY', 'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 
                           'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL']
    incorrect_resnames=[]
    for resname, record_types in zip(non_protein_ag.residues.resnames,
                                     non_protein_ag.residues.record_types):
        if resname not in compatible_resnames and np.all(record_types != 'HETATM'):
            incorrect_resnames.append(f'{resname}')
    if len(incorrect_resnames) > 0:
        raise NameError(f"""Your system still contains resname that are not recognisable by PROPKA: {', '.join(incorrect_resnames)}
If it is an amino acid, please use custom_resname_correction argument to add the corrections for the resname(s) identified above.
If it is a ligand, please use hetatm_resname argument to convert the record type of the ligand to HETATM.""")
    else:
        logger.info('Resname and record type check passed.\n')
        
def check_terminal_oxygens(universe):
    """
    Check if the terminal oxygens in the universe are named 'O' and 'OXT
    Parameters
    ----------
    universe : MDAnalysis.Universe
        The MDAnalysis universe containing the system to be checked.
    Raises
    ------
    NameError
        If no terminal oxygens are found with the names 'O' and 'OXT',
        an exception will be raised.
    """
    terminals = []
    for index, name in zip(universe.residues.resindices, universe.residues.names):
        if np.isin('O', name).any() and np.isin('OXT', name).any():
            ag = universe.select_atoms(f'resindex {index}')
            terminals.append(f'{ag.residues.resnames[0]}{ag.residues.resids[0]}')
    if len(terminals) > 0:
        logger.info(f"""Terminal oxygen check passed, involving: 
{', '.join(terminals)}\n""")
    else:
        raise NameError('No terminal oxygen named O and OXT, please either modify from your topology_file or via custom_terminal_oxygens argument')
