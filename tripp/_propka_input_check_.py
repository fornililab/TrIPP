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
from tripp._correction_dictionary_ import corrected_amino_acids

logger = logging.getLogger(__name__)
    
def check_resname_HETATM(non_protein_ag):
    """
    Checks if the resnames in the non-protein atoms universe are compatible with PROPKA.
    Parameters
    ----------
    non_protein_ag : MDAnalysis.AtomGroup
        The AtomGroup (MDAnalysis selection syntax) for residues with resnames not found in either the built-in or user-defined conversion dictionaries.
        
    Raises
    ------
    NameError
        If there are resnames that are not recognized by PROPKA, and their record type
        is not 'HETATM', an exception will be raised.
    """
    compatible_resnames = np.unique(list(corrected_amino_acids.values()))
    incorrect_resnames=[]
    for resname, record_types in zip(non_protein_ag.residues.resnames,
                                     non_protein_ag.residues.record_types):
        if resname not in compatible_resnames and np.all(record_types != 'HETATM'):
            incorrect_resnames.append(f'{resname}')
    if len(incorrect_resnames) > 0:
        raise NameError(f"""Your system still contains resname(s) not recognised by PROPKA: {', '.join(incorrect_resnames)}
For amino acids, please use the custom_resname_correction argument to indicate valid PROPKA resnames for the residues indicated above.
For ligands, please use the hetatm_resname argument to convert the record type of the ligand to HETATM.""")
    else:
        logger.info('Resname and record type check passed.\n')
        
def check_terminal_oxygens(universe):
    """
    Checks if the C-terminal oxygen atoms in the universe are named 'O' and 'OXT
    Parameters
    ----------
    universe : MDAnalysis.Universe
        The MDAnalysis universe containing the system to be checked.
    Raises
    ------
    NameError
        If no C-terminal oxygen atoms are found with the names 'O' and 'OXT',
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
        raise NameError('No terminal oxygen atom named O and OXT was found, please either modify your topology_file or use the custom_terminal_oxygens argument')