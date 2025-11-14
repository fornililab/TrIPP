"""
    This file is part of the TrIPP software
    (https://github.com/fornililab/TrIPP).
    Copyright (c) Christos Matsingos, Ka Fu Man and Arianna Fornili.

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
    
def check_resname_HETATM(non_protein_ag, predictor):
    """
    Checks if the resnames in the non-protein atoms universe are compatible with PROPKA.
    Parameters
    ----------
    non_protein_ag : MDAnalysis.AtomGroup
        The AtomGroup (MDAnalysis selection syntax) for residues with resnames not found in either the built-in or user-defined conversion dictionaries.
    predictor : str
        The pKa predictor to be used.
    Raises
    ------
    NameError
        If there are resnames that are not recognized by PROPKA, and their record type
        is not 'HETATM', an exception will be raised.
    """
    if predictor == 'propka':
        predictor_name = predictor.upper()
        compatible_resnames = np.unique(list(corrected_amino_acids.values()))
    elif predictor in ['pkai', 'pkai+']:
        predictor_name = predictor[0] + predictor[1:].upper()
        compatible_resnames = np.unique(list(corrected_amino_acids.values())).tolist() + ['NTR','CTR'] # Add NTR and CTR for N- and C-termini as compatible resnames for pKAI and pKAI+
        
    incorrect_resnames=[]
    for resname, record_types in zip(non_protein_ag.residues.resnames,
                                     non_protein_ag.residues.record_types):
        if resname not in compatible_resnames and np.all(record_types != 'HETATM'):
            incorrect_resnames.append(f'{resname}')
    if len(incorrect_resnames) > 0:
        raise NameError(f"""Your system still contains resname(s) not recognised by {predictor_name}: {', '.join(incorrect_resnames)}
For amino acids, please use the custom_resname_correction argument to indicate valid {predictor_name} resnames for the residues indicated above.
For ligands, please use the hetatm_resname argument to convert the record type of the ligand to HETATM.""")
    else:
        logger.info('Resname and record type check passed.\n')
        
def check_terminal_oxygens(universe, correct_terminal_oxygens):
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
        if np.isin(correct_terminal_oxygens[0], name).any() and np.isin(correct_terminal_oxygens[1], name).any():
            ag = universe.select_atoms(f'resindex {index}')
            terminals.append(f'{ag.residues.resnames[0]}{ag.residues.resids[0]}')
    if len(terminals) > 0:
        logger.info(f"""Terminal oxygen check passed, involving: 
{', '.join(terminals)}\n""")
    else:
        raise NameError(f'No terminal oxygen atom named {correct_terminal_oxygens[0]} and {correct_terminal_oxygens[1]} was found, please either modify your topology_file or use the custom_terminal_oxygens argument')