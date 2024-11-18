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
def input_check(universe):
    
    """ 
    Function that checks all the naming of residues has been modified to 
    be PROPKA compatible. If not, user will have to add the correction 
    manually in _correction_dictionary_.py. If the system contains 
    non-protein atoms and user named its record type to 'ATOM', PROPKA 
    would not recognise them. We recommend changing them to 'HETATM'.
    """ 
    
    non_protein_ag = universe.select_atoms('not protein')
    
    check_terminal_oxygens(universe)
    
    incorrect_resnames = check_resname_HETATM(non_protein_ag)
    if len(incorrect_resnames) > 0:
        raise NameError(f"""Your system still contain resname not recognisable by PROPKA: {', '.join(incorrect_resnames)}
If it is an amino acid, please use custom_resname_correction argument to add the corrections for the resname(s) identified above.
If it is a ligand, please use hetatm_resid argument to convert the record type of the ligand to HETATM.""")
    else:
        logger.info('PROPKA compatibility checked.')
    
def check_resname_HETATM(non_protein_ag):
    compatible_resnames = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 
                           'GLY', 'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 
                           'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL']
    incorrect_resnames=[]
    for resname, resid, record_types in zip(non_protein_ag.residues.resnames,
                                            non_protein_ag.residues.resids,
                                            non_protein_ag.residues.record_types):
        if resname not in compatible_resnames and np.all(record_types == 'HETATM'):
            logger.info(f'Your system contains non protein residues {resname}{resid}, but the record type is HETATM so proceeding.')
        else:
            incorrect_resnames.append(f'{resname}{resid}')
    return incorrect_resnames

def check_terminal_oxygens(universe):
    last_residue_names = universe.residues.names[-1]
    if 'O' not in last_residue_names or 'OXT' not in last_residue_names:
        raise NameError('Terminal oxygens is not named O and OXT, please either modify from your topology_file or via argument terminal_oxygens')
            
    
#     non_protein_ag = universe.select_atoms('not protein')
#     if len(non_protein_ag) > 0:
#         incorrect_record_types = []
#         for record_type,resid in zip(non_protein_ag.residues.record_types,non_protein_ag.residues.resids):
#             if 'ATOM' in record_type:
#                 incorrect_record_types.append(resid)
        
#         if len(incorrect_record_types) >0:
#             print(f"""The record type of the following resids contains \'ATOM\':\n
# {incorrect_record_types}\n
# These will not be recognised by PROPKA, TrIPP will proceed without any modifications.\n
# We recommend changing your topology file of their record type to \'HETATM\' if you want to take them into account during pKa prediction.""")
#         else:
#             print('The record type for all non-protein atoms are \'HETATM\'.')