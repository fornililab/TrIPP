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

def input_check(universe):
    
    """ 
    Function that checks all the naming of residues has been modified to 
    be PROPKA compatible. If not, user will have to add the correction 
    manually in _correction_dictionary_.py. If the system contains 
    non-protein atoms and user named its record type to 'ATOM', PROPKA 
    would not recognise them. We recommend changing them to 'HETATM'.
    """ 
    
    compatible_resnames = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 
                           'GLY', 'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 
                           'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL']
    protein_ag = universe.select_atoms('protein')
    incorrect_resnames=[]
    for resname,resid in zip(protein_ag.residues.resnames,protein_ag.residues.resids):
        if resname not in compatible_resnames:
            incorrect_resnames.append(f'{resname}{resid}')
    
    if len(incorrect_resnames) > 0:
        print(f"""Amino acids in your system still contain resname not recognisable by PROPKA:\n
{incorrect_resnames}\n
Please modify the _correction_dictionary_.py to add the corrections for the resnames identified above.""")
    else:
        print('All resnames can be recognised by PROPKA.')
    
    non_protein_ag = universe.select_atoms('not protein')
    if len(non_protein_ag) > 0:
        incorrect_record_types = []
        for record_type,resid in zip(non_protein_ag.residues.record_types,non_protein_ag.residues.resids):
            if 'ATOM' in record_type:
                incorrect_record_types.append(resid)
        
        if len(incorrect_record_types) >0:
            print(f"""The record type of the following resids contains \'ATOM\':\n
{incorrect_record_types}\n
These will not be recognised by PROPKA, TrIPP will proceed without any modifications.\n
We recommend changing your topology file of their record type to \'HETATM\' if you want to take them into account during pKa prediction.""")
        else:
            print('The record type for all non-protein atoms are \'HETATM\'.')
    print('PROPKA prediction proceeding...')