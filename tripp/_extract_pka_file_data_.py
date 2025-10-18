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

import numpy as np 

def propka_pka_buriedness_data_parser(molecule, time):
    """
    Parses the pKa and buriedness data from the PROPKA output files.
    Parameters
    ----------
    molecule: propka.molecule.Molecule object
        The PROPKA molecule object containing the pKa calculation results.
    time: int
        The time corresponding to the pKa calculation.
    Returns
    -------
    data: dict
        A dictionary containing residue identifiers, pKa values and buriedness
        values for the titratable groups in the molecule.
    """
    groups = molecule.conformations[molecule.conformation_names[0]].get_titratable_groups()

    residue_identifier_list = []
    pka_list = []
    buriedness_list = []
    for group in groups:
        if group.atom.name == 'N':
            residue_name = 'N+'
        elif group.atom.name == 'OXT':
            residue_name = 'C-'
        else:
            residue_name = group.atom.res_name
        residue_id = str(group.atom.res_num)
        residue_identifier_list.append(residue_name + residue_id + ':' + group.atom.chain_id)
        pka_list.append(round(group.pka_value,2))
        buriedness_list.append(round(group.buried*100))
    data = {time: {'residue_identifier_list':np.array(residue_identifier_list),
                  'pka_list':np.array(pka_list),
                  'buriedness_list':np.array(buriedness_list)}}
    return data

def pypka_pka_data_parser(tit_result, time):
    processed_line = []
    for line in str(tit_result).split('\n'):
        tmp = []
        if 'Not In Range' in line or 'Predicted Isoelectric Point' in line:
            continue
        for element in line.split(' '):
            if element != '':
                tmp.append(element)
        if tmp:
            processed_line.append(tmp)
    result_arr = np.array(processed_line[1:]) # col0:chain, col1:resid, col2:resname, col3:pka
    result_arr[:,2] = list(map(lambda x: x.replace('NTR','N+').replace('CTR','C-'), result_arr[:,2])) #Convert NTR and CTR to N+ and C-
    residue_identifier_list = np.char.array(result_arr[:,2]) + np.char.array(result_arr[:,1]) + ':' + np.char.array(result_arr[:,0]) 
    pka_arr = result_arr[:,3]
    data = {time: {'residue_identifier_list':np.array(residue_identifier_list),
                  'pka_list':pka_arr,
                  'buriedness_list':None}}
    return data