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
import propka
def propka_pka_buriedness_data_parser(molecule, time):
    """
    Parses the pKa and buriedness data from the PROPKA molecule object.
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
        residue_name = group.residue_type
        if not propka.group.is_protein_group(group.parameters,group.atom): # Use PROPKA's internal function to check if the group is a protein group
            continue                                                       # Skip non-protein groups
        residue_id = str(group.atom.res_num)
        residue_identifier_list.append(residue_name + residue_id + ':' + group.atom.chain_id)
        pka_list.append(round(group.pka_value,2))
        buriedness_list.append(round(group.buried*100))
    data = {time: {'residue_identifier_list':np.array(residue_identifier_list),
                  'pka_list':np.array(pka_list),
                  'buriedness_list':np.array(buriedness_list)}}
    return data

def pkai_pka_data_parser(result, time):
    """
    Parses the pKa result from the pkai/pkai+ output files.
    Parameters
    ----------
    result: list
        The pKa result from pkai/pkai+ as a list of lists.
        col0:chain, col1:resid, col2:resname, col3:pka
    time: int
        The time corresponding to the pKa calculation.
    Returns
    -------
    data: dict
        A dictionary containing residue identifiers, pKa values and buriedness
        values for the titratable groups in the molecule.
    """
    result_arr = np.array(result) # col0:chain, col1:resid, col2:resname, col3:pka
    residue_identifier_list = np.char.array(result_arr[:,2]) + np.char.array(result_arr[:,1]) + ':' + np.char.array(result_arr[:,0])
    pka_arr = result_arr[:,3]
    data = {time: {'residue_identifier_list':np.array(residue_identifier_list),
                  'pka_list':pka_arr,
                  'buriedness_list':None}}
    return data