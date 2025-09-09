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

def extract_pka_buriedness_data(file, time):
    """
    Extracts pKa and buriedness (buried ratio) values from a PROPKA output file.
    Parameters
    ----------
    file: str
        The path to the PROPKA temporary output file.
    time: str
        The time of the frame processed by PROPKA, used as a key in the returned dictionary.
    Returns
    -------
    data: dict
        A nested dictionary where each time point is mapped to a dictionary containing the following keys:
        - 'residue_identifier_list': an array of residue identifiers in the format 'RESID:CHAIN'.
        - 'pka_list': an array of pKa values corresponding to the residues.
        - 'buriedness_list': an array of buriedness values corresponding to the residues.
        Elements in the three arrays are in the same order (element i from each array refers to the same residue)
    """
    compatible_resnames = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 
                           'GLY', 'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 
                           'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL', 'N+', 'C-'] # N+ and C- added for N- and C-termini in the PROPKA file
    pkafile = open(file, 'r')
    residue_identifier_list = []
    pka_list = [] 
    buriedness_list = []
    for line in pkafile: 
        line_processed = line.rstrip() 
        line_list = line_processed.strip().split()
        # identify the line containing the pka, only considering lines with amino acids (compatible_resnames).
        if len(line_list) > 15 and line_list[0] in compatible_resnames:
            residue_identifier_list.append(line_list[0]+line_list[1]+':'+line_list[2])
            pka_list.append(line_list[3])
            buriedness_list.append(line_list[4])

    data = {time: {'residue_identifier_list':np.array(residue_identifier_list),
                  'pka_list':np.array(pka_list),
                  'buriedness_list':np.array(buriedness_list)}}
    return data