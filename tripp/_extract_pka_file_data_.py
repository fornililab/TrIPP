import numpy as np 

def extract_pka_buriedness_data(file, time):
    """
    Extracts pKa and buriedness data from a PROPKA output file.
    Parameters
    ----------
    file: str
        The path to the PROPKA temporary output file.
    time: str
        The time point for which the data is extracted, used as a key in the returned dictionary.
    Returns
    -------
    data: dict
        A dictionary containing the following keys:
        - 'residue_identifier_list': A list of residue identifiers in the format 'RESID:CHAIN'.
        - 'pka_list': A list of pKa values corresponding to the residues.
        - 'buriedness_list': A list of buriedness values corresponding to the residues.
        The three lists are the same order so that later can be used to match the pKa and
        buriedness values to the residue identifiers.
    """
    compatible_resnames = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 
                           'GLY', 'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 
                           'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL', 'N+', 'C-'] # N+ and C- added for terminal in propka file
    pkafile = open(file, 'r')
    residue_identifier_list = []
    pka_list = [] 
    buriedness_list = []
    for line in pkafile: 
        line_processed = line.rstrip() 
        line_list = line_processed.strip().split()
        # identify the line containing the pka, only considering line with amino acids (compatible_resnames).
        if len(line_list) > 15 and line_list[0] in compatible_resnames:
            residue_identifier_list.append(line_list[0]+line_list[1]+':'+line_list[2])
            pka_list.append(line_list[3])
            buriedness_list.append(line_list[4])

    data = {time: {'residue_identifier_list':np.array(residue_identifier_list),
                  'pka_list':np.array(pka_list),
                  'buriedness_list':np.array(buriedness_list)}}
    return data