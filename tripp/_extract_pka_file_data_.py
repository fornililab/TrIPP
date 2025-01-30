import numpy as np 

def extract_pka_buriedness_data(file, chains, time):
    compatible_resnames = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 
                           'GLY', 'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 
                           'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL', 'N+', 'C-'] # N+ and C- added for terminal in propka file
    chains = list(chains)
    pkafile = open(file, 'r')
    resname_list = []
    pka_list = [] 
    buriedness_list = []
    chain_list = []
    for line in pkafile: 
        line_processed = line.rstrip() 
        line_list = line_processed.strip().split()
        # make sure each chain is looped through and save to different file
        for chain in chains:
            # identify the line containing the pka, only considering line with amino acids (compatible_resnames).
            if len(line_list) > 15 and line_list[2] == chain and line_list[0] in compatible_resnames:
                resname_list.append(line_list[0]+line_list[1])
                pka_list.append(line_list[3])
                buriedness_list.append(line_list[4])
                chain_list.append(line_list[2])
    data = {}
    for chain in chains:
        data[chain] = {'Time [ps]':time,
                    'resname_list':np.array(resname_list)[np.array(chain_list) == chain],
                    'pka_list':np.array(pka_list)[np.array(chain_list) == chain],
                    'buriedness_list':np.array(buriedness_list)[np.array(chain_list) == chain]}
    return data