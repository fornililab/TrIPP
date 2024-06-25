import pandas as pd 
from tripp._detect_disulphide_bonds_ import detect_disulphide_bonds 
import os
def sort_pka_df(cores, topology_file, output_file, extract_surface_data, mutation, disulphide_bond_detection, universe, chain):
        
    if type(mutation) == int: 
        out = f'{output_file}_{mutation}' 
    elif type(mutation) == list: 
        out = f'{output_file}_{"_".join(map(str, mutation))}' 
    else: 
        out = output_file 
    if os.path.isfile(f'{out}_pka.csv'):
        os.remove(f'{out}_pka.csv')
    if os.path.isfile(f'{out}_surf.csv'):
        os.remove(f'{out}_surf.csv')
        
    for core in range(cores):
        with open(f'temp_pka_worker{core}.csv','r') as input, open(f'{out}_pka.csv','a') as output:
            for line in input:
                output.write(line)
        os.remove(f'temp_pka_worker{core}.csv')
        if extract_surface_data == True:
            with open(f'temp_surf_worker{core}.csv','r') as input,  open(f'{out}_surf.csv','a') as output:
                for line in input:
                    output.write(line)
            os.remove(f'temp_surf_worker{core}.csv')
            
    if disulphide_bond_detection == True:
        df_pka = pd.read_csv(f'{out}_pka.csv') 
        disulphide_cysteines = detect_disulphide_bonds(topology_file) 
        df_pka.drop(disulphide_cysteines, axis=1, inplace=True)
        df_pka.to_csv(f'{out}_pka.csv', index=False)
        if extract_surface_data == True:
            df_surf = pd.read_csv(f'{out}_surf.csv')
            df_surf.drop(disulphide_cysteines, axis=1, inplace=True)
            df_surf.to_csv(f'{out}_surf.csv', index=False)