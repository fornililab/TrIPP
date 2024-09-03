import pandas as pd 
from tripp._detect_disulphide_bonds_ import detect_disulphide_bonds 
import os
def sort_pka_df(cores, topology_file, output_directory, output_prefix, extract_buridness_data, mutation, disulphide_bond_detection):
        
    if os.path.isfile(f'{output_directory}/{output_prefix}_pka.csv'):
        os.remove(f'{output_directory}/{output_prefix}_pka.csv')
    if os.path.isfile(f'{output_directory}/{output_prefix}_buridness.csv'):
        os.remove(f'{output_directory}/{output_prefix}_buridness.csv')
        
    for core in range(cores):
        with open(f'{output_directory}/.temp_pka_worker{core}.csv','r') as input, open(f'{output_directory}/{output_prefix}_pka.csv','a') as output:
            for line in input:
                output.write(line)
        os.remove(f'{output_directory}/.temp_pka_worker{core}.csv')
        
        if extract_buridness_data == True:
            with open(f'{output_directory}/.temp_buridness_worker{core}.csv','r') as input, open(f'{output_directory}/{output_prefix}_buridness.csv','a') as output:
                for line in input:
                    output.write(line)
            os.remove(f'{output_directory}/.temp_buridness_worker{core}.csv')
            
    if disulphide_bond_detection == True:
        df_pka = pd.read_csv(f'{output_directory}/{output_prefix}_pka.csv')
        disulphide_cysteines = detect_disulphide_bonds(topology_file)
        df_pka.drop(disulphide_cysteines, axis=1, inplace=True)
        df_pka.to_csv(f'{output_directory}/{output_prefix}_pka.csv', index=False)
        
        if extract_buridness_data == True:
            df_buridness = pd.read_csv(f'{output_directory}/{output_prefix}_buridness.csv')
            df_buridness.drop(disulphide_cysteines, axis=1, inplace=True)
            df_buridness.to_csv(f'{output_directory}/{output_prefix}_buridness.csv', index=False)