import pandas as pd 
from tripp._detect_disulphide_bonds_ import detect_disulphide_bonds 
from tripp._missing_data_recovery_ import missing_data_recovery 

def sort_pka_df(topology_file, output_file, extract_surface_data, mutation, disulphide_bond_detection, universe, chain):
        
    if type(mutation) == int: 
        out = f'{output_file}_{mutation}' 
    elif type(mutation) == list: 
        out = f'{output_file}_{"_".join(map(str, mutation))}' 
    else: 
        out = output_file 
    
    missing_data_recovery(pka_file=f'{out}_pka.csv', surf_file=None, universe=universe, mutation=mutation, chain=chain) 
    df_pka = pd.read_csv(f'{out}_pka.csv') 
    df_pka = df_pka.sort_values('Time [ps]') 
    
    if disulphide_bond_detection == True: 
        disulphide_cysteines = detect_disulphide_bonds(topology_file) 
        df_pka.drop(disulphide_cysteines, axis=1, inplace=True) 
    
    df_pka.to_csv(f'{out}_pka.csv', index=False)

    if extract_surface_data == True: 
        missing_data_recovery(pka_file=None, surf_file=f'{out}_surf.csv', universe=universe, mutation=mutation, chain=chain) 
        df_surf = pd.read_csv(f'{out}_surf.csv') 
        df_surf = df_surf.sort_values('Time [ps]')
        if disulphide_bond_detection == True: 
            df_surf.drop(disulphide_cysteines, axis=1, inplace=True) 
        df_surf.to_csv(f'{out}_surf.csv', index=False)