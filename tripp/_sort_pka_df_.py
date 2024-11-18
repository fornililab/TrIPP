import pandas as pd 
import os
import glob
   

def sort_pka_df(output_directory, output_prefix, pka_dfs, disulphide_cysteines_list):
    sorted_pka_df = pd.concat(pka_dfs).astype(float).sort_values('Time [ps]')
    if len(disulphide_cysteines_list) > 0:
        sorted_pka_df.drop(disulphide_cysteines_list, axis=1, inplace=True)
    sorted_pka_df.to_csv(f'{output_directory}/{output_prefix}_pka.csv', index=False)
    [os.remove(tmp_pka_csv) for tmp_pka_csv in glob.glob(f'{output_directory}/.temp_pka_worker*')]
    
def sort_buriedness_df(output_directory, output_prefix, buriedness_dfs, disulphide_cysteines_list):
    sorted_buriedness_df = pd.concat(buriedness_dfs).astype(float).sort_values('Time [ps]')
    if len(disulphide_cysteines_list) > 0:
        sorted_buriedness_df.drop(disulphide_cysteines_list, axis=1, inplace=True)
    sorted_buriedness_df.to_csv(f'{output_directory}/{output_prefix}_buriedness.csv', index=False)
    [os.remove(tmp_buriedness_csv) for tmp_buriedness_csv in glob.glob(f'{output_directory}/.temp_buriedness_worker*')]