import pandas as pd 
import logging
import numpy as np

logger = logging.getLogger(__name__)

def output_df(output_directory, output_prefix, data, disulphide_cysteines_list, extract_buriedness_data):
    for chain in data[0][0].keys():
        pka_full_data = []
        buriedness_full_data = []
        time = []
        for data_worker in data:
            for pka_buriedness_t_data in data_worker:
                pka_full_data.append(pka_buriedness_t_data[chain]['pka_list'])
                buriedness_full_data.append(pka_buriedness_t_data[chain]['buriedness_list'])
                time.append(pka_buriedness_t_data[chain]['Time [ps]'])
        pka_df = pd.DataFrame(pka_full_data, columns=pka_buriedness_t_data[chain]['resname_list'])
        pka_df.insert(0,'Time [ps]',time)

        try:
            pka_df.astype(float).sort_values('Time [ps]', inplace=True)
        except ValueError:
            contains_star = pka_df.applymap(lambda x: '*' in x if isinstance(x, str) else False)
            coupled_residues = ", ".join(list(pka_df.loc[:, contains_star.any()].columns))
            logging.info(f"""PROPKA detected that your chain {chain} has the following coupled residues:
{coupled_residues}
run -d in the optargs to consider alternative state.
""")
            pka_df = pka_df.applymap(lambda x: x.replace('*','') if isinstance(x, str) else x)
            pka_df.astype(float).sort_values('Time [ps]', inplace=True)
            
        if extract_buriedness_data:
            buriedness_df = pd.DataFrame(buriedness_full_data, columns=pka_buriedness_t_data[chain]['resname_list'])
            buriedness_df.insert(0,'Time [ps]',time)
            buriedness_df['Time [ps]'] = buriedness_df['Time [ps]'].astype(float)
            buriedness_df.sort_values('Time [ps]', inplace=True)

        if len(disulphide_cysteines_list) > 0:
            pka_df.drop(disulphide_cysteines_list, axis=1, inplace=True)
            if extract_buriedness_data:
                buriedness_df.drop(disulphide_cysteines_list, axis=1, inplace=True)
        
        pka_df.to_csv(f'{output_directory}/{output_prefix}_pka_chain{chain}.csv', index=False)
        
        if extract_buriedness_data:
            buriedness_df.to_csv(f'{output_directory}/{output_prefix}_buriedness_chain{chain}.csv', index=False)