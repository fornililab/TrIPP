import pandas as pd
import logging

logger = logging.getLogger(__name__)

def log_header():
    header = """
888888888888         88  88888888ba   88888888ba   
     88              88  88      "8b  88      "8b  
     88              88  88      ,8P  88      ,8P  
     88  8b,dPPYba,  88  88aaaaaa8P'  88aaaaaa8P'  
     88  88P'   "Y8  88  88""""""'    88""""""'    
     88  88          88  88           88           
     88  88          88  88           88           
     88  88          88  88           88           


The Trajectory Iterative pKa Predictor (TrIPP) 
Written by: Christos Matsingos, Ka Fu Man, and Arianna Fornili 

If you are using TrIPP, please cite: 

-----------------------------------------------------------------
"""
    return header

def pka_statistics_table(df):
    pka_statistics = []
    for residue, pKaValues in df.items():
        if 'Time [ps]' == residue:
            continue
        pka_statistics.append([residue,
                            "{:.2f}".format(pKaValues.mean()),
                            "{:.2f}".format(pKaValues.median()),
                            "{:.2f}".format(pKaValues.std())])
    pka_statistics_table = pd.DataFrame(pka_statistics,
                                        columns=['Residue',
                                                'Mean',
                                                'Median',
                                                'Standard_Deviation'])
    tmp = pka_statistics_table.to_string(index=False).split('\n')
    tmp = [','.join(ele.split()) for ele in tmp]
    tmp = [element.replace('_',' ') for element in tmp]
    pka_statistics_table_in_str = '\r\n'.join(tmp)
    return pka_statistics_table_in_str

def trajectory_log(output_directory,
                   output_prefix, 
                   extract_buriedness_data,
                   mutation_selection, 
                   disulphide_cysteines_list,
                   optargs,
                   cores,
                   trajectory_slices,
                   start,
                   end):
    if len(disulphide_cysteines_list) > 0:
        disulphide_bond_detection = True
    else:
        disulphide_bond_detection = False
        
    logger.info(f"""-----------------------------------------------------------------                

Start time: {start}
End time: {end}

PARAMETERS:
Output directory: {output_directory}
Output prefix: {output_prefix}
Number of cores: {cores}
Trajectory slices:
{trajectory_slices}
Mutation: {mutation_selection}
Extract buriedness: {extract_buriedness_data}
Remove cysteines from pKa or buriedness CSV: {disulphide_bond_detection}
List of cysteines removed: {disulphide_cysteines_list}
PropKa optional arguments: {optargs}

-----------------------------------------------------------------
""")
    df = pd.read_csv(f'{output_directory}/{output_prefix}_pka.csv')
    df.drop(columns=['Time [ps]'],inplace=True)
    chains = set([x.split(':')[-1] for x in df.columns])
    for chain in chains:
        logger.info(f"""pKa Statistics for chain {chain}:
{pka_statistics_table(df)}

-----------------------------------------------------------------
""")