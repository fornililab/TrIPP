import pandas as pd
import logging

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

def pka_statistics_table(output_directory,
                        output_prefix):
    df = pd.read_csv(f'{output_directory}/{output_prefix}_pka.csv')
    pka_statistics = []
    for residue,pKaValues in df.items():
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
                   chain, 
                   mutation_arg, 
                   disulphide_cysteines_list,
                   optargs,
                   cores,
                   trajectory_slices,
                   start,
                   end):
    logger = logging.getLogger(__name__)
    if len(disulphide_cysteines_list) > 0:
        disulphide_bond_detection = True
    else:
        disulphide_bond_detection = False
    logger.info(f"""
-----------------------------------------------------------------                

Start time: {start}
End time: {end}

Parmeters selected:

Output directory: {output_directory}
Output prefix: {output_prefix}
Number of cores: {cores}
Trajectory slices:
{trajectory_slices}
Chain: {chain}
Mutation: {mutation_arg}
Extract buriedness: {extract_buriedness_data}
Remove cysteines from pKa or buriedness CSV: {disulphide_bond_detection}
{disulphide_cysteines_list}
PropKa optional arguments: {optargs}

-----------------------------------------------------------------

pKa Statistics:
{pka_statistics_table(output_directory,output_prefix)}

-----------------------------------------------------------------
""")