import pandas as pd

def trajectory_log(output_directory,
                   output_prefix, 
                   extract_buridness_data,
                   chain, 
                   mutation, 
                   disulphide_bond_detection,
                   optargs,
                   cores,
                   trajectory_slices,
                   start,
                   end):
    
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

-----------------------------------------------------------------""" 

    
    with open(f'{output_directory}/{output_prefix}.log','w') as output:
        output.write(f"""{header}

Start time: {start}
End time: {end}

Parmeters selected:

Output directory: {output_directory}
Output prefix: {output_prefix}
Number of cores:{cores}
Trajectory slices: {trajectory_slices}
Chain: {chain}
Mutation: {mutation}
Extract buridness: {extract_buridness_data}
Remove cysteines from pKa calulations: {disulphide_bond_detection}
PropKa optional arguments: {optargs}
-----------------------------------------------------------------

pKa Statistics:
{pka_statistics_table(output_directory,output_prefix)}

""")

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
                                                'Standard Deviation'])
    tmp = pka_statistics_table.to_string(index=False).split('\n')
    tmp = [','.join(ele.split()) for ele in tmp]
    tmp = [element.replace('_',' ') for element in tmp]
    pka_statistics_table_in_str = '\r\n'.join(tmp)
    return pka_statistics_table_in_str