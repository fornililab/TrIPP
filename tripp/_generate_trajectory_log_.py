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
    """
    Generates a formatted string representation of pKa statistics from a DataFrame.
    Parameters
    ----------
    df: pd.DataFrame
        A DataFrame containing pKa values for different residues at different time points.
        Columns (labelled with residue identifiers) correspond to residues and 
        and rows to time points.
    Returns
    -------
    pka_statistics_table_in_str: str
        A string representation of the pKa statistics table.
    """
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
                   disulphide_cys_col,
                   optargs,
                   cores,
                   trajectory_slices,
                   start,
                   end):
    """
    Logs the parameters and results of the pKa calculation.
    Parameters
    ----------
    output_directory : str
        The directory where the output files are saved.
    output_prefix : str
        The prefix for the output files.
    extract_buriedness_data : bool
        Whether to extract buriedness (buried ratio) data. 
    mutation_selection : str
        Selection string for pseudo-mutations.
    disulphide_cys_col : list | None
        List of disulphide bonded cysteines.
    optargs : dict
        Optional arguments for the PROPKA calculation.
    cores : int
        The number of cores to use for the calculation.
    trajectory_slices : list
        A list of the trajectory slices that are processed in parallel.
    start : str
        Wall-clock start time.
    end : str
        Wall-clock end time.
    """
    if disulphide_cys_col:
        save_disulphide_pka = False
    else:
        save_disulphide_pka = True
        
    logger.info(f"""-----------------------------------------------------------------                

Start time: {start}
End time: {end}

PARAMETERS:
Output directory: {output_directory}
Output prefix: {output_prefix}
Number of cores: {cores}
Trajectory slices: {trajectory_slices}
Pseudo-mutations: {mutation_selection}
Extract buriedness: {extract_buriedness_data}
Save disulphide bond cysteines in csv: {save_disulphide_pka}
List of cysteines removed: {disulphide_cys_col}
PROPKA optional arguments: {optargs}

-----------------------------------------------------------------
""")
    df = pd.read_csv(f'{output_directory}/{output_prefix}_pka.csv')
    df.drop(columns=['Time [ps]'],inplace=True)
    chains = set([x.split(':')[-1] for x in df.columns])
    for chain in chains:
        df_chain = df[df.columns[df.columns.str.split(':').str[-1] == chain]]
        logger.info(f"""pKa Statistics for chain {chain}:
{pka_statistics_table(df_chain)}

-----------------------------------------------------------------
""")