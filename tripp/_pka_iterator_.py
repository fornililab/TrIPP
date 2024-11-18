import MDAnalysis as mda 
from tripp._edit_pdb_ import mutate 
from propka import run 
from tripp._extract_pka_file_data_ import extract_pka_file_data 
import os
import pandas as pd


def pka_iterator(trajectory_slice, universe,
                 output_directory, mutation, chain,
                 extract_buriedness_data, optargs=[]):
    """
    Function to run propka.run.single on the distributed trajectory slice.
    
    Parameters
    ----------
    trajectory_slices: list of int
        Inherited from the Trajectory class initialisation, where trajectory
        slicing is performed.
    universe: MDAnalysis Universe object
        Inherited from the Trajectory class initialisation, which is the
        corrected universe
        
    """
    
    pid = os.getpid()
    
    temp_name = f'{output_directory}/.temp_{pid}'

    start = trajectory_slice[0]
    end = trajectory_slice[1]

    # Note if you run two pka_iterator process in the same folder, this would 
    if os.path.isfile(f'{output_directory}/.temp_pka_worker_{pid}.csv'):
        os.remove(f'{output_directory}/.temp_pka_worker_{pid}.csv')

    if os.path.isfile(f'{output_directory}/.temp_buriedness_worker_{pid}.csv'):
        os.remove(f'{output_directory}/.temp_buriedness_worker_{pid}.csv')

    cwd = os.getcwd()
    pka_data = []
    buriedness_data = []
    for ts in universe.trajectory[start:end]:
        with mda.Writer(f'{temp_name}.pdb') as w:
            w.write(universe)
        if mutation is not None:
            mutate(temp_name, mutation)

        os.chdir(output_directory)
        run.single(f'.temp_{pid}.pdb', optargs=optargs)
        os.chdir(cwd)

        time = ts.time
        # Write pKa csv
        pka, buriedness = extract_pka_file_data(f'{temp_name}.pka', chain=chain, time=time)
        pka_data.append(pka[1])
        # Write buriedness data if True
        if extract_buriedness_data:
            buriedness_data.append(buriedness[1])

        os.remove(f'{temp_name}.pdb') 
        os.remove(f'{temp_name}.pka')
        
    pka_df = pd.DataFrame(pka_data, columns = pka[0])
    # pka_df.to_csv(f'{output_directory}/.temp_pka_worker_{pid}.csv', sep=',', index=False)
    
    buriedness_df = None
    if extract_buriedness_data:
        buriedness_df = pd.DataFrame(buriedness_data, columns=buriedness[0])
        # buriedness_df.to_csv(f'{output_directory}/.temp_buriedness_worker_{pid}.csv', sep=',', index=False)
    return pka_df, buriedness_df
