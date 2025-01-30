import MDAnalysis as mda 
from tripp._edit_pdb_ import mutate 
from propka import run 
from tripp._extract_pka_file_data_ import extract_pka_buriedness_data 
import os
import pandas as pd
import logging
import io

def pka_iterator(trajectory_slice, universe,
                 output_directory, mutation_selections, chain,
                 optargs=[]):
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
    # Redirect warning from propka.group to a file
    logger = logging.getLogger('propka.group')
    log_capture_string = io.StringIO()
    handler = logging.StreamHandler(log_capture_string)
    logger.addHandler(handler)
    log_contents = None
    
    pid = os.getpid()
    
    temp_name = f'{output_directory}/.temp_{pid}'

    start = trajectory_slice[0]
    end = trajectory_slice[1]
    
    cwd = os.getcwd()
    data = []
    for ts in universe.trajectory[start:end]:
        if mutation_selections is not None:
            mutate(universe, mutation_selections, temp_name)
        else:
            with mda.Writer(f'{temp_name}.pdb') as w:
                w.write(universe)
        os.chdir(output_directory)
        run.single(f'.temp_{pid}.pdb', optargs=optargs)
        
        if log_capture_string.getvalue():
            log_contents = ("-----------------------------------------------------------------\n"+
                            "\n"
                            f"PROPKA warning occured on frame {ts.frame}:\n" + 
                            log_capture_string.getvalue())
        os.chdir(cwd)

        time = ts.time
        # Write pKa csv
        data_dictionary = extract_pka_buriedness_data(f'{temp_name}.pka', chains=chain, time=time)
        data.append(data_dictionary)

        os.remove(f'{temp_name}.pdb') 
        os.remove(f'{temp_name}.pka')
    
    logger.removeHandler(handler)
    log_capture_string.close()
    
    return data, log_contents
