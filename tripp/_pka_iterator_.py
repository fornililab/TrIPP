import MDAnalysis as mda 
from tripp._edit_pdb_ import mutate 
from propka import run 
from tripp._extract_pka_file_data_ import extract_pka_buriedness_data 
import os
import glob
import logging
import io

def pka_iterator(trajectory_slice, universe,
                 output_directory, mutation_selections,
                 optargs=[]):
    """
    Function to run propka.run.single on the distributed trajectory slice.
    
    Parameters
    ----------
    trajectory_slices: list of int
        Trajectory slices from the Trajectory class initialisation.
    universe: MDAnalysis.universe object
        Modified MDAnalysis universe from the Trajectory class initialisation.
    output_directory: str
        Directory to write the PROPKA output files to.
    mutation_selections: str
        Selection string in MDAnalysis format (only for pseudo-mutations) 
    optargs: list of str, default=[]
        PROPKA predictions can be run with optional arguments
        (see https://propka.readthedocs.io/en/latest/command.html).
        For example, if optargs is set to `["-k"]`, propka will run with the -k flag
        (protons from the input file are kept).
    """
    # Redirect warning from propka.group to a file
    logger = logging.getLogger('propka')
    logger.propagate = False 
    logger.setLevel(logging.WARNING)
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
            log_contents = (f"PROPKA warning raised for frame {ts.frame}:\n" + 
                            log_capture_string.getvalue()+
                            '\n')
        os.chdir(cwd)

        time = ts.time
        # Extract pKa and buriedness data from the generated .pka file
        # and append it to the data list.
        temp_pka_file = glob.glob(f'{temp_name}*.pka')[0] # .pka could have different suffixes ie: _alt_state.pka
        data_dictionary = extract_pka_buriedness_data(temp_pka_file, time=time)
        data.append(data_dictionary)

        os.remove(f'{temp_name}.pdb') 
        os.remove(temp_pka_file)
    
    logger.removeHandler(handler)
    log_capture_string.close()
    
    return data, log_contents
