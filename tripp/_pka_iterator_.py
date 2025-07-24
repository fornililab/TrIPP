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
        Inherited from the Trajectory class initialisation, where trajectory
        slicing is performed.
    universe: MDAnalysis.universe object
        Inherited from the Trajectory class initialisation, which is the
        corrected universe
    output_directory: str
        Directory to write the propka output files to.
    mutation_selections: list of str
        List of selection strings for the atoms to mutate.
    optargs: list of str, default=[]
        PROPKA prediction can be run with optional arguments as indicated
        in their documentation. Each flag is string separated by a comma.
        For example, `["-k","--pH=7.2"]` keeps the proton in your topology
        and setting pH 7.2 for ie: stability calculations, respectively.
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
            log_contents = (f"PROPKA warning occured on frame {ts.frame}:\n" + 
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
