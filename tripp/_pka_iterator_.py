"""
    This file is part of the TrIPP software
    (https://github.com/fornililab/TrIPP).
    Copyright (c) Christos Matsingos, Ka Fu Man and Arianna Fornili.

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, version 3.

    This program is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program. If not, see <http://www.gnu.org/licenses/>.
"""

import MDAnalysis as mda 
from tripp._edit_pdb_ import mutate 
from propka import run 
from tripp._extract_pka_file_data_ import propka_pka_buriedness_data_parser, pkai_pka_data_parser
import os
import logging
import io
import sys

def write_temp_pdb(trajectory_slice, universe, temp_name):
    """
    Function to write a temporary pdb file for the trajectory slice.
    
    Parameters
    ----------
    trajectory_slices: list of int
        Trajectory slices from the Trajectory class initialisation.
    universe: MDAnalysis.universe object
        Modified MDAnalysis universe from the Trajectory class initialisation.
    temp_name: str
        Name of the temporary pdb file to be written.
    """
    start = trajectory_slice[0]
    end = trajectory_slice[1]
    
    for ts in universe.trajectory[start:end]:
        with mda.Writer(f'{temp_name}.pdb') as w:
            w.write(universe)
            
def propka_predictor(pdb_file, optargs, frame):
    # Redirect warning from propka to a file
    logger = logging.getLogger('propka')
    logger.propagate = False 
    logger.setLevel(logging.WARNING)
    log_capture_string = io.StringIO()
    handler = logging.StreamHandler(log_capture_string)
    logger.addHandler(handler)
    molecule = run.single(pdb_file, optargs=optargs, write_pka=False)
    log_contents = None
    if log_capture_string.getvalue():
        log_contents = (f"PROPKA warning raised for frame {frame}:\n" + 
                        log_capture_string.getvalue()+
                        '\n')
    logger.removeHandler(handler)
    log_capture_string.close()
    
    return molecule, log_contents

def pkai_predictor(pdb_file, predictor, frame):
    from pkai.pKAI import pKAI
    # Redirect warning from pKAI to a file
    logger = logging.getLogger('pKAI')
    logger.propagate = False
    logger.setLevel(logging.WARNING)
    log_capture_string = io.StringIO()
    handler = logging.StreamHandler(log_capture_string)
    logger.addHandler(handler)
    # Redirect stdout to capture Titration output
    text_trap = io.StringIO()
    sys.stdout = text_trap
    
    if predictor.lower() == 'pkai':
        model_name = 'pKAI'
    elif predictor.lower() == 'pkai+':
        model_name = 'pKAI+'
    pka_result = pKAI(pdb=pdb_file, model_name=model_name, device='cpu', threads=1)
    log_contents = None
    if log_capture_string.getvalue():
        log_contents = (f"pKAI warning raised for frame {frame}:\n" +
                        log_capture_string.getvalue() +
                        '\n')
    logger.removeHandler(handler)
    log_capture_string.close()
    sys.stdout = sys.__stdout__

    return pka_result, log_contents

def pka_iterator(trajectory_slice, universe,
                 output_directory, mutation_selections,
                 predictor='propka', optargs=[]):
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
    predictor: str, default='propka'
        The pKa predictor to use. Options are 'propka', 'pkai' and 'pkai+'.
    optargs: list of str, default=[]
        PROPKA predictions can be run with optional arguments
        (see https://propka.readthedocs.io/en/latest/command.html).
        For example, if optargs is set to `["-k"]`, propka will run with the -k flag
        (protons from the input file are kept).
        Only valid for PROPKA predictor, this parameter will be ignored if
        pKAI/pKAI+ predictor is used.
    """
    
    pid = os.getpid()
    
    temp_name = f'{output_directory}/.temp_{pid}'

    start = trajectory_slice[0]
    end = trajectory_slice[1]
    
    cwd = os.getcwd()
    data = []
    for ts in universe.trajectory[start:end]:
        time = ts.time
        frame = ts.frame
        if mutation_selections is not None:
            mutate(universe, mutation_selections, temp_name)
        else:
            with mda.Writer(f'{temp_name}.pdb') as w:
                w.write(universe)
        os.chdir(output_directory)
        temp_pdb_file = f'.temp_{pid}.pdb'
        if predictor == 'propka':
            molecule, log_contents = propka_predictor(temp_pdb_file, optargs, frame)
            data_dictionary = propka_pka_buriedness_data_parser(molecule, time=time)
            os.chdir(cwd)
            os.remove(f'{temp_name}.pdb')
        elif predictor.lower() == 'pkai' or predictor.lower() == 'pkai+':
            result, log_contents = pkai_predictor(temp_pdb_file, predictor, frame)
            data_dictionary = pkai_pka_data_parser(result, time=time)
            os.chdir(cwd)
            os.remove(f'{temp_name}.pdb')
        data.append(data_dictionary)
                        
    return data, log_contents