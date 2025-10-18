#!/home/ivan/bin/miniconda3/envs/tripp/bin/python
from tripp import Trajectory
from tripp.analysis import calculate_difference_to_model
from tripp.analysis import PCProjectionScreening  
from tripp._create_mda_universe_ import create_mda_universe
import pandas as pd 
import matplotlib.pyplot as plt 
import seaborn as sns
import numpy as np
import re
from itertools import combinations
from tripp._model_pka_values_ import model_pka_values
from tripp import Visualization  # PyMOL installation required, user needs to specify the path to the PyMOL executable (see Part 2)
from tripp import Clustering

# Start of user-defined variables -------------------------

topology = '1AKI_minstr.pdb'  # name of the reference PDB file used to define the topology of the system
input_dir = 'data'            # name of the folder containing the input files
traj_prefix_list = ['1AKI_md1_300ns'] # list of input trajectories (without the extension)
traj_format = 'xtc'           # extension of trajectory input files (only xtc files have been tested)
output_dir = 'output'         # name of the folder where all output files and folders will be saved

# End of user-defined variables -------------------------

ntraj = len(traj_prefix_list)
topology_file = input_dir + '/' + topology
trajectory_files = [input_dir + '/' + x + '.' + traj_format for x in traj_prefix_list]

pka_files = [output_dir + '/' + x + '/' + x + '_pka.csv' for x in traj_prefix_list]
buriedness_files = [output_dir + '/' + x + '/' + x + '_buriedness.csv' for x in traj_prefix_list]
trajectory_dict = {x: input_dir + '/' + x + '.' + traj_format for x in traj_prefix_list}

output_directories = [output_dir + '/' + x for x in traj_prefix_list]
output_prefix = traj_prefix_list

for idx in range(ntraj):
    print(f'Processing {trajectory_files[idx]}')
    TrIPP_Traj = Trajectory(topology_file = topology_file,
                            trajectory_file = trajectory_files[idx],
                            output_directory = output_directories[idx],
                            output_prefix = output_prefix[idx],
                            hetatm_resname = None, # str or list of str of PDB residue names for non-protein molecules 
                                                   # that we want to be taken into account by PROPKA in the pKa calculation.
                                                   # Their record type will be set to 'HETATM'. e.g.: 'ADP' or ['ADP','PI2']
                            custom_terminal_oxygens = None,  # names of terminal oxygen atoms when different from ['O', 'OXT'], 
                                                             # provided as a list, e.g.: ['OC1', 'OC2']
                            custom_resname_correction = None,   # dictionary of custom protein residue names not included in
                                                                # the hard-coded TrIPP dictionary (tripp._correction_dictionary_.py). 
                                                                # Can be given as e.g. {'XXX':'ASP'}, where 'XXX' is the residue
                                                                # name in the PDB file and 'ASP' is the corresponding
                                                                # PROPKA name.
                            cpu_core_number = 16)   # number of CPU cores. Set to -1 to use all available cores.
    
    TrIPP_Traj.run(extract_buriedness_data = True, mutation_selections = None,save_disulphide_pka = False, optargs = [])
    print('---------------------------------------------------------------------------------')
