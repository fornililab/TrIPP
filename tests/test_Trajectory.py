from tripp.Trajectory import Trajectory
import os 
import glob
import pandas as pd
import pytest
import numpy as np

def get_data_path(filename):
    return os.path.join(os.path.dirname(__file__), "data", filename)
def check_dir(path):
    if not os.path.isdir(path):
        # Make directory if not present
        os.makedirs(path) 
    else:
        # Remove files named .temp* in the directory before proceeding.
        [os.remove(file) for file in glob.glob(f'{path}/.temp*')]
    return path

class TestTrajectory:
    def parm1():
        topology_file = get_data_path('myosin_elc_test.pdb')
        trajectory_file = get_data_path('myosin_elc_test.xtc')
        output_directory = check_dir('test_output/myosin_elc_test')
        output_prefix = 'myosin_elc_test'
        cpu_core_number = -1
        hetatm_resname = ['ADP','PI2','MG']
        custom_terminal_oxygens = None
        custom_resname_correction = None
        parameters = (topology_file, trajectory_file, output_directory, output_prefix, cpu_core_number, hetatm_resname, custom_terminal_oxygens, custom_resname_correction)
        return parameters
    @pytest.mark.parametrize("parameters", [parm1()])
    def test_Trajectory_mysoin_elc(self, parameters):
        TrIPP_Traj = Trajectory(*parameters)
        TrIPP_Traj.run(extract_buriedness_data=True,
                       mutation_selections=None,
                       save_disulphide_pka=False,
                       optargs=[])
        assert os.path.isfile(f'{parameters[2]}/{parameters[3]}_pka.csv')
        df = pd.read_csv(f'{parameters[2]}/{parameters[3]}_pka.csv')
        chains = np.unique(np.vstack(df.drop(columns=['Time [ps]']).columns.str.split(':'))[:,1])
        assert (chains == np.unique(TrIPP_Traj.corrected_universe.atoms.chainIDs)).all()
        df = pd.read_csv(f'{parameters[2]}/{parameters[3]}_buriedness.csv')
        chains = np.unique(np.vstack(df.drop(columns=['Time [ps]']).columns.str.split(':'))[:,1])
        assert (chains == np.unique(TrIPP_Traj.corrected_universe.atoms.chainIDs)).all()
        
        
        