from tripp.Visualization import Visualization
import pytest
import numpy as np
import os
import glob

def check_dir(path):
    if not os.path.isdir(path):
        # Make directory if not present
        os.makedirs(path) 
    else:
        # Remove files named .temp* in the directory before proceeding.
        [os.remove(file) for file in glob.glob(f'{path}/*')]
    return path

PyMOL_path = False
while not PyMOL_path:
    PyMOL_path = input('Please provide path to PyMOL executable to test Visualization class or type "skip":')
    
def fix_parm():
    topology_file = 'data/lyso.pdb'
    pka_file = 'reference_output/lyso_test_default/lyso_test_default_pka.csv'
    correlation_file = 'reference_output/lyso_test_PearsonCorrelation.csv'
    start = 0
    end = -1
    parameters = (topology_file, pka_file, correlation_file, start, end)
    return parameters
@pytest.fixture(params=[fix_parm()])
def return_param(request):
    return request.param
class TestVisualization:
    def parm1():
        topology_file = 'data/lyso.pdb'
        pka_file = 'reference_output/lyso_test_default/lyso_test_default_pka.csv'
        correlation_file = None
        start = 0
        end = -1
        parameters = (topology_file, pka_file, correlation_file, start, end)
        return parameters
    def parm2():
        topology_file = 'data/lyso.pdb'
        pka_file = 'reference_output/lyso_test_default/lyso_test_default_pka.csv'
        correlation_file = None
        start = 100
        end = 1000
        parameters = (topology_file, pka_file, correlation_file, start, end)
        return parameters
    def parm3():
        topology_file = 'data/lyso.pdb'
        pka_file = 'reference_output/lyso_test_default/lyso_test_default_pka.csv'
        correlation_file = 'reference_output/lyso_test_PearsonCorrelation.csv'
        start = 0
        end = -1
        parameters = (topology_file, pka_file, correlation_file, start, end)
        return parameters
        
    @pytest.mark.parametrize("parameters", [parm1(), parm2(), parm3()])
    def test_Visualization(self, parameters):
        TrIPP_Vis = Visualization(*parameters)
        pka_df = TrIPP_Vis.pka_values
        assert len(pka_df['Time [ps]'].to_list()) == sum((pka_df['Time [ps]'] >= TrIPP_Vis.start) & (pka_df['Time [ps]'] <= TrIPP_Vis.end))
        if parameters[2]:
            assert len(TrIPP_Vis.correlation_file['Residue']) == len(TrIPP_Vis.pka_values.columns) - 1 # -1 because of Time [ps] column
        
    def gen_pse_parm1():
        output_directory = check_dir('test_output/TestVisualization')
        output_prefix = 'lyso_gen_pse_parm1'
        coloring_method = 'mean'
        lower_limit = 0
        upper_limit = 14
        correlation_threshold = None
        color_palette = 'red_white_blue'
        parameters = (output_directory, output_prefix, coloring_method, lower_limit, upper_limit, correlation_threshold, color_palette)
        return parameters
    def gen_pse_parm2():
        output_directory = check_dir('test_output/TestVisualization')
        output_prefix = 'lyso_gen_pse_parm2'
        coloring_method = 'difference_to_model_value'
        lower_limit = 0
        upper_limit = 14
        correlation_threshold = None
        color_palette = 'red_white_blue'
        parameters = (output_directory, output_prefix, coloring_method, lower_limit, upper_limit, correlation_threshold, color_palette)
        return parameters
    def gen_pse_parm3():
        output_directory = check_dir('test_output/TestVisualization')
        output_prefix = 'lyso_gen_pse_parm3'
        coloring_method = 'correlation'
        lower_limit = -1
        upper_limit = 1
        correlation_threshold = 0.5
        color_palette = 'red_white_blue'
        parameters = (output_directory, output_prefix, coloring_method, lower_limit, upper_limit, correlation_threshold, color_palette)
        return parameters
    @pytest.mark.parametrize("parameters", [gen_pse_parm1(), gen_pse_parm2(), gen_pse_parm3()])
    def test_gen_pse(self, return_param, parameters):
        TrIPP_Vis = Visualization(*list(return_param))
        TrIPP_Vis.gen_pse(PyMOL_path, *parameters)
        if parameters[2] == 'mean':
            assert os.path.isfile(f'{parameters[0]}/{parameters[1]}_mean.pdb'), 'PDB (mean) file was not created'
            assert os.path.isfile(f'{parameters[0]}/{parameters[1]}_mean.pse'), 'PyMOL session (mean) file was not created'
        elif parameters[2] == 'difference_to_model_value':
            assert os.path.isfile(f'{parameters[0]}/{parameters[1]}_difference_to_model_value.pdb'), 'PDB (difference to model value) file was not created'
            assert os.path.isfile(f'{parameters[0]}/{parameters[1]}_difference_to_model_value.pse'), 'PyMOL session (difference to model value) file was not created'
        elif parameters[2] == 'correlation':
            assert os.path.isfile(f'{parameters[0]}/{parameters[1]}_correlation.pdb'), 'PDB (correlation) file was not created'
            assert os.path.isfile(f'{parameters[0]}/{parameters[1]}_correlation.pse'), 'PyMOL session (correlation) file was not created'