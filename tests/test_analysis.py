from tripp.analysis import calculate_difference_to_model, PCProjectionScreening
import os
import pandas as pd
import glob
import pytest

def check_dir(path):
    if not os.path.isdir(path):
        # Make directory if not present
        os.makedirs(path) 
    else:
        # Remove files named .temp* in the directory before proceeding.
        [os.remove(file) for file in glob.glob(f'{path}/*')]
    return path

class TestAnalysis:
    def test_calculate_difference_to_model(self):
        input_directory = 'reference_output/lyso_test_default'
        output_prefix = 'lyso_test_default'
        calculate_difference_to_model(input_directory, output_prefix)
        assert os.path.isfile(f'{input_directory}/{output_prefix}_difference_to_model.csv')
    def parm1():
        output_directory = check_dir('test_output/PCProjectionScreening')
        output_prefix = 'lyso_test_default_pka_projection_correlation'
        pka_file = 'reference_output/lyso_test_default/lyso_test_default_pka.csv'
        projection_file = 'data/lyso_test_pc1.csv'
        method = 'Pearson'
        start = 0
        end = -1
        header = None
        return (output_directory, output_prefix, pka_file, projection_file, method, start, end, header)
    def parm2():
        output_directory = check_dir('test_output/PCProjectionScreening')
        output_prefix = 'lyso_test_default_pka_projection_correlation'
        pka_file = 'reference_output/lyso_test_default/lyso_test_default_pka.csv'
        projection_file = 'data/lyso_test_pc1.csv'
        method = 'Spearman'
        start = 0
        end = -1
        header = None
        return (output_directory, output_prefix, pka_file, projection_file, method, start, end, header)
    def parm3():
        output_directory = check_dir('test_output/PCProjectionScreening')
        output_prefix = 'lyso_test_default_pka_projection_correlation'
        pka_file = 'reference_output/lyso_test_default/lyso_test_default_pka.csv'
        projection_file = 'data/lyso_test_pc1.csv'
        method = 'Pearson'
        start = 1000
        end = 2000
        header = None
        return (output_directory, output_prefix, pka_file, projection_file, method, start, end, header)
    
    @pytest.mark.parametrize("output_directory, output_prefix, pka_file, projection_file, method, start, end, header", [
    parm1(), parm2(), parm3()])
    def test_PCProjectionScreening_success(self, output_directory, output_prefix, pka_file, projection_file, method, start, end, header):
        PCProjectionScreening(output_directory=output_directory,
                              output_prefix=output_prefix,
                              pka_file=pka_file,
                              projection_file=projection_file,
                              method=method,
                              start=start,
                              end=end,
                              header=header)
        assert os.path.isfile(f'{output_directory}/{output_prefix}.csv')
        df = pd.read_csv(f'{output_directory}/{output_prefix}.csv')
        pka_df = pd.read_csv(pka_file)
        assert 'Residue' in df.columns
        assert 'Correlation' in df.columns
        assert 'p-value' in df.columns
        assert df['Residue'].to_list() == pka_df.columns.drop('Time [ps]').to_list()
        
    def parm4():
        output_directory = check_dir('test_output/PCProjectionScreening')
        output_prefix = 'lyso_test_default_pka_projection_correlation'
        pka_file = 'reference_output/lyso_test_default/lyso_test_default_pka.csv'
        projection_file = 'data/lyso_test_pc1.csv'
        method = 'apsgh'
        start = 1000
        end = 2000
        header = None
        return (output_directory, output_prefix, pka_file, projection_file, method, start, end, header)
    def parm5():
        output_directory = check_dir('test_output/PCProjectionScreening')
        output_prefix = 'lyso_test_default_pka_projection_correlation'
        pka_file = 'reference_output/lyso_test_default/lyso_test_default_pka.csv'
        projection_file = 'data/lyso_test_pc1.csv'
        method = 'Pearson'
        start = 1000
        end = -5
        header = None
        return (output_directory, output_prefix, pka_file, projection_file, method, start, end, header)
    def parm6():
        output_directory = check_dir('test_output/PCProjectionScreening')
        output_prefix = 'lyso_test_default_pka_projection_correlation'
        pka_file = 'reference_output/lyso_test_default/lyso_test_default_pka.csv'
        projection_file = 'data/lyso_test_pc1.csv'
        method = 'Pearson'
        start = 6
        end = -5
        header = None
        return (output_directory, output_prefix, pka_file, projection_file, method, start, end, header)
    @pytest.mark.parametrize("output_directory, output_prefix, pka_file, projection_file, method, start, end, header", [
    parm4(), parm5(), parm6()])
    def test_PCProjectionScreening_failure(self, output_directory, output_prefix, pka_file, projection_file, method, start, end, header):
        with pytest.raises(Exception):
            PCProjectionScreening(output_directory=output_directory,
                              output_prefix=output_prefix,
                              pka_file=pka_file,
                              projection_file=projection_file,
                              method=method,
                              start=start,
                              end=end,
                              header=header)