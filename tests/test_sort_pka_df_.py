from tripp._sort_pka_df_ import output_df
import pandas as pd
import os 
import glob 

data = ([{0.0:{'residue_identifier_list': ['ASP1:A', 'LYS2:A', 'CYS3:A'], 'pka_list': [3.71, 10.79, 99.99], 'buriedness_list': [5.0, 2.0, 99.99]}},
        {1.0:{'residue_identifier_list': ['ASP1:A', 'LYS2:A', 'CYS3:A'], 'pka_list': [3.72, 10.78, 99.99], 'buriedness_list': [5.1, 2.1, 99.99]}},
        {2.0:{'residue_identifier_list': ['ASP1:A', 'LYS2:A', 'CYS3:A'], 'pka_list': [3.73, 10.77, 99.99], 'buriedness_list': [5.2, 2.2, 99.99]}}],
        [{5.0:{'residue_identifier_list': ['ASP1:A', 'LYS2:A', 'CYS3:A'], 'pka_list': [3.74, 10.76, 99.99], 'buriedness_list': [5.3, 2.3, 99.99]}},
        {6.0:{'residue_identifier_list': ['ASP1:A', 'LYS2:A', 'CYS3:A'], 'pka_list': [3.75, 10.75, 99.99], 'buriedness_list': [5.4, 2.4, 99.99]}},
        {7.0:{'residue_identifier_list': ['ASP1:A', 'LYS2:A', 'CYS3:A'], 'pka_list': [3.76, 10.74, 99.99], 'buriedness_list': [5.5, 2.5, 99.99]}}],
        [{3.0:{'residue_identifier_list': ['ASP1:A', 'LYS2:A', 'CYS3:A'], 'pka_list': [3.77, 10.73, 99.99], 'buriedness_list': [5.6, 2.6, 99.99]}},
        {4.0:{'residue_identifier_list': ['ASP1:A', 'LYS2:A', 'CYS3:A'], 'pka_list': [3.78, 10.72, 99.99], 'buriedness_list': [5.7, 2.7, 99.99]}}])

def check_dir(path):
    if not os.path.isdir(path):
        # Make directory if not present
        os.makedirs(path) 
    else:
        # Remove files named .temp* in the directory before proceeding.
        [os.remove(file) for file in glob.glob(f'{path}/*')]
    return path

class Test_sort_pka_df:
    def test_output_df_no_cys(self):
        output_directory = check_dir('test_output/output_df_no_cys')
        disulphide_cys_col = output_df(
            output_directory=output_directory,
            output_prefix='test_no_cys',
            data=data,
            save_disulphide_pka=False,
            extract_buriedness_data=True
        )
        assert os.path.isfile(f'{output_directory}/test_no_cys_pka.csv')
        assert os.path.isfile(f'{output_directory}/test_no_cys_buriedness.csv')
        pka_df = pd.read_csv(f'{output_directory}/test_no_cys_pka.csv')
        buriedness_df = pd.read_csv(f'{output_directory}/test_no_cys_buriedness.csv')
        assert pka_df['Time [ps]'].tolist() == [0.0,1.0,2.0,3.0,4.0,5.0,6.0,7.0]
        assert not any('CYS' in col for col in pka_df.columns.tolist())
        assert buriedness_df['Time [ps]'].tolist() == [0.0,1.0,2.0,3.0,4.0,5.0,6.0,7.0]
        assert not any('CYS' in col for col in buriedness_df.columns.tolist())
        assert disulphide_cys_col == ['CYS3:A']
        
    def test_output_df_with_cys(self):
        output_directory = check_dir('test_output/output_df_with_cys')
        disulphide_cys_col = output_df(
            output_directory=output_directory,
            output_prefix='test_with_cys',
            data=data,
            save_disulphide_pka=True,
            extract_buriedness_data=True
        )
        assert os.path.isfile(f'{output_directory}/test_with_cys_pka.csv')
        assert os.path.isfile(f'{output_directory}/test_with_cys_buriedness.csv')
        pka_df = pd.read_csv(f'{output_directory}/test_with_cys_pka.csv')
        buriedness_df = pd.read_csv(f'{output_directory}/test_with_cys_buriedness.csv')
        assert pka_df['Time [ps]'].tolist() == [0.0,1.0,2.0,3.0,4.0,5.0,6.0,7.0]
        assert 'CYS3:A' in pka_df.columns.tolist()
        assert buriedness_df['Time [ps]'].tolist() == [0.0,1.0,2.0,3.0,4.0,5.0,6.0,7.0]
        assert 'CYS3:A' in buriedness_df.columns.tolist()
        assert disulphide_cys_col is None