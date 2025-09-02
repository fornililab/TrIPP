from tripp._pca_ import pca
from tripp._create_clustering_matrix_ import create_clustering_matrix
import numpy as np
import os
import glob
import pandas as pd
import pytest

def get_data_path(filename):
    return os.path.join(os.path.dirname(__file__), "data", filename)
def check_dir(path):
    if not os.path.isdir(path):
        # Make directory if not present
        os.makedirs(path) 
    else:
        # Remove files named .temp* in the directory before proceeding.
        [os.remove(file) for file in glob.glob(f'{path}/*')]
    return path
def make_pka_or_buriedness_df(file, trajectory_file):
        """
        Parser from Clustering.py
        """
        if isinstance(file, str):
            df = pd.read_csv(file, index_col="Time [ps]")
            df["Trajectories"] = "Unnamed trajectory"

        elif isinstance(file, list):
            df_list = []
            for index, f in enumerate(file):
                df = pd.read_csv(f, index_col="Time [ps]")
                df["Trajectories"] = list(trajectory_file.keys())[index]
                df_list.append(df)
            df = pd.concat(df_list)

        elif file is None:
            df = None
        return df 
    
class TestPCA:
    def parm():
        topology_file = get_data_path('lyso.pdb')
        trajectory_file = {
            "traj1": get_data_path('lyso_test.xtc'),
            "traj2": get_data_path('lyso_test.xtc')}
        pka_file = ['reference_output/lyso_test_default/lyso_test_default_pka.csv',
                    'reference_output/lyso_test_default/lyso_test_default_pka.csv']
        pka_df = make_pka_or_buriedness_df(pka_file, trajectory_file)
        selections = ['chainID A and resid 35', 'chainID A and resid 52', 'chainID A and resid 48']
        include_distances = True
        buriedness_file = ['reference_output/lyso_test_default/lyso_test_default_buriedness.csv',
                           'reference_output/lyso_test_default/lyso_test_default_buriedness.csv']
        buriedness_df = make_pka_or_buriedness_df(buriedness_file, trajectory_file)
        include_buriedness = True
        return (topology_file, trajectory_file, pka_df, selections, include_distances, buriedness_df, include_buriedness)
    @pytest.mark.parametrize("topology_file, trajectory_file, pka_df, selections, include_distances, buriedness_df, include_buriedness", [
     parm()])
    def test_pca(self, topology_file, trajectory_file, pka_df, selections, include_distances, buriedness_df, include_buriedness):
        clustering_matrix, times, frames, trajectory_names, trajectory_dict = create_clustering_matrix(topology_file=topology_file,
                                 trajectory_file=trajectory_file,
                                 pka_df=pka_df,
                                 selections=selections,
                                 include_distances=include_distances,
                                 buriedness_df=buriedness_df,
                                 include_buriedness=include_buriedness)
        n_components, cummulative_variance, clustering_matrix_transformed = pca(clustering_matrix)
        assert cummulative_variance >= 90
        assert clustering_matrix_transformed.shape == (clustering_matrix.shape[0],n_components)