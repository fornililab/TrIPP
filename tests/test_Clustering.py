from tripp.Clustering import Clustering
import pandas as pd
import numpy as np
import os
import pytest
import glob

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

def parm1():
    topology_file = get_data_path('lyso.pdb')
    trajectory_file = {
        "traj1": get_data_path('lyso_test.xtc'),
        "traj2": get_data_path('lyso_test.xtc')}
    pka_file = ['reference_output/lyso_test_default/lyso_test_default_pka.csv',
                'reference_output/lyso_test_default/lyso_test_default_pka.csv']
    selections = ['chainID A and resid 35', 'chainID A and resid 52', 'chainID A and resid 48']
    include_distances = True
    buriedness_file = ['reference_output/lyso_test_default/lyso_test_default_buriedness.csv',
                        'reference_output/lyso_test_default/lyso_test_default_buriedness.csv']
    include_buriedness = True
    dimensionality_reduction = True
    output_directory = None
    output_prefix = None
    return (topology_file, trajectory_file, pka_file, selections, output_directory, output_prefix, buriedness_file, include_distances, include_buriedness, dimensionality_reduction)

@pytest.fixture(params=[parm1()])
def return_param(request):
    return request.param
class TestClustering():
    def test_Clustering(self, return_param):
        output_directory = check_dir('test_output/TestClustering')
        output_prefix = 'lyso_test_clustering'
        params = list(return_param)
        params[4] = output_directory
        params[5] = output_prefix
        Clust = Clustering(*params)
        
        selections = params[3]
        include_distances = params[7]
        include_buriedness = params[8]
        feature_count = len(selections)
        if include_distances:
            feature_count += (len(selections)*(len(selections)-1))/2
        if include_buriedness:
            feature_count += len(selections)
        assert (Clust.times.flatten() == np.array(Clust.pka_df.index)).all()
        assert (Clust.frames.flatten() == np.hstack([np.arange(len(Clust.trajectory_dict[i].trajectory)) for i in Clust.trajectory_dict.keys()])).all()
        assert (Clust.trajectory_names.flatten() == np.hstack([np.repeat(i, len(Clust.trajectory_dict[i].trajectory)) for i in Clust.trajectory_dict.keys()])).all()
        assert Clust.cummulative_variance >= 90
        assert Clust.clustering_matrix.shape == (Clust.clustering_matrix.shape[0],Clust.n_components)
    
    def test_kmedoids_auto(self, return_param):
        output_directory = check_dir('test_output/TestClustering/TestKMedoidsAuto')
        output_prefix = 'lyso_test_clustering'
        params = list(return_param)
        params[4] = output_directory
        params[5] = output_prefix
        Clust = Clustering(*params)
        Clust.kmedoids(automatic=True)
        assert os.path.exists(f'{output_directory}/{output_prefix}_KMedoids_cluster.csv')
        df = pd.read_csv(f'{output_directory}/{output_prefix}_KMedoids_cluster.csv')
        assert len(glob.glob(f'{output_directory}/*_KMedoids_C*.pdb')) == len(df['Clusters'].unique())

    def test_greedy_auto(self, return_param):
        output_directory = check_dir('test_output/TestClustering/TestGreedyAuto')
        output_prefix = 'lyso_test_clustering'
        params = list(return_param)
        params[4] = output_directory
        params[5] = output_prefix
        Clust = Clustering(*params)
        Clust.greedy(automatic=True)
        assert os.path.exists(f'{output_directory}/{output_prefix}_greedy_cluster.csv')
        df = pd.read_csv(f'{output_directory}/{output_prefix}_greedy_cluster.csv')
        assert len(glob.glob(f'{output_directory}/*_greedy_C*.pdb')) == len(df['Clusters'].unique())

    def test_dbscan_auto(self, return_param):
        output_directory = check_dir('test_output/TestClustering/TestDBSCANAuto')
        output_prefix = 'lyso_test_clustering'
        params = list(return_param)
        params[4] = output_directory
        params[5] = output_prefix
        Clust = Clustering(*params)
        Clust.dbscan(automatic=True)
        assert os.path.exists(f'{output_directory}/{output_prefix}_DBSCAN_cluster.csv')
        df = pd.read_csv(f'{output_directory}/{output_prefix}_DBSCAN_cluster.csv')
        if -1 in df['Clusters'].values:
            assert len(glob.glob(f'{output_directory}/*_DBSCAN_C*.pdb')) == len(df['Clusters'].unique()) - 1
        else:
            assert len(glob.glob(f'{output_directory}/*_DBSCAN_C*.pdb')) == len(df['Clusters'].unique())

    def test_hdbscan_auto(self, return_param):
        output_directory = check_dir('test_output/TestClustering/TestHDBSCANAuto')
        output_prefix = 'lyso_test_clustering'
        params = list(return_param)
        params[4] = output_directory
        params[5] = output_prefix
        Clust = Clustering(*params)
        Clust.hdbscan(automatic=True,min_cluster_size_range=(5,50,5), min_samples_range=(5,60,1))
        assert os.path.exists(f'{output_directory}/{output_prefix}_HDBSCAN_cluster.csv')
        df = pd.read_csv(f'{output_directory}/{output_prefix}_HDBSCAN_cluster.csv')
        if -1 in df['Clusters'].values:
            assert len(glob.glob(f'{output_directory}/*_HDBSCAN_C*.pdb')) == len(df['Clusters'].unique()) - 1
        else:
            assert len(glob.glob(f'{output_directory}/*_HDBSCAN_C*.pdb')) == len(df['Clusters'].unique())
    