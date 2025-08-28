from tripp import Trajectory
from tripp import Clustering
from tripp import Visualization
from tripp.analysis import PCProjectionScreening  
import numpy as np
import pytest
import os
import shutil

# Clean up test_output directory if it exists
test_output_dir = 'test_output'
if os.path.isdir(test_output_dir):
    # Remove all contents and the directory itself
    shutil.rmtree(test_output_dir)

def get_data_path(filename):
    return os.path.join(os.path.dirname(__file__), "data", filename)
    
topology_file = get_data_path('lyso.pdb')
trajectory_file = get_data_path('lyso_test.xtc')
pka_file = 'test_output/lyso_test_default/lyso_test_default_pka.csv'

PyMOL_path = False
while not PyMOL_path:
    PyMOL_path = input('Please provide path to PyMOL executable to test Visualization class or type "skip":')

class TestInstallation:
    def compare_output(file1, file2):
        with open(file1, 'r') as f1, open(file2, 'r') as f2:
            fileone = f1.readlines()
            filetwo = f2.readlines()
            for line in filetwo:
                if line not in fileone:
                    raise AssertionError('Two files are not identical, installation failed.')
        return True
    def tripp_traj_default(output_directory, output_prefix):
        TrIPP_Traj = Trajectory(topology_file=topology_file,  
                                    trajectory_file=trajectory_file,  
                                    output_directory = output_directory,  
                                    output_prefix = output_prefix,  
                                    hetatm_resname=None,  
                                    custom_terminal_oxygens=None, 
                                    custom_resname_correction=None,  
                                    cpu_core_number=-1)
        return TrIPP_Traj

    
    def test_Trajectory_lyso_default(self):
        output_directory = 'test_output/lyso_test_default'
        output_prefix = 'lyso_test_default'
        TrIPP_Traj = TestInstallation.tripp_traj_default(output_directory, output_prefix)
        TrIPP_Traj.run(extract_buriedness_data=True,
                       mutation_selections=None,
                       save_disulphide_pka=False,
                       optargs=[])
        assert TestInstallation.compare_output('test_output/lyso_test_default/lyso_test_default_pka.csv','reference_output/lyso_test_default/lyso_test_default_pka.csv')
        assert TestInstallation.compare_output('test_output/lyso_test_default/lyso_test_default_buriedness.csv','reference_output/lyso_test_default/lyso_test_default_buriedness.csv')
    
    def test_Trajectory_lyso_mutation(self):
        output_directory = 'test_output/lyso_test_44'
        output_prefix = 'lyso_test_44'
        TrIPP_Traj = TestInstallation.tripp_traj_default(output_directory, output_prefix)
        TrIPP_Traj.run(extract_buriedness_data=True,
                       mutation_selections='chainID A and resid 44',
                       save_disulphide_pka=False,
                       optargs=[])
        assert TestInstallation.compare_output('test_output/lyso_test_44/lyso_test_44_pka.csv','reference_output/lyso_test_44/lyso_test_44_pka.csv')
        assert TestInstallation.compare_output('test_output/lyso_test_44/lyso_test_44_buriedness.csv','reference_output/lyso_test_44/lyso_test_44_buriedness.csv')
        
    def test_Clustering_lyso(self):
        output_directory = 'test_output'
        TrIPP_Clust = Clustering(topology_file=topology_file,
                                 trajectory_file=trajectory_file,
                                 pka_file=pka_file,
                                 buriedness_file=None,
                                 selections=['chainID A and resid 35', 'chainID A and resid 52'],
                                 output_directory=output_directory,
                                 output_prefix='lyso_test',
                                 include_distances=False,
                                 include_buriedness=False,
                                 dimensionality_reduction=False)
        ref_clustering_matrix = np.loadtxt('reference_output/lyso_test_clustering_matrix.csv', delimiter=',')
        assert np.all(ref_clustering_matrix == np.round(TrIPP_Clust.clustering_matrix,3))
    
    def test_kmedoids_lyso(self):
        output_directory = 'test_output'
        TrIPP_Clust = Clustering(topology_file=topology_file,
                            trajectory_file=trajectory_file,
                            pka_file=pka_file,
                            buriedness_file=None,
                            selections=['chainID A and resid 35', 'chainID A and resid 52'],
                            output_directory=output_directory,
                            output_prefix='lyso_test',
                            include_distances=False,
                            include_buriedness=False,
                            dimensionality_reduction=False)
        TrIPP_Clust.kmedoids(n_clusters=2, random_state=1)
        
        assert TestInstallation.compare_output('test_output/lyso_test_KMedoids_cluster.csv','reference_output/lyso_test_KMedoids_cluster.csv')
    
    @pytest.mark.skipif(PyMOL_path.lower() == 'skip', reason='Skipping Visualization test')
    def test_Visualization_lyso(self):
        output_directory = 'test_output'
        
        TrIPP_Vis = Visualization(topology_file=topology_file, 
                                  pka_file=pka_file) 
        TrIPP_Vis.gen_pse(pymol_path=PyMOL_path,
                          output_directory=output_directory,
                          output_prefix='lyso_test',
                          coloring_method='mean',
                          lower_limit=0,
                          upper_limit=14,
                          color_palette='red_white_blue')
        
        assert TestInstallation.compare_output('test_output/lyso_test_mean.pdb','reference_output/lyso_test_mean.pdb')
        assert os.path.isfile('test_output/lyso_test_mean.pse'), 'PyMOL session file was not created, please check your PyMOL path.'

    def test_PCProjectionScreening(self):
        output_directory = 'test_output'
        PCProjectionScreening(output_directory=output_directory,
                              output_prefix='lyso_test_PearsonCorrelation',
                              pka_file=pka_file,
                              projection_file='data/lyso_test_pc1.csv',
                              method='Pearson')
        assert TestInstallation.compare_output('test_output/lyso_test_PearsonCorrelation.csv','reference_output/lyso_test_PearsonCorrelation.csv')
    
    @pytest.mark.skipif(PyMOL_path.lower() == 'skip', reason='Skipping Visualization test')
    def test_PCProjectionScreening_Visualization(self):
        output_directory = 'test_output'
        Corr_Vis = Visualization(topology_file=topology_file,
                             correlation_file='test_output/lyso_test_PearsonCorrelation.csv')

        Corr_Vis.gen_pse(pymol_path=PyMOL_path, 
                        output_directory=output_directory,
                        output_prefix='lyso_test_PearsonCorrelation',
                        coloring_method='correlation',
                        lower_limit=-1,
                        upper_limit=1,
                        correlation_threshold=0.5,
                        color_palette='red_white_blue')
        
        assert TestInstallation.compare_output('test_output/lyso_test_PearsonCorrelation_correlation.pdb','reference_output/lyso_test_PearsonCorrelation_correlation.pdb')
        assert os.path.isfile('test_output/lyso_test_PearsonCorrelation_correlation.pse'), 'PyMOL session file was not created, please check your PyMOL path.'