from tripp import Trajectory
from tripp import Visualization
from tripp import Clustering
from tripp.analysis import PCProjectionScreening  
import numpy as np
import os
import shutil

def compare_output(file1, file2):
    with open(file1, 'r') as f1, open(file2, 'r') as f2:
        fileone = f1.readlines()
        filetwo = f2.readlines()
        for line in filetwo:
            if line not in fileone:
                raise AssertionError('one of your')
    return True

if __name__ == '__main__':
    if os.path.isdir('test_output'):
        shutil.rmtree('test_output')
    
    # Trajectory class test
    topology_file = 'files/test.pdb'
    trajectory_file = 'files/test.xtc'
    output_directory = 'test_output/MD1'
    output_prefix = 'MD1'
    TrIPP_Traj = Trajectory(topology_file=topology_file,  
                                trajectory_file=trajectory_file,  
                                output_directory = output_directory ,  
                                output_prefix = output_prefix,  
                                hetatm_resname=None,  
                                custom_terminal_oxygens=None, 
                                custom_resname_correction=None,  
                                cpu_core_number=-1)
    TrIPP_Traj.run(extract_buriedness_data=True,
                    chain='A',
                    mutation_selections=None,
                    disulphide_bond_detection=True,
                    optargs=[])

    output_directory = 'test_output/MD1_44'
    output_prefix = 'MD1_44'
    TrIPP_Traj = Trajectory(topology_file=topology_file,  
                                trajectory_file=trajectory_file,  
                                output_directory = output_directory ,  
                                output_prefix = output_prefix,  
                                hetatm_resname=None,  
                                custom_terminal_oxygens=None, 
                                custom_resname_correction=None,  
                                cpu_core_number=-1)
    TrIPP_Traj.run(extract_buriedness_data=True,
                    chain='A',
                    mutation_selections='chainID A and resid 44',
                    disulphide_bond_detection=True,
                    optargs=[])

    try:
        assert compare_output('test_output/MD1/MD1_pka_chainA.csv','reference_output/MD1/MD1_pka_chainA.csv')
        assert compare_output('test_output/MD1/MD1_buriedness_chainA.csv','reference_output/MD1/MD1_buriedness_chainA.csv')
        assert compare_output('test_output/MD1_44/MD1_44_pka_chainA.csv','reference_output/MD1_44/MD1_44_pka_chainA.csv')
        assert compare_output('test_output/MD1_44/MD1_44_buriedness_chainA.csv','reference_output/MD1_44/MD1_44_buriedness_chainA.csv')
        print('Trajectory class test passed')
    except AssertionError as e:
        raise e

    # Clustering class test
    output_directory = 'test_output'
    pka_file = 'test_output/MD1/MD1_pka_chainA.csv'
    TrIPP_Clust = Clustering(topology_file=topology_file,
                            trajectory_file=trajectory_file,
                            pka_file=pka_file, 
                            buriedness_file=None,
                            selections=['chainID A and resid 35','chainID A and resid 52'], 
                            output_directory=output_directory,
                            output_prefix='All_md', 
                            include_distances=False,  
                           include_buriedness=False,
                            dimensionality_reduction=False)


    np.savetxt('test_output/All_md_clustering_matrix.csv',TrIPP_Clust.clustering_matrix,delimiter=',',fmt='%.3f')


    try:
        assert compare_output('test_output/All_md_clustering_matrix.csv','reference_output/All_md_clustering_matrix.csv')
        print('Clustering class test passed')
    except AssertionError as e:
        raise e

    # Visualization class test
    PyMOL_path = None
    while PyMOL_path is None:
        PyMOL_path = input('Please provide path to PyMOL executable to test Visualization class:')
    TrIPP_Vis = Visualization(topology_file=topology_file, 
                              pka_file=pka_file) 
    
    TrIPP_Vis.gen_pse(pymol_path=PyMOL_path,
                    output_directory=output_directory,
                    output_prefix='All_md_pka',
                    chain = 'A',
                    coloring_method='mean',  
                    lower_limit=0,  
                    upper_limit=14,  
                    color_palette='red_white_blue')  #


    try:
        assert compare_output('test_output/All_md_pka_mean.pdb','reference_output/All_md_pka_mean.pdb')
        print('Visualization class test passed')
    except AssertionError as e:
        raise e
    except FileNotFoundError:
        if not os.path.isfile('test_output/All_md_pka_mean.pdb') or not os.path.isfile('test_output/All_md_pka_mean.pse'):
            raise TypeError('Path to PyMOL executable is incorrect. Visualization class test failed')

    # Projection and pKa correlation scanning
    PCProjectionScreening(output_directory=output_directory,
                          output_prefix='MD1_pc1_projection_pKa_PearsonCorrelation',
                          pka_file='test_output/MD1/MD1_pka_chainA.csv',
                          projection_file='files/test_pc1.csv',
                          method='Pearson')
    Corr_Vis = Visualization(topology_file=topology_file,
                             correlation_file='test_output/MD1_pc1_projection_pKa_PearsonCorrelation.csv')

    Corr_Vis.gen_pse(pymol_path=PyMOL_path, 
                     output_directory=output_directory,
                     output_prefix='MD1_pc1_projection_pKa_PearsonCorrelation',
                     coloring_method='correlation',
                     lower_limit=-1,
                     upper_limit=1,
                     correlation_threshold=0.5,
                     color_palette='red_white_blue')
    try:
        assert compare_output('test_output/MD1_pc1_projection_pKa_PearsonCorrelation.csv','reference_output/MD1_pc1_projection_pKa_PearsonCorrelation.csv')
        assert compare_output('test_output/MD1_pc1_projection_pKa_PearsonCorrelation_correlation.pdb','reference_output/MD1_pc1_projection_pKa_PearsonCorrelation_correlation.pdb')
        print('Correlation scanning test passed')
    except AssertionError as e:
        raise e
