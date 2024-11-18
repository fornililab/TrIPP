from tripp import Trajectory
from tripp import Visualization
from tripp import Clustering
import numpy as np
import os

def compare_output(file1, file2):
    with open(file1, 'r') as f1, open(file2, 'r') as f2:
        fileone = f1.readlines()
        filetwo = f2.readlines()
        for line in filetwo:
            if line not in fileone:
                raise AssertionError('one of your')
    return True

if __name__ == '__main__':
    # Trajectory class test
    topology_file = 'files/1AKI_minstr.pdb'
    trajectory_file = 'files/1AKI_md1_300ns.xtc'
    output_directory = 'test_output/MD1'
    output_prefix = 'MD1'
    TrIPP_Traj = Trajectory(topology_file=topology_file,  
                                trajectory_file=trajectory_file,  
                                output_directory = output_directory ,  
                                output_prefix = output_prefix,  
                                hetatm_resid=None,  
                                custom_terminal_oxygens=None, 
                                custom_resname_correction=None,  
                                cpu_core_number=-1)
    TrIPP_Traj.run(extract_buriedness_data=True,
                    chain='A',
                    mutation=None,
                    disulphide_bond_detection=True,
                    optargs=[])

    topology_file = 'files/1AKI_minstr.pdb'
    trajectory_file = 'files/1AKI_md1_300ns.xtc'
    output_directory = 'test_output/MD1_44'
    output_prefix = 'MD1_44'
    TrIPP_Traj = Trajectory(topology_file=topology_file,  
                                trajectory_file=trajectory_file,  
                                output_directory = output_directory ,  
                                output_prefix = output_prefix,  
                                hetatm_resid=None,  
                                custom_terminal_oxygens=None, 
                                custom_resname_correction=None,  
                                cpu_core_number=-1)
    TrIPP_Traj.run(extract_buriedness_data=True,
                    chain='A',
                    mutation=44,
                    disulphide_bond_detection=True,
                    optargs=[])

    try:
        assert compare_output('test_output/MD1/MD1_pka.csv','reference_output/MD1/MD1_pka.csv')
        assert compare_output('test_output/MD1/MD1_buriedness.csv','reference_output/MD1/MD1_buriedness.csv')
        assert compare_output('test_output/MD1_44/MD1_44_pka.csv','reference_output/MD1_44/MD1_44_pka.csv')
        assert compare_output('test_output/MD1_44/MD1_44_buriedness.csv','reference_output/MD1_44/MD1_44_buriedness.csv')
        print('Trajectory class test passed')
    except AssertionError as e:
        raise e

    # Clustering class test
    output_directory = 'test_output'
    pka_file = 'test_output/MD1/MD1_pka.csv'
    TrIPP_Clust = Clustering(topology_file=topology_file,
                            trajectory_file=trajectory_file,
                            pka_file=pka_file, 
                            buriedness_file=None,
                            residues=[35, 52], 
                            output_directory=output_directory,
                            output_prefix='All_md', 
                            include_distances=False,  
                            include_buriedness=False,
                            dimensionality_reduction=False)


    np.savetxt('All_md_clustering_matrix.csv',TrIPP_Clust.clustering_matrix,delimiter=',',fmt='%.3f')


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
                            pka_file=pka_file) # Can be str or list of str, with the path to the pKa file(s)

    TrIPP_Vis.gen_pse(pymol_path=PyMOL_path,
                    output_directory=output_directory,
                    output_prefix='All_md_pka',
                    chain = 'A',
                    coloring_method='mean',  # Output PyMOL session with the mean pKa 
                    lower_limit=0,  # Minimum of the colour bar
                    upper_limit=14,  # Maximum of the colour bar
                    color_palette='red_white_blue')  #


    try:
        assert compare_output('test_output/All_md_pka_mean.pdb','reference_output/All_md_pka_mean.pdb')
        print('Visualization class test passed')
    except AssertionError as e:
        raise e
    except FileNotFoundError:
        if not os.path.isfile('test_output/All_md_pka_mean.pdb') or not os.path.isfile('test_output/All_md_pka_mean.pse'):
            raise TypeError('Path to PyMOL executable is incorrect. Visualization class test failed')
