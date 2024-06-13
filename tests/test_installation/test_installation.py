from tripp import Trajectory

file_directory = '../../tutorial/files' 
output_directory = '../../tutorial/output' 
test_installation = './'

topology_file = f'{file_directory}/1AKI_clean_min_cg_pbc_p.pdb' #topology, can be any format readable by MDAnalysis

for md in range(1,4): 
    trajectory_file = f'{file_directory}/1AKI_clean_eq_md{md}_T300_300ns_p_fit_skip100.xtc' #trajectory, can be any format readable by MDAnalysis
    TrIPP_Traj = Trajectory(trajectory_file=trajectory_file, 
                            topology_file=topology_file, 
                            cpu_core_number=-1) #Setting the argument cpu_number to -1 allows for the use of all available CPUs 
    
    TrIPP_Traj.run(output_file = f'{test_installation}/1AKI_MD{md}_local_file', 
                   extract_surface_data=True, #If this is set to False, only pKa data will be extracted, else both pKa and buridness data will be extracted 
                   chain='A', 
                   mutation=None, 
                   disulphide_bond_detection=True) #Automatically detects disulphide bonds and removes the cysteines from pKa calculation 