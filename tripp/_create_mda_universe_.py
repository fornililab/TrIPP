import MDAnalysis as mda 
import numpy as np 

def create_mda_universe(topology_file, trajectory_file): 

    universe = mda.Universe(topology_file, trajectory_file) 
    topology = universe._topology 
    
    #Check if chainID is empty or not, if so default chain A for the whole system.
    if '' in topology.chainIDs.values:
        topology.chainIDs.values = np.full(len(topology.chainIDs.values),'A',dtype=str) 
        print('Your topology file contains no chain identity. Will add chain A for your whole system by default') 
    
    return universe 
