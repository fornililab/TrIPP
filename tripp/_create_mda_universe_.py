import MDAnalysis as mda 
import numpy as np 

def create_mda_universe(topology_file, trajectory_file): 

    if trajectory_file == None: 
        universe = mda.Universe(topology_file) 
    
    else: 
        universe = mda.Universe(topology_file, trajectory_file) 
    
    topology = universe._topology 
    
    #Check if chainID is empty or not, if so default chain A for the whole system.
    if '' in topology.chainIDs.values:
        topology.chainIDs.values = np.full(len(topology.chainIDs.values),'A',dtype=str) 
        print('Your topology file contains no chain identity. Will add chain A for your whole system by default') 
    
    if hasattr(universe.atoms, 'formalcharges') == False: 
        universe.add_TopologyAttr('formalcharges', np.full(len(topology.chainIDs.values),0,dtype=float))
        print('Your topology file contains no formal charges. Will set formal charges to 0 for your whole system by default') 
    
    return universe 