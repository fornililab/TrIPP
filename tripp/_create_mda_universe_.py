import MDAnalysis as mda 
import numpy as np 
from tripp._correction_dictionary_ import corrected_amino_acids, corrected_atom_names
from tripp._propka_input_check_ import input_check
def create_mda_universe(topology_file, trajectory_file): 

    if trajectory_file == None: 
        universe = mda.Universe(topology_file) 
    
    else: 
        universe = mda.Universe(topology_file, trajectory_file) 
    
    topology = universe._topology 
    
    #Check if chainID is empty or not, if so default chain A for the whole system.
    if '' in topology.chainIDs.values:
        topology.chainIDs.values = np.full(len(topology.chainIDs.values),'A',dtype=str) 
        print('Your topology file contains no chain identity for at least one atom. Will add chain A for your whole system by default') 
    
    if hasattr(universe.atoms, 'formalcharges') == False: 
        universe.add_TopologyAttr('formalcharges', np.full(len(topology.chainIDs.values),0,dtype=float))
        print('Your topology file contains no formal charges. Will set formal charges to 0 for your whole system by default') 
    
    return universe 

def create_propka_compatible_universe(universe):
    propka_resnames = [corrected_amino_acids[resname] if resname in corrected_amino_acids else resname for resname in universe.residues.resnames]
    propka_atom_names = [corrected_atom_names[name] if name in corrected_atom_names else name for names in universe.residues.names for name in names]
    
    universe.residues.resnames = propka_resnames
    universe.atoms.names = propka_atom_names
    
    input_check(universe)
    return universe    
    