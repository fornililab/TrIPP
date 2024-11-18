import MDAnalysis as mda 
import numpy as np 
from tripp._correction_dictionary_ import corrected_amino_acids, corrected_atom_names
from tripp._propka_input_check_ import input_check
import logging

logger = logging.getLogger(__name__)
def create_mda_universe(topology_file, trajectory_file):
    
    if trajectory_file is None:
        universe = mda.Universe(topology_file)
    
    else: 
        universe = mda.Universe(topology_file, trajectory_file) 
    
    topology = universe._topology
    
    # Check if chainID is empty or not, if so default chain A 
    # for the whole system.
    if '' in topology.chainIDs.values:
        topology.chainIDs.values = np.full(len(topology.chainIDs.values),
                                           'A',
                                           dtype=str) 
        logger.info('Your topology file contains no chain identity for at least one atom. Will add chain A for your whole system by default')
    
    if hasattr(universe.atoms, 'formalcharges') is False:
        universe.add_TopologyAttr('formalcharges', np.full(len(topology.chainIDs.values), 0, dtype=float))
        logger.info('Your topology file contains no formal charges. Will set formal charges to 0 for your whole system by default')
    
    added_universe = universe.copy()
    return added_universe

def create_propka_compatible_universe(universe, hetatm_resid, terminal_oxygens, custom_resname_correction):
    corrected_universe = universe.copy()
    # When user provided custom_resname_correction, append them to the dictionary temporarily.
    if custom_resname_correction is not None:
        for current_resname, correct_resname in custom_resname_correction.items():
            corrected_amino_acids[current_resname] = correct_resname
            
    # Correct resname to propka compatible with predefined list.
    corrected_resnames = []
    corrected_resnames_trace = []
    for resname, resid in zip(corrected_universe.residues.resnames, corrected_universe.residues.resids):
        if resname in corrected_amino_acids:
            corrected_resnames_trace.append(f'{resname}->{corrected_amino_acids[resname]}')
            resname = corrected_amino_acids[resname]
        corrected_resnames.append(resname)
        
    corrected_universe.residues.resnames = corrected_resnames
    
    if len(corrected_resnames_trace) > 0:
        corrected_resnames_trace = '\n'.join(set(corrected_resnames_trace))
        logger.info(f"""The resnames of the following residues have been modified to be compatible with PROPKA:
{corrected_resnames_trace}""")
    
    # When Hetatm resid is provided, for example the system contain ligand, it will be corrected to HETATM.
    if hetatm_resid is not None:
        correction_ag = corrected_universe.select_atoms(f'resid {hetatm_resid}')
        correction_ag.atoms.record_types = np.full(len(correction_ag.atoms.record_types),'HETATM')
        logger.info(f'The record type of residue {correction_ag.residues.resnames[0]}{correction_ag.residues.resids[0]} is modified to HETATM')
    
    # When terminal_oxygens is provided by user as list of string,
    # the terminal oxygens will be renamed according to the dictionary
    # in _correction_dictionary_.corrected_atom_names.
    if terminal_oxygens is None:
        last_residue_names = corrected_universe.residues[-1].atoms.names
        correct_terminal_oxygens = ['O', 'OXT']
        corrected_terminal_oxygens_trace = []
        for current_name, correct_name in corrected_atom_names.items():
            if current_name in last_residue_names:
                last_residue_names[last_residue_names == current_name] = correct_name
                corrected_terminal_oxygens_trace.append(f'{current_name} -> {correct_name}')
        corrected_universe.residues[-1].atoms.names = last_residue_names
        corrected_terminal_oxygens_trace = '\n'.join(corrected_terminal_oxygens_trace)
        logger.info(f"""Terminal oxygens will be modified with our predefined dictionary:
{corrected_terminal_oxygens_trace}""")
    # Check if argument is a list
    elif not isinstance(terminal_oxygens,list):
        raise TypeError(f'terminal_oxygens argument must be a list, not {type(terminal_oxygens)}')
    # Check if argument is a list of string
    elif not all(isinstance(element,str) for element in terminal_oxygens):
        raise TypeError(f'terminal_oxygens argument must be a list of strings, at least one of your element is not a string')
    # Check if arguments are of length 2
    elif len(terminal_oxygens) != 2:
        raise ValueError(f'Length of terminal_oxygens is not 2, but {len(terminal_oxygens)}')

    # When terminal_oxygens is provided by user as list of string,
    # the terminal oxygens names will be renamed to O and OXT.
    elif isinstance(terminal_oxygens, list) and len(terminal_oxygens) == 2:
        last_residue_names = corrected_universe.residues[-1].atoms.names
        correct_terminal_oxygens = ['O', 'OXT']
        
        corrected_terminal_oxygens_trace = []
        for current_name, correct_name in zip(terminal_oxygens, correct_terminal_oxygens):
            if current_name in last_residue_names:
                last_residue_names[last_residue_names == current_name] = correct_name
                corrected_terminal_oxygens_trace.append(f'{current_name} -> {correct_name}')
                
        corrected_universe.residues[-1].atoms.names = last_residue_names
        
        corrected_terminal_oxygens_trace = '\n'.join(corrected_terminal_oxygens_trace)
        logger.info(f"""Terminal oxygens list supplied and are modified:
{corrected_terminal_oxygens_trace}""")
    
    # Check the corrected_universe for PROPKA compatibility.
    input_check(corrected_universe)
    return corrected_universe