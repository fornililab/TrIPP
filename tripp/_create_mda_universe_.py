import MDAnalysis as mda 
import numpy as np 
from tripp._correction_dictionary_ import corrected_amino_acids, corrected_atom_names
from tripp._propka_input_check_ import check_resname_HETATM, check_terminal_oxygens
import logging
import pandas as pd

logger = logging.getLogger(__name__)
def create_mda_universe(topology_file, trajectory_file):
    """
    Create an MDAnalysis Universe object from a topology file and an 
    optional trajectory file. Checks for chain identity and formal 
    charges, and adds default values if necessary. Note that the 
    check for chain identity is not exhaustive.
    
    Parameters
    ----------
    topology_file : str
        Path to the topology file.
    trajectory_file : str or None
        Path to the trajectory file. If None, only the topology is loaded.
    Return
    -------
    universe : MDAnalysis.universe
        MDAnalysis Universe object containing the topology and trajectory 
        data.
    Raises
    ------
    Exception
        If the topology file contains residues with residue IDs but the
        chain IDs are empty, or if the duplicated residue IDs have the
        same chain ID.
    """
    if trajectory_file is None:
        universe = mda.Universe(topology_file)
    
    else: 
        universe = mda.Universe(topology_file, trajectory_file) 
    
    topology = universe._topology
    _, counts = np.unique(universe.residues.resids, return_counts=True)
    # Check if chainID of any residue is empty or not, 
    # if so default chain A for the whole system.
    if '' in topology.chainIDs.values and (counts == 1).all():
        topology.chainIDs.values = np.full(len(topology.chainIDs.values),
                                           'A',
                                           dtype=str) 
        logger.info(f"""Your topology file contains no chain identity for at least one atom. 
Will add chain A for your whole system by default""")
    elif '' in topology.chainIDs.values and (counts > 1).any():
        raise Exception('Your system contain duplicated resid and at least one atom does not have a chain ID, please assign chain IDs before running TrIPP.')
    
    elif len(set(topology.chainIDs.values)) == 1 and (counts > 1).any():
        raise Exception('Your system contain duplicated resid, but all atoms in your system have the same chain ID, please re-assign the chain IDs before running TrIPP')
        
    if hasattr(universe.atoms, 'formalcharges') is False:
        universe.add_TopologyAttr('formalcharges', np.full(len(topology.chainIDs.values), 0, dtype=float))
        logger.info("""Your topology file contains no formal charges. 
Will set formal charges to 0 for your whole system by default.""")
    
    added_universe = universe.copy()
    return added_universe

def create_propka_compatible_universe(universe,
                                      hetatm_resname,
                                      custom_terminal_oxygens,
                                      custom_resname_correction):
    """
    Converting the loaded universe from create_mda_universe function 
    to PROPKA compatible format, including the resname of standard amino 
    acids and atom names of terminal oxygens. The record type of a resname 
    can also be modified with user's input in the case of a ligand. Checking 
    will be performed for resname, record type and terminal oxygens. See 
    _propka_input_check_ for more information.
    
    Parameters
    ----------
    universe: MDAnalysis.universe
        MDAnalysis universe that has been processed through create_mda_universe.
    hetatm_resname: str, list, default=None
        Resname of the hetatm that need to be modified. See log for info
        if encountered error.
    custom_terminal_oxygens: list, default=None
        PROPKA only recognizes O and OXT for the terminal oxygens. Either 
        provide a list of length 2 containing the atom name of your terminal
        oxygens, or None for correction according to our pre-defined
        dictionary (See _correction_dictionary_.py for the pre-defined
        name)
    custom_resname_correction: dict, default=None
        Our pre-defined dictionary might not contain all resname to be
        corrected. (See _correction_dictionary_.py for the pre-defined
        resname) Use this parameter to add a custom resname correction.
        ie: {'ASPH':'ASP'}
    Returns
    -------
    corrected_universe: MDAnalysis.universe
        MDAnalysis universe that is compatible with PROPKA, including
        20 standard amino acid resnames, O and OXT for C-terminal(s), 
        HETATM for non-protein atoms (when hetatm_resname is used).
    """
    corrected_universe = universe.copy()
    # When user provided custom_resname_correction, append them to the dictionary temporarily.
    if custom_resname_correction is not None:
        for current_resname, correct_resname in custom_resname_correction.items():
            corrected_amino_acids[current_resname] = correct_resname
            
    # Correct resname to propka compatible with predefined list.
    corrected_resnames = []
    corrected_resnames_trace = []
    for resname, residx in zip(corrected_universe.residues.resnames,
                              corrected_universe.residues.resindices):
        if resname in corrected_amino_acids:
            corrected_resnames_trace.append(f'{resname}->{corrected_amino_acids[resname]}')
            resname = corrected_amino_acids[resname]
        # Special care taken for MSE -> MET, where the SE atomname needs to be modified as well.
        elif resname == 'MSE':
            corrected_resnames_trace.append(f'{resname}->MET (SE->SD)')
            resname = 'MET'
            atom_names = corrected_universe.select_atoms(f'resindex {residx}').atoms.names
            atom_types = corrected_universe.select_atoms(f'resindex {residx}').atoms.types
            corrected_universe.select_atoms(f'resindex {residx}').atoms.names = np.char.replace(atom_names.astype(str), 'SE', 'SD')
            corrected_universe.select_atoms(f'resindex {residx}').atoms.types = np.char.replace(atom_types.astype(str), 'SE', 'S')
            
        corrected_resnames.append(resname)
        
    corrected_universe.residues.resnames = corrected_resnames
    
    if len(corrected_resnames_trace) > 0:
        corrected_resnames_trace = '\n'.join(set(corrected_resnames_trace))
        logger.info(f"""The resnames of the following residues have been modified to be compatible with PROPKA:
{corrected_resnames_trace}""")
    
    # When hetatm_resname is provided, for example if the system contains ligand, 
    # the record type will be corrected to HETATM.
    if hetatm_resname is not None:
        if isinstance(hetatm_resname, list):
            hetatm_resname = ' '.join([str(i) for i in hetatm_resname])
        correction_ag = corrected_universe.select_atoms(f'resname {hetatm_resname}')
        correction_ag.atoms.record_types = np.full(len(correction_ag.atoms.record_types),'HETATM')
        corrected_hetatm = []
        for resid, resname in zip(correction_ag.residues.resids, correction_ag.residues.resnames):
            corrected_hetatm.append(f'{resname}{resid}')
        logger.info(f"""The record type of the following residues has been modified to HETATM:
{', '.join(corrected_hetatm)}""")
    
    check_resname_HETATM(corrected_universe.select_atoms('not protein'))
    
    # When custom_terminal_oxygens is not provided by user as list of string,
    # the terminal oxygens will be renamed according to the dictionary
    # in _correction_dictionary_.corrected_atom_names.
    if custom_terminal_oxygens is None:
        corrected_universe, corrected_terminal_oxygens_trace = modifiy_terminal_oxygens(corrected_universe, corrected_atom_names)
        logger.info(f"""Terminal oxygens will be modified with our predefined dictionary:
{corrected_terminal_oxygens_trace}""")
        
    # Check if argument is a list
    elif not isinstance(custom_terminal_oxygens,list):
        raise TypeError(f'custom_terminal_oxygens argument must be a list, not {type(custom_terminal_oxygens)}')
    
    # Check if argument is a list of string
    elif not all(isinstance(element,str) for element in custom_terminal_oxygens):
        raise TypeError(f'custom_terminal_oxygens argument must be a list of strings, at least one of your element is not a string')
    
    # Check if argument is a list of length 2
    elif len(custom_terminal_oxygens) != 2:
        raise ValueError(f'Length of custom_terminal_oxygens is not 2, but {len(custom_terminal_oxygens)}')

    # When custom_terminal_oxygens is provided by user as a list of string,
    # the terminal oxygens names will be renamed to O and OXT.
    elif isinstance(custom_terminal_oxygens, list) and len(custom_terminal_oxygens) == 2:
        correct_terminal_oxygens = ['O', 'OXT']
        to_be_corrected_dict = {i:j for i,j in zip(custom_terminal_oxygens, correct_terminal_oxygens)}
        corrected_universe, corrected_terminal_oxygens_trace = modifiy_terminal_oxygens(corrected_universe, to_be_corrected_dict)
        
        logger.info(f"""custom_terminal_oxygens supplied and are modified:
{corrected_terminal_oxygens_trace}""")
    
    check_terminal_oxygens(corrected_universe)
        
    return corrected_universe

def modifiy_terminal_oxygens(corrected_universe, to_be_corrected_dict):
    """
    Function to modify the terminal oxygens in the corrected_universe 
    based on to_be_corrected_dict.
    
    Parameters
    ----------
    corrected_universe: MDAnalysis.universe
        Universe that is to be corrected.
    to_be_corrected_dict: dict
        Dictionary for which the keys are the atom names of the terminal oxygens
        in your topology, and the value is the atom name of O or OXT.
        Provided through custom_terminal_oxygens/corrected_atom_names.
    Return
    ----------
    corrected_universe: MDAnalysis.universe
        Universe that the atom names of terminal oxygens have been corrected.
    corrected_terminal_oxygens_trace: str
        Information of the processed terminal oxygens.
    """
    ag = corrected_universe.select_atoms('protein')
    corrected_terminal_oxygens_trace = []
    for index, name in zip(ag.residues.resindices, ag.residues.names):
        if np.isin([i.strip() for i in list(to_be_corrected_dict.keys())], name).any():
            corrected_name = pd.Series(name).replace(to_be_corrected_dict).to_numpy()
            ag_to_correct = corrected_universe.select_atoms(f'resindex {index}')
            before_correction = ag_to_correct.atoms.names.copy()
            ag_to_correct.atoms.names = corrected_name
            after_correction = ag_to_correct.atoms.names.copy()
            pnames = np.setdiff1d(before_correction, after_correction)
            anames = np.setdiff1d(after_correction, before_correction)
            resid = ag_to_correct.residues.resids[0]
            resname = ag_to_correct.residues.resnames[0]
            corrected_terminal_oxygens_trace.append(f"{resname}{resid}: {', '.join([f'{pname} -> {aname.strip()}' for pname,aname in zip(pnames,anames)])}")
    corrected_terminal_oxygens_trace = '\n'.join(corrected_terminal_oxygens_trace)
    return corrected_universe, corrected_terminal_oxygens_trace

