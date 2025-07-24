"""
    @release_date  : $release_date
    @version       : $release_version
    @author        : Christos Matsingos, Ka Fu Man 
    
    This file is part of the TrIPP software
    (https://github.com/fornililab/TrIPP).
    Copyright (c) 2024 Christos Matsingos, Ka Fu Man and Arianna Fornili.

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, version 3.

    This program is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program. If not, see <http://www.gnu.org/licenses/>.
"""

def mutate(universe, mutation_selections, temp_name): 
    """ 
    Function that deletes all atoms of a residue except for a methyl group. 
    The residue is then renamed to alanine. 
    Parameters
    ----------
    universe: MDAnalysis.universe
        The MDAnalysis universe with the trajectory and topology loaded. 
        This universe should be corrected by create_propka_compatible_universe 
        function.
    mutation_selections: str
        A selection string in MDAnalysis selection algebra to select the 
        residue for which the mutation is to be performed.
    temp_name: str
        The name of the temporary PDB file to be created with the mutated residue.
    Exceptions
    ----------
    Exception
        If the residue type is GLY, an exception is raised because GLY cannot be mutated
        due to the absence of the CB atom.
    """ 
    replace_name = ' '.join(['N','HN','H','CA','HA','CB','O','C'])
    mutation_ag = universe.select_atoms(mutation_selections)
    if mutation_ag.residues.resnames == 'GLY':
        raise Exception('GLY cannot be mutated because of the CB not present.')
    mutation_ag.residues.resnames = 'ALA'
    # Selecting all but not the mutation_selection, and also the mutation_selection but only those of replace_name.
    mutation_ag = universe.select_atoms(f"(all and not ({mutation_selections})) or ({mutation_selections} and name {replace_name})")
    mutation_ag.write(f'{temp_name}.pdb')
