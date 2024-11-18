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


def edit_pdb(pdbfile): 

    """ 
    Function that edits pdb file naming of residues so that it is compatible 
    with PROPKA. 
    """ 
    
    from tripp._correction_dictionary_ import corrected_amino_acids

    with open(pdbfile, 'r') as file: 
        data = file.read() 
        for incorrect_resname, corrected_resname in corrected_amino_acids.items(): 
            data = data.replace(incorrect_resname, corrected_resname)            
    with open(pdbfile, 'w') as file:
        file.write(data)
    
def mutate(temp_name, mutation): 


    """ 
    Function that deletes all atoms of a residue except for a methyl group. The 
    residue is then rename to alanine. 
    """ 

    
    if type(mutation) == int:
        mutation = [mutation] 
    else: 
        mutation = mutation 
    residue_type = [] 
    lines = [] 
    changed_lines = [] 
    deleted_lines = [] 
    with open(f'{temp_name}.pdb', 'r') as f: 
        lines = f.readlines() 
    for index, item in enumerate(lines): 
        line = item.split() 
        if line[0] == 'ATOM' and int(line[5]) in mutation:  # int(line[5]) is the resid column in the .temp_*.pdb
            if line[2] == 'N' or line[2] == 'HN' or line[2] == 'CA' or line[2] == 'HA' or line[2] == 'CB' or line[2] == 'O' or line[2]  == 'C': 
                changed_lines.append(index) 
            else: 
                deleted_lines.append(index) 
            residue_type.append(line[3]) 

    residue_type_list = list(set(residue_type)) 
    
    with open(f'{temp_name}.pdb', 'w') as fp:
        for number, line in enumerate(lines): 
            if number not in deleted_lines and number in changed_lines:
                for a in residue_type_list: 
                    line = line.replace(a, 'ALA') 
                fp.write(line) 
            elif number not in deleted_lines and number not in changed_lines:
                fp.write(line)

def find_mutation(mutations, corrected_universe):
    if type(mutations) == int:
        mutations = [mutations]
    selected_mutations = []
    for mutation in mutations:
        resname = corrected_universe.select_atoms(f'resid {mutation}').residues.resnames[0]
        resid = corrected_universe.select_atoms(f'resid {mutation}').residues.resids[0]
        selected_mutations.append(f'{resname}{resid}')
    selected_mutations = ', '.join(set(selected_mutations))
    mutation_arg = f'{mutations} \n{selected_mutations}'
    return mutation_arg
