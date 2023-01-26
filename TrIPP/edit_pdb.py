def edit_pdb(pdbfile): 
    with open(pdbfile, 'r') as file: 
        data = file.read() 
        corrected_amino_acids = { 
            'HSD' : 'HIS', 
            'HSE' : 'HIS', 
            'HSP' : 'HIS', 
            'LSN' : 'LYS', 
            'GLUP' : 'GLU ', 
            'ASPP' : 'ASP ' 
            }
        for correction in corrected_amino_acids: 
            data = data.replace(correction, corrected_amino_acids[correction])
    with open(pdbfile, 'w') as file: 
        file.write(data) 

def mutate(temp_name, mutation): 
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
        if line[0] == 'ATOM' and int(line[5]) in mutation: 
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