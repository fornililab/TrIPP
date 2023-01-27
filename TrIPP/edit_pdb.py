def edit_pdb(pdbfile): 
    with open(pdbfile, 'r') as file: 
        data = file.read() 
        corrected_amino_acids = { 
            'ALAD' : 'ALA', 
            'ARGN' : 'ARG ', 
            'ASN1' : 'ASN ', 
            'ASPH' : 'ASP ', 
            'ASPP' : 'ASP ',  
            'CALA' : 'ALA ', 
            'CARG' : 'ARG ', 
            'CASF' : 'ASP ', 
            'CASN' : 'ASN ', 
            'CASP' : 'ASP ', 
            'CCYS' : 'CYS ', 
            'CCYX' : 'CYS ', 
            'CGLN' : 'GLN ', 
            'CGLU' : 'GLU ', 
            'CGLY' : 'GLY ', 
            'CHID' : 'HIS ', 
            'CHIE' : 'HIS ', 
            'CHIP' : 'HIS ', 
            'CILE' : 'ILE ', 
            'CLEU' : 'LEU ', 
            'CLYS' : 'LYS ', 
            'CMET' : 'MET ', 
            'CPHE' : 'PHE ', 
            'CPRO' : 'PRO ', 
            'CSER' : 'SER ', 
            'CTHR' : 'THR ', 
            'CTRP' : 'TRP ', 
            'CTYR' : 'TYR ', 
            'CVAL' : 'VAL ', 
            'CYS1' : 'CYS ', 
            'CYS2' : 'CYS ', 
            'CYSH' : 'CYS ', 
            'GLUH' : 'GLU ', 
            'GLUP' : 'GLU ', 
            'HIS1' : 'HIS ', 
            'HIS2' : 'HIS ', 
            'HISA' : 'HIS ', 
            'HISB' : 'HIS ', 
            'HISD' : 'HIS ', 
            'HISE' : 'HIS ', 
            'HISH' : 'HIS ', 
            'LYSH' : 'LYS ', 
            'NALA' : 'ALA ', 
            'NARG' : 'ARG ', 
            'NASN' : 'ASN ', 
            'NASP' : 'ASP ', 
            'NCYS' : 'CYS ', 
            'NCYX' : 'CYS ', 
            'NGLN' : 'GLY ', 
            'NGLU' : 'GLU ', 
            'NGLY' : 'GLY ', 
            'NHID' : 'HIS ', 
            'NHIE' : 'HIS ', 
            'NHIP' : 'HIS ', 
            'NILE' : 'ILE ', 
            'NLEU' : 'LEU ', 
            'NLYS' : 'LYS ', 
            'NMET' : 'MET ', 
            'NPHE' : 'PHE ', 
            'NPRO' : 'PRO ', 
            'NSER' : 'SER ', 
            'NTHR' : 'THR ', 
            'NTRP' : 'TRP ', 
            'NTYR' : 'TYR ', 
            'NVAL' : 'VAL ', 
            'PGLU' : 'GLU ',
            'ASF' : 'ASP', 
            'ASH' : 'ASP', 
            'CYM' : 'CYS', 
            'CYN' : 'CYS', 
            'CYX' : 'CYS', 
            'GLH' : 'GLU', 
            'HID' : 'HIS', 
            'HIE' : 'HIS', 
            'HIP' : 'HIS', 
            'HSD' : 'HIS', 
            'HSE' : 'HIS', 
            'HSP' : 'HIS', 
            'LYN' : 'LYS', 
            'LSN' : 'LYS', 
            'MSE' : 'MET'
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