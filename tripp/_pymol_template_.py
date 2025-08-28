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


import subprocess
import os
def gen_pymol_template(tempfactors_topology_file,
                       pymol_path, 
                       pse_output_filename, 
                       values_df, 
                       lower_limit, 
                       upper_limit, 
                       color_palette): 

    """
    Function to colour-map residue properties on a reference structure using PyMOL.
    The values of the parameters are determined by gen_pse()
    from Visualization class.

    Parameters:
    tempfactors_topology_file: str
    PDB file with the input topology and tempfactors assigned.
    
    pymol_path: str
    Path to PyMOL. The script will spawn a subprocess shell to run a 
    python script in PyMOL.

    pse_output_filename: str
    Name of the PyMOL session file (pse) produced as output. The name is 
    a combination of pse_output_prefix and the colouring method.

    values_df: Pandas DataFrame
    DataFrame with two columns, one for the residue identifier and one for
    the value to be mapped (pKa, pKa difference or correlation)

    lower_limit: int or float
    Lower bound of the color scale in the PyMOL session. Residues with values
    below this threshold are assigned the gradient’s minimum color.

    upper_limit: int or float
    Upper bound of the color scale in the PyMOL session. Residues with values
    above this threshold are assigned the gradient’s maximum color.

    color_palette: str
    Colour palette. The default is set to 'red_white_blue'. See PyMOL spectrum
    for allowed colour palettes (a three-colour palette is recommended).
    """
    with open('.pymol_template.py','a') as output:
        output.write(f"""cmd.load('{tempfactors_topology_file}', 'protein_str')
cmd.show("cartoon", 'protein_str')
cmd.color("white", "protein_str")\n""")
    names = []
    residue_identifiers, values = (columns for _,columns in values_df.items()) 
    for residue_identifier, value in zip(residue_identifiers,values):
        residue = residue_identifier.split(':')[0]
        chain = residue_identifier.split(':')[-1]
        rounded_value = round(value,2)
        if 'N+' in residue:
            resid = residue[2:]
            name = f'NTR{resid}_{chain}'
            selection = f'(bb. and not elem O and not elem C and byres protein_str and chain {chain} and resi {resid}) extend 1 and not elem C'
            label_sel = f'{name} and bb. and elem N'
        elif 'C-' in residue:
            resid = residue[2:]
            name = f'CTR{resid}_{chain}'
            selection = f'((bb. and byres protein_str and chain {chain} and resi {resid}) and elem C and not name CA) extend 1 and not name CA'
            label_sel = f'{name} and bb. and elem C and not name CA'
        else:
            name = f'{residue}_{chain}'
            resid = residue[3:]
            selection = f'((byres protein_str and chain {chain} and resi {resid})&(sc.|(n. CA|n. N&r. PRO))) and not name H1+H2+H3 and not (protein_str and chain {chain} and resi {resid} & name C extend 1 &! name CA)'
            label_sel = f'{name} and name CB'
        names.append(name)
        with open('.pymol_template.py', 'a') as output:
            output.write(f"""cmd.create('{name}', '{selection}') 
cmd.show('licorice', '{name}') 
cmd.spectrum('b','{color_palette}','{name}',{lower_limit},{upper_limit})
cmd.label('{label_sel}','{rounded_value}')\n""")
    
    sorted_residues = ' '.join(sorted(names, key=lambda x: (x[-1], int(x[3:-2]))))
    with open('.pymol_template.py', 'a') as output:
        output.write(f"""cmd.order('{sorted_residues}')
cmd.ramp_new('colorbar', 'none', [{lower_limit}, ({lower_limit} + {upper_limit})/2, {upper_limit}], {color_palette.split('_')})
cmd.set('label_size','-2')
cmd.set('label_position','(1.2,1.2,1.2)')
cmd.orient('protein_str')
cmd.bg_color('white')
cmd.set('orthoscopic')
cmd.set('depth_cue',0)
cmd.save('{pse_output_filename}')
cmd.quit()\n""")
    subprocess.run([f'{pymol_path} -c .pymol_template.py'],shell=True)
    os.remove('.pymol_template.py')