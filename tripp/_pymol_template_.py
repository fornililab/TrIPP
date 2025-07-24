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
    Function that can visualize the values of residues using PyMOL.
    No parameters need to be adjusted as they are determined from the gen_pse()
    in Visualization class.

    Parameters:
    tempfactors_topology_file: str
    The input of the topology_file pdb file which has tempfactors assigned.
    
    pymol_path: str
    The path to the PyMOL software needs to be specified. The script will
    spawn a subprocess shell to run a python script in PyMOL. Preventing
    packaging issue.

    pse_output_filename: str
    The output name for pse session, automatically combined with the
    pse_output_prefix and the coloring_method in gen_pse().

    values_df: Pandas dataframe
    Pandas dartaframe containing one column of the name of residue with resid,
    and another column of the value (pKa or correlation).

    lower limit: int or float
    Determines lower limit used to colour the reisdues in the PyMOL session.
    Any value below the limit is coloured using the lowest end of the color
    gradient used.

    upper limit: int or float
    Determines upper limit used to colour the reisdues in the PyMOL session.
    Any value above the limit is coloured using the highest end of the color
    gradient used.

    color_palette: str
    color palettes used to color the residues in the PyMOL session according
    to the pKa value. The default is set to 'red_white_blue'. See PyMOL spectrum
    for allowed color palettes. Three colors palette is suggested.
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