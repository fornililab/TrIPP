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
def visualize_pka(tempfactors_structure, pymol_path, pse_output_filename, pka_value_summary, lower_limit, upper_limit, color_palette): 
    
    """ 
    Function that can visualize the pka values of residues using PyMOL.
    No parameters need to be adjusted as they are determined from the color_pka()
    in Visualization class.
    
    Parameters:
    tempfactors_structure: str
    The input of the structure pdb file which has tempfactors assigned.
    
    pymol_path: str
    The path to the PyMOL software needs to be specified. The script will 
    spawn a subprocess shell to run a python script in PyMOL. Preventing
    packaging issue.
    
    pse_output_filename: str
    The output name for pse session, automatically combined with the 
    pse_output_prefix and the coloring_method in color_pka().
    
    pka_value_summary: Pandas Series
    Pandas series containing one column of the name of residue with resid,
    and another column of the predicted pKa.
    
    lower limit: int or float
    Determines lower limit used to colour the reisdues in the PyMOL session. Any 
    value below the limit is coloured using the lowest end of the color gradient 
    used. 
    
    upper limit: int or float
    Determines upper limit used to colour the reisdues in the PyMOL session. Any 
    value above the limit is coloured using the highest end of the color gradient 
    used. 

    color_palette: str
    color palettes used to color the residues in the PyMOL session according to 
    the pKa value. The default is set to 'red_white_blue'. See PyMOL spectrum for
    allowed color palettes. Three colors palette is suggested.
    """
    with open('.pymol_template.py','a') as output:
        output.write(f"""cmd.load('{tempfactors_structure}', 'protein_str')
cmd.show("cartoon", 'protein_str')
cmd.color("white", "protein_str")\n""")
    Nterm_atoms = 'N+H1+H2+H3'
    Cterm_atoms = 'C+OC1+OC2+OT1+OT2'
    names = []
    for residue,predicted_pka in pka_value_summary.items():
        rounded_predicted_pka = round(predicted_pka,2)
        if 'N+' in residue:
            name = 'NTR'
            resid = residue[2:]
            selection = f'resi {resid} and name {Nterm_atoms}'
        elif 'C-' in residue:
            name = 'CTR'
            resid = residue[2:]
            selection = f'resi {resid} and name {Cterm_atoms}'
        else:
            name = residue
            names.append(name)
            resid = residue[3:]
            selection = f'resi {resid}'
        with open('.pymol_template.py','a') as output:
            output.write(f"""cmd.select('{name}', '{selection}') 
cmd.show('licorice', '{name}') 
cmd.spectrum('b','{color_palette}','{name}',{lower_limit},{upper_limit})
cmd.label('{name} and name CB','{rounded_predicted_pka}')\n""")
    
    sorted_residues = ' '.join(['NTR']+sorted(names, key=lambda x: int(x[3:]))+['CTR'])
    with open('.pymol_template.py','a') as output:
        output.write(f"""cmd.order('{sorted_residues}') 
cmd.ramp_new('colorbar', 'none', [{lower_limit}, ({lower_limit} + {upper_limit})/2, {upper_limit}], {color_palette.split('_')})
cmd.set('label_size','-2')
cmd.set('label_position','(1.2,1.2,1.2)') 
cmd.save('{pse_output_filename}')
cmd.quit()\n""")
    subprocess.run([f'{pymol_path} .pymol_template.py'],shell=True)
    os.remove('.pymol_template.py')