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


import pymol 
from pymol import cmd 
import MDAnalysis as mda 

def visualise_pka(structure, color_dictionary, lower_limit, upper_limit, color_palette): 
    


    """ 
    Function that can visualise the pka vlues of residues using PyMOL. The structure 
    and color_dictionary are used as input. The color_dicitonary contains ranges of 
    pKa values associated with lists of residues in those ranges and hues. The function 
    creates a PyMOL session file. 
    """ 
    

    u = mda.Universe(structure) 
    protein = u.select_atoms('protein') 
    nter = protein[0].resid 
    cter = protein[-1].resid 

    cmd.load(structure, 'protein_str')

    cmd.show("cartoon", 'protein_str')
    cmd.color("white", "protein_str") 

    nt = f'(resi {nter} and (name N+H1+H2+H3))' 
    ct = f'(resi {cter} and (name C+O+OT1+OT2+OXT))' 
    names = [] 

    for entry in color_dictionary: 
        if len(color_dictionary[entry]['residues']) > 0: 
            cmd.set_color(entry, color_dictionary[entry]['color']) 
            for residue in color_dictionary[entry]['residues']: 
                name = residue 
                if name == 'NTR': 
                    selection = nt 
                elif name == 'CTR': 
                    selection = ct 
                else: 
                    selection = f'resi {residue[3:]} and not ({ct} or {nt})' 
                    names.append(name) 
                
                cmd.select(name, selection) 
                cmd.show('licorice', name) 
                cmd.color(entry, name) 
    
    def sort_by_number(residue):
        index = int(residue[3:])
        return index

    sorted_residues = ['NTR']+sorted(names, key=sort_by_number)+['CTR']  

    sorted_res = ' '.join(sorted_residues) 

    cmd.order(sorted_res) 
    cmd.ramp_new('colorbar', 'none', [lower_limit, (lower_limit + upper_limit)/2, upper_limit], color_palette.split('_')) 
    cmd.save(structure.replace('pdb', 'pse')) 