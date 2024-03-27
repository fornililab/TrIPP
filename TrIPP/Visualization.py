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


import pandas as pd 
from tripp.color_palettes import color_palettes 
from tripp.visualise_pka import visualise_pka 
from tripp.model_pka_values import model_pka_values 

class Visualization: 


    """ 
    This class allows for the visualisation of pKa values using PyMOL. The 
    class is called using a strucutre of the protein and a CSV file generated 
    by the Trajectory class containing all pKa values. The method color_pka 
    can be called which generated a PyMOL session file. The arguments of 
    color_pka are the following: 
    
    coloring_method: 'mean' or 'difference_to_model_value' 
    used to determine how the color of each residue is produced. Can be 'mean', 
    where the mean pKa value accross all frames is used or 
    'difference_to_model_value' where the mean pKa value is calculated 
    and the difference to the model value of the amino acid in solution is 
    used. 
    
    lower limit: int or float 
    determines lower limit used to colour the reisdues in the PyMOL session. Any 
    value below the limit is coloured using the lowest end of the color gradient 
    used. 
    
    upper limit: int or float 
    determines upper limit used to colour the reisdues in the PyMOL session. Any 
    value above the limit is coloured using the highest end of the color gradient 
    used. 

    color_palette: 'red_white_blue', 'blue_white_red', 'red_white_green', and 'green_white_red' 
    color palettes used to color the residues in the PyMOL session according to 
    pKa value. The defaults is set to 'red_white_blue'. 
    """


    def __init__(self, structure, pka_file):

        self.structure = structure
        self.pka_file = pka_file
        

    def color_pka(self, coloring_method, lower_limit, upper_limit, color_palette='red_white_blue'): 
        
        #load pKa values and remove time column 
        pka_values = pd.read_csv(self.pka_file) 
        del pka_values['Time [ps]'] 
        
        #calculation of values depending on colouring method 
        if coloring_method == 'mean': 
            pka_values_summary = pka_values.mean(axis=0) 

        elif coloring_method == 'difference_to_model_value': 
            pka_values_mean = pka_values.mean(axis=0) 
            for residue, value in pka_values_mean.items(): 
                if 'N+' in residue: 
                    pka_values_mean[residue] = pka_values_mean[residue]-model_pka_values['NTR'] 
                elif 'C-' in residue: 
                    pka_values_mean[residue] = pka_values_mean[residue]-model_pka_values['CTR'] 
                else: 
                    pka_values_mean[residue] = pka_values_mean[residue]-model_pka_values[residue[0:3]] 
            pka_values_summary = pka_values_mean 
        
        #create a color_dictionary where the residues and hues of colour are associated to pKa intervals 
        color_dictionary = {} 
        color_intervals = [] 
        color_step = lower_limit 
        color_palette_pka = color_palettes[color_palette] 
        for i in range(len(color_palette_pka)): 
            color_intervals.append([round(color_step, 2), round(color_step+(upper_limit-lower_limit)/len(color_palette_pka), 2)]) 
            color_step+=round((upper_limit-lower_limit)/len(color_palette_pka), 2) 
        for index, item in enumerate(color_intervals): 
            color_dictionary['pka_group_'+str(index+1)] = {'pka interval' : color_intervals[index], 'color' : color_palette_pka[index], 'residues' : []}  
        
        for residue, value in pka_values_summary.items(): 
            if 'N+' in residue: 
                resid = 'NTR' 
            elif 'C-' in residue: 
                resid = 'CTR' 
            else: 
                resid = residue
            if value < color_dictionary['pka_group_1']['pka interval'][0]: 
                color_dictionary['pka_group_1']['residues'].append(resid) 
            elif value > color_dictionary['pka_group_51']['pka interval'][1]: 
                color_dictionary['pka_group_51']['residues'].append(resid) 
            else: 
                for entry in color_dictionary: 
                    if value >= color_dictionary[entry]['pka interval'][0] and value < color_dictionary[entry]['pka interval'][1]: 
                        color_dictionary[entry]['residues'].append(resid) 
        
        visualise_pka(self.structure, color_dictionary, lower_limit, upper_limit, color_palette) 

 