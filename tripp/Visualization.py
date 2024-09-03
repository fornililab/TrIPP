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


import MDAnalysis as mda
import numpy as np
import pandas as pd 
from tripp._visualize_pka_ import visualize_pka 
from tripp._model_pka_values_ import model_pka_values 
from tripp._create_mda_universe_ import create_mda_universe 

class Visualization: 


    """ 
    This class allows for the visualization of pKa values using PyMOL. The 
    class is called using a strucutre of the protein and a CSV file generated 
    by the Trajectory class containing all pKa values.
    
    """


    def __init__(self, structure, pka_file):

        self.structure = structure
        self.pka_file = pka_file 
        self.u = create_mda_universe(topology_file=self.structure, trajectory_file=None) 
        
        if type(self.pka_file) == list: 
            pka_file_l = [] 
            for file in pka_file: 
                pka_file_l.append(pd.read_csv(file)) 
            self.pka_values = pd.concat(pka_file_l) 
        
        else: 
            self.pka_values = pd.read_csv(self.pka_file) 
        

    def color_pka(self, pymol_path, output, coloring_method='mean', lower_limit=0, upper_limit=14, color_palette='red_white_blue'): 
        """
        The method color_pka can be called which generate a PyMOL session file.
        
        Parameters:
        
        pymol_path: str
        Path to the PyMOL software needs to be specified. The script will spawn
        a subprocess shell to run a python script in PyMOL. Preventing packaging 
        issue.
        
        output: str
        The output prefix for the PyMOL .pse file. The prefix will be combined 
        with the coloring_method ('mean' or 'difference_to_model_value') to give
        the pse_output_filename.
        
        coloring_method: str, int, float default 'mean'
        To determine how the color of each residue is produced. Can be 'mean', 
        where the mean pKa value accross all frames is used or 
        'difference_to_model_value' where the mean pKa value is calculated 
        and the difference to the model value of the amino acid in solution is 
        used. A specific timestep can also be selected for visualisation by setting the 
        coloring_method to the timestep in question. 
        
        lower limit: int or float, default 0
        Determines lower limit used to colour the reisdues in the PyMOL session. Any 
        value below the limit is coloured using the lowest end of the color gradient 
        used. 
        
        upper limit: int or float, default 14
        Determines upper limit used to colour the reisdues in the PyMOL session. Any 
        value above the limit is coloured using the highest end of the color gradient 
        used. 

        color_palette: str, default 'red_white_blue'
        color palettes used to color the residues in the PyMOL session according to 
        the pKa value. The default is set to 'red_white_blue'. See PyMOL spectrum for
        allowed color palettes. Three colors palette is suggested.
        """

        #calculation of values depending on colouring method 
        if coloring_method == 'mean': 
            del self.pka_values['Time [ps]'] 
            pka_values_summary = self.pka_values.mean(axis=0) 
            tempfactors_output_structure= f"{output}_mean.pdb"

        elif coloring_method == 'difference_to_model_value': 
            del self.pka_values['Time [ps]'] 
            pka_values_mean = self.pka_values.mean(axis=0) 
            for residue, value in pka_values_mean.items(): 
                if 'N+' in residue: 
                    pka_values_mean[residue] = pka_values_mean[residue]-model_pka_values['NTR'] 
                elif 'C-' in residue: 
                    pka_values_mean[residue] = pka_values_mean[residue]-model_pka_values['CTR'] 
                else: 
                    pka_values_mean[residue] = pka_values_mean[residue]-model_pka_values[residue[0:3]] 
            pka_values_summary = pka_values_mean
            tempfactors_output_structure = f"{output}_difference_to_model_value.pdb" 

        elif type(coloring_method) == int or type(coloring_method) == float: 

            if type(self.pka_file) == list and len(self.pka_file) != 1: 
                print('Please select only one CSV file when selecting a specific time step.') 
                
            else: 
                self.pka_values = self.pka_values[self.pka_values['Time [ps]'] == coloring_method] 
                del self.pka_values['Time [ps]'] 
                pka_values_summary = self.pka_values.mean(axis=0) 
                tempfactors_output_structure = f"{output}_time{coloring_method}.pdb" 
        
        # GROMACS atom naming scheme, other naming scheme will not be valid, user may 
        # need to add it by themselves for Nterm and Cterm.
        Nterm_atoms = 'N H1 H2 H3'
        Cterm_atoms = 'C O OC1 OC2 OT1 OT2 OXT'
        # Looping the pka_values_summary which contains one column of the name for 
        # residue and resid, and the other column of predicted pKa. The tempfactor 
        # (previously known as bfactor) of individual residue is assigned according 
        # to the predicted pKa from pka_values_summary. The structure with the 
        # tempfactor is written as pdb and a PyMOL session is generated as pse.
        for residue, predicted_pka in pka_values_summary.items():
            if 'N+' in residue or 'C-' in residue:
                resid = int(residue[2:])
            else:
                resid = int(residue[3:])
            rounded_predicted_pka = round(predicted_pka,2)
            if 'N+' in residue:
                ag = self.u.select_atoms(f'resid {resid} and name {Nterm_atoms}')
            elif 'C-' in residue:
                ag = self.u.select_atoms(f'resid {resid} and name {Cterm_atoms}')
            elif resid == self.u.residues.resids[0]:
                ag = self.u.select_atoms(f'resid {resid} and not name {Nterm_atoms}')
            elif resid == self.u.residues.resids[-1]:
                ag = self.u.select_atoms(f'resid {resid} and not name {Cterm_atoms}')
            else:
                ag = self.u.select_atoms(f'resid {resid}')
            ag.tempfactors = np.full(ag.tempfactors.shape,rounded_predicted_pka)
        ag = self.u.select_atoms('all')
        ag.write(tempfactors_output_structure)
        pse_output_filename = tempfactors_output_structure.replace('.pdb', '.pse')
        print(tempfactors_output_structure)
        visualize_pka(tempfactors_output_structure, pymol_path, pse_output_filename, pka_values_summary, lower_limit, upper_limit, color_palette) 

        