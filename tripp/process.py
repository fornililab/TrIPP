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
from tripp._model_pka_values_ import model_pka_values 
import numpy as np 

def calculate_difference_to_model(output_file): 

    df_pka = pd.read_csv(f'{output_file}_pka.csv') 

    df_dif = pd.DataFrame() 
    df_dif['Time [ps]'] = df_pka['Time [ps]'] 
    del df_pka['Time [ps]'] 
    residues = np.array(df_pka.columns) 
    for residue in residues: 
        if 'N+' in residue: 
            residue_modified = residue.replace('N+', 'NTR') 
        elif 'C-' in residue: 
            residue_modified = residue.replace('C-', 'CTR') 
        else: 
            residue_modified = residue 
        
        amino_acid = residue_modified[:3] 
        model_value = model_pka_values[amino_acid] 
        residue_data = np.around(df_pka[residue].to_numpy() - model_value, decimals=3) 
        df_dif[residue] = residue_data 
    
    df_dif.to_csv(f'{output_file}_difference_to_model.csv', index=False) 