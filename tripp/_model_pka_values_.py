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


""" 
Dictionary that contains the model pKa values of amino acids in solution used in 
PROPKA (https://github.com/jensengroup/propka/blob/27e0ad2f8d653404f57c78f1cae2ec32cb5adb68/propka/propka.cfg#L5-L13)
""" 


model_pka_values = {'ASP' : 3.8, 
                    'GLU' : 4.5, 
                    'CTR' : 3.2, 
                    'HIS' : 6.5, 
                    'NTR' : 8, 
                    'CYS' : 9, 
                    'TYR' : 10, 
                    'LYS' : 10.5, 
                    'ARG' : 12.5} 