"""
    This file is part of the TrIPP software
    (https://github.com/fornililab/TrIPP).
    Copyright (c) Christos Matsingos, Ka Fu Man and Arianna Fornili.

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
and pKAI (https://github.com/bayer-science-for-a-better-life/pKAI/blob/11cd7973c936dd6fa1b38654e7d0ae30bec55bfc/pKAI/protein.py#L4-L13)
""" 

model_propka_values = {'ASP' : 3.8, 
                        'GLU' : 4.5, 
                        'CTR' : 3.2, 
                        'HIS' : 6.5, 
                        'NTR' : 8, 
                        'CYS' : 9, 
                        'TYR' : 10, 
                        'LYS' : 10.5, 
                        'ARG' : 12.5} 

model_pkai_values= {"ASP": 3.79,
                    "CTR": 2.90,
                    "CYS": 8.67,
                    "GLU": 4.20,
                    "HIS": 6.74,
                    "LYS": 10.46,
                    "NTR": 7.99,
                    "TYR": 9.59}