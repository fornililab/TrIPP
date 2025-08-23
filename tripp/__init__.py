"""
    @release_date  : $release_date
    @version       : $release_version
    @author        : Christos Matsingos
    
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

from .Trajectory import Trajectory #Addition of code so that it is not necessary to write from tripp.Trajectory import Trajectory but rather from tripp import Trajectory 
from .Visualization import Visualization 
from .Clustering import Clustering 