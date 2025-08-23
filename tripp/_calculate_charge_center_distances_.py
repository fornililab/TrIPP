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

import numpy as np


def calculate_charge_center_distances(positions):
    """
    Function that calculates the distances between the charge centers
    of n residues. The function takes as input a n x 3 matrix
    containing the coordinates of the charge centers
    and calculates n(n-1)/2 distances.
    
    Parameters
    ----------
    positions : np.ndarray
        An n x 3 numpy array containing the coordinates of the charge centers.
    Returns
    -------
    distance_array : np.ndarray
        A 1D numpy array containing the distances between each pair of charge centers
    """

    distances = []
    for i in range(len(positions)):
        center_i = positions[i]
        for j in range(i + 1, len(positions)):
            center_j = positions[j]
            sqr_sum = np.sum(np.square(center_i - center_j))
            d = np.sqrt(sqr_sum)
            distances.append(d)

    distance_array = np.array(distances)

    return distance_array
