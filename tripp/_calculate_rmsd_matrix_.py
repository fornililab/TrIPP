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
from tqdm import tqdm

def calculate_rmsd_matrix(clustering_matrix, frames):
    """
    Function that calculates the pairwise Euclidean distance matrix for the greedy clustering algorithm.
    
    Parameters
    ----------
    clustering_matrix : np.ndarray
        The clustering matrix (feature matrix) created by the create_clustering_matrix function.
    frames : list
        A list of frames corresponding to the points to be clustered.
    Returns
    -------
    rmsd_array : np.ndarray
        A 2D numpy array containing the Euclidean distance values between each pair of points.
    """

    def calculate_rmsd(point1, point2):

        sqr_sum = np.sum(np.square(point1 - point2))
        rmsd = np.sqrt(sqr_sum)

        return rmsd

    rmsd_array = np.zeros((len(frames), len(frames)))
    print('Building distance matrix...')
    for i in tqdm(range(len(frames))):
        point1 = clustering_matrix[i]
        for j in range(len(frames)):
            point2 = clustering_matrix[j]
            rmsd_array[i, j] = calculate_rmsd(point1, point2)

    return rmsd_array
