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

from sklearn_extra.cluster import KMedoids
import numpy as np


def kmedoids_clustering(
    n_clusters,
    metric,
    method,
    init,
    max_iter,
    random_state,
    clustering_matrix,
    frames,
    trajectory_names,
):
    """
    Function to run KMedoids clustering.
    """

    kmedoids_clustering = KMedoids(
        n_clusters=n_clusters,
        metric=metric,
        method=method,
        init=init,
        max_iter=max_iter,
        random_state=random_state,
    ).fit(clustering_matrix)

    labels = kmedoids_clustering.labels_
    medoid_indices = kmedoids_clustering.medoid_indices_
    cluster_centers = np.ravel(frames[medoid_indices])
    cluster_centers_trajectories = np.ravel(trajectory_names[medoid_indices])

    return labels, cluster_centers, medoid_indices, cluster_centers_trajectories
