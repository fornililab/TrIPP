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
    Function to run KMedoids clustering from sklearn_extra.
    
    Parameters
    ----------
    n_clusters: int
        The number of clusters to form.
    metric: str
        The metric to use for distance computation.
    method: str
        The method to use for clustering.
    init: str
        The initialization method for the medoids.
    max_iter: int
        The maximum number of iterations for the clustering algorithm.
    random_state: int
        Random seed for reproducibility.
    clustering_matrix: np.ndarray
        The input data for clustering.
    frames: np.ndarray
        The frames corresponding to the clustering data.
    trajectory_names: np.ndarray
        The names of the trajectories corresponding to the clustering data.
    Returns
    -------
    labels: np.ndarray
        An array of cluster labels for each point in the clustering matrix.
    cluster_centers: np.ndarray
        The coordinates of the cluster centers (medoids).
    medoid_indices: np.ndarray
        The indices of the medoids in the clustering matrix.
    cluster_centers_trajectories: np.ndarray
        A list of indices of the cluster centroids mapped onto the trajectory.
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
