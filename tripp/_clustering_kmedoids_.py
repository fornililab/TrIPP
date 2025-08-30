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
    Function to run K-Medoids clustering from sklearn_extra.
    Standard K-Medoids parameters can be found in the sklearn_extra documentation:
    https://scikit-learn-extra.readthedocs.io/en/stable/generated/sklearn_extra.cluster.KMedoids.html
    
    Parameters
    ----------
    n_clusters: int
        The number of output clusters.
    metric: str
        The metric to use for distance computation.
    method: str
        The algorithm used for clustering.
    init: str
        The medoid initialization method.
    max_iter: int
        The maximum number of iterations for the clustering algorithm.
    random_state: int
        Random seed used in the medoid initialisation if init = 'random'.
    clustering_matrix: np.ndarray
        The clustering matrix (feature matrix) created by the create_clustering_matrix function.
    frames: np.ndarray
        A list of frames corresponding to the points to be clustered.
    trajectory_names: np.ndarray
        An array of names of the trajectories to be clustered.
    Returns
    -------
    labels: np.ndarray
        An array of cluster labels for each point in the feature matrix.
    cluster_centers: list
        Indices of the cluster centers (medoids) within their individual trajectories (local indices). 
        Indices start from 0 for each trajectory.
    cluster_center_indices: list
        Indices of the cluster centers (medoids) within the full feature matrix (global indices).
    cluster_centers_trajectories: list
        Trajectory names of the cluster centers (medoids).
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
    medoid_indices = list(kmedoids_clustering.medoid_indices_)
    cluster_centers = list(np.ravel(frames[medoid_indices]))
    cluster_centers_trajectories = list(np.ravel(trajectory_names[medoid_indices]))

    return labels, cluster_centers, medoid_indices, cluster_centers_trajectories
