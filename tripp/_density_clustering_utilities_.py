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

def find_smallest_distance(cluster, frames, labels, trajectory_names, clustering_matrix):
    """Find the point in the cluster that is closest to the centroid of the cluster.
    Parameters
    ----------
    cluster : list
        A list of unique cluster labels.
    frames : list
        A list of frames corresponding to feature matrix points.
    labels : np.ndarray
        An array of cluster labels for each point in the feature matrix.
    trajectory_names : np.ndarray
        A 1D array containing the name of each frame in the feature matrix.
    clustering_matrix : np.ndarray
        The clustering matrix (feature matrix) created by the create_clustering_matrix function.
    Returns
    -------
    center : int
        The index of the centroid frame in its individual trajectory.
    center_index : int
        The index of the centroid frame in the full feature matrix (i.e., global index).
    center_trajectory : str
        The name of the trajectory to which the centroid frame belongs.
    """
    # The centroid is defined as the point which has the smallest
    # distance from the mean of all points in the cluster.
    cluster_frames = frames[labels == cluster]
    cluster_trajectories = trajectory_names[labels == cluster]
    cluster_members = clustering_matrix[labels == cluster]
    cluster_mean = np.mean(cluster_members, axis=0)
    dist_mean = []
    for item in cluster_members:
        dist_mean.append(np.sqrt(np.sum(np.square(item - cluster_mean))))
    smallest_index = np.argmin(np.array(dist_mean))
    center = cluster_frames[smallest_index][0]
    center_trajectory = cluster_trajectories[smallest_index]
    center_index_frame_mask = np.ravel(frames == center)
    center_index_trajectory_mask = trajectory_names == center_trajectory
    center_index_mask = center_index_frame_mask & center_index_trajectory_mask
    center_index = np.where(center_index_mask)[0][0]
    return center, center_index, center_trajectory

def find_density_clustering_centroids(labels, frames, trajectory_names, clustering_matrix):
    """Find the centroids of the density clusters.
    Parameters
    ----------
    labels : np.ndarray
        An array of cluster labels for each point in the feature matrix.
    frames : list
        A list of frames corresponding to feature matrix points.
    trajectory_names : np.ndarray
        A 1D array containing the name of each frame in the feature matrix.
    clustering_matrix : np.ndarray
        The clustering matrix (feature matrix) created by the create_clustering_matrix function.
    Returns
    -------
    cluster_centers: list
        The indices of the cluster centers within their respective individual trajectories.
    cluster_center_indices: list
        The indices of the cluster centers in the full feature matrix (i.e., global indices).
    cluster_centers_trajectories: list
         A list of trajectory names of the cluster centers.
    """
    
    # A list of clusters is made to iterate over.
    clusters = list(set(labels))

    if -1 in clusters:
        # The cluster label -1 is removed since it is the label
        # for outliers.
        clusters.remove(-1)

    cluster_centers = []
    cluster_center_indices = []
    cluster_centers_trajectories = []

    for cluster in clusters:
        center, center_index, center_trajectory = find_smallest_distance(cluster, frames, labels, trajectory_names, clustering_matrix)
        cluster_centers.append(center)
        cluster_center_indices.append(center_index)
        cluster_centers_trajectories.append(center_trajectory)

    return cluster_centers, cluster_center_indices, cluster_centers_trajectories