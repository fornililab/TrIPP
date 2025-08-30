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

def sort_clusters(labels, cluster_centers, cluster_center_indices, cluster_centers_trajectories): 
    """
    Sorts clusters based on their population and reassigns labels.
    Parameters
    ----------
    labels : np.ndarray
        An array of cluster labels for each point in the clustering matrix                 (feature matrix).
    cluster_centers: list
        Indices of the cluster centers within their individual trajectories (local indices). 
        Indices start from 0 for each trajectory.
    cluster_center_indices: list
        Indices of the cluster centers within the full feature matrix (global indices).
    cluster_centers_trajectories: list
        Trajectory names of the cluster centers.
    """
    clusters = np.unique(labels) 
    clusters = clusters[clusters != -1] 
    labels_without_outliers = labels[labels != -1]
    
    frequency = np.bincount(labels_without_outliers, minlength=len(clusters))  
    
    new_labels = np.zeros_like(labels, dtype=int) 
    new_labels[labels == -1] = -1 
    
    for i in range(len(clusters)): 
        index_biggest = np.argmax(frequency) 
        biggest_cluster = clusters[index_biggest] 
        labels_mask = labels == biggest_cluster 
        new_labels[labels_mask] = i 
        frequency[index_biggest] = 0 
        
    
    new_cluster_centers = np.zeros_like(clusters, dtype=int) 
    new_cluster_center_indices = np.zeros_like(clusters, dtype=int) 
    new_cluster_centers_trajectories = np.zeros_like(clusters, dtype=object) 

    for i in range(len(clusters)): 
        cluster_index = cluster_center_indices[i] 
        new_label = new_labels[cluster_index] 
        new_cluster_center_indices[new_label] = cluster_center_indices[i] 
        new_cluster_centers[new_label] = cluster_centers[i] 
        new_cluster_centers_trajectories[new_label] = cluster_centers_trajectories[i] 
    
    return new_labels, new_cluster_centers, new_cluster_center_indices, new_cluster_centers_trajectories 