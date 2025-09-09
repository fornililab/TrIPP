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

import numpy as np 

def determine_cluster_population_ratio(labels, max_cluster_population, clustering_method): 
    """
    Calculates the relative population (fraction of frames) for each cluster and 
    compares the maximum value to a given threshold
    Parameters
    ----------
    labels: np.ndarray
        Array where each element is the cluster ID assigned to a data point
    max_cluster_population: float
        The maximum relative population.
    clustering_method: str
        The clustering method used, e.g., 'DBSCAN'.
    Returns
    -------
    bool:
        True if the maximum relative population of any cluster is greater than 
        or equal to `max_cluster_population` (default 95%) , False otherwise.
    """
    unique_labels = list(set(labels)) 
    cluster_population = np.zeros(len(unique_labels)) 

    labels = np.ravel(labels) 

    if clustering_method == 'DBSCAN' or clustering_method == 'HDBSCAN': 
        labels = labels[labels!=-1] 

    for label in unique_labels: 
        cluster_labels = labels[labels==label] 
        cluster_elements = len(cluster_labels) 
        cluster_population[label] = cluster_elements/len(labels) 
    
    if np.max(cluster_population) >= max_cluster_population: 
        return True 
    
    else: 
        return False 