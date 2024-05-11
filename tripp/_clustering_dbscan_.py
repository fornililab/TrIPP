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

from sklearn.cluster import DBSCAN 
import numpy as np 

def dbscan_clustering_(eps, min_samples, metric, metric_params, algorithm, leaf_size, p, n_jobs, clustering_matrix, frames, find_centroid): 

    """
    Function to run DBSCAN clustering. 
    """

    dbscan_clustering = DBSCAN(eps=eps, min_samples=min_samples, metric=metric, metric_params=metric_params, algorithm=algorithm, leaf_size=leaf_size, p=p, n_jobs=n_jobs) 
    dbscan_clustering.fit(clustering_matrix) 
    labels = dbscan_clustering.labels_ 

    def find_dbscan_centroid(): 
        
        clusters = list(set(labels)) #A list of clusters is made to iterate over. 
        
        if -1 in clusters: 
            clusters.remove(-1) #The cluster label -1 is removed since it is the label for outliers. 
        
        cluster_centers = [] 
        cluster_center_indices = [] 

        def find_smallest_distance(cluster): 
            #The centroid is defined as the point which has the smallest distance from the mean of all points in the cluster. 
            cluster_mean = np.mean(clustering_matrix, axis=0) 
            cluster_frames = frames[labels==cluster] 
            cluster_members = clustering_matrix[labels==cluster] 
            dist_mean = [] 
            for item in cluster_members: 
                dist_mean.append(np.sqrt(np.sum(np.square(item-cluster_mean)))) 
            smallest_index = np.argmin(dist_mean) 
            center = cluster_frames[smallest_index] 
            center_index = np.where(frames==center)[0][0]
            return center, center_index  
        
        for cluster in clusters: 
            center, center_index = find_smallest_distance(cluster) 
            cluster_centers.append(center) 
            cluster_center_indices.append(center_index) 
        
        return cluster_centers, cluster_center_indices 
    
    if find_centroid == True: 
        cluster_centroids, cluster_centroid_indices = find_dbscan_centroid() 
        return labels, cluster_centroids, cluster_centroid_indices 
    
    elif find_centroid == False: 
        return labels 
