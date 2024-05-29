import numpy as np 

def determine_cluster_population_ratio(labels, max_cluster_population, clustering_method): 

    unique_labels = list(set(labels)) 
    cluster_population = np.zeros(len(unique_labels)) 

    labels = np.ravel(labels) 

    if clustering_method == 'DBSCAN': 
        labels = labels[labels!=-1] 

    for label in unique_labels: 
        cluster_labels = labels[labels==label] 
        cluster_elements = len(cluster_labels) 
        cluster_population[label] = cluster_elements/len(labels) 
    
    if np.max(cluster_population) >= max_cluster_population: 
        return True 
    
    else: 
        return False 