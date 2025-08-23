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