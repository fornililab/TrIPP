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
from tripp._generate_trajectory_log_ import log_header

def generate_clustering_summary(trajectory_file, 
                                topology_file, 
                                pka_file, 
                                selections, 
                                include_distances, 
                                include_buriedness, 
                                clustering_method, 
                                automatic, 
                                silhouette_scores, 
                                n_components, 
                                cummulative_variance, 
                                buriedness_file):
    """
    Generates a summary of the clustering process.
    Parameters
    ----------
    trajectory_file: str or dict
        The path to the trajectory file or a dictionary with trajectory names as keys and file paths as values.
    topology_file: str
        The path to the topology file.
    pka_file: str or list
        The path to the pKa file or a list of paths if multiple pKa files are used.
    selections: list
        A list of residue selections used for clustering.
    include_distances: bool
        Whether distances between charge centers were included in the clustering.
    include_buriedness: bool
        Whether buriedness (buried ratio) values were included in the clustering.
    clustering_method: str
        The clustering method used, e.g., 'DBSCAN'.
    automatic: bool
        Whether automatic clustering (automatic scanning of clustering parameters) was performed.
    silhouette_scores: pd.DataFrame
        A DataFrame containing silhouette scores for different clustering parameters.
    n_components: int or None
        The number of principal components used for dimensionality reduction, or None if no reduction was done.
    cummulative_variance: float
        The cumulative variance explained by the principal components used.
    buriedness_file: str or list
        The path to the buriedness file or a list of paths if multiple buriedness files are used.
    Returns
    -------
    summary: str
        A formatted string summarizing the clustering process, including file information, residue selections, 
        whether distances and buriedness were included, clustering method, dimensionality reduction details, 
        automatic clustering parameters, and the best silhouette score parameters.
    """
    #information on files 
    if type(trajectory_file) == str: 
        trajectory_name = 'Unnnamed Trajectory' 
        if include_buriedness == False: 
            trajectory_file_summary= f'{trajectory_name} \nTrajectory file: {trajectory_file} \npKa file: {pka_file} \n\n' 
        elif include_buriedness == True: 
            trajectory_file_summary= f'{trajectory_name} \nTrajectory file: {trajectory_file} \npKa file: {pka_file}\nBuriedness file: {buriedness_file} \n\n' 

    else: 
        if include_buriedness == False: 
            trajectory_file_summary = '' 
            for trajectory_index, trajectory_name in enumerate(trajectory_file.keys()): 
                trajectory_file_summary+=f'{trajectory_name} \nTrajectory file: {trajectory_file[trajectory_name]} \npKa file: {pka_file[trajectory_index]} \n\n' 
        elif include_buriedness == True: 
            trajectory_file_summary = '' 
            for trajectory_index, trajectory_name in enumerate(trajectory_file.keys()): 
                trajectory_file_summary+=f'{trajectory_name} \nTrajectory file: {trajectory_file[trajectory_name]} \npKa file: {pka_file[trajectory_index]}\nBuriedness file: {buriedness_file[trajectory_index]} \n\n' 

    file_summary = f"""
Topology file: {topology_file} 

Trajectories: 

{trajectory_file_summary}
""" 
    
    #information on what residues were used for the clustering 
    residue_summary = f'Clustering was performed using {", ".join(selections[:-1])} and {selections[-1]}.' 

    #information on whether distances between charge centers were used for the clustering 
    if include_distances == True: 
        n = len(selections) 
        num_distances = int((n*(n-1))/2)
        if num_distances == 1: 
            include_distances_summary = f'In total {num_distances} distance between charge centers was included in the clustering.' 
        else: 
            include_distances_summary = f'In total {num_distances} distances between charge centers were included in the clustering.' 
    
    elif include_distances == False: 
        include_distances_summary = f'No distances between charge centers were included in the clustering.' 
    
    #information on whether buriedness measures were used for the clustering 
    if include_buriedness== True: 
        n = len(selections) 
        if n == 1: 
            include_buriedness_summary = f'In total {n} buriedness measure was included in the clustering.' 
        else: 
            include_buriedness_summary = f'In total {n} buriedness measures were included in the clustering.' 
    
    elif include_buriedness == False: 
        include_buriedness_summary = f'No buriedness measures were included in the clustering.' 
    
    #information on clustering method 
    clustering_method_summary = f'The {clustering_method} method was used for the clustering.' 

    #information on dimensionality reduction 
    if n_components == None: 
        dimensionality_reduction_summary = 'No dimensionality reduction was performed.' 
    
    else: 
        dimensionality_reduction_summary = f'Dimensionality reduction was performed using PCA.\nIn total {n_components} principal components were used with a cummulative variance of {cummulative_variance}%.' 
    
    #information on the silhouette score 
    best_index = silhouette_scores['Average silhouette score'].idxmax() 
    best_params = '' 
    for params in silhouette_scores.columns: 
        best_params+=f'{params}: {silhouette_scores[params].iloc[best_index]}\n' 
    
    best_params_summary = f'Clustering was performed using the following parameters: \n{best_params} '
    
    #information on all parameters tested 
    if automatic == True: 
        automatic_clustering_summary = f'Automatic clustering was selected. \nThe following parameters were tested: \n\n{silhouette_scores.to_markdown(index=False)}' 
    
    elif automatic == False: 
        automatic_clustering_summary = 'No automatic clustering was performed' 
    
    summary = f"""{log_header()} 

{file_summary} 
{residue_summary} 
{include_distances_summary} 
{include_buriedness_summary} 
{clustering_method_summary} 
{dimensionality_reduction_summary} 

{automatic_clustering_summary} 

{best_params_summary} 

-----------------------------------------------------------------

Cluster details: 
""" 
    
    return summary 