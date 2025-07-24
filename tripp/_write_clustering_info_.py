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

from tripp._sort_clusters_ import sort_clusters 
import pandas as pd 
import MDAnalysis as mda 

def write_clustering_info(summary, trajectory_dict, pka_df, times, frames, trajectory_names, labels, cluster_centers, cluster_center_indices, cluster_centers_trajectories, output_directory, output_prefix, clustering_method): 
    """
    Write clustering information to a log file and save cluster centers as PDB files.
    Parameters
    ----------
    summary : str
        A summary string for the clustering parameters coming from generate_clustering_summary
        function.
    trajectory_dict : dict
        A dictionary where keys are trajectory names and values are MDAnalysis Universe objects.
    pka_df : pd.DataFrame
        A DataFrame containing pKa values with rows indexed by time and columns by residue identifier.
    times : np.ndarray
        An array of times corresponding to the clustering points.
    frames : np.ndarray
        An array of frames corresponding to the clustering points.
    trajectory_names : np.ndarray
        A 1D numpy array containing the names of the trajectories corresponding to the clustering data.
    labels : np.ndarray
        An array of cluster labels for each point in the clustering matrix.
    cluster_centers : list
        A list of cluster centers.
    cluster_center_indices: list
        A list of indices of the cluster centers mapped onto the clustering matrix. 
    cluster_centers_trajectories : list
        A list of indices of the cluster centers mapped onto the trajectory.
    output_directory : str
        The path to the directory where output files will be saved.
    output_prefix : str
        A prefix for the clustering output files.
    clustering_method : str
        The clustering method used (e.g., 'DBSCAN', 'KMedoids').
    Returns
    -------
    labels : np.ndarray
        An array of cluster labels for each point in the clustering matrix sorted by cluster population.
    cluster_center_indices: list
        A list of indices of the cluster centers mapped onto the clustering matrix sorted by cluster population.
    """
    labels, cluster_centers, cluster_center_indices, cluster_centers_trajectories = sort_clusters(labels=labels, cluster_centers=cluster_centers, cluster_center_indices=cluster_center_indices, cluster_centers_trajectories=cluster_centers_trajectories)

    def write_clustering_log():
        df = pd.DataFrame({'Times' : times.flatten(), 'Frames' : frames.flatten(), 'Labels' : labels.flatten(), 'Trajectories' : trajectory_names.flatten()}) 
        cluster_data = {} 
        for i in range(len(cluster_centers)): 
            df_i = df[df['Labels']==i]
            cluster_data[f'Cluster {i}'] = {'Population' : f'{round((len(df_i)/len(df))*100,2)}%', 
                                            'Centroid frame' : str(cluster_centers[i]), 
                                            'Centroid time' : str(times.flatten()[cluster_center_indices[i]]), 
                                            'Centroid trajectory' : str(trajectory_names.flatten()[cluster_center_indices[i]])} 
            for traj in trajectory_dict.keys(): 
                cluster_data[f'Cluster {i}'][traj] = {'Times' : ', '.join(map(str, df_i[df_i['Trajectories']==traj]['Times'])), 
                                                      'Frames' : ', '.join(map(str, df_i[df_i['Trajectories']==traj]['Frames']))} 
        
        if clustering_method == 'DBSCAN': 
            df_i = df[df['Labels']==-1]
            cluster_data['Cluster -1 (Outliers)'] = {'Population' : f'{round((len(df_i)/len(df))*100,2)}%'} 
            for traj in trajectory_dict.keys(): 
                cluster_data['Cluster -1 (Outliers)'][traj] = {'Times' : ', '.join(map(str, df_i[df_i['Trajectories']==traj]['Times'])), 
                                                               'Frames' : ', '.join(map(str, df_i[df_i['Trajectories']==traj]['Frames']))} 
        
        with open(f'{output_directory}/{output_prefix}_{clustering_method}.log', 'w') as l: 
            l.write(summary) 
            for i in range(len(cluster_centers)): 
                l.write(f'\nCluster {i}\n\n') 
                for key in cluster_data[f'Cluster {i}'].keys(): 
                    if type(cluster_data[f'Cluster {i}'][key]) == dict: 
                        l.write(f'{key}\n\n'+'Times'+'\n'+cluster_data[f'Cluster {i}'][key]['Times']+'\n\n') 
                        l.write(f'Frames'+'\n'+cluster_data[f'Cluster {i}'][key]['Frames']+'\n\n') 
                    else: 
                        l.write(f'{key}\t\t'+cluster_data[f'Cluster {i}'][key]+'\n\n') 
            
            if clustering_method == 'DBSCAN': 
                l.write('\nCluster -1 (Outliers)\n\n') 
                for key in cluster_data[f'Cluster -1 (Outliers)'].keys(): 
                    if type(cluster_data['Cluster -1 (Outliers)'][key]) == dict: 
                        l.write(f'{key}\n\n'+'Times'+'\n'+cluster_data['Cluster -1 (Outliers)'][key]['Times']+'\n\n') 
                        l.write(f'Frames'+'\n'+cluster_data['Cluster -1 (Outliers)'][key]['Frames']+'\n\n') 
                    else: 
                        l.write(f'{key}\t\t'+cluster_data['Cluster -1 (Outliers)'][key]+'\n\n') 
    
    write_clustering_log()  
    
    def write_cluster_centers(): 
        for index, cluster_center_frame in enumerate(cluster_centers): 
            universe = trajectory_dict[cluster_centers_trajectories[index]] 
            for ts in universe.trajectory: 
                if ts.frame == cluster_center_frame: 
                    pdb_file = f'{output_directory}/{output_prefix}_{clustering_method}_C{index}.pdb' 
                    with mda.Writer(pdb_file) as w: 
                        w.write(universe) 
    
    write_cluster_centers() 

    def write_new_dataframe(): 
        pka_df['Clusters'] = labels 
        pka_df.to_csv(f'{output_directory}/{output_prefix}_{clustering_method}_cluster.csv') 
    
    write_new_dataframe() 
    return labels, cluster_center_indices
