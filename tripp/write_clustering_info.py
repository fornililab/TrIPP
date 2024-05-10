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

from tripp.sort_clusters import sort_clusters 
import pandas as pd 
import MDAnalysis as mda 

def write_clustering_info(universe, pka_file, times, frames, labels, cluster_centers, cluster_indices, log_file, clustering_method, sil_score): 
    
    labels, cluster_centers, cluster_indices = sort_clusters(labels=labels, cluster_indices=cluster_indices, cluster_centers=cluster_centers, clustering_method=clustering_method) 
    
    def write_clustering_log(): 
        df = pd.DataFrame({'Times' : times.flatten(), 'Frames' : frames.flatten(), 'Labels' : labels.flatten()}) 
        cluster_data = {} 
        for i in range(len(cluster_centers)): 
            df_i = df[df['Labels']==i]
            cluster_data[f'Cluster {i}'] = {'Population' : f'{round((len(df_i)/len(df))*100,2)}%', 
                                            'Centroid frame' : str(cluster_centers[i]), 
                                            'Centroid time' : str(times.flatten()[cluster_indices[i]]), 
                                            'Times' : ', '.join(map(str, df_i['Times'])), 
                                            'Frames' : ', '.join(map(str, df_i['Frames']))} 
        
        if clustering_method == 'DBSCAN': 
            df_i = df[df['Labels']==-1]
            cluster_data['Cluster -1 (Outliers)'] = {'Population' : f'{round((len(df_i)/len(df))*100,2)}%', 
                                                     'Times' : ', '.join(map(str, df_i['Times'])), 
                                                     'Frames' : ', '.join(map(str, df_i['Frames']))} 
        
        header = """
The Trajectory Iterative pKa Predictor (TrIPP) 
Written by: Christos Matsingos, Ka Fu Man, and Arianna Fornili 

If you are using TrIPP, please cite: 

-----------------------------------------------------------------
"""

        summary = f"""
Summary of clustering results: 

Clustering done using the {clustering_method} method. 
A total of {len(set(labels))} clusters were generated. 
Calculated average silhouette score: {sil_score}. 

"""
        
        with open(f'{log_file}_{clustering_method}.log', 'w') as l: 
            l.write(header)
            l.write(summary) 
            for i in range(len(cluster_centers)): 
                l.write(f'\nCluster {i}\n\n') 
                for key in cluster_data[f'Cluster {i}'].keys(): 
                    l.write(f'{key}\t\t'+cluster_data[f'Cluster {i}'][key]+'\n\n') 
            
            if clustering_method == 'DBSCAN': 
                l.write('\nCluster -1 (Outliers)\n\n') 
                for key in cluster_data[f'Cluster -1 (Outliers)'].keys(): 
                    l.write(f'{key}\t\t'+cluster_data['Cluster -1 (Outliers)'][key]+'\n\n') 
    
    write_clustering_log() 
    
    def write_cluster_centers(): 
        for index, cluster_index in enumerate(cluster_indices): 
            for ts in universe.trajectory: 
                if ts.frame == cluster_index: 
                    pdb_file = f'{log_file}_C{index}.pdb' 
                    with mda.Writer(pdb_file) as w: 
                        w.write(universe) 
    
    write_cluster_centers() 

    def write_new_dataframe(): 
        pka_df = pd.read_csv(pka_file) 
        pka_df['Clusters'] = labels 
        pka_df.to_csv(pka_file.replace('.csv', '_cluster.csv'), index=False) 
    
    write_new_dataframe() 
