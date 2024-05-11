from tripp import Clustering 
import pandas as pd 
import matplotlib.pyplot as plt 
import seaborn as sns 

#Clustering('/Volumes/chris_drive/Simulations/GPR68/MD/Inactive/M0498/Mutations/I12A/MD1/GPR68I_I12A_M0498_POPC_eq_md_450ns_p_fit_skip100.xtc', '/Volumes/chris_drive/Simulations/GPR68/MD/Inactive/M0498/Mutations/I12A/MD1/GPR68I_I12A_M0498_POPC_min_sd.pdb', '/Users/christos/Downloads/1_pka.csv', [100, 145, 146], '/Users/christos/Downloads/log').kmedoids(automatic=True, n_clusters=2, metric='euclidean', method='alternate', init='heuristic', max_iter=300, random_state=None)   
#Clustering('/Volumes/chris_drive/Simulations/GPR68/MD/Inactive/M0498/Mutations/I12A/MD1/GPR68I_I12A_M0498_POPC_eq_md_450ns_p_fit_skip100.xtc', '/Volumes/chris_drive/Simulations/GPR68/MD/Inactive/M0498/Mutations/I12A/MD1/GPR68I_I12A_M0498_POPC_min_sd.pdb', '/Users/christos/Downloads/1_pka.csv', [100, 145, 146], '/Users/christos/Downloads/log').gromos(automatic=True, max_clusters=20)   
Clustering('/Volumes/chris_drive/Simulations/GPR68/MD/Inactive/M0498/Mutations/I12A/MD1/GPR68I_I12A_M0498_POPC_eq_md_450ns_p_fit_skip100.xtc', '/Volumes/chris_drive/Simulations/GPR68/MD/Inactive/M0498/Mutations/I12A/MD1/GPR68I_I12A_M0498_POPC_min_sd.pdb', '/Users/christos/Downloads/1_pka.csv', [100, 145, 146], '/Users/christos/Downloads/log').dbscan(automatic=True)    


df = pd.read_csv('/Users/christos/Downloads/1_pka_cluster.csv') 

fig, ax = plt.subplots(1,3) 

sns.scatterplot(data=df, x='GLU100', y='LYS145', hue='Clusters', ax=ax[0], s=1, palette='tab10') 
sns.scatterplot(data=df, x='LYS145', y='GLU146', hue='Clusters', ax=ax[1], s=1, palette='tab10') 
sns.scatterplot(data=df, x='GLU146', y='GLU100', hue='Clusters', ax=ax[2], s=1, palette='tab10') 

plt.show() 