from tripp._clustering_dbscan_ import dbscan_clustering
import numpy as np

class TestGreedy:
    def test_greedy_clustering(self):
        frames = np.array([i for i in range(0,31)] * 3).reshape(-1,1)
        trajectory_names = np.array(['traj1'] * 31 + ['traj2'] * 31 + ['traj3'] * 31)
        known_cluster_data = np.load('data/known_cluster_data.npy')
        labels, cluster_centers, cluster_center_indices, cluster_centers_trajectories = dbscan_clustering(eps=2, min_samples=20, metric='euclidean', metric_params=None, algorithm='auto', leaf_size = 30, p = None, n_jobs = None, clustering_matrix=known_cluster_data, frames = frames, find_centroid=True, trajectory_names=trajectory_names)
        assert len(set(labels)) == 3
        assert cluster_centers == [4, 22, 30]
        assert cluster_center_indices == [4, 53, 92]
        assert cluster_centers_trajectories == ['traj1', 'traj2', 'traj3']
        assert np.allclose(known_cluster_data[cluster_center_indices], [[-0.19279722, -0.1965899], [5.09255669, 4.92882875], [-5., 5.]], atol=1e-6)