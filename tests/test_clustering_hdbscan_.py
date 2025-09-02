from tripp._clustering_hdbscan_ import hdbscan_clustering
import numpy as np

class TestGreedy:
    def test_greedy_clustering(self):
        frames = np.array([i for i in range(0,31)] * 3).reshape(-1,1)
        trajectory_names = np.array(['traj1'] * 31 + ['traj2'] * 31 + ['traj3'] * 31)
        known_cluster_data = np.load('data/known_cluster_data.npy')
        labels, cluster_centers, cluster_center_indices, cluster_centers_trajectories = hdbscan_clustering(min_cluster_size=20, min_samples=5, clustering_matrix=known_cluster_data,frames=frames,find_centroid=True,trajectory_names=trajectory_names)
        assert len(set(labels)) == 3
        assert cluster_centers == [22, 4, 30]
        assert cluster_center_indices == [53, 4, 92]
        assert cluster_centers_trajectories == ['traj2', 'traj1', 'traj3']
        assert np.allclose(known_cluster_data[cluster_center_indices], [[5.09255669, 4.92882875], [-0.19279722, -0.1965899], [-5., 5.]], atol=1e-6)