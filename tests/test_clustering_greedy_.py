from tripp._clustering_greedy_ import greedy_clustering
from tripp._calculate_rmsd_matrix_ import calculate_rmsd_matrix
import numpy as np

class TestGreedy:
    def test_greedy_clustering(self):
        frames = [i for i in range(0,31)] * 3
        trajectory_names = ['traj1'] * 31 + ['traj2'] * 31 + ['traj3'] * 31
        known_cluster_data = np.load('data/known_cluster_data.npy')
        rmsd_matrix = calculate_rmsd_matrix(clustering_matrix=known_cluster_data, frames=np.array(frames))
        labels, cluster_centers, cluster_center_indices, cluster_centers_trajectories = greedy_clustering(cutoff=3, rmsd_matrix=rmsd_matrix, frames=np.array(frames), trajectory_names=np.array(trajectory_names))
        assert len(set(labels)) == 3
        assert cluster_centers == [2, 7, 0]
        assert cluster_center_indices == [2, 38, 62]
        assert cluster_centers_trajectories == ['traj1', 'traj2', 'traj3']
        assert np.allclose(known_cluster_data[cluster_center_indices], [[-0.07111913, 0.54394042], [5.20969185, 5.50370631], [-5.17115188, 4.58314601]], atol=1e-6)