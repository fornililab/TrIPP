from tripp._clustering_kmedoids_ import kmedoids_clustering
import numpy as np

# centroids = np.array([[0, 0], [5, 5], [-5, 5]])
# cluster_sizes = [30, 30, 30]
# data = []
# for i, centroid in enumerate(centroids):
#     num_points = cluster_sizes[i]
#     cluster_data = np.random.normal(loc=centroid, scale=1.0, size=(num_points, len(centroid)))
#     data.append(cluster_data)
#     data.append(centroids[i])

# known_cluster_data = np.vstack(data)
# np.save('data/known_cluster_data.npy', known_cluster_data)
class TestKMedoids:
    def test_kmedoids_clustering(self):
        frames = [i for i in range(0,31)] * 3
        trajectory_names = ['traj1'] * 31 + ['traj2'] * 31 + ['traj3'] * 31
        known_cluster_data = np.load('data/known_cluster_data.npy')
        labels, cluster_centers, cluster_center_indices, cluster_centers_trajectories = kmedoids_clustering(n_clusters=3, clustering_matrix=known_cluster_data, frames=np.array(frames), trajectory_names=np.array(trajectory_names), metric='euclidean', method='alternate', init='k-medoids++', max_iter=300, random_state=20)
        assert len(set(labels)) == 3
        assert cluster_centers == [30, 30, 30]
        assert cluster_center_indices == [92, 61, 30]
        assert cluster_centers_trajectories == ['traj3', 'traj2', 'traj1']
        assert (known_cluster_data[cluster_center_indices] == [[-5, 5], [5, 5], [0, 0]]).all()