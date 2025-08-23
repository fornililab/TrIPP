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

from sklearn.decomposition import PCA

def pca(clustering_matrix): 
    """
    Performs PCA on the feature matrix to reduce dimensionality. 
    Parameters
    ----------
    clustering_matrix: np.ndarray
        The clustering (feature) matrix to be transformed.
    Returns
    -------
    n_components: int
        The number of principal components (PCs) used for dimensionality reduction. 
        It is determined as the minimum number of principal components that explain
        at least 90% of the total variance. 
    cummulative_variance: float
        The cumulative variance explained by the selected principal components, expressed as a percentage.
    clustering_matrix_transformed: np.ndarray
        The transformed clustering matrix after PCA.
    """
    pca_class = PCA() 
    pca_class.fit(clustering_matrix) 
    var_ratio = pca_class.explained_variance_ratio_ 
    n_components = 0
    cummulative_variance = 0
    for item in var_ratio: 
        n_components += 1
        cummulative_variance += item 
        if cummulative_variance >= 0.9: 
            break 
    
    cummulative_variance = round(cummulative_variance*100, 2) 
    
    pca_class = PCA(n_components=n_components) 
    clustering_matrix_transformed = pca_class.fit_transform(clustering_matrix) 

    return n_components, cummulative_variance, clustering_matrix_transformed 