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

import numpy as np
from tripp._calculate_charge_center_distances_ import calculate_charge_center_distances
from tripp._determine_charge_center_ import determine_charge_center
from scipy.stats import zscore
from tripp._create_mda_universe_ import create_mda_universe


def create_clustering_matrix(trajectory_file, topology_file, pka_df,
                             selections, include_distances,
                             buriedness_df, include_buriedness):
    """
    Function that generates the clustering_matrix. The clustering_matrix
    is a 2D numpy array, where columns correspond to normalized pKa
    values, extracted from pka_df, and normalized interresidue
    distances, calculated using the calculate_charge_center_distances
    function. Normalization is done using the z-score.
    
    Parameters
    ----------
    trajectory_file: str or dict
        When str, it is the path of the file containing the trajectory. The same
        formats permited by MDAnalysis can be used. When dict, the clustering is done using
        multiple trajectories specified in a dictionary, where the name of each trajectory
        is used as a key and the path as an object: {'MD1' : 'file1', 'MD2' : 'file2', ...}.
        The same topology file is used for all trajectories.
    topology_file: str
        Path to the topology file. The same formats allowed by MDAnalysis can be used.
    pka_df: pd.DataFrame
        A DataFrame containing pKa values indexed by time, residue identifier, and
        trajectory name. This DataFrame is processed by make_pka_or_buriendess_df function
        in Clustering class.
    selections: list of str
        List of selections (MDAnalysis selection algebra) for which the clustering will be done.
        The residues have a pKa value assigned to them by PROPKA. Note if your system is
        multichain, your topology and selection must include chain identity.
    include_distances: bool
        If True, the relative positions (as distances between the charge centers) are used as
        additional features for the clustering alongside the pKa values.
    buriedness_df: pd.DataFrame
        A DataFrame containing buried ratio (buriedness) values indexed by time, residue identifier,
        and trajectory name. This DataFrame is processed by make_pka_or_buriendess_df function
        in Clustering class.
    include_buriedness: bool
        If True, buriedness values will be included in the clustering matrix.
        Default is False.
    Returns
    -------
    clustering_matrix: np.ndarray
        A 2D numpy array containing the normalized pKa values and distances.
    times: np.ndarray
        A 1D numpy array containing the times corresponding to the clustering data.
    frames: np.ndarray
        A 1D numpy array containing the frames corresponding to the clustering data.
    trajectory_names: np.ndarray
        A 1D numpy array containing the names of the trajectories corresponding to the clustering data.
    trajectory_dict: dict
        A dictionary where keys are trajectory names and values are MDAnalysis Universe objects.
    """

    if isinstance(trajectory_file, str):
        trajectory_dict = {
            "Unnamed trajectory": create_mda_universe(
                topology_file=topology_file,
                trajectory_file=trajectory_file
                )
            }

    elif isinstance(trajectory_file, dict):
        trajectory_dict = {}
        for key in trajectory_file.keys():
            traj = trajectory_file[key]
            trajectory_dict[key] = create_mda_universe(
                topology_file=topology_file, trajectory_file=traj
            )

    for index, key in enumerate(trajectory_dict.keys()):
        df_traj_pka = pka_df[pka_df["Trajectories"] == key]
        if include_buriedness:
            df_traj_buriedness = buriedness_df[buriedness_df["Trajectories"] == key]
        elif not include_buriedness:
            df_traj_buriedness = None
        if index == 0:
            clustering_matrix, times, frames, trajectory_names = (
                extract_pka_distances_buriedness(key,
                                                 trajectory_dict[key],
                                                 df_traj_pka,
                                                 df_traj_buriedness,
                                                 selections, 
                                                 include_distances,
                                                 include_buriedness)
            )

        else:
            partial_clustering_matrix, partial_times, partial_frames, partial_trajectory_names = (
                extract_pka_distances_buriedness(key, 
                                                 trajectory_dict[key],
                                                 df_traj_pka,
                                                 df_traj_buriedness,
                                                 selections,
                                                 include_distances,
                                                 include_buriedness))
            clustering_matrix = np.concatenate((clustering_matrix, partial_clustering_matrix), axis=0)
            times = np.concatenate((times, partial_times), axis=0)
            frames = np.concatenate((frames, partial_frames), axis=0)
            trajectory_names = np.concatenate((trajectory_names, partial_trajectory_names), axis=0)

    return clustering_matrix, times, frames, trajectory_names, trajectory_dict


def extract_pka_distances_buriedness(trajectory_name, universe,
                                     df_traj_pka, df_traj_buriedness,
                                     selections, include_distances,
                                     include_buriedness):

        pka = []
        times = []
        frames = []
        trajectory_names = []
        d = []
        b = []

        for ts in universe.trajectory:
            times.append([ts.time])
            frames.append([ts.frame])
            trajectory_names.append(trajectory_name)
            positions = []
            pkas = []
            buriedness = []
            for selection in selections:
                charge_center, residue_identifier = determine_charge_center(
                                    universe, selection
                                )
                
                pka_residue = df_traj_pka.loc[ts.time, residue_identifier]
                pkas.append(pka_residue)

                if include_distances is True:
                    positions.append(charge_center)

                if include_buriedness is True:
                    buriedness_residue = df_traj_buriedness.loc[
                        ts.time, residue_identifier
                    ]
                    buriedness.append(buriedness_residue)

            pka.append(np.array(pkas))

            if include_distances is True:
                distances = calculate_charge_center_distances(positions)
                d.append(distances)

            if include_buriedness is True:
                b.append(buriedness)

        pka = np.array(pka)

        if include_distances is True and include_buriedness is True:
            d = np.array(d)
            b = np.array(b)
            clustering_matrix = np.concatenate(
                (zscore(d, axis=None), zscore(b, axis=None), zscore(pka, axis=None)),
                axis=1,
            )

        elif include_distances is True and include_buriedness is False:
            d = np.array(d)
            clustering_matrix = np.concatenate(
                (zscore(d, axis=None), zscore(pka, axis=None)), axis=1
            )

        elif include_distances is False and include_buriedness is True:
            b = np.array(b)
            clustering_matrix = np.concatenate(
                (zscore(b, axis=None), zscore(pka, axis=None)), axis=1
            )

        elif include_distances is False and include_buriedness is False:
            clustering_matrix = zscore(pka, axis=None)

        times = np.array(times)
        frames = np.array(frames)
        trajectory_names = np.array(trajectory_names)

        return clustering_matrix, times, frames, trajectory_names
