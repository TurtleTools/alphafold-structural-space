from portein import get_best_transformation, apply_transformation, find_size
import prody as pd
from geometricus import MomentInvariants, SplitType
from scipy.signal import resample
import numpy as np
import matplotlib.pyplot as plt


def get_coords_topic_scores(coords, topic_id, h_matrix_norm, shapemer_to_index):
    def shapemer_to_topic_value(s_string):
        if s_string in shapemer_to_index:
            return h_matrix_norm[topic_id][shapemer_to_index[s_string]]
        else:
            return 0

    resolution_kmer = 4
    resolution_radius = 6
    weights = np.zeros(coords.shape[0])
    protein_invariants = MomentInvariants.from_coordinates("protein_id",
                                                           coords,
                                                           None,
                                                           split_type=SplitType.KMER_CUT)

    def get_similarity(x1, x2, gamma=0.03):
        return np.exp(-gamma * np.sum((coords[x1] - coords[x2]) ** 2, axis=-1))

    shapemers = (np.log1p(protein_invariants.moments) * resolution_kmer).astype(int)
    for i, x in enumerate(shapemers):
        weight = shapemer_to_topic_value(f"k{x[0]}i{x[1]}i{x[2]}i{x[3]}")
        for index in range(i, i + 16):
            weights[i] += weight * get_similarity(i, index)
    protein_invariants = MomentInvariants.from_coordinates(
        "protein_id",
        coords,
        None,
        split_type=SplitType.RADIUS,
        split_size=10)
    shapemers = (np.log1p(protein_invariants.moments) * resolution_radius).astype(int)
    for i, x in enumerate(shapemers):
        weight = shapemer_to_topic_value(f"r{x[0]}i{x[1]}i{x[2]}i{x[3]}")
        for index in protein_invariants.split_indices[i]:
            weights[index] += weight * get_similarity(i, index)
    return weights

def get_protein_topic_scores(path, topic_id, h_matrix_norm, shapemer_to_index, matplotlib=True):
    pdb = pd.parsePDB(str(path))
    pdb_alpha = pdb.select("protein and calpha")
    opacities = pdb_alpha.getBetas() / 100
    coords = pdb_alpha.getCoords()
    weights = get_coords_topic_scores(coords, topic_id, h_matrix_norm, shapemer_to_index)
    if matplotlib:
        coords = apply_transformation(coords, get_best_transformation(coords))
        return coords, weights, opacities
    else:
        matrix = get_best_transformation(coords)
        pdb = pd.applyTransformation(pd.Transformation(matrix), pdb)
        for i, res in enumerate(pdb.iterResidues()):
            res.setBetas([weights[i]] * len(res))
        return pdb


def plot_protein(coords, weights, opacities, max_value, upsample_rate=3):
    coords = resample(coords[:, :2], upsample_rate * coords.shape[0])
    weights = np.repeat(weights, upsample_rate)
    opacities = np.repeat(opacities, upsample_rate)
    colors = [plt.cm.coolwarm(int(256 * (x / max_value))) for x in weights]
    fig, ax = plt.subplots(figsize=find_size(coords, height=5, width=None))
    for i in range(coords.shape[0] - upsample_rate):
        ax.plot(coords[:, 0][i:i + 2], coords[:, 1][i:i + 2],
                lw=2, color=colors[i], alpha=opacities[i])
    plt.axis("off")
    return fig, ax
