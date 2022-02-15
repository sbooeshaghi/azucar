from .utils import read_str_list, read_markers_ec, nd, get_marker_centroids
from .gmm import ImprovedGaussianMixture
from scipy.stats import entropy
from sklearn.neighbors import KDTree
from collections import defaultdict
import numpy as np
from scipy.io import mmread

# @title mx assign
def mx_assign(
    matrix_fn: str,
    barcodes_fn: str,
    genes_fn: str,
    markers_ec_fn: str,
    groups_fn: str,
    out_assignments_fn: str,
):
    # mx assign assumes matrix has been filtered and that genes are ordered by their numbering in matrix_ec
    groups = []
    read_str_list(groups_fn, groups)

    genes = []
    read_str_list(genes_fn, genes)

    barcodes = []
    read_str_list(barcodes_fn, barcodes)

    markers_ec = defaultdict(list)
    read_markers_ec(markers_ec_fn, markers_ec)

    # read in gene count matrix
    G = mmread(matrix_fn).toarray()

    n_clusters = len(markers_ec.keys())
    n_samples, n_features = G.shape
    # print(n_clusters, *G.shape, sep=", ")

    # get the centroids for the existing data

    data = {"X": {"raw": G, "raw_log1p": np.log1p(G), "norm": G / nd(G.sum())[:, None]}}

    method = "raw"

    # have to initialize the clusters by first mean centering alternative is to zscore the means

    X_mean = data["X"][method].mean(0)
    X_init = data["X"][method] - X_mean
    centroids_init = get_marker_centroids(X_init, markers_ec, "max")

    tree = KDTree(centroids_init, metric="euclidean")
    nearest_dist, nearest_ind = tree.query(X_init, k=1)

    # assign cells to clusters
    p = 1
    resp = np.ones((n_samples, n_clusters)) * (1 - p) / (n_clusters - 1)
    resp[np.arange(n_samples), nearest_ind.flatten()] = p

    # initialize params
    nk = resp.sum(axis=0) + 10 * np.finfo(resp.dtype).eps

    # then once we have the means, add the previously subtracted means back
    means_init = np.dot(resp.T, X_init) / nk[:, np.newaxis]
    means_init += X_mean

    # alternative to uniform weights is nk / n_samples (using the new assignments)
    uniform_weights = np.array([1.0 / n_clusters] * n_clusters)

    # alternative is to compute precisions by first doing M-step to get gaus params
    identity_precisions = np.repeat(
        np.array([np.eye(data["X"][method].shape[1])]), n_clusters, 0
    )

    gmm_params = {
        "n_components": n_clusters,
        "means_init": None,  # to be added
        "weights_init": None,  # to be added
        "precisions_init": None,  # to be added
        "random_state": 0,
        "reg_covar": 1e-8,
        "verbose": 2,
        "n_init": 1,
        "max_iter": 1000,
        "tol": 1e-4,
        "init_params": "random",
        "covariance_type": "full",
    }

    params = {
        **gmm_params,
        "means_init": means_init,  # centroids,
        "weights_init": uniform_weights,
        "precisions_init": identity_precisions,
    }

    # actually do the GMM
    gmm = ImprovedGaussianMixture(**params)

    # print(method, json.dumps(params, indent=4, default=str))

    labels = gmm.fit_predict(data["X"][method], B=markers_ec)
    means = gmm.means_
    prob = gmm.predict_proba(data["X"][method])
    ent = entropy(prob, axis=1)

    # make df
    df = pd.DataFrame(G, columns=[f"{i}_ipf" for i in genes])

    df["label_id"] = labels
    df["label"] = df["label_id"].map({i: groups[i] for i in range(len(groups))})
    df["ent"] = ent
    df.index = barcodes

    for idx, p in enumerate(prob.T):
        df[f"mahalanobis_{idx}"] = p

    df.to_csv(out_assignments_fn, index=False, sep="\t")
