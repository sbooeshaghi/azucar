from collections import defaultdict
import numpy as np


def nd(arr):
    return np.asarray(arr).reshape(-1)


def read_markers(
    fname,
    markers_ec=defaultdict(list),
    celltype=defaultdict(),
    marker_genes=defaultdict(),
):
    with open(fname, "r") as f:
        for idx, line in enumerate(f.readlines()):
            ct, genes = line.strip().split("\t")
            celltype[ct] = idx

            # two things
            # 1. make marker_genes list
            # 2. make markers_ec
            for g in genes.split(","):
                gidx = len(marker_genes)

                # check if the gene has been added already
                if g in marker_genes.keys():  # gene repeated
                    gidx = marker_genes[g]
                else:
                    marker_genes[g] = gidx

                # for the cell type index, add the marker gene index
                markers_ec[celltype[ct]].append(marker_genes[g])

            # sort the marker genes
            markers_ec[celltype[ct]] = sorted(markers_ec[celltype[ct]])


def read_int_list(fname, lst=[]):
    with open(fname) as f:
        for idx, i in enumerate(f.readlines()):
            lst.append(int(i.strip()))


def read_str_list(fname, lst=list):
    with open(fname, "r") as f:
        for idx, line in enumerate(f.readlines()):
            lst.append(line.strip())


def write_list(fname, lst=list):
    with open(fname, "w") as f:
        for idx, ele in enumerate(lst):
            f.write(f"{ele}\n")


def write_markers(fname, markers):
    with open(fname, "w") as f:
        for k, v in markers.items():
            f.write(f"{k}\t")
            n = len(v)
            for idx, i in enumerate(v):
                f.write(f"{i}")
                if idx < n - 1:
                    f.write(",")
            f.write("\n")


def map_dict_list_keys_values(dct, k_lst, v_lst, nd=defaultdict(list)):
    for k, v in dct.items():
        nd[k_lst[k]] = [v_lst[i] for i in v]


def map_dict_list_keys(dct, k_lst, v_lst, nd=defaultdict(list)):
    for k, v in dct.items():
        nd[k_lst[k]] = v


def map_dict_list_values(dct, k_lst, v_lst, nd=defaultdict(list)):
    for k, v in dct.items():
        nd[k] = [v_lst[i] for i in v]


def read_markers_ec(fname, markers_ec=defaultdict(list)):
    with open(fname, "r") as f:
        for idx, line in enumerate(f.readlines()):
            ct_id, gene_ids = line.strip().split("\t")
            markers_ec[int(ct_id)] = [int(i) for i in gene_ids.split(",")]


def read_markers_str(fname, markers=defaultdict(list)):
    with open(fname, "r") as f:
        for idx, line in enumerate(f.readlines()):
            ct_id, gene_ids = line.strip().split("\t")
            markers[ct_id] = [i for i in gene_ids.split(",")]


def drop_markers(markers_ec, drop_ids):
    if len(drop_ids) == 0:
        return

    for k, v in markers_ec.items():

        gidx = len(v) - 1
        while gidx > -1:
            mg = markers_ec[k][gidx]

            if mg in drop_ids:
                markers_ec[k].pop(gidx)
            else:
                to_sub = 0
                for d in drop_ids:
                    if d < mg:
                        to_sub += 1
                markers_ec[k][gidx] -= to_sub
            gidx -= 1


# testing data
# drop_ids = set([2, 3, 34, 42])
# truth = {0: [0, 1, 2, 3, 4, 5, 6, 7, 8],
#         1: [7, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18],
#         2: [19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29],
#         3: [0, 3,  6,  7, 24, 29, 30, 31, 32, 33],
#         4: [0, 4,  5,  6,  7, 18, 23, 24, 30, 31, 34],
#         5: [2, 22, 23, 24, 30, 35, 36, 37, 38, 39],
#         6: [0, 3,  4,  6,  7, 24, 30, 31, 32, 40]}
# drop_markers(markers_ec, set(drop_genes))
# markers_ec == truth


def get_marker_centroids(X, markers_ec, method="mean"):
    n_clusters = len(list(markers_ec.keys()))
    _, n_features = X.shape

    marker_centroids = np.ones((n_clusters, n_features)) * 1e-5

    for k, v in markers_ec.items():
        submx = X[:, v]
        if method == "max":
            repl = submx.max(0)
        else:
            repl = submx.mean(0)
        marker_centroids[k][v] = repl
    return marker_centroids


def get_centroids(X, z):
    clusters = np.sort(np.unique(z))

    (n_clusters,) = clusters.shape
    _, n_features = X.shape

    centroids = np.ones((n_clusters, n_features))
    for i, g in enumerate(clusters):
        centroids[i] = X[np.where(z == g)[0]].mean(0)
    return centroids
