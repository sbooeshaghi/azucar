import numpy as np
from ipfn import ipfn

from .utils import nd
from scipy.io import mmread, mmwrite
from scipy.sparse import csr_matrix
from sklearn.preprocessing import normalize


def do_ipf(mtx, axis_uniform=1):
    rows, cols = mtx.shape

    n = 1
    if axis_uniform == 1:
        aggregates = [
            np.ones(rows) * n / rows,  # rows, each cell uniform
            nd(mtx.sum(0)) / mtx.sum(),  # columns, each tag proportional
        ]
    elif axis_uniform == 0:
        aggregates = [
            nd(mtx.sum(1)) / mtx.sum(),  # rows, each cell proportional
            np.ones(cols) * n / cols,  # columns, each tag uniform
        ]
    elif axis_uniform == -1:
        aggregates = [
            np.ones(rows) * n / rows,  # rows, each cell uniform,
            np.ones(cols) * n / cols,  # columns, each tag uniform
        ]
    dimensions = [[0], [1]]

    IPF = ipfn.ipfn(mtx, aggregates, dimensions, max_iteration=100_000)
    m = IPF.iteration()
    return m


def mx_norm(matrix_fn, out_matrix_fn, how="ipf", target_sum=None):
    mtx = mmread(matrix_fn).toarray()

    if how == "ipf":
        mtx_ipf = do_ipf(
            mtx.copy(), axis_uniform=1
        )  # 0: uniform cols, 1: uniform rows, -1: both uniform
        s = mtx.sum()
        if target_sum:
            s = target_sum
        mmwrite(out_matrix_fn, csr_matrix(mtx_ipf * s))
    elif how == "log1p":
        log1p = np.log1p(mtx)
        mtx_log1p = csr_matrix(log1p)
        mmwrite(out_matrix_fn, mtx_log1p)
    elif how == "sf":
        sf = mtx.sum(1).flatten()
        mtx_sf = mtx / (sf / sf.mean())[:, None]
        if target_sum:
            mtx_sf = normalize(mtx, norm="l1") * target_sum
        mmwrite(out_matrix_fn, csr_matrix(mtx_sf))


def read_str_list(fname, lst=list):
    with open(fname, "r") as f:
        for idx, line in enumerate(f.readlines()):
            lst.append(line.strip())
