import numpy as np
from .utils import read_str_list, write_list
from scipy.io import mmread, mmwrite
from scipy.sparse import csr_matrix


def sanitize_mtx(mtx: np.ndarray):
    cell_count_mask = mtx.sum(1) > 0  # count for each cell
    gene_count_mask = mtx.sum(0) > 0  # count for each gene

    genes_detected_mask = (mtx > 0).sum(1) > 0  # n genes per cell
    cells_detected_mask = (mtx > 0).sum(0) > 0  # n cells per gene
    row_mask = np.logical_and(cell_count_mask, genes_detected_mask)
    col_mask = np.logical_and(gene_count_mask, cells_detected_mask)

    return (row_mask, col_mask)


def mx_sanitize(
    matrix_fn: str,
    barcodes_fn: str,
    genes_fn: str,
    out_matrix_fn: str,
    out_barcodes_fn: str,
    out_genes_fn: str,
):
    # the barcode names
    barcodes = []
    read_str_list(barcodes_fn, barcodes)

    # the gene names
    genes = []
    read_str_list(genes_fn, genes)

    # read in matrix and select columns and write back to disc
    M = mmread(matrix_fn).toarray()  # noqa

    # sanitize gene count matrix (remove cells / genes) and remove genes
    # from  marker_ec
    row_mask, col_mask = sanitize_mtx(M)
    barcodes = np.array(barcodes)[row_mask]
    genes = np.array(genes)[col_mask]

    mtx = M[row_mask][:, col_mask].astype(int)

    # write everything
    mmwrite(out_matrix_fn, csr_matrix(mtx))
    write_list(out_barcodes_fn, barcodes)
    write_list(out_genes_fn, genes)
