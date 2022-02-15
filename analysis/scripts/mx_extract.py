from scipy.io import mmread, mmwrite
from scipy.sparse import csr_matrix
from .utils import write_list, read_int_list, read_str_list


def mx_extract(matrix_fn, md_fn, select_fn, out_matrix_fn, out_md_fn, axis=1):
    M = mmread(matrix_fn).toarray()

    # column indices to select from gene matrix
    sel = []
    read_int_list(select_fn, sel)

    # read in axis metadata (in this case its the gene list)
    md = []
    read_str_list(md_fn, md)

    # # markers.ec, maps groups (indices) to marker genes (indices)
    # markers_ec = defaultdict(list)
    # read_markers_ec(markers_ec_fn, markers_ec)

    # # the gene names
    # genes = []
    # read_str_list(genes_fn, genes)

    # TODO MAKE SURE DROP GENES USED
    # drop_genes = np.arange(M.shape[1])[~col_mask]
    # drop_markers(markers_ec, set(drop_genes))

    # mtx = M[row_mask][:,col_mask].astype(int)
    # mtx_ipf = do_ipf(mtx.copy())

    # mmwrite(matrix_e_fn, csr_matrix(mtx_ipf[:,sel]))
    mmwrite(out_matrix_fn, csr_matrix(np.take(M, sel, axis=axis)))
    write_list(out_md_fn, np.take(md, sel))
