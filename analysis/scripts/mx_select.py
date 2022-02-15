from collections import defaultdict
from .utils import read_markers, write_list


def read_genes(genes_fname, genes=defaultdict()):
    with open(genes_fname) as f:
        for idx, line in enumerate(f.readlines()):
            gene = line.strip()
            genes[gene] = idx


def sel_genes(genes, marker_genes, sel=[]):
    mg_inv = {v: k for k, v in marker_genes.items()}
    for idx in range(len(mg_inv)):
        # this maps the marker gene name index to the gene index
        # in order of the marker_genes file
        sel.append(genes[mg_inv[idx]])


def mx_select(markers_fname, genes_fname, out_select_fn):
    # select should be extensible to axis and genes -> md (metadata)
    markers_ec = defaultdict(list)
    celltypes = defaultdict()
    marker_genes = defaultdict()
    # this is duplicated from index, not ideal but w/e maybe ok
    # ideally would want to give it markers.ec
    read_markers(markers_fname, markers_ec, celltypes, marker_genes)

    genes = defaultdict()
    read_genes(genes_fname, genes)

    sel = []
    sel_genes(genes, marker_genes, sel)
    write_list(out_select_fn, sel)
