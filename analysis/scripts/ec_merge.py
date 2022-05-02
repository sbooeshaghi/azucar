from collections import defaultdict
from .utils import read_markers_str, write_markers


def merge(ecs1, ecs2):
    m = defaultdict(set)

    uniq = set(list(ecs1.keys()) + list(ecs2.keys()))
    for k in uniq:
        m[k].update(ecs2.get(k, []))
        m[k].update(ecs1.get(k, []))

    # convert from set to list
    m = dict(zip(m.keys(), map(list, m.values())))
    return m


def ec_merge(ec_a_fn, ec_b_fn, ec_merge_fn):
    ec_a = defaultdict(list)
    read_markers_str(ec_a_fn, ec_a)

    ec_b = defaultdict(list)
    read_markers_str(ec_b_fn, ec_b)

    m = merge(ec_a, ec_b)

    write_markers(ec_merge_fn, m)
