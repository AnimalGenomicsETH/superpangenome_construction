#!/usr/bin/env python
import argparse
import gzip
from dataclasses import dataclass


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-g", "--graph", help="graph to calculate statistics", required=True)
    parser.add_argument("-o", "--output", help="where to put the statistics", required=True)
    parser.add_argument("-r", "--ref", help="reference path", required=True, nargs='+')
    parser.add_argument("-t", "--grtype", help="type of graph", choices=["cactus", "pggb", "minigraph"], required=True)
    return parser.parse_args()


def node_calc(node_cum):
    total_node = 0
    len_node = 0
    ref_node = 0
    ref_node_len = 0
    non_ref_node = 0
    non_ref_node_len = 0

    for key, value in node_cum.items():
        id = value.id
        nodelen = value.nodelen
        rtype = value.rtype

        total_node += 1
        len_node += nodelen
        if rtype:
            ref_node += 1
            ref_node_len += nodelen
        else:
            non_ref_node += 1
            non_ref_node_len += nodelen
    return total_node, len_node, ref_node, ref_node_len, non_ref_node, non_ref_node_len


def edge_calc(edge_cum, node_cum):
    total_rr = 0
    total_rn = 0
    total_nn = 0
    total_edges = 0

    for comp in edge_cum:
        parent, child = comp
        parstat = node_cum[parent].rtype
        childstat = node_cum[child].rtype

        total_edges += 1
        if parstat and childstat:
            total_rr += 1
        elif (parstat and not childstat) or (not parstat and childstat):
            total_rn += 1
        elif not parstat and not childstat:
            total_nn += 1
    return total_edges, total_rr, total_rn, total_nn


@dataclass
class node_stats:
    id: str
    nodelen: int
    rtype: int


if __name__ == "__main__":
    args = parse_args()
    graph = args.graph
    grtype = args.grtype
    ref = args.ref
    # path in cactus is doubled from contig name
    if grtype == "cactus":
        ref = []
        #ref = args.ref + "." + args.ref
        for comp in args.ref:
            ref.append(str(comp) + "." + str(comp))

    output_file = args.output

    outfile = open(output_file, "w")

    node_cum = dict()
    #node_stats = namedtuple("node_stats", ["id", "nodelen", "rtype"])

    edge_cum = list()

    if graph[-2:] == "gz":
        input_file = gzip.open(graph, "rt")
    else:
        input_file = open(graph, "r")

    for line in input_file:
        line_comp = line.strip().split()
        if line.startswith("S"):
            node_id = line_comp[1]
            node_len = len(line_comp[2])
            # rtype 0 as nonref
            node_cum[node_id] = node_stats(id=node_id, nodelen=node_len, rtype=0)
            # in minigraph we use the rank to determine ref path
            if grtype == "minigraph":
                node_cum[node_id].rtype = 1 if int(line_comp[-1].split(":")[-1]) == 0 else 0
        elif line.startswith("L"):
            parent_node = line_comp[1]
            child_node = line_comp[3]
            edge_cum.append([parent_node, child_node])
        elif line.startswith("P"):
            if line_comp[1] in ref:
                # split all node in path
                for comp in line_comp[2].split(","):
                    node_comp = comp[:-1]
                    node_cum[node_comp].rtype = 1
    input_file.close()

    # total_node, len_node, ref_node, ref_node_len, non_ref_node, non_ref_node_len = node_calc(node_stats)
    print(graph, grtype, "nodes", *node_calc(node_cum), file=outfile)
    print(graph, grtype, "edges", *edge_calc(edge_cum, node_cum), file=outfile)
    outfile.close()
