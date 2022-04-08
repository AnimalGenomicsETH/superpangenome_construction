#!/usr/bin/env python
import argparse
import gzip
from dataclasses import dataclass
from collections import defaultdict
import re


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-g", "--graph", help="graph to calculate statistics", required=True)
    parser.add_argument("-o", "--output", help="where to put the statistics", required=True)
    parser.add_argument("-r", "--ref", help="reference path", required=True, nargs='+')
    # TODO: need to work on this, to generalize, now using regex from contig name 
    parser.add_argument("-n", "--nonref", help="non-reference path sample", nargs='+')
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
        nodeid = value.nodeid
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


def small_node_calc(node_cum):
    total_node = 0
    len_node = 0
    ref_node = 0
    ref_node_len = 0
    non_ref_node = 0
    non_ref_node_len = 0

    for key, value in node_cum.items():
        if value.nodelen < 50:
            nodeid = value.nodeid
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

def core_flexible_calc(node_count,breed_count,node_cum):
    node_tally = defaultdict(int)
    nodelen_tally = defaultdict(int)
    # max_path = 0
    # node count is node id, no path
    for key,values in node_count.items():
        node_tally[int(values)] += 1
        nodelen = node_cum[key].nodelen
        nodelen_tally[int(values)] += nodelen
        # if int(values) > max_path:
        #     max_path = int(values)

    # for node count
    core_genome = node_tally[breed_count]
    private_genome = node_tally[1]
    flexible_genome = sum([value for key,value in node_tally.items() if key > 1 and key < breed_count ])
    all_node = sum(node_tally.values())

    # for node len
    core_len = nodelen_tally[breed_count]
    private_len = nodelen_tally[1]
    flexible_len = sum([value for key,value in nodelen_tally.items() if key > 1 and key < breed_count ])
    all_len = sum(nodelen_tally.values())
    

    return core_genome, private_genome, flexible_genome, all_node, core_len, private_len, flexible_len, all_len

@dataclass
class node_stats:
    nodeid: str
    nodelen: int
    rtype: int
    nodeseq: str


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
            # this is little hack
            # TODO: formalize ref name
            #ref_cac = str(comp) + "_UCD." + str(comp)
            ref_cac = str(comp) + "." +  str(comp)
            print(ref_cac)
            ref.append(ref_cac)

    output_file = args.output

    outfile = open(output_file, "w")

    node_cum = dict()
    #node_stats = namedtuple("node_stats", ["id", "nodelen", "rtype"])

    edge_cum = list()
    node_count = defaultdict(int)
    all_breed = set()

    if graph[-2:] == "gz":
        input_file = gzip.open(graph, "rt")
    else:
        input_file = open(graph, "r")

    for line in input_file:
        line_comp = line.strip().split()
        if line.startswith("S"):
            node_id = line_comp[1]
            node_len = len(line_comp[2])
            node_seq = line_comp[2]
            # rtype 0 as nonref
            node_cum[node_id] = node_stats(nodeid=node_id, nodelen=node_len, nodeseq=node_seq, rtype=0)
            # in minigraph we use the rank to determine ref path
            # the minigraph has been converted as compatible pggb graph
            #if grtype == "minigraph":
               # node_cum[node_id].rtype = 1 if int(line_comp[-1].split(":")[-1]) == 0 else 0
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

            ## core-flexible analysis
            ## need to deal with the non-ref path 
            # use this regex
            # re.findall(pattern=r"[A-Za-z]+",string="1_Brahman")
            path_info = line_comp[1] 
            # add breed info 
            node_breed = re.findall(pattern=r"[A-Za-z]+",string=path_info)[0]
            all_breed.add(node_breed)
            if not re.search(pattern=re.compile('minigraph',re.IGNORECASE),string=path_info):
                for comp in line_comp[2].split(","):
                     node_comp = comp[:-1]
                     node_count[node_comp] += 1

    input_file.close()

    for key,values in node_cum.items():
        if not values.rtype:
            if values.nodelen >= 50:
                print(f">{values.nodeid}")
                print(values.nodeseq)


    # total_node, len_node, ref_node, ref_node_len, non_ref_node, non_ref_node_len = node_calc(node_stats)
    # print(graph, grtype, "nodes", *node_calc(node_cum), file=outfile)
    # print(graph, grtype, "edges", *edge_calc(edge_cum, node_cum), file=outfile)
    # print(graph, grtype, "small_nodes", *small_node_calc(node_cum), file=outfile)
    # print(graph, grtype, "core-flexible_analysis", *core_flexible_calc(node_count, len(all_breed), node_cum), file=outfile)

    # make the output become more informative
    node_info = ["total_nodes", "len_nodes", "ref_nodes", "ref_nodes_len", "non_ref_nodes", "non_ref_nodes_len"]
    edge_info = ["total_edges", "Ref-Ref_edges", "Ref-Nonref_edges", "Nonref-Nonref_edges"]
    flexible_info = ["core_nodes", "private_nodes", "flexible_nodes",
                     "all_nodes", "core_len", "private_len", "flexible_len", "all_len"]
    node_info_small = [x + "_small" for x in node_info]

    for nodeinf, nodestat in zip(node_info, node_calc(node_cum)):
        print(nodeinf, nodestat, file=outfile)

    for edgeinf, edgestat in zip(edge_info, edge_calc(edge_cum,node_cum)):
        print(edgeinf, edgestat, file=outfile)
    
    # work for core-flexible for minigraph
    #for flexcomp, flexstat in zip(flexible_info, core_flexible_calc(node_count, len(all_breed), node_cum)):
        #print(flexcomp, flexstat, file=outfile)

    for node_info_small, nodestat in zip(node_info_small, small_node_calc(node_cum)):
        print(node_info_small, nodestat, file=outfile)

    outfile.close()
