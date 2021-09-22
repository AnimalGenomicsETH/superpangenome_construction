#!/usr/bin/env python
"""
collect statistics from bam file

"""

import argparse
import pysam
import re



def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input",
                        help="input alignment")
    parser.add_argument("-o", "--output",
                        help="output stat")
    return(parser.parse_args())


if __name__ == "__main__":
    args = parse_args()
    input_file = args.input
    output_file = args.output
    output_stream = open(output_file, "w")
    samfile = pysam.AlignmentFile(input_file, "rb", threads=10)
    # dict to save the stat to count

    statcount = {"totmap": 0, "unmap": 0, "mapq0": 0, "mapq10": 0, "mapq60": 0,
                 "al99": 0, "alp": 0, "clip": 0, "totun": 0, "ppair": 0}

    for read in samfile:
        statcount["totmap"] += 1
        if read.is_unmapped:
            statcount["unmap"] += 1
        else:
            if read.mapping_quality == 0:
                statcount["mapq0"] += 1
            if read.mapping_quality > 10:
                statcount["mapq10"] += 1
            if read.mapping_quality == 60:
                statcount["mapq60"] += 1
            if read.is_proper_pair:
                statcount["ppair"] += 1
            edis = read.get_tag("NM")
            identity = (read.query_alignment_length - edis) / read.query_alignment_length
            if identity >= 0.99:
                statcount["al99"] += 1
            if re.search(r"S|H", read.cigarstring):
                statcount["clip"] += 1
            if not re.search(r"S|H", read.cigarstring) and identity == 1:
                statcount["alp"] += 1

    stat = ["totmap", "unmap", "mapq0", "mapq10", "mapq60", "al99",
            "alp", "clip", "ppair"]
    # map_linear/map_stat_{rep}.tsv

    rep = input_file.split("/")[-1].split("_")[-1].split(".")[0]
    print("linear", "linear", rep, " ".join(str(statcount[x]) for x in stat), file=output_stream)
    output_stream.close()
