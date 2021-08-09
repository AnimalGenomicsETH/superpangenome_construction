#!/usr/bin/env python
import argparse

def parse_args():
        parser = argparse.ArgumentParser(description = __doc__,
                                     formatter_class = argparse.RawDescriptionHelpFormatter)
        parser.add_argument("-i","--input",help="input assembly file")
        parser.add_argument("-s","--suffix",help="suffix to add in the contig")
        parser.add_argument("-o","--output",help="output file")
        return parser.parse_args()


if __name__ == "__main__":
    args=parse_args()
    assemb=args.input
    suffix=args.suffix
    outname=args.output

    with open(assemb) as infile, open(outname,"w") as outfile:
        for line in infile:
            token=line.strip().split()
            if line.startswith(">"):
                print(token[0]+"_"+suffix,file=outfile)
            else:
                outfile.write(line)
                     


        