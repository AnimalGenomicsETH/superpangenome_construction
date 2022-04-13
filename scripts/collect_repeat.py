#!/usr/bin/env python

import argparse
import re

def extract_repeats(fname):
    elements = ('Total interspersed repeats', 'Small RNA', 'Satellites', 'Simple repeats', 'Low complexity','DNA transposons','Unclassified','Retroelements','LINEs','SINEs','LTR elements')
    repeats = dict()

    with open(fname,'r') as file_:
        for line in file_:
            if 'file name' in line:
                repeats['chrm'] = int(re.findall(r"[0-9]+",line)[0])
                # repeats['primary'] = str(repeats['chrm']) in chrms 
            elif 'total length' in line:
                repeats['length'] = int(line.split('(')[1].split()[0])
            elif 'GC level' in line:
                repeats['GC'] = float(line.split()[2]) * repeats['length']
            elif 'bases masked' in line:
                repeats['masked'] = int(line.split()[2])
                repeats['prop_masked'] = re.findall(r"[0-9]+\.[0-9]+",line)[0]
            for element in elements:
                if element in line:
                    repeats[element] = int(line.split('bp')[0].split()[-1])
    return repeats

def main():
    parser = argparse.ArgumentParser(description='masker.')
    parser.add_argument('--input', nargs="+", required=True)
    parser.add_argument('--prog', type=str, required=True)
    args = parser.parse_args()
    elements = ('Total interspersed repeats', 'Small RNA', 'Satellites', 'Simple repeats', 'Low complexity','DNA transposons','Unclassified','Retroelements','LINEs','SINEs','LTR elements')
    
    for fname in args.input:
        repeats=extract_repeats(fname)
        for key,value in repeats.items():
            if key in elements:
                print(args.prog,repeats["chrm"],key,value,repeats['length'],repeats['masked'],sep="\t")

if __name__ == '__main__':
    main()
