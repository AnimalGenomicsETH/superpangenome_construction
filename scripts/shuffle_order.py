#!/usr/bin/env python

import random
import argparse

def parse_args():
    parser = argparse.ArgumentParser(description = __doc__,
                    formatter_class = argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-s","--shuffle")
    parser.add_argument("-p","--putend")
    parser.add_argument("-b","--backbone")
    return parser.parse_args()

def random_shuffle_breeds(breed_list,shuffle_size=10):
    """

    Randomly shuffled breeds input while maintaining backbone genome
    """

    breed_to_shuffle=breed_list[1:]
    shuffled_breeds = []

    for i in range(0,shuffle_size):
        random.shuffle(breed_to_shuffle)
        shuffled_breeds.append([breed_list[0],*breed_to_shuffle])

    return shuffled_breeds


def put_assemb_end(breed_list):
    """

    Put assembly at the end of order
    """

    breed_to_shuffle=breed_list[1:]
    breed_comb=[]
    for x in range(0,len(breed_to_shuffle)):
        breed_list_end = [breed_list[0]]
        end_breed=breed_to_shuffle[x]
        for comp in breed_to_shuffle:
            if comp != end_breed:
                breed_list_end.append(comp)
        breed_list_end.append(end_breed)
        breed_comb.append(breed_list_end)
    return breed_comb

def random_backbone(mash_distance):
    """

    Given mash distance, randomly use backbone and 
    order the assemblies according to the distance to the 
    choice of backbone
    """
    import pandas as pd
    
    breed_df = pd.read_csv(mash_distance,sep="\t",index_col=False,names=["breed1","breed2","reldist","absdist"])
    breed_df["breed1"] = breed_df.breed1.str.extract(r"([A-Z][a-zA-Z]+[^_])")
    breed_df["breed2"] = breed_df.breed2.str.extract(r"([A-Z][a-zA-Z]+[^_])")
    breed_list=breed_df.breed1.unique()
    breed_order = []
    for breed in breed_list:
        extract_df = breed_df[breed_df.breed2==breed].sort_values("reldist")
        breed_order.append(extract_df.breed1.tolist())
    return(breed_order)


if __name__ == "__main__":
    args=parse_args()
    
    if args.shuffle:
        breed_list=open(args.shuffle).readlines()[0].strip().split()
        for comp in random_shuffle_breeds(breed_list):
            print(*comp)
    if args.putend:
        breed_list=open(args.putend).readlines()[0].strip().split()
        for comp in put_assemb_end(breed_list):
            print(*comp)
    if args.backbone:
        for comp in random_backbone(args.backbone):
            print(*comp)
