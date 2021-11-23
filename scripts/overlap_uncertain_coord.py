#!/usr/bin/env python
"""

Given two interval files, overlap with tolerated x bp differences

"""
import sys

first_file=sys.argv[1]
second_file=sys.argv[2]
tolerated = int(sys.argv[3])
first_graph=sys.argv[4]
second_graph=sys.argv[5]


in1=open(first_file,"r")
in2=open(second_file,"r")

max_pos=0

in1.readline()
in2.readline()

first_comp=in1.readline().strip().split()
second_comp=in2.readline().strip().split()
first_coord=int(first_comp[1])
second_coord=int(second_comp[1])
coord_sel=first_coord
comb_coord=[]

while True:
    if first_coord <= second_coord:
        
        if abs(first_coord-second_coord) <= tolerated:
            #overlap.append([first_comp,second_comp])
            #print([first_comp,second_comp])
            if coord_sel == first_coord:
                #comb_coord.append(second_coord)
                comb_coord.extend([first_coord,first_graph])
            else:
                print(*comb_coord)
                coord_sel = first_coord
                comb_coord=[first_coord,first_graph,second_coord,second_graph]
            #print(first_coord,second_coord)

        try:
            first_comp =  in1.readline().strip().split()
            first_coord = int(first_comp[1])
        except:
            break
    
    if first_coord > second_coord:

        if abs(first_coord-second_coord) <= tolerated:
            #overlap.append([first_comp,second_comp])
            #print([first_comp,second_comp])
            #print(first_coord,second_coord)
            if coord_sel == first_coord:
                comb_coord.extend([second_coord,second_graph])
            else:
                print(*comb_coord)
                coord_sel = first_coord
                comb_coord=[first_coord,first_graph,second_coord,second_graph]

        try:
            second_comp = in2.readline().strip().split()
            second_coord = int(second_comp[1])
        except:
            break

in1.close()
in2.close()
      

      



