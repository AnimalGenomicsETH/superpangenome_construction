import numpy as np
from scipy.cluster import hierarchy
from scipy.spatial.distance import squareform
import sys

names, vals = [], []

with open(sys.argv[1],'r') as fin:
    for i,line in enumerate(fin):
        parts = line.rstrip().split()
        names.append(parts[0].split('/')[1].split('.')[0])
        vals.append(parts[1:]+[0])

Q = np.asarray([np.pad(a, (0, len(vals) - len(a)), 'constant', constant_values=0) for a in vals],dtype=float)
Z = hierarchy.linkage(squareform((Q+Q.T)),method='average',optimal_ordering=True)

def get_newick(node, parent_dist, leaf_names, newick='') -> str:
    if node.is_leaf():
        return f'{leaf_names[node.id]}:{parent_dist - node.dist}{newick}'
    else:
        if len(newick) > 0:
            newick = f'):{parent_dist - node.dist}{newick}'
        else:
            newick = ');'
        newick = get_newick(node.get_left(), node.dist, leaf_names, newick=newick)
        newick = get_newick(node.get_right(), node.dist, leaf_names, newick=f',{newick}')
        newick = f'({newick}'
        return newick

tree = hierarchy.to_tree(Z, False)

print(get_newick(tree, tree.dist, names))
print('\n'.join([f'{N} {P}' for (N,P) in zip(names,sys.argv[2:])]))
