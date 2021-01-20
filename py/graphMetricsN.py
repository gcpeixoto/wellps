#!/usr/bin/env python3
'''
    Call package Networkx to compute betweeness centrality from
    edge list passed by WELLPS. 

    This program is not compiled to be directly executable because
    is not necessary yet. However. This should be improved in future 
    version.

    From MATLAB, IT MUST be executed with Python >= 3.5
    
        !python3 py/betwByNetworkx.py edge_file outfile

    @gustavo
'''

import os, sys

# get Anaconda Python site-packages dir
if 'CONDA_PREFIX' in os.environ.keys():
    sys.path.append(os.path.join(os.environ['CONDA_PREFIX'],'lib/python3.7/site-packages'))
else:
    sys.path.append(os.path.join(os.environ['HOME'],'anaconda3/lib/python3.7/site-packages'))

    
import networkx as nx

# create graph
G = nx.Graph()

# read edge file 
with open(sys.argv[1]) as f:
    e = [tuple(map(int, i.split(' '))) for i in f]

# fill graph with edges
G.add_edges_from(e)
    
# compute degree centrality
deg = nx.degree_centrality(G)
print('---> WARNING: Networkx computes degree centrality as a fraction.\n \
                          Do not expect to have integer values as in SNAP!\n')

# compute betweeness centrality
bet = nx.betweenness_centrality(G,normalized=False,endpoints=True)

# compute closeness centrality
cln = nx.closeness_centrality(G)


# save to file
with open(sys.argv[2], 'w') as fw:
    fw.write('#Node\tBetweeness\tCloseness\tDegree\n')    
    for n in bet.keys():
        fw.write(f'{n}\t{bet[n]}\t{cln[n]}\t{deg[n]}\n')  
