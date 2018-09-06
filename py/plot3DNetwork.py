#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  6 12:24:16 2018

@author: gustavo
"""

import sys
sys.path.append('/anaconda/lib/python3.6/site-packages/')

from scipy.io import loadmat

# igraph was best installed with: 
# conda install -c conda-forge python-igraph
import igraph as ig

import plotly.plotly as py
import plotly.graph_objs as go


f = '../mat/DRT_geometric_ln_13_Metrics.mat'
f2 = '../mat/DRT_geometric_ln_13.mat'
nc = 1 # component number
nc -= 1 # address

# load .mat file and make its handling MATLAB-like 
# see: https://docs.scipy.org/doc/scipy/reference/tutorial/io.html
S = loadmat(f,squeeze_me=True,struct_as_record=False)
S2 = loadmat(f2,squeeze_me=True,struct_as_record=False)
s = S['Maux'] 
s2 = S2['drtSt'] 

# get cluster adjacency matrix
adj = s.adjMatrix[nc]

# cast to int
adj = adj.astype(int)

# source and target nodes
srcNodes = adj[:,0].todense()
tgtNodes = adj[:,1].todense()
nnodes = adj.shape[0]

# edges
edges=[(srcNodes[k,0], tgtNodes[k,0]) for k in range(nnodes)]


# create graph object
G = ig.Graph(edges, directed=False)

coords = s2.compVoxelCoords[nc]

Xn = [coords[k][0] for k in range(nnodes)]# x-coordinates of nodes
Yn = [coords[k][1] for k in range(nnodes)]# y-coordinates
Zn = [coords[k][2] for k in range(nnodes)]# z-coordinates

Xe=[]
Ye=[]
Ze=[]
for e in edges:
    Xe+=[coords[e[0]][0],coords[e[1]][0], None]# x-coordinates of edge ends
    Ye+=[coords[e[0]][1],coords[e[1]][1], None]  
    Ze+=[coords[e[0]][2],coords[e[1]][2], None]  
    
    
trace1=go.Scatter3d(x=Xe,
               y=Ye,
               z=Ze,
               mode='lines',
               line=dict(color='rgb(125,125,125)', width=1),
               hoverinfo='none'
               )

trace2=go.Scatter3d(x=Xn,
               y=Yn,
               z=Zn,
               mode='markers',
               name='actors',
               marker=dict(symbol='circle',
                             size=6,
                             #color=group,
                             colorscale='Viridis',
                             line=dict(color='rgb(50,50,50)', width=0.5)
                             ),
               #text=labels,
               hoverinfo='text'
               )

axis=dict(showbackground=False,
          showline=False,
          zeroline=False,
          showgrid=False,
          showticklabels=False,
          title=''
          )

layout = go.Layout(
         title="test",
         width=1000,
         height=1000,
         showlegend=False,
         scene=dict(
             xaxis=dict(axis),
             yaxis=dict(axis),
             zaxis=dict(axis),
        ),
     margin=dict(
        t=100
    ),
    hovermode='closest' )
     
     
data=[trace1, trace2]
fig = go.Figure(data=data, layout=layout)

py.iplot(fig, filename='graph')
