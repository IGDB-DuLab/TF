# -*- coding: utf-8 -*-
"""
Created on Wed May 27 10:51:32 2020

@author: ZZG
"""

import os
import pandas as pd
import umap
import umap.plot

# UMAP

'''
from : https://umap-learn.readthedocs.io/en/latest/reproducibility.html
UMAP is a stochastic algorithm – 
it makes use of randomness both to speed up approximation steps, and to aid in solving hard optimization problems. 
This means that different runs of UMAP can produce different results. UMAP is relatively stable – 
thus the variance between runs should ideally be relatively small – but different runs may have variations none the less. 
To ensure that results can be reproduced exactly UMAP allows the user to set a random seed state. 

'''

file_path = r'./data'
# select TF ( expression cell percent great than zero and less than 0.9 )
input_matrix = pd.read_csv(os.path.join(file_path, r'UMAP_matrix.txt'), sep='\t', index_col=(0,1) )
label = input_matrix.index.get_level_values(1)


# MAP parameter
i = 10
j = 0.8

mapper = umap.UMAP(n_neighbors=i,min_dist=j, n_components=2, metric='euclidean').fit(input_matrix.values)
u = mapper.fit_transform(input_matrix)

# output UMAP position
u_pd= pd.DataFrame(u)
u_pd.columns = ['UMAP_x','UMAP_y']
u_pd.index = input_matrix.index

# output 
# u_pd : two dimension position of UMAP




# louvain


from sklearn.neighbors import kneighbors_graph
import community as community_louvain
import networkx as nx

# parameter
n_neighbors = 5
graph_k = pd.DataFrame(kneighbors_graph(u_pd, n_neighbors=n_neighbors, mode='connectivity').toarray())
graph_k.columns = list(u_pd.index.get_level_values(0))
graph_k.index = list(u_pd.index.get_level_values(0))



G = nx.from_pandas_adjacency(graph_k)
partition = community_louvain.best_partition(G)


# output cluster
louvain_cluster = []
for key, value in partition.items():
    louvain_cluster.append([key, value])
louvain_cluster_pd = pd.DataFrame(louvain_cluster)
louvain_cluster_pd.columns = ['cell', 'cluster']

# output 
# louvain_cluster_pd : sub-cluster

