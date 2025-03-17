### import packages ###
import pandas as pd
import numpy as np
import scipy
import networkx as nx
import matplotlib as plt
import leidenalg as la
import igraph as ig
import os

### set paths and options ###
# Path to edgelist tsv.gz or tsv file

path_el = "Your_edgelist_file"

### load as pandas dataframe ###
dataframe_el = pd.read_csv(path_el,
                           sep='[;\t]',
                           engine='python',
                           compression="infer",
                           keep_default_na=False,
                           na_values='nan')

dataframe_el["Edge"] = dataframe_el["Source"]+"-"+dataframe_el["Target"]

# remove all entries with interaction==ptmod #
dataframe_el = dataframe_el[dataframe_el.Interaction != "ptmod"]

# reformat #
edgelist_graph = list(dataframe_el["Source"]+" "+dataframe_el["Target"])

# read as graph #
graph = nx.parse_edgelist(edgelist_graph, create_using=nx.DiGraph())

# convert to igraph #
ig_graph = ig.Graph.from_networkx(graph)

# node list #
nodes_list = np.asarray(ig_graph.vs["_nx_name"])

# find leiden partition #
leiden = la.find_partition(ig_graph,
                           partition_type=la.ModularityVertexPartition,
                           seed=0,
                           n_iterations=100)

# show sizes and numbers of clusters #
for part, i in zip(list(leiden), range(0,len(list(leiden)))):
    print(i, len(part))

# create dict to match nodes to their community #
partitions = dict()
for partition, nodes in enumerate(leiden):
    mapping = dict(zip(nodes_list[nodes], [partition]*len(nodes)))
    partitions.update(mapping)
    
# match edges to their communities or set them to -1 if they are between communities #
def community(source,target,partitions):
    community = partitions[source] if partitions[source] == partitions[target] else -1
    return community

dataframe_el["Community"] = dataframe_el.apply(lambda row: community(row["Source"],row["Target"],partitions),axis=1)

# show number of edges per community #
dataframe_el["Community"].value_counts()

# show number of interactions per cluster #
summary = dataframe_el.groupby(["Community","Interaction"]).size().unstack().replace(np.nan,0).astype(int)
summary["n_edges"] = summary.sum(axis=1)
summary["percent_known"] = (1-(summary["unknown"]/summary["n_edges"]))*100
