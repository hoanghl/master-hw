#%%
################################################################################
# DATA16001: Network Analysis (2024)
# Homework 1
# Boilerplate code for Problem 1
# Last Updated: January 15, 2024
################################################################################

import numpy as np
import networkx as nx
import matplotlib.pyplot as plt

#%%


def generate_G_nm(n:int, m:int,seed:int=42) -> nx.Graph:
    """Return an instance of Erdos-Renyi graph

    Args:
        n (int): number of nodes
        m (int): number of edges

    Returns:
        nx.Graph: NetworkX Graph object
    """
    ############################################################################
    # TODO: Your code here!
    # Note: You can use in-built NetworkX functions for generation
    ############################################################################
    return graph

#%%
def generate_B_nm(n:int, m:int,seed:int=42) -> nx.Graph:
    """Return an instance of the B(n,m) graph

    Args:
        n (int): number of nodes
        m (int): number of edges

    Returns:
        nx.Graph: NetworkX Graph object
    """
    ############################################################################
    # TODO: Your code here!
    # NOTE: Add edges based on the pseudocode

    
    ############################################################################
    assert G.number_of_edges() == m
    assert G.numer_of_nodes() == n
    return G
#%%
def load_real_world_network(name="air_traffic") -> nx.Graph:
    """
    Read from file and return ``name`` graph.

    INPUT:
    - ``name`` -- name of real-world network

    OUTPUT:
    - NetworkX graph object

    """
    ############################################################################
    # TODO: Your code here!
    # NOTE: You will have to provide path to data file. 
    # Graph is stored in edgelist format. 
    # Each line of the file has two numbers separated by a space and represents 
    # an edge between the two nodes.
    ############################################################################
    G.name = name
    return G
#%%

def compute_clustering_coefficient(graph:nx.Graph)->float:
    """
    Compute average clustering coefficient of ``graph``.

    INPUT:
    - ``graph`` -- NetworkX graph object

    OUTPUT:
    - average clustering coefficient of ``graph`` (type: float)

    """
    ############################################################################
    # TODO: Your code here!
    # NOTE: You must implement your own method here. 
    # You may ony use the in-built NetworkX method to verify your answer.
    ############################################################################
    return C

#%%
def get_degree_distribution(graph:nx.Graph):
    """
    Return degree values and cumulative count of number of nodes for each degree value.

    INPUT:
    - ``graph`` -- NetworkX graph object

    OUTPUT:
    - ``x_graph`` -- degree values of nodes in ``graph`` in sorted order (type: list)
    - ``y_graph`` -- number of nodes of degree `d` for all degree values `d` (type: list)
    """
    ############################################################################
    # TODO: Your code here!
    ############################################################################
    return x_graph, y_graph

#%%
def plot_degree_distributions(graph:nx.Graph,G_nm:nx.Graph,B_nm:nx.Graph):
    """
    Draw degree distribution plot of the real world graph and the two model instances
    """
    ############################################################################
    # TODO: Your code here!
    # NOTE: Use the `get_degree_distribution` function for each graph
    ############################################################################

#%% 
def print_clustering_coefficient(graph:nx.Graph,G_nm:nx.Graph,B_nm:nx.Graph):
    """
    Draw degree distribution plot of the real world graph and the two model instances
    """
    print(f"Clustering coefficient of {graph.name} = {compute_clustering_coefficient(graph):.4f}")
    print(f"Clustering coefficient of G({n},{m}) = {compute_clustering_coefficient(G_nm):.4f}")
    print(f"Clustering coefficient of B({n},{m}) = {compute_clustering_coefficient(B_nm):.4f}")

#%%
def print_num_nodes_edges(graph:nx.Graph):
    """
    Print number of nodes and edges of a given graphs.
    """
    print(f"{graph.name} Graph: Number of Nodes = {graph.number_of_nodes()}, Number of Edges = {graph.number_of_edges()}")

#%%
if __name__ == '__main__':
    # load real-world graph
    air_traffic_graph = load_real_world_network(name="air_traffic")
    n,m = air_traffic_graph.number_of_nodes(),air_traffic_graph.number_of_edges()
    # print number of edges
    print_num_nodes_edges(air_traffic_graph)

    # generate instance of G(n,m)
    air_traffic_G_nm = generate_G_nm(n,m)
    # generate instance of B(n,m)
    air_traffic_B_nm = generate_B_nm(n,m)
    # print average clustering coefficient
    print_clustering_coefficient(air_traffic_graph,air_traffic_G_nm,air_traffic_B_nm)
    # plot degree distribution
    plot_degree_distributions(air_traffic_graph,air_traffic_G_nm,air_traffic_B_nm)


    acad_collab_graph = load_real_world_network(name="academic_collaboration")
    n,m = acad_collab_graph.number_of_nodes(),acad_collab_graph.number_of_edges()
    print_num_nodes_edges(acad_collab_graph)
    
    # generate instance of G(n,m)
    acad_collab_G_nm = generate_G_nm(n,m)
    # generate instance of B(n,m)
    acad_collab_B_nm = generate_B_nm(n,m)
    # print average clustering coefficient
    print_clustering_coefficient(acad_collab_graph,acad_collab_G_nm,acad_collab_B_nm)
    # plot degree distribution
    plot_degree_distributions(acad_collab_graph,acad_collab_G_nm,acad_collab_B_nm)



# %%
