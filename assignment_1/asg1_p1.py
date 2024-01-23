# %%
################################################################################
# DATA16001: Network Analysis (2024)
# Homework 1
# Boilerplate code for Problem 1
# Last Updated: January 15, 2024
################################################################################

import random

import matplotlib.pyplot as plt
import networkx as nx
import pandas as pd
import numpy as np
import seaborn as sns

# %%


def generate_G_nm(n: int, m: int, seed: int = 42) -> nx.Graph:
    """Return an instance of Erdos-Renyi graph

    Args:
        n (int): number of nodes
        m (int): number of edges

    Returns:
        nx.Graph: NetworkX Graph object
    """
    ############################################################################
    # TODO: Your code here!
    # NOTE: You can use in-built NetworkX functions for generation
    # FIXME: HoangLe [Jan-18]: Done
    ############################################################################
    graph = nx.gnm_random_graph(n=n, m=m, seed=seed)

    return graph


# %%
def generate_B_nm(n: int, m: int, seed: int = 42) -> nx.Graph:
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
    # FIXME: HoangLe [Jan-19]: Done
    ############################################################################
    k = max(1, np.ceil(m/n).astype(np.int32))

    G = nx.Graph()

    for i in range(1, k+2):
        G.add_node(i)

    for i in range(2, k+2):
        G.add_edge(1, i)

    for i in range(k+2, n+1):
        deg = pd.Series(dict(G.degree))
        prob = (deg / sum(deg)).values
        subset_U = np.random.choice(deg.index, k, replace=False, p=prob)
        
        G.add_node(i)

        for u in subset_U:
            G.add_edge(i, u)

    while G.number_of_edges() > m:
        G.remove_edge(*random.choice(list(G.edges)))

    assert G.number_of_edges() == m
    assert G.number_of_nodes() == n
    return G


# %%
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
    # FIXME: HoangLe [Jan-18]: Done
    # Graph is stored in edgelist format.
    # Each line of the file has two numbers separated by a space and represents
    # an edge between the two nodes.
    ############################################################################
    filename = f"{name}.edgelist"

    G = nx.read_edgelist(filename)

    G.name = name
    return G


# %%


def compute_clustering_coefficient(graph: nx.Graph) -> float:
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
    # FIXME: HoangLe [Jan-18]: Done
    # You may ony use the in-built NetworkX method to verify your answer.
    ############################################################################
    def _get_cluster_coef(g: nx.Graph, node):
        neighbors = nx.neighbors(g, node)
        g_s = g.subgraph(neighbors).copy()
        g_s.remove_edges_from(nx.selfloop_edges(g_s))
        e_v = nx.number_of_edges(g_s)
        d_v = nx.degree(g, node)
        c_v = 2 * e_v / (d_v * (d_v - 1)) if d_v > 1 else 0

        return c_v
    
    g = graph.copy()
    g.remove_edges_from(nx.selfloop_edges(g))

    local_cluster_coef = []
    for node in g.nodes:
        c_v = _get_cluster_coef(g, node)

        local_cluster_coef.append(c_v)

    C = sum(local_cluster_coef) / len(local_cluster_coef)

    assert C == nx.average_clustering(graph)

    return C


# %%
def get_degree_distribution(graph: nx.Graph):
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
    # FIXME: HoangLe [Jan-18]: Done
    ############################################################################
    x_graph = sorted((d for n, d in graph.degree()))
    degree_counts = pd.Series(x_graph).value_counts()
    y_graph = [
        degree_counts[i] if i in degree_counts else 0
        for i in range(degree_counts.index.max() + 1)
    ]

    return x_graph, y_graph


# %%
def plot_degree_distributions(graph: nx.Graph, G_nm: nx.Graph, B_nm: nx.Graph):
    """
    Draw degree distribution plot of the real world graph and the two model instances
    """
    ############################################################################
    # TODO: Your code here!
    # NOTE: Use the `get_degree_distribution` function for each graph
    # FIXME: HoangLe [Jan-18]: Done
    ############################################################################
    def _get_deg_seq(g: nx.Graph):
        return sorted((d for n, d in g.degree()))

    fig = plt.figure(figsize=(15, 8))

    graphs_info = [
        {"graph": graph, "name": "real_world"},
        {"graph": G_nm, "name": "G"},
        {"graph": B_nm, "name": "B"},
    ]
    for idx, info in enumerate(graphs_info):
        ax = fig.add_subplot(1, 3, idx + 1)
        deg_seq = _get_deg_seq(info['graph'])
        deg_seq = pd.Series(deg_seq, dtype=np.int32)
        deg_seq = deg_seq[deg_seq > 0]

        sns.histplot(deg_seq, ax=ax)
        ax.set_title(info['name'])

        ax.set_yscale('log')
        ax.set_xscale('log')

    plt.suptitle(f"Different graphs inferred from {graph.name}")

    plt.show()


# %%
def print_clustering_coefficient(graph: nx.Graph, G_nm: nx.Graph, B_nm: nx.Graph):
    """
    Draw degree distribution plot of the real world graph and the two model instances
    """
    print(
        f"Clustering coefficient of {graph.name} = {compute_clustering_coefficient(graph):.4f}"
    )
    print(
        f"Clustering coefficient of G({n},{m}) = {compute_clustering_coefficient(G_nm):.4f}"
    )
    print(
        f"Clustering coefficient of B({n},{m}) = {compute_clustering_coefficient(B_nm):.4f}"
    )


# %%
def print_num_nodes_edges(graph: nx.Graph):
    """
    Print number of nodes and edges of a given graphs.
    """
    print(
        f"{graph.name} Graph: Number of Nodes = {graph.number_of_nodes()}, Number of Edges = {graph.number_of_edges()}"
    )


# %%
if __name__ == "__main__":
    # load real-world graph
    air_traffic_graph = load_real_world_network(name="air_traffic")
    n, m = air_traffic_graph.number_of_nodes(), air_traffic_graph.number_of_edges()
    # print number of edges
    print_num_nodes_edges(air_traffic_graph)

    # generate instance of G(n,m)
    air_traffic_G_nm = generate_G_nm(n, m)
    # generate instance of B(n,m)
    air_traffic_B_nm = generate_B_nm(n, m)
    # print average clustering coefficient
    print_clustering_coefficient(air_traffic_graph, air_traffic_G_nm, air_traffic_B_nm)
    # plot degree distribution
    plot_degree_distributions(air_traffic_graph, air_traffic_G_nm, air_traffic_B_nm)

    acad_collab_graph = load_real_world_network(name="academic_collaboration")
    n, m = acad_collab_graph.number_of_nodes(), acad_collab_graph.number_of_edges()
    print_num_nodes_edges(acad_collab_graph)

    # generate instance of G(n,m)
    acad_collab_G_nm = generate_G_nm(n, m)
    # generate instance of B(n,m)
    acad_collab_B_nm = generate_B_nm(n, m)
    # print average clustering coefficient
    print_clustering_coefficient(acad_collab_graph, acad_collab_G_nm, acad_collab_B_nm)
    # plot degree distribution
    plot_degree_distributions(acad_collab_graph, acad_collab_G_nm, acad_collab_B_nm)


# %%
