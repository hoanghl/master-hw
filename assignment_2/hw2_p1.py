################################################################################
# DATA16001: Network Analysis (2024)
# Homework 2
# Boilerplate code for Exercise 1
# Last Updated: January 15, 2024
################################################################################

from itertools import combinations

import networkx as nx
import matplotlib.pyplot as plt
from pathlib import Path

def get_strong_weak_links(G, threshold):
    """
    A function that returns the list of weak and strong links in the graph based on the threshold of edge weight
    
    Args:
        G : NetworkX graph object
        threshold (float) : weight threshold for strong ties
    """
    ###############################################################################
    # TODO: your code here
    # FIXME: HoangLe [Jan-24]: Done
    ###############################################################################
    
    # Hint: get weights using nx.get_edge_attributes 
    strong_links = []
    weak_links = []

    for edge, w in nx.get_edge_attributes(G, 'weight').items():
        if w >= threshold:
            strong_links.append(edge)
        else:
            weak_links.append(edge)
    return weak_links, strong_links



def draw_network_with_tie_strength(G, weak_links, strong_links):
    """
    Plots the network and links as described in exercise handout
    Args:
        G : NetworkX graph object 
        weak_links : list of tuples (u,v) of weak edges 
        strong_links : list of tuples (u,v) of strong edges

    """
    ############################################################################
    # TODO: Your code here!
    # FIXME: HoangLe [Jan-24]: Done
    ############################################################################
    # get layout
    layout = nx.circular_layout(G)

    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_subplot(111)

    # draw nodes
    nx.draw_networkx_nodes(G, pos=layout, node_size=2000, node_color='#d78cf5', ax=ax)

    # draw node labels 
    nx.draw_networkx_labels(G, pos=layout, font_size=14, ax=ax)

    # draw strong links
    nx.draw_networkx_edges(G, pos=layout, edgelist=strong_links, width=2, style='solid', edge_color='red')

    # draw weak links 
    nx.draw_networkx_edges(G, pos=layout, edgelist=weak_links, width=2, style='dashed', edge_color='blue')

    plt.title("Barn Swallow Contact Network")
    plt.show()

def check_stcp(G,weak_links,strong_links):
    """
    Checks the Strong Triadic Closure Property of all nodes in the graph


    Args:
        G : NetworkX grpah object
        weak_links : list of tuples (u,v) of weak edges
        strong_links : list of tuples (u,v) of strong edges

    Returns:
        dict{node:boolean}: a dict of nodes with True or False for satisfying STPC  
    """

    ############################################################################
    # TODO: Your code here!
    ############################################################################
    stcp_validity = {}
    # Note: G is undirected, so edge (u,v) is same as (v,u).
    # Hint: use itertools.combinations to get all pairs of two edges of nodes 

    def are_nodes_connected(g: nx.Graph, u, v):
        return v in g.neighbors(u)

    for node in G.nodes:
        neighbors_strong = [
            u for u in G.neighbors(node)
            if ((node, u) in strong_links) or ((u, node) in strong_links)
        ]

        if len(neighbors_strong) < 2:
            is_violated = False
        else:
            is_violated = False
            for (u, v) in combinations(neighbors_strong, 2):
                if not are_nodes_connected(G, u, v):
                    is_violated = True

                    break

        stcp_validity[node] = not is_violated

    return stcp_validity 


if __name__ == '__main__':
    network_file = Path(__file__).parent / "aves-barn-swallow-contact-network.edges"
    # keep the network file in a subfolder called data
    if not network_file.exists():
        raise FileNotFoundError(f"Cannot find network file at {network_file.resolve()}. Please check path.")

    aves = nx.read_weighted_edgelist(network_file,nodetype=int)

    # Choose threshold from exercise 
    weak_links, strong_links = get_strong_weak_links(aves,threshold=3)

    draw_network_with_tie_strength(aves,weak_links,strong_links)

    stcp_validity = check_stcp(aves,weak_links,strong_links)
    if sum(stcp_validity.values()) == len(stcp_validity):
        print("All nodes satisfy STCP property")
    else:
        stcp_violaters = [n for n,validity in stcp_validity.items() if not validity]
        print(f"Nodes {stcp_violaters} violate Strong Triadic Closure Property")
