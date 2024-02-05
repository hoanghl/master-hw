################################################################################
# DATA16001: Network Analysis (2024)
# Homework 2
# Boilerplate code for Exercise 2
# Last Updated: January 15, 2024
################################################################################
from collections import Counter, defaultdict, deque
from itertools import combinations
from pathlib import Path

import matplotlib.pyplot as plt
import networkx as nx
import pandas as pd


def load_signed_network(path):
    """Reads from file at path and return graph

    Args:
        path : path location of the csv edge list file

    Returns:
        G : NetworkX Directed Graph object
    """
    ############################################################################
    # TODO: Your code here!
    # NOTE: Graph is directed and stored in edgelist format.
    # FIXME: HoangLe [Jan-29]: Done
    # Each line of the file has 4 numbers. For example
    # Source Target Rating Time
    # 1 2 4      1289241911.72836
    # The time attribute is NOT needed.
    # The rating attribute is needed and is the weight of the edge (u,v)
    ############################################################################

    # Hint: load csv into pandas first. Then use load dataframe into networkx
    # use edge_attr to choose columns for edge weights
    df = pd.read_csv(path)

    G = nx.from_pandas_edgelist(
        df,
        source="Source",
        target="Target",
        edge_attr="Rating",
        create_using=nx.DiGraph,
    )
    assert nx.is_directed(G)
    return G


def get_asymm_edges_diffs(G):
    """Find the asymmetric edges in a directed graph. Also find the absolute difference in rating in asymmetric edges.

    An asymmetric edge is defined as a pair of nodes u, v where w_{u,v} != w_{v,u}.

    Args:
        G : NetworkX Directed Graph

    Returns:
        asymmetric_edges : list of tuples (u,v) where w_{u,v} != w_{v,u}. Only keep one copy (u,v) for an asymmetric edge
        absolute_diffs : dict of counts of absolute difference of rating in asymmetric edges
    """
    ############################################################################
    # TODO: Your code here!
    # FIXME: HoangLe [Jan-29]: Partially done
    ############################################################################
    assert nx.is_directed(G)
    asymmetric_edges = []
    # Note: process an asymmetric edge (u,v) ONLY ONCE and store one entry
    # Note: both (u,v) and (v,u) need to be present in graph to qualify as asymmetric
    # Hint: keep absolute differences in a list and then use `collections.Counter` to get dict of counts
    weight_attr = "Rating"

    edges_asym = set()
    list_abs_diff = []
    for u, v in G.edges:
        if (
            (v, u) in G.edges
            and G.edges[(u, v)][weight_attr] != G.edges[(v, u)][weight_attr]
            and (v, u) not in edges_asym
        ):
            edges_asym.add((u, v))
            diff = abs(G.edges[(u, v)][weight_attr] - G.edges[(v, u)][weight_attr])
            list_abs_diff.append(diff)

    asymmetric_edges = list(edges_asym)
    absolute_diffs = dict(Counter(list_abs_diff))

    return asymmetric_edges, absolute_diffs


def plot_absolute_diffs(absolute_diffs):
    fig, ax = plt.subplots()
    rects = ax.bar(absolute_diffs.keys(), absolute_diffs.values())
    ax.set_xticks(list(absolute_diffs.keys()))
    ax.set_yscale("log")
    ax.set_xlabel("Absolute Difference in ratings")
    ax.set_ylabel("Number of asymmetrical edges")
    # Hide the right and top spines
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)

    def autolabel(rects, ax):
        """Attach a text label above each bar in *rects*, displaying its height."""
        for rect in rects:
            height = rect.get_height()
            ax.annotate(
                "{}".format(height),
                xy=(rect.get_x() + rect.get_width() / 2, height),
                xytext=(0, 3),  # 3 points vertical offset
                textcoords="offset points",
                ha="center",
                va="bottom",
            )

    autolabel(rects, ax)
    fig.tight_layout()
    plt.savefig("absolute_diffs.png", dpi=200)


def convert_undirected(G, asymmetric_edges):
    """Remove both directions of asymmetric edges and convert to undirected

    Args:
        G : NetworkX Directed Graph object
        asymmetric_edges : list of edges (u,v) where w_{u,v} != w_{v,u}

    Returns:
        G_und : NetowrkX Undirected Graph object
    """
    assert nx.is_directed(G)
    #######################################################
    # TODO: Your code here!
    # FIXME: HoangLe [Jan-29]: Done
    #######################################################
    # Note: remove both (u,v) and (v,u) edges from G
    # Then convert to undirected graph
    asymmetric_edges = set(asymmetric_edges)

    G_und = nx.Graph()
    for edge in G.edges:
        if edge in asymmetric_edges or (edge[1], edge[0]) in asymmetric_edges:
            continue
        weight = G.edges[edge]["Rating"]
        if weight > 0:
            weight = 1
        elif weight == 0:
            weight = 0
        else:
            weight = -1

        G_und.add_edge(*edge, weight=weight)
    assert not nx.is_directed(G_und)
    return G_und


def get_supernode_subgraphs(G):
    """A function that returns the supernodes.
    The supernodes are the connected components of the original graphs only considering positive edges between nodes.


    Args:
        G : NetworkX Undirected Graph

    Returns:
        supernodes : list of subgraphs induced by nodes in the supernodes]
        negative_edges : set of (u,v) for each negative edge in G
    """

    #######################################################
    # TODO: Your code here!
    # FIXME: HoangLe [Jan-26]: Done
    #######################################################

    # get negative edges in the graph
    negative_edges = set()

    for edge, w in nx.get_edge_attributes(G, "weight").items():
        if w < 0:
            negative_edges.add(edge)

    # Removing negative edges from the copy
    positive_graph = G.copy()
    positive_graph.remove_edges_from(negative_edges)

    # supernodes are the connected components in the positive graph with edges from the original network
    supernodes = [G.subgraph(c) for c in nx.connected_components(positive_graph)]

    # Make sure return subgraphs and not just nodes of connected components
    assert all(isinstance(obj, nx.Graph) for obj in supernodes)
    return supernodes, negative_edges


def get_reduced_graph(supernodes, negative_edges):
    """Method creates a reduced graph from supernodes and negative edges from original graph

    Args:
        supernodes : list of supernodes subgraphs
        negative_edges : set of (u,v) of negative edges

    Returns:
        reduced_graph: graph where nodes correspond to supernodes and edges connecting supernodes
    """
    supernode_labels = ["s" + str(i) for i in range(len(supernodes))]

    mapping = {}
    for i, s in enumerate(supernodes):
        mapping.update(dict.fromkeys(s.nodes, supernode_labels[i]))

    # remap negative edges (u,v) in original graph to (s_u,s_v)
    # s_u, s_v are labels for supernodes that u and v belong to respectively

    # create a set so that duplicate edges are counted only once
    remapped_negative_edges = set(
        map(lambda e: tuple(map(mapping.get, e)), negative_edges)
    )
    # verify there aren't any negative edges inside the same supernode
    # for eg, (s0,s0)
    intra_negative_edges = set((u, v) for u, v in remapped_negative_edges if u == v)
    assert len(intra_negative_edges) == 0

    # create reduced graph
    reduced_graph = nx.empty_graph()
    # add supernodes labels
    reduced_graph.add_nodes_from(supernode_labels)
    # add supernode edges
    reduced_graph.add_edges_from(remapped_negative_edges)

    return reduced_graph


def bfs_check(G: nx.Graph, src):
    """Function that performs Breadth First Search on the given graph from the specified source.
    Returns None if Odd Length Cycle is found in the graph.
    Otherwise return levels of nodes.

    Args:
        G : NetworkX undirected graph Object
        src : Node to begin BFS

    Returns:
        level: A dict mapping each node to a level corresponding to traversal depth
    """

    visited = dict.fromkeys(G.nodes, False)
    level = dict.fromkeys(G.nodes, 0)
    visited[src] = True

    q = deque()
    q.append(src)

    #######################################################
    # TODO: Your code here!
    # FIXME: HoangLe [Jan-26]: Done
    #######################################################
    # Begin your BFS loop here
    # Note: update the level dict during traversal
    # Check for odd length cycles
    cur_level = 0
    while len(q) > 0:
        # 1. For each node being currently in queue, put to queue
        while len(q) > 0:
            node = q.popleft()
            if level[node] == cur_level:
                visited[node] = True

                for node_adj in G.adj[node].keys():
                    if visited[node_adj] is False:
                        q.append(node_adj)
                        level[node_adj] = cur_level + 1
            elif level[node] > cur_level:
                # If reaches here, it means all nodes of the current level have been looped
                q.appendleft(node)
                cur_level += 1
                break
            else:
                raise ValueError(
                    f"cur_level is greater than node_level = {level[node]}"
                )

        # 2. Check every pairs of nodes within the queue whether there is an edge between a pair
        for u, v in combinations(q, 2):
            if (u, v) in G.edges:
                # logger.info(f"Found Odd Length Cycle at: {(u, v)}")
                # Found Odd Length Cycle at edge (u, v)
                return None

    return level


def check_balance(G):
    """Method that checks the balance of the graph G using the supernode algorithm

    Args:
        G : NetworkX Undirected signed graph

    Returns:
        bipartite_mapping: dict of two keys corresponding to disjoin sets of the nodes of G if balanced. Returns None if unbalanced.
    """

    # balance only defined for directed graphs
    assert not nx.is_directed(G)

    # get supernodes of graph
    supernodes, negative_edges = get_supernode_subgraphs(G)

    # verify if supernodes are valid
    for s in supernodes:
        weights = set(w for _, _, w in s.edges(data="weight"))
        if -1 in weights:
            print("Supernode has internal negative edges")
            return None

    print(f"There are {len(supernodes)} supernodes and all are valid.")

    # obtain the reduced graph
    reduced_graph = get_reduced_graph(supernodes, negative_edges)

    levels = bfs_check(reduced_graph, "s0")

    if levels is None:
        print("Odd cycle found in reduced graph")
        return None

    # graph is balanced
    # get bipartitpe partition of the nodes of the original graph
    bipartite_mapping = defaultdict(list)
    for k, v in levels.items():
        # all supernode labels are s1,s2, ...,
        supernode_id = int(k[1:])
        nodes = list(supernodes[supernode_id].nodes)
        if v % 2 == 0:
            bipartite_mapping["X"].extend(nodes)
        else:
            bipartite_mapping["Y"].extend(nodes)

    return bipartite_mapping


if __name__ == "__main__":
    # keep the data files in a subfolder called data

    simple_file = Path(__file__).parent / "simple.edges"

    simple_graph = nx.read_weighted_edgelist(simple_file, nodetype=int)

    if check_balance(simple_graph) is None:
        print("Simple Graph is Unbalanced\n")
    else:
        raise ValueError("ERROR: Simple Graph should be unbalanced. Recheck code")

    # modify the edge (2,4) as positive and check for balance
    simple_graph_modified = simple_graph.copy()
    simple_graph_modified[2][4]["weight"] = 1

    bipartite_mapping = check_balance(simple_graph_modified)
    if bipartite_mapping is None:
        raise ValueError(
            "ERROR: Modified Simple Graph should be balanced. Recheck code"
        )
    else:
        print("Bipartite mapping of modified simple graph nodes")
        print("Group X : ", bipartite_mapping["X"])
        print("Group Y: ", bipartite_mapping["Y"])
        print()

    bitcoin_file = Path(__file__).parent / "soc-sign-bitcoinotc.csv"
    if not bitcoin_file.exists():
        raise FileNotFoundError(
            f"Cannot find network file at {bitcoin_file.resolve()}. Please check path"
        )
    bitcoin_directed = load_signed_network(bitcoin_file)
    asymmetric_edges, absolute_diffs = get_asymm_edges_diffs(bitcoin_directed)

    print(
        f"There are {len(asymmetric_edges)} pairs of nodes that have asymmetric edges"
    )
    plot_absolute_diffs(absolute_diffs)

    # convert directed network to undirected
    bitcoin_undirected = convert_undirected(bitcoin_directed, asymmetric_edges)

    if check_balance(bitcoin_undirected):
        print("Bitcoin Network is balanced")
    else:
        print("Bitcoin Network is not balanced")
