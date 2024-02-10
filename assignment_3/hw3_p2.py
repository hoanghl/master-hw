################################################################################
# DATA16001: Network Analysis (2024)
# Homework 2
# Boilerplate code for Exercise 2
# Last Updated: Jan 15, 2024
################################################################################

from pathlib import Path
from collections import defaultdict
import random

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy
import seaborn as sns


import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import sys
from typing import Dict, Any, Final, List

#! NOTE: INSTALL epydemic using `pip install epydemic`
from epydemic import CompartmentedModel, SynchronousDynamics, StochasticDynamics
import copy
import pandas as pd
import seaborn as sns

plt.style.use('seaborn-v0_8')
plt.rcParams.update({'font.size': 8})

class SIR_custom(CompartmentedModel):
    """The Susceptible-Infected-Removed model.
    Susceptible nodes are infected by infected neighbours, and recover to
    removed."""

    # Model parameters
    #: Parameter for probability of infection on contact.
    P_INFECT: Final[str] = "pInfect"
    #: Parameter for probability of removal (recovery).
    P_REMOVE: Final[str] = "pRecover"
    P_SEEDS: Final[str] = "pSeeds"  # : Parameter for initial infected nodes
    # Possible dynamics states of a node for SIR dynamics
    #: Compartment for nodes susceptible to infection.
    SUSCEPTIBLE: Final[str] = "S"
    INFECTED: Final[str] = "I"  #: Compartment for nodes infected.
    #: Compartment for nodes recovered/removed.
    REMOVED: Final[str] = "R"

    # Locus containing the edges at which dynamics can occur
    SI: Final[str] = "SI"  #: Edge able to transmit infection.

    def __init__(self):
        super().__init__()

    def build(self, params: Dict[str, Any]):
        """Build the SIR model.

        :param params: the model parameters"""
        super().build(params)

        # collect parameters which are passed to us in a dictionary
        pInfect = params[self.P_INFECT]
        pRemove = params[self.P_REMOVE]

        # create compartments or states
        self.addCompartment(self.SUSCEPTIBLE, 0.0)
        self.addCompartment(self.INFECTED, 0.0)
        self.addCompartment(self.REMOVED, 0.0)

        # create and track the loci where events will occur
        self.trackEdgesBetweenCompartments(
            self.SUSCEPTIBLE, self.INFECTED, name=self.SI
        )
        self.trackNodesInCompartment(self.INFECTED)

        # define the events that will occur
        self.addEventPerElement(self.SI, pInfect, self.infect)
        self.addEventPerElement(self.INFECTED, pRemove, self.remove)

    def setUp(self, params: Dict[str, Any]):
        # initialise all nodes to an empty compartment
        # (so we can assume all nodes have a compartment attribute)
        g = self.network()
        for n in g.nodes():
            g.nodes[n][self.COMPARTMENT] = None
        # mark edges as unoccupied
        for _, _, data in g.edges(data=True):
            data[self.OCCUPIED] = False
        # Go through all nodes
        for node in g.nodes():
            # if node is in seed set mark as INFECTED
            if node in params[self.P_SEEDS]:
                self.changeCompartment(node, SIR_custom.INFECTED)
            # Otherwise mark as SUSCEPTIBLE
            else:
                self.changeCompartment(node, SIR_custom.SUSCEPTIBLE)

    def infect(self, t: float, e: Any):
        """Perform an infection event. This changes the compartment of
        the susceptible-end node to :attr:`INFECTED`. It also marks the edge
        traversed as occupied.

        :param t: the simulation time
        :param e: the edge transmitting the infection, susceptible-infected"""
        (n, _) = e
        self.changeCompartment(n, self.INFECTED)
        self.markOccupied(e, t)

    def remove(self, t: float, n: Any):
        """Perform a removal event. This changes the compartment of
        the node to :attr:`REMOVED`.

        :param t: the simulation time (unused)
        :param n: the node"""
        self.changeCompartment(n, self.REMOVED)


def getGraph(file: Path):
    """Load the graph from a file

    Args:
        file (Path): Path to the edgelist file

    Returns:
        networxX.Graph: The graph
    """
    ###############################################################################
    # TODO: your code here
    # FIXME: HoangLe [Feb-06]: Done
    ###############################################################################
    G = nx.read_edgelist(file)
    return G


def SIRsimulation(G: nx.Graph, infected_seeds: List[int], p: float, r: float) -> Dict:
    """Runs one simulation of the SIR model on given graph
    and

    Args:
        G (nx.Graph): graph to run simulation
        infected_seeds (List[int]): Nodes to be initially infected
        p (float): The probability of infection spreading to neighbors
        r (float): The probability of recovery after infection

    Returns:
        Dict: results of nodes in different states of S,I and R
    """
    param = dict()
    param[SIR_custom.P_INFECT] = p  # infection probability
    param[SIR_custom.P_REMOVE] = r  # probability that node gets recovered
    # set the nodes to be initial infected seed nodes
    param[SIR_custom.P_SEEDS] = infected_seeds
    # create model
    m = SIR_custom()
    # create experiment
    e = StochasticDynamics(m, G)
    # set the experiment parameters and run the experiments
    rc = e.set(param).run()
    # gather the experiment results
    results = rc["results"]
    return results


def calculate_average_spread(
    candidate: int,
    seed_set: List[int],
    G: nx.Graph,
    p: float,
    r: float = 1,
    sigma: int = 5,
) -> int:
    """Function to compute the infection spread after adding
    candidate node to existing seed set

    Args:
        candidate (int): The node to add to seed set
        seed_set (List[int]): The existing set of initially infected nodes
        G (nx.Graph): The graph
        p (float): the probability of infection
        r (float): The probability of recovery. Defaults to 1.
        sigma (int): Number of simulations to get average spread. Defaults to 5.

    Returns:
        int: The spread from SIR experiment
    """
    ###############################################################################
    # TODO: your code here
    # FIXME: HoangLe [Feb-10]: Done
    # hint: use the SIRsimulation to run multiple simulations
    ###############################################################################
    new_seed = seed_set + [candidate]
    spreads = []
    for _ in range(sigma):
        result = SIRsimulation(G, new_seed, p, r)
        spreads.append(result['I'] + result['R'])

    spread = int(sum(spreads)*1.0 / len(spreads))

    return spread

def greedy(
    p: float,
    G: nx.Graph,
    candidate_set: List[int],
    k: int,
    r: float,
    verbose: bool = False,
    sigma: int = 5,
):
    """Greedy algorithm to select seed set of k
    with maximum spread from candidate set.

    Args:
        p (float): The probability of infection
        G (nx.Graph): The graph
        candidate_set (List[int]): The set of candidates for seed nodes
        k (int): The number of seed nodes required
        r (float): The probability of recovery
        verbose (bool): Whether to print debug statements
        sigma (int): Number of simulations to get average spread


    Returns:
        (List[int],List[Int]): The list of seeds and the spread when adding that candidate
    """
    assert len(candidate_set) > k
    spreads = []
    seeds = []
    ###############################################################################
    # TODO: your code here
    # FIXME: HoangLe [Feb-10]: Done
    # Implement the Greedy algorithm. Use the calculate_spread function
    # keep seeds and record average spread at each iteration
    ###############################################################################
    for _ in range(k):
        # Find node maximizing the expected spread
        max_spread = -10
        optimal_node = -10
        for c in candidate_set:
            spread = calculate_average_spread(c, seeds, G, p, r, sigma)
            if spread > max_spread:
                max_spread = spread
                optimal_node = c

        # Update seed set
        seeds.append(optimal_node)
        spreads.append(max_spread)

        # Discard optimal node from candidate seed
        candidate_set.remove(optimal_node)

    return seeds, spreads

def experiment(
    G: nx.Graph, 
    list_p: List[float],
    r: float,
    k: int, 
    num_init_cand: int = 7,
    num_cand_set: int = 3,
    sigma: int = 5,
    verbose: bool = False
) -> dict[float, list]:
    """Run experiment

    Args:
        G (nx.Graph): Graph to do experiment with SIR
        list_p (List[float]): List of infection probabilities
        r (float): Recovery probability
        k (int): Number of seeding nodes
        num_init_cand (int, optional): Number of candidates picked randomly from G. Defaults to 7.
        num_cand_set (int, optional): Number of candidate sets to do experiment. Defaults to 3.
        sigma (int, optional): Num. simulations to approximate expected spread. Defaults to 5.
        verbose (bool, optional): Flag for debugging purpose. Defaults to False.
    
    Returns:
        dict[float, list]: results
    """

    # Generate experiment configurations
    configurations = []
    for p in list_p:
        for _ in range(num_cand_set):
            candidates = random.sample(list(G.nodes), num_init_cand)
            configurations.append({
                'p' : p,
                'candidates': candidates
            })

    # Do experiments with 6 different configurations
    results = defaultdict(list)
    for conf in configurations:
        p, candidates = conf['p'], conf['candidates'] 
        _, spreads = greedy(p, G, candidates, k, r, verbose, sigma)

        results[p].append(spreads)

    return results

def visualize(results: dict):
    """Plot line charts

    Args:
        results (dict): results from experiments
    """

    # Convert to dataframe
    records = []

    for p, list_spreads in results.items():
        for nth_set, spreads in enumerate(list_spreads):
            for k, spread in enumerate(spreads):
                records.append({
                    'p': p,
                    'nth_set': nth_set,
                    'k': k+1,
                    'spread': spread
                })

    df = pd.DataFrame.from_records(records)

    # Plot
    fig = plt.figure(figsize=(10, 5))

    for idx, p in enumerate(df['p'].unique()):
        ax = fig.add_subplot(1, 2, idx + 1)
        df_p = df[df['p'] == p]
        sns.lineplot(df_p, x='k', y='spread', hue='nth_set', style="nth_set", markers=True, ax=ax)
        ax.set_title(f"p = {p}")

    fig.tight_layout()
    plt.savefig("as3_p2.png", dpi=200)

if __name__ == "__main__":
    path = Path("h3_graph_data.edgelist")
    G = getGraph(path)

    # Do experiments
    list_p = [0.1, 0.5]
    results = experiment(G, list_p, r=1, k=5, num_init_cand=7, num_cand_set=3, sigma=5)
    
    # Visualize
    visualize(results)

