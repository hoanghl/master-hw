import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import seaborn as sns
from tqdm import trange

plt.style.use("seaborn-v0_8")
plt.rcParams.update({"font.size": 8})

if __name__ == "__main__":
    # Plot the figure c-isolated probability
    step = 0.05
    n = 500
    cs = np.linspace(0, 10, num=int(10 / 0.05) + 1)
    ps = cs / (n - 1)
    fig = plt.figure(figsize=(5, 5))

    ax = fig.add_subplot(111)
    ax.plot(cs, (1 - ps) ** (n - 1))
    ax.set_xlabel("c")
    ax.set_ylabel("isolated probability")
    ax.set_title("Different c values and probability of node being isolated")

    plt.show()

    # Plot other
    c_connected = []
    for _ in trange(100):
        for c in cs:
            p = c / (n - 1)
            g = nx.erdos_renyi_graph(n, p)

            if nx.is_connected(g):
                # print(f"is_connected at c = {c}")
                c_connected.append(c)
                break
    c_connected = np.array(c_connected)
    p_connected = c_connected / (n - 1)

    fig = plt.figure(figsize=(8, 4))
    ax = fig.add_subplot(111)
    sns.histplot(c_connected, ax=ax, bins=30)
    ax.set_title("Histogram of the smallest c values with which G(n, p) is connected")
    plt.show()

    fig = plt.figure(figsize=(8, 4))
    ax = fig.add_subplot(111)
    sns.histplot(p_connected, ax=ax, bins=30)
    ax.vlines([np.log(n) / n], ymin=0, ymax=8, linestyles="dashed", colors=["red"])
    ax.set_title("Histogram of the smallest p values with which G(n, p) is connected.")
    plt.show()
