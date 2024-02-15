################################################################################
# DATA16001: Network Analysis (2023)
# Homework 4
# Boilerplate code for Exercise 2
# Last Updated: Feb 10, 2023
################################################################################
import pickle
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
import scipy
import seaborn as sns
from loguru import logger
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import accuracy_score
from sklearn.model_selection import train_test_split
from tqdm import trange

plt.style.use("seaborn-v0_8")
plt.rcParams.update({"font.size": 8})

path_figs = Path("figures")

logger.remove(0)
logger.add(sys.stderr, level="SUCCESS")


def load_dataset(dataset_number: int):
    """Load the matrices from the files

    Args:
        dataset_number (int): The dataset number

    Returns:
        A,X,Y: The adjacency, original features and labels for the dataset
    """
    ###############################################################################
    # TODO: your code here
    # FIXME: HoangLe [Feb-13]: Done
    ###############################################################################
    path_adj = f"hw4_p2_data/graph_{dataset_number}_A.pkl"
    path_feat = f"hw4_p2_data/graph_{dataset_number}_X.pkl"
    path_tgt = f"hw4_p2_data/graph_{dataset_number}_Y.pkl"

    with open(path_adj, "rb") as f_adj, open(path_feat, "rb") as f_feat, open(
        path_tgt, "rb"
    ) as f_tgt:
        A, X, Y = pickle.load(f_adj), pickle.load(f_feat), pickle.load(f_tgt)

    return A, X, Y


def compute_laplacian_embeddings(A: np.ndarray):
    """Computes the Laplacian eigenvector embeddings of a given adjacency matrix.

    Args:
        A (np.ndarray): The adjacency matrix of a graph

    Returns:
        np.ndarray: The 2nd and 3rd smallest magnitude eigenvectors
    """
    ###############################################################################
    # TODO: your code here
    # Hint: Use scipy.sparse.linalg.eigsh
    # FIXME: HoangLe [Feb-15]: Done
    ###############################################################################
    D = np.diag(A.sum(axis=1))
    L = (D - A).astype(np.float32)

    _, eigenvec = scipy.sparse.linalg.eigsh(L, k=3, which="SM")

    return eigenvec[:, 1:]


def plot_scatter(
    A: np.ndarray, X: np.ndarray, Y: np.ndarray, dataset_number: int
) -> None:
    """Make the scatter plots

    Args:
        A (np.ndarray): The adjacency matrix
        X (np.ndarray): The original feature matrix
        Y (np.ndarray): The label vector
        dataset_number (int): the dataset number
    """
    ###############################################################################
    # TODO: your code here
    # Hint: Use the compute_laplacian_embeddings to get the embeddings
    # FIXME: HoangLe [Feb-15]: Done
    ###############################################################################
    eigenvec = compute_laplacian_embeddings(A)
    df_coordinates = pd.DataFrame(
        {
            "point": range(1, len(X) + 1),
            "x_coor_eigen": eigenvec[:, 0],
            "y_coor_eigen": eigenvec[:, 1],
            "x_coor_feat": X[:, 0],
            "y_coor_feat": X[:, 1],
            "tgt": Y,
        }
    )

    fig = plt.figure(figsize=(10, 5))

    ax = fig.add_subplot(121)
    sns.scatterplot(df_coordinates, x="x_coor_feat", y="y_coor_feat", hue="tgt", ax=ax)
    ax.set_title("Coordinates from node features")

    ax = fig.add_subplot(122)
    sns.scatterplot(
        df_coordinates, x="x_coor_eigen", y="y_coor_eigen", hue="tgt", ax=ax
    )
    ax.set_title("Coordinates from 2nd, 3rd eigenvalues")

    fig.suptitle(f"Plots for dataset number {dataset_number}", fontsize=16)

    path_fig = path_figs / f"{dataset_number}" / f"hw4_p2_scatter_{dataset_number}.png"
    fig.savefig(path_fig, bbox_inches="tight")


def plot_layouts(A: np.ndarray, dataset_number: int) -> None:
    """Make the NetworkX layout plots

    Args:
        A (np.ndarray): The adjacency matrix
        dataset_number (int): The dataset number
    """
    ###############################################################################
    # TODO: your code here
    # Hint: Use the networkx layout functions
    # FIXME: HoangLe [Feb-15]: Done
    ###############################################################################
    G = nx.from_numpy_array(A)

    fig = plt.figure(figsize=(10, 5))

    ax = fig.add_subplot(121)
    pos = nx.spectral_layout(G)
    nx.draw(G, pos=pos)
    nx.draw(G, pos=pos, node_color="violet", ax=ax)
    ax.set_title("Graph with spectral layout")

    ax = fig.add_subplot(122)
    pos = nx.spring_layout(G)
    nx.draw(G, pos=pos)
    nx.draw(G, pos=pos, node_color="violet", ax=ax)
    ax.set_title("Graph with spring layout")

    fig.suptitle(
        f"Graph with different layouts of dataset number {dataset_number}", fontsize=16
    )

    path_fig = path_figs / f"{dataset_number}" / f"hw4_p2_layout_{dataset_number}.png"
    fig.savefig(path_fig, bbox_inches="tight")


def fit_predict(
    Xtrain_: pd.DataFrame,
    ytrain: pd.DataFrame,
    Xtest_: pd.DataFrame,
    ytest: pd.DataFrame,
) -> float:
    """Train and predict with test

    Args:
        Xtrain_ (pd.DataFrame): training feature
        ytrain (pd.DataFrame): training target
        Xtest_ (pd.DataFrame): testing feature
        ytest (pd.DataFrame): testing target

    Returns:
        float: accuracy on test set
    """
    ###############################################################################
    # TODO: your code here
    # Hint: Use the scikit learn library
    # Note: Split the features and labels into 40% for training and 60% for testing
    # and keep the splits consistent. For example, using a seed
    # FIXME: HoangLe [Feb-15]: Done
    ###############################################################################
    Xtrain, Xtest = Xtrain_.to_numpy(), Xtest_.to_numpy()
    if len(Xtrain.shape) == 1:
        Xtrain = Xtrain.reshape(-1, 1)
        Xtest = Xtest.reshape(-1, 1)

    # fit model using training features and labels
    clf = RandomForestClassifier()
    clf.fit(Xtrain, ytrain)

    # Compute accuracy using test labels
    acc = accuracy_score(ytest, clf.predict(Xtest))

    return acc


def plot_accs(
    A: np.ndarray,
    X: np.ndarray,
    Y: np.ndarray,
    dataset_number: int,
    train_test_ratio: float = 0.4,
):
    """Make the accuracy bar plots as per assignment

    Args:
        A (np.ndarray): the adjacency matrix
        X (np.ndarray): the original feature matrix
        Y (np.ndarray): the label vector
        dataset_number (int): the dataset number
    """

    ###############################################################################
    # TODO: your code here
    # Note: Gather results for all feature matrices mentioned in assignment
    ###############################################################################
    # Get eigenvec

    logger.info("Get eigenvectors")

    eigenvec = compute_laplacian_embeddings(A)

    # Create dataframe of features
    features = pd.DataFrame(
        {
            "u2": eigenvec[:, 0],
            "u3": eigenvec[:, 1],
            "feat1": X[:, 0],
            "feat2": X[:, 1],
        }
    )

    # Split train-test
    logger.info("Split train-test")

    Xtrain, Xtest, ytrain, ytest = train_test_split(
        features, Y, test_size=train_test_ratio
    )

    # Start training/testing with different scenarios
    scenarios = {
        "u2": ["u2"],
        "u3": ["u3"],
        "u2,u3": ["u2", "u3"],
        "feat1,feat2": ["feat1", "feat2"],
        "u2,u3,feat1,feat2": ["u2", "u3", "feat1", "feat2"],
    }
    results = []

    for name, scenario in scenarios.items():
        logger.info(f"Train model: {name}")

        Xtrain_, Xtest_ = Xtrain[scenario], Xtest[scenario]
        acc = fit_predict(Xtrain_, ytrain, Xtest_, ytest)

        results.append({"feature": name, "accuracy": acc})

    # Plot
    logger.info("Plot barchart")

    df_result = pd.DataFrame.from_records(results)
    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot(111)
    sns.barplot(df_result, x="feature", y="accuracy", ax=ax)
    ax.set_title(
        f"Accuracy of RF models trained with different inputs of dataset {dataset_number}"
    )

    path_fig = path_figs / f"{dataset_number}" / f"hw4_p2_accuracy_{dataset_number}.png"
    fig.savefig(path_fig, bbox_inches="tight")


if __name__ == "__main__":
    for dataset_number in trange(1, 4, desc="Dataset"):
        path_fig = path_figs / f"{dataset_number}"
        path_fig.mkdir(exist_ok=True, parents=True)

        # load the matrices
        A, X, Y = load_dataset(dataset_number)
        # make the scatter plots
        plot_scatter(A, X, Y, dataset_number)
        # make the layout plots
        plot_layouts(A, dataset_number)
        # make the accuracy plots
        plot_accs(A, X, Y, dataset_number)
