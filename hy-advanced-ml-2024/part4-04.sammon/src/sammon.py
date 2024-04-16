import sys
import numpy as np
from sklearn.decomposition import PCA
from sklearn.metrics import euclidean_distances


def sammon(D, P, eta, tau, itercnt):
    """
    Sammon projection

    Parameters
    ----------
    D : an array of size (n, n)
        distance matrix between n data points, must be symmetric matrix
    P : an array of size (n, k)
        initial coordinates for the n data points in k-dimensional space
    eta : real
        parameter regulating the gradient step size (deterministic step)
    tau : real
        parameter regulating the gradient step size (deterministic step)
    itercnt : int
        number of iterations

    Returns
    -------
    P : an array of size (n, k)
        final coordinates for the n data points
    """

    cnt = D.shape[0]

    # place your code here

    return P


def main(argv):
    X = np.loadtxt(argv[1])
    k = int(argv[2])
    eta = float(argv[3])
    tau = float(argv[4])
    itercnt = int(argv[5])

    D = euclidean_distances(X)

    pca = PCA(n_components=k)
    P = pca.fit_transform(X)
    print('PCA:')
    print(P)

    print('Sammon:')
    print(sammon(D, P, eta, tau, itercnt))



# This allows the script to be used as a module and as a standalone program
if __name__ == "__main__": 
    if len(sys.argv) != 6:
        print('usage: python %s filename proj_dim eta tau itercnt' % sys.argv[0])
    else:
        main(sys.argv)
