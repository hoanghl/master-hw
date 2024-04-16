import sys
import numpy as np


def error(X, W, H, reg):
    """ 
    Regularized error

    Parameters
    ----------
    X : an array of size (n, m)
        input data matrix, misssing values are marked as NaN
    W : an array of size (n, k)
        left factor
    H : an array of size (k, m)
        right factor
    reg : real
        weight for regularization terms of W and H

    Returns
    -------
    err : real
        regularized error
    """

    err = 0
    # place your code here
    return err


def solve(X, W, reg):
    """
    Finds H solving X ~ WH with regularization

    Parameters
    ----------
    X : an array of size (n, m)
        input data matrix, misssing values are marked as NaN
    W : an array of size (n, k)
        left factor
    reg : real
        weight for regularization term

    Returns
    -------
    H : an array of size (k, m)
        right factor
    """

    m = X.shape[1]
    k = W.shape[1]
    H = np.zeros((k, m))

    # place your code here

    return H


def als(X, W, reg, itercnt):
    """
    Alternative least squares algorithm

    Parameters
    ----------
    X : an array of size (n, m)
        input data matrix, misssing values are marked as NaN
    W : an array of size (n, k)
        initial left factor
    reg : real
        weight for regularization terms of W and H
    itercnt : int
        number of iterations

    Returns
    -------
    W : an array of size (n, k)
        final left factor
    H : an array of size (k, m)
        final right factor
    err : an array of size itercnt
        ith element: error after i + 1 iterations. Single iteration = optimize
        H given W and optimize W given H.
    """

    err = np.zeros(itercnt)

    # place your code here

    return W, H, err



def main(argv):
    np.random.seed(2022) # Forces random to be repeatable. Remove when random seed is desired. 

    X = np.genfromtxt(argv[1])
    k = int(argv[2])
    reg = float(argv[3])
    itercnt = int(argv[4])

    n, m = X.shape

    W = np.random.random((n, k))

    W, H, err = als(X, W, reg, itercnt)
    print('W')
    print(W)
    print('H')
    print(H)
    print('errors')
    print(err)
    print('product matrix')
    print(W @ H)



# This allows the script to be used as a module and as a standalone program
if __name__ == "__main__": 
    if len(sys.argv) != 5:
        print('usage: python %s filename number_of_factors regularizer iteration_count' % sys.argv[0])
    else:
        main(sys.argv)
