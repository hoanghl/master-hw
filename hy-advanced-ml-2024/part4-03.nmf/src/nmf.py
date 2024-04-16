import sys
import numpy as np


def nmf(X, W, H, itercnt):
    """ 
    NMF algorithm

    Parameters
    ----------
    X : an array of size (n, m)
        input data matrix
    W : an array of size (n, k)
        initial left factor, must be non-negative
    H : an array of size (k, m)
        initial right factor, must be non-negative
    itercnt : int
        number of iterations	

    Returns
    -------
    W : an array of size (n, k)
        final left factor, must be non-negative
    H : an array of size (k, m)
        final right factor, must be non-negative
    err: an array of size itercnt+1
        ith element correspongs to the error after the ith iteration. The 0th
        entry is the error with the initial weights.
    """

    err = np.zeros(itercnt + 1)

    # place your code here

    return W, H, err


def main(argv):
    np.random.seed(2022) # Forces random to be repeatable. Remove when random seed is desired. 

    X = np.loadtxt(argv[1])
    k = int(argv[2])
    itercnt = int(argv[3])

    n, m = X.shape

    W = np.random.random((n, k))
    H = np.random.random((k, m))

    W, H, err = nmf(X, W, H, itercnt)
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
    if len(sys.argv) != 4:
        print('usage: python %s filename number_of_factors iteration_count' % sys.argv[0])
    else:
        main(sys.argv)
