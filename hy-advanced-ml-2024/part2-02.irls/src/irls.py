import sys

import numpy as np


def sigmoid(z):
    return 1 / (1 + np.exp(-z))


def get_misclass(y_hat, y):
    return np.mean(y_hat != y)


def nll(y, X, w):
    N = y.shape[0]

    result = 0
    for i in range(N):
        result += np.log(1 + np.exp(-y[i] * w.T @ X[i]))

    return result


def get_g(X, w, y):
    mu = sigmoid(X @ w)
    y_ = (y + 1) / 2

    g = X.T @ (mu - y_)

    return g


def get_inv_H(X, w):
    mu = sigmoid(X @ w)
    S = np.diag(mu * (1 - mu))
    H = X.T @ S @ X
    H_inv = np.linalg.inv(H)

    return H_inv


def pred(w, X):
    prob = X @ w
    y_pred = (sigmoid(prob) >= 0.5).astype(np.int32) * 2 - 1

    return y_pred


def irls(X, labels, w, itercnt):
    """
    IRLS algorithm

    Parameters
    ----------
    X : an array of size (n, k)
            training input data for the classifier
    labels : an array of size n
            training labels for the classifier, elements must be 0/1
    w : an array of size k
            initial weights
    itercnt : int
            number of iterations

    Returns
    -------
    w : an array of size k
            weights after itercnt iterations
    err: an array of size itercnt + 1
            ith element correspongs to the error (objective function minimized in
            logistic regression) after the ith iteration. The 0th entry is the
            error with the initial weights.
    misclass: an array of size itercnt + 1
            ith element correspongs to the misclassification proportion after the
            ith iteration. The 0th entry is the misclassification proportion with
            the initial weights.
    """

    err = np.zeros(itercnt + 1)
    misclass = np.zeros(itercnt + 1)

    y = labels * 2 - 1  # label 0 -> y = -1, label 1 -> y = 1

    err[0] = nll(y, X, w)
    misclass[0] = get_misclass(pred(w, X), y)

    for i in range(itercnt):
        w = w - get_inv_H(X, w) @ get_g(X, w, y)

        err[i + 1] = nll(y, X, w)
        y_hat = pred(w, X)
        misclass[i + 1] = get_misclass(y_hat, y)

    return w, err, misclass


def main(argv):
    D = np.loadtxt(argv[1])
    labels = D[
        :, 0
    ].copy()  # copy is needed, otherwise next line will mess up the splice
    D[
        :, 0
    ] = 1  # replace the label column of D with constant, now the first feature gives us the bias term

    itercnt = int(argv[2])

    w = np.zeros(D.shape[1])
    w, err, misclass = irls(D, labels, w, itercnt)

    print("weights:")
    print(w)
    print("error:")
    print(err)
    print("misclassification rate:")
    print(misclass)


# This allows the script to be used as a module and as a standalone program
if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("usage: python %s filename" % sys.argv[0])
    else:
        main(sys.argv)
