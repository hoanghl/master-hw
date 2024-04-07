import sys

import numpy as np


def nb(X, labels):
    """
    Computes the weight vector w and and bias b corresponding
    to the Naive Bayes classifier with Bernoulli components.

    Parameters
    ----------
    X : an array of size (n, k)
            training input data for the classifier, elements must be 0/1
    labels : an array of size n
            training labels for the classifier, elements must be 0/1

    Returns
    -------
    w : an array of size k
            weights corresponding to the classifier
    bias: real number
            bias term corresponding to the classifier
    """

    cnt, k = X.shape
    w = np.zeros(k)
    b = 0

    # Estimate the parameters in NB
    theta1 = {j: 0 for j in range(k)}
    theta2 = {j: 0 for j in range(k)}
    pi = 0

    for i in range(cnt):
        # Count pi
        if labels[i] == 1:
            pi += 1

        # Count theta
        theta = theta1 if labels[i] == 1 else theta2
        for j in range(k):
            if X[i, j] == 1:
                theta[j] += 1
    for j in range(k):
        theta1[j], theta2[j] = theta1[j] / pi, theta2[j] / (cnt - pi)
    pi = pi / cnt

    # Estimate w and b
    s = 0
    for j in range(k):
        w[j] = (
            np.log(theta1[j])
            - np.log(1 - theta1[j])
            + np.log(1 - theta2[j])
            - np.log(theta2[j])
        )

        s += np.log(1 - theta1[j]) - np.log(1 - theta2[j])
    b = np.log(pi) - np.log(1 - pi) + s

    return w, b


def main(argv):
    D = np.loadtxt(argv[1])
    X = D[:, 1:]
    labels = D[:, 0]
    print(nb(X, labels))


# This allows the script to be used as a module and as a standalone program
if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("usage: python %s filename" % sys.argv[0])
    else:
        main(sys.argv)
