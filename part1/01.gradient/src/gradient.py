import sys

import numpy as np


class Gradient:
    def __init__(self, mode, eta, tau, k):
        """
        Initializes gradient class.

        Parameters
        ----------
        mode : string
                mode for determining gradient size, should be either determ or adagrad
        eta : real
                parameter for step size, used by determ and adagrad
        tau : real
                parameter for step size, used by determ and adagrad
        k : real
                the dimension of variable space, used to initialize
                self.grads_squared, a sum of squared gradients, used by adagrad
        """

        self.mode = mode  # either 'determ' or 'adagrad'
        self.eta = eta
        self.tau = tau
        self.grads_squared = np.zeros(k)  # used by adagrad
        self.round = 0  # used by determ

    def step_size(self):
        """
        Computes current step size using the current members.
        determ uses eta, tau, and round while adagrad uses
        eta, tau, and grads_squared.
        """

        if self.mode == "determ":
            return self.eta / (1 + self.eta * self.tau * self.round)
        elif self.mode == "adagrad":
            # place your code here
            return self.
        return 1

    def update_step(self, grad):
        """
        Updates the round and grads_squared members.

        Parameters
        ----------
        grad : an array of size k
                the current gradient, the size must match to grads_squared
        """

        # place your code here

    def gradient(self, X, y, w):
        """
        Returns the gradient d(||Xw - y||^2 / n) / dw

        Parameters
        ----------
        X : an array of size (n, k)
                input matrix of real numbers
        y : an array of size n
                targets
        w : an array of size k
                weights
        Returns
        -------
        g : an array of size k
                the gradient for the weights
        """

        # place your code here
        return None

    def epoch(self, X, y, w, M):
        """
        Performs gradient descent with mini-batches. Input targets are split
        in mini-batches of M consecutive data points (the last batch is smaller
        than M if M doesn't divide n) and gradient descent is applied to each
        batch. The order in which gradient is applied must match the batch
        order.


        Parameters
        ----------
        X : an array of size (n, k)
                input matrix of real numbers
        y : an array of size n
                targets
        w : an array of size k
                initial weights
        M : int
                size of the minibatch, note that M doesn't necessarily divide n

        Returns
        -------
        trajectory : an array of size (r, k)
                trajectory[i, :] is the weight after the (i + 1)th minibatch, for
                example trajectory[0, :] is the weight after the first minibatch.
                r is the total number of minibatches, r = ceil(n / M).
        """

        cnt = X.shape[0]
        trajectory = np.zeros((int(np.ceil(cnt / M)), len(w)))
        # place your code here
        return trajectory

    def train(self, X, y, w, M, epochcnt, shuffle=True):
        """
        Performs gradient descent with mini-batches with several epochs.

        Parameters
        ----------
        X : an array of size (n, k)
                input matrix of real numbers
        y : an array of size n
                targets
        w : an array of size k
                initial weights
        M : int
                size of the minibatch, note that M doesn't necessarily divide n
        epochcnt : int
                number of epochs, iterations over the whole data
        shuffle : bool, optional
                If true, data is permuted before each epoch, default value: True.

        Returns
        -------
        trajectory : an array of size (epochcnt*r, k)
                complete history of trajectories for the weights, concatenation
                of trajectories produced by of self.epoch. Here, r is the
                number of minibatches in one epoch.
        """

        cnt = X.shape[0]
        Xr = X
        yr = y
        trajectory = np.zeros((0, len(w)))
        for i in range(epochcnt):
            if shuffle:
                # Permute the data points
                ind = np.random.permutation(cnt)
                Xr = X[ind, :]
                yr = y[ind]
            trajectory = np.r_[trajectory, self.epoch(Xr, yr, w, M)]
            w = trajectory[-1, :]
        return trajectory


def main(argv):
    D = np.loadtxt(argv[1])
    y = D[:, 0]
    D = D[:, 1:]

    mode = argv[2]
    eta = float(argv[3])
    tau = float(argv[4])
    batchsize = int(argv[5])
    itercnt = int(argv[6])

    w = np.zeros(D.shape[1])

    grad = Gradient(mode, eta, tau, len(w))
    print(grad.train(D, y, w, batchsize, itercnt))


# This allows the script to be used as a module and as a standalone program
if __name__ == "__main__":
    if len(sys.argv) != 7:
        print(
            "usage: python %s filename mode eta tau batch_size epochcnt" % sys.argv[0]
        )
    else:
        main(sys.argv)
