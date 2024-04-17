import sys
import numpy as np


def logsumrows(X):
    """
    Computes the sums of rows of log-numbers

    Parameters
    ----------
    X : an array of size (n, m)
        matrix of log-numbers

    Returns
    -------
    s : an array of size n
        ith element is the sum of the ith row
    """

    # place your code here
    m = X.shape[1]

    M = -np.max(X, axis=-1)
    M_ = np.repeat(M[:, None], m, -1)
    out = np.log(np.sum(np.exp(X + M_), axis=-1)) - M
    return out


def computeparameters(R, X):
    """
    Computes the optimal parameters for the Gaussian mixture model
    
    Parameters
    ----------
    X : an array of size (n, m)
        input data matrix
    R : an array of size (n, k)
        responsibilities: R[i, j] = probability of data point i belonging to
        jth component.

    Returns
    -------
    prior : an array of size k
        prior probabilities of the components
    mu : an array of size (k, m)
        optimal means: mu[i, :] is the optimal mean of the ith component
    C : an array of size (k, m, m)
        optimal covariances: C[i, :, :] is the optimal covariance matrix of the
        ith component
    """

    k = R.shape[1]
    cnt, dim = X.shape

    prior = np.zeros(k)
    mu = np.zeros((k, dim))
    C = np.zeros((k, dim, dim))

    # place your code here
    for j in range(k):
        mu_comp = 0
        for n in range(cnt):
            mu_comp += R[n, j] * X[n]
        mu[j] = 1.0 / np.sum(R[:, j]) * mu_comp

        s_comp = 0
        for n in range(cnt):
            s_comp += R[n, j] * (X[n][:, None] - mu[j][:, None]) @ (X[n][:, None] - mu[j][:, None]).T
        C[j] = s_comp / np.sum(R[:, j])

        prior[j] = 1.0 / cnt * np.sum(R[:, j])
    

    return prior, mu, C


def computeparametersdiagonal(R, X):
    """
    Computes the optimal parameters for the Gaussian mixture model with
    diagonal covariance matrices.
    
    Parameters
    ----------
    X : an array of size (n, m)
        input data matrix
    R : an array of size (n, k)
        responsibilities: R[i, j] = probability of data point i belonging to
        jth component.

    Returns
    -------
    prior : an array of size k
        prior probabilities of the components
    mu : an array of size (k, m)
        optimal means: mu[i, :] is the optimal mean of the ith component
    C : an array of size (k, m, m)
        optimal covariances: C[i, :, :] is the optimal covariance diagonal
        matrix of the ith component 
    """

    def get_Ai(i, n):
        A = np.zeros((n, n))
        A[i, i] = 1
        return A

    k = R.shape[1]
    cnt, dim = X.shape

    prior = np.zeros(k)
    mu = np.zeros((k, dim))
    C = np.zeros((k, dim, dim))

    # place your code here
    for j in range(k):
        mu_comp = 0
        for n in range(cnt):
            mu_comp += R[n, j] * X[n]
        mu[j] = 1.0 / np.sum(R[:, j]) * mu_comp

        s_comp = 0
        for n in range(cnt):
            s_comp += R[n, j] * (X[n][:, None] - mu[j][:, None]) @ (X[n][:, None] - mu[j][:, None]).T
        C[j] = np.diag(np.diag(s_comp / np.sum(R[:, j])))

        prior[j] = 1.0 / cnt * np.sum(R[:, j])

    return prior, mu, C


def computeparameterssame(R, X):
    """
    Computes the optimal parameters for the Gaussian mixture model with
    equal covariance matrices.
    
    Parameters
    ----------
    X : an array of size (n, m)
        input data matrix
    R : an array of size (n, k)
        responsibilities: R[i, j] = probability of data point i belonging to
        jth component.

    Returns
    -------
    prior : an array of size k
        prior probabilities of the components
    mu : an array of size (k, m)
        optimal means: mu[i, :] is the optimal mean of the ith component
    C : an array of size (k, m, m)
        optimal covariances: C[i, :, :] is the optimal covariance matrix of the
        ith components, the covariance matrices must be the same
    """

    k = R.shape[1]
    cnt, dim = X.shape

    prior = np.zeros(k)
    mu = np.zeros((k, dim))
    C = np.zeros((k, dim, dim))

    # place your code here
    for j in range(k):
        mu_comp = 0
        for n in range(cnt):
            mu_comp += R[n, j] * X[n]
        mu[j] = 1.0 / np.sum(R[:, j]) * mu_comp

        prior[j] = 1.0 / cnt * np.sum(R[:, j])

    for j in range(k):
        for n in range(cnt):
            C += R[n, j] * (X[n][:, None] - mu[j][:, None]) @ (X[n][:, None] - mu[j][:, None]).T
    C *= 1.0 / np.sum(R)
    return prior, mu, C


def computeparametersspherical(R, X):
    """
    Computes the optimal parameters for the Gaussian mixture model with
    equal diagonal spherical covariance matrices.
    
    Parameters
    ----------
    X : an array of size (n, m)
        input data matrix
    R : an array of size (n, k)
        responsibilities: R[i, j] = probability of data point i belonging to
        jth component.

    Returns
    -------
    prior : an array of size k
        prior probabilities of the components
    mu : an array of size (k, m)
        optimal means: mu[i, :] is the optimal mean of the ith component
    C : an array of size (k, m, m)
        Optimal covariances: C[i, :, :] is the optimal covariance diagonal
        matrix of the ith component. The numbers on the diagonals must be equal.
    """

    k = R.shape[1]
    cnt, dim = X.shape

    prior = np.zeros(k)
    mu = np.zeros((k, dim))
    C = np.zeros((k, dim, dim))

    # place your code here
    for j in range(k):
        mu_comp = 0
        for n in range(cnt):
            mu_comp += R[n, j] * X[n]
        mu[j] = 1.0 / np.sum(R[:, j]) * mu_comp

        prior[j] = 1.0 / cnt * np.sum(R[:, j])

    sigma = 0
    for j in range(k):
        for n in range(cnt):
            sigma += R[n, j] * (X[n][:, None] - mu[j][:, None]).T @ (X[n][:, None] - mu[j][:, None])
    sigma /= (dim * np.sum(R))
    C = sigma * np.identity(dim)
    C = C[None, :].repeat(k, axis=0)
    return prior, mu, C


def computeresponsibilities(X, prior, mu, C):
    """
    Computes responsibilities: R[i, j] = probability of data point i belonging
    to jth component.

    Parameters
    ----------
    X : an array of size (n, m)
        input data matrix
    prior : an array of size k
        prior probabilities of the components
    mu : an array of size (k, m)
        mu[i, :] is the mean of the ith component
    C : an array of size (k, m, m)
        C[i, :, :] is the covariance matrix of the ith component
    
    Returns
    -------
    R : an array of size (n, k)
        responsibilities: R[i, j] = probability of data point i belonging to
        jth component.
    """

    k = prior.shape[0]
    cnt, m = X.shape
    R = np.zeros((cnt, k))

    # place your code here
    L = np.zeros((cnt, k))

    for i in range(cnt):
        for j in range(k):
            L[i, j] = -m/2.0 * np.log(2 * np.pi) \
                - 1.0/2 * np.log(np.linalg.det(C[j])) \
                - 1.0/2 * (X[i] - mu[j]).T @ np.linalg.inv(C[j]) @ (X[i] - mu[j]) \
                + np.log(prior[j])
            
    R = L - logsumrows(L)[:, None]
    R = np.exp(R)

    return R


def em(X, R, itercnt, stats):
    """
    EM algorithm: computes model parameters given the responsibilities and
    computes the responsibilities given the parameters.  Repeats itercnt times.

    Parameters
    ----------
    X : an array of size (n, m)
        input data matrix
    R : an array of size (n, k)
        initial responsibilities: R[i, j] = probability of data point i
        belonging to jth component.
    itercnt : int
        number of iterations
    stats : function
        Function for computing the model parameters given the responsibilities,
        for example, computeparameters

    Returns
    -------
    R : an array of size (n, k)
        final responsibilities: R[i, j] = probability of data point i
        belonging to jth component.
    prior : an array of size k
        prior probabilities of the components
    mu : an array of size (k, m)
        mu[i, :] is the mean of the ith component
    C : an array of size (k, m, m)
        C[i, :, :] is the covariance matrix of the ith component

    """
    n, k = R.shape
    m = X.shape[1]

    # place your code here
    for _ in range(itercnt):
        prior, mu, C = stats(R, X)

        R = computeresponsibilities(X, prior, mu, C)

    return R, prior, mu, C


def main(argv):
    np.random.seed(2022) # Forces random to be repeatable. Remove when random seed is desired. 

    X = np.loadtxt(argv[1])
    k = int(argv[2])
    mode = argv[3]
    itercnt = int(argv[4])

    n, m = X.shape
    R = np.random.rand(X.shape[0], k)
    R = R / R.sum(axis=1)[:,np.newaxis]
    print(R)

    if mode == 'normal':
        R, prior, mu, C = em(X, R, itercnt, computeparameters)
    elif mode == 'diag':
        R, prior, mu, C = em(X, R, itercnt, computeparametersdiagonal)
    elif mode == 'same':
        R, prior, mu, C = em(X, R, itercnt, computeparameterssame)
    elif mode == 'sphere':
        R, prior, mu, C = em(X, R, itercnt, computeparameterssame)
    else:
        print("Mode %s unrecognized" % mode)
        return

    print('R')
    print(R)
    print('priors')
    print(prior)
    print('means')
    print(mu)
    print('covariance matrices')
    print(C)



# This allows the script to be used as a module and as a standalone program
if __name__ == "__main__": 
    if len(sys.argv) != 5:
        print('usage: python %s filename number_of_factors iteration_count' % sys.argv[0])
    else:
        main(sys.argv)
