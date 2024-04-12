import sys
import numpy as np

class LDA:
    def fit(self, X, y, w):
        """
        Fits weights (self.w) and bias (self.b) for the variant of LDA classifier.

        Parameters
        ----------
        X : an array of size (n, k)
            training input data for the classifier
        y : an array of size n
            training labels for the classifier, elements must be 1 or -1
        w : an array of size n
            weights for the data points
        """

        C = np.linalg.inv(np.cov(X.T, aweights=w) + 0.001*np.eye(X.shape[1])) # add diagonal matrix to deal with singular matrices
        ind0 = y == -1
        ind1 = y == 1
        m0 = np.sum(X[ind0, :] * w[ind0, np.newaxis], axis=0) / np.sum(w[ind0])
        m1 = np.sum(X[ind1, :] * w[ind1, np.newaxis], axis=0) / np.sum(w[ind1])

        self.w = C @ (m1 - m0)
        self.b = self.find_threshold(X, y, w)
        if self.b == None:
            self.w = -self.w
            self.b = self.find_threshold(X, y, w)

    
    def find_threshold(self, X, y, w):
        """ 
        Helper function for finding the optimal bias given the weight.

        Parameters
        ----------
        X : an array of size (n, k)
            training input data for the classifier
        y : an array of size n
            training labels for the classifier, elements must be 1 or -1
        w : an array of size n
            weights for the data points

        Returns
        -------
        bias: int
            optimal bias
        """

        score = X @ self.w
        ind = np.argsort(score)

        err = np.cumsum(y[ind]*w[ind])
        i = np.argmin(err)
        if i == len(err) - 1:
            return None
        return (score[ind[i]] + score[ind[i + 1]]) / 2


    def predict(self, X):
        """ 
        Classifies data points in X.

        Parameters
        ----------
        X : an array of size (n, k)
            input data for the classifier

        Returns
        -------
        y : an array of size n
            prediction for the n data points, either -1 or 1
        """

        return np.sign(X @ self.w - self.b) 



def adaboost(X, y, itercnt):
    """ 
    Adaboost algorithm

    Parameters
    ----------
    X : an array of size (n, k)
        training input data for the classifier
    y : an array of size n
        training labels for the classifier, elements must be 1 or -1
    itercnt : int
        number of iterations

    Returns
    -------
    output : an array of size n
        ith element = weighted total vote for the (i + 1)th data point
    err_individual : an array of size itercnt
        ith element = weighted misclassification error of the (i + 1)th
        individual classifier
    err_ensemble : an array of size itercnt
        ith element = misclassification error of the ensemble after i + 1
        iterations
    err_exponential : an array of size itercnt
        ith element = exponential loss of the ensemble after i + 1 iterations
    """

    cnt, k = X.shape

    err_individual = np.zeros(itercnt)
    err_ensemble = np.zeros(itercnt)
    err_exponential = np.zeros(itercnt)
    output = np.zeros(cnt)

    # place your code here
    classifier = LDA()

    weights = np.ones((cnt, )) / cnt
    ensemble_pred = 0
    for i in range(itercnt):
        classifier.fit(X, y, weights)
        y_pred = classifier.predict(X)

        epsilon_m = 1.0 / np.sum(weights) * np.sum(weights * (y != y_pred).astype(np.int32))
        err_individual[i] = epsilon_m

        beta_m = 1.0 / 2 * np.log((1 - epsilon_m) / epsilon_m)

        weights = weights * np.exp(-beta_m * y * y_pred)

        ensemble_pred += beta_m * y_pred
        err_ensemble[i] = (np.sign(ensemble_pred) != y).astype(np.int32).sum() * 1.0 / len(y)

        exponential = 0
        for j in range(cnt):
            exponential += np.exp(-y[j] * ensemble_pred[j])
        err_exponential[i] = 1.0 / cnt * exponential

        output += beta_m * y_pred
    
    return output, err_individual, err_ensemble, err_exponential


def main(argv):
    D = np.loadtxt(argv[1])
    labels = D[:,0].copy() # copy is needed, otherwise next line will mess up the splice
    D[:,0] = 1 # replace the label column of D with constant, now the first feature gives us the bias term

    itercnt = int(argv[2])

    output, err_individual, err_ensemble, err_exponential = adaboost(D, labels, itercnt)

    print('individual error:')
    print(err_individual)
    print('ensemble error:')
    print(err_ensemble)
    print('exponential error:')
    print(err_exponential)

    # uncomment these if you want to save the output in a file
    #np.savetxt('output.txt', output)
    #np.savetxt('err_individual.txt', err_individual)
    #np.savetxt('err_ensemble.txt', err_ensemble)
    #np.savetxt('err_exponential.txt', err_exponential)


# This allows the script to be used as a module and as a standalone program
if __name__ == "__main__": 
    if len(sys.argv) != 3:
        print('usage: python %s filename' % sys.argv[0])
    else:
        main(sys.argv)
