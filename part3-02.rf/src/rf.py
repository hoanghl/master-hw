import sys
import numpy as np

def entropy(C):
	""" 
	Computes the entropy of each column of matrix C.

	Parameters
	----------
	C : an array of size (n, k)
		input matrix, must be positive, columns do not need to sum to 1.

	Returns
	-------
	e : an array of size k
		e[i] = entropy of the ith column (normalized to sum to 1)
	"""

	with np.errstate(divide='ignore'): # suppresses warning for log(0)
		L = np.log(C)
		N = np.log(np.sum(C, axis=0, keepdims=True))
	L[np.isinf(L)] = 0 # avoids 0*-inf
	N[np.isinf(N)] = 0 # avoids 0*-inf
	return -np.sum(C*(L - N), axis=0)



def entropy_gain(X, y):
	""" 
	Computes the entropy gain in labels when split is done with ith feature.

	Parameters
	----------
	X : an array of size (n, k)
		input matrix, entries must be 0/1
	y : an array of size n
		y, entries must be +/-1

	Returns
	-------
	gain : an array of size k
		gain[i] =  entropy gain in labels when split is done with ith feature
	"""

	cnt, k = X.shape
	labels = (y + 1) / 2 # remap labels to 0/1 for easier handling
	m = np.sum(labels)
	
	P1 = np.zeros((2, k))
	P0 = np.zeros((2, k))

	cnt1 = np.sum(X, axis=0) # count(feature = 1)
	cnt0 = cnt - cnt1        # count(feature = 0)

	P1[1,:] = labels @ X      # count(feature = 1, labels = 1)
	P1[0,:] = cnt1 - P1[1,:]  # count(feature = 1, labels = 0)
	P0[1,:] = m - P1[1,:]     # count(feature = 0, labels = 1)
	P0[0,:] = cnt0 - P0[1,:]  # count(feature = 0, labels = 0)

	Q = np.array([cnt - m, m])
	
	return entropy(Q) - entropy(P1) - entropy(P0) 



def fit(X, y, featurecnt):
	""" 
	Constructs a decision tree: finds an optimal split among the candidates,
	and recurses on the split data.

	Parameters
	----------
	X : an array of size (n, k)
		input matrix, entries must be 0/1
	y : an array of size n
		labels, entries must be +/-1
	featurecnt: int
		number of features sampled, from which the best split is selected

	Returns
	-------
	r : a node representing the root to the decision tree
		r is either 0/1 if it is a leaf, or a triple (f, n1, n2), where f is
		the feature used for splitting, n1 is the root of the left branch and
		n2 is the root of the right branch. 
	"""

	# Compute the gain for every feature
	g = entropy_gain(X, y)

	# List of features for which split is possible 
	# Note: this list can be shorter than featurecnt
	feasible = np.where(g > 0.00001)[0]

	# No candidates left, make a leaf with majority label
	if len(feasible) == 0:
		return int(np.sign(0.5 + np.sum(y))) # 0.5 trick makes sure that we never land on 0

	# Select best feature
	# Change the code to use sampling
	ind = feasible[np.argmax(g[feasible])]

	# place your code here
	
	split = X[:,ind] == 1
	return (ind,
		fit(X[~split, :], y[~split], featurecnt),
		fit(X[split, :], y[split], featurecnt))



def predict_sample(T, x):
	""" 
	Predicts the output using the tree T for input x. Traverses the tree
	recursively until the leaf is found.

	Parameters
	----------
	T : tree
		constructed using fit
	x : an array of size k
		input data

	Returns
	-------
	y : -1 or 1
		predicted label using T for x
	"""

	if type(T) is tuple:
		# not a leaf, find the correct branch and recurse
		ind = int(x[T[0]])
		return predict_sample(T[1 + ind], x)
	return T # leaf, return the prediction


def predict(T, X):
	""" 
	Predicts the output using the tree T for input matrix X
	using predict_sample.

	Parameters
	----------
	T : tree
		constructed using fit
	X : an array of size (n, k)
		input matrix, entries must be 0/1

	Returns
	-------
	p : an array of size n
		predictions, entries are either +1 or -1
	"""

	cnt = X.shape[0]
	p = np.zeros(cnt)
	for i in range(cnt):
		p[i] = predict_sample(T, X[i, :])
	return p
		
	


def rf(Xtrain, ltrain, treecnt, samplecnt, featurecnt, Xtest, ltest):
	""" 
	Randomized forest classifier.  Trains the classifier using the training
	data and computes the misclassification rate using the testing data.

	Parameters
	----------
	Xtrain : an array of size (n, k)
		training input matrix, entries must be 0/1
	ltrain : an array of size n
		training labels, entries must be 1 or -1
	treecnt : int
		number of trees
	samplecnt : int
		number of data points sampled for each tree
	featurecnt: int
		number of features sampled, from which the best split is selected
	Xtest : an array of size (m, k)
		test input matrix, entries must be 0/1
	ltest : an array of size m
		test labels, entries must be 1 or -1

	Returns
	-------
	p : an array of size m
		p[i] = vote total for the ith data point
	misclass : an array of size treecnt
		misclass[i] = misclassification rate after (i + 1)th iteration
	"""

	cnt, k = Xtrain.shape

	p = np.zeros(Xtest.shape[0])
	misclass = np.zeros(treecnt)

	# place your code here
	
	return p, misclass



def main(argv):
	np.random.seed(2022) # Forces random to be repeatable. Remove when random seed is desired. 

	D = np.loadtxt(argv[1])
	ltrain = D[:,0]
	Xtrain = D[:,1:] 

	D = np.loadtxt(argv[2])
	ltest = D[:,0]
	Xtest = D[:,1:] 

	treecnt = int(argv[3])

	cnt, k = Xtrain.shape

	p, misclass = rf(Xtrain, ltrain, treecnt, cnt, int(np.sqrt(k)), Xtest, ltest)

	print('Votes for testing data: negative -> -1 class, positive -> 1 class')
	print(p)
	print('Misclassification rate')
	print(misclass)



# This makes sure the main function is not called immediatedly
# when TMC imports this module
if __name__ == "__main__": 
	if len(sys.argv) != 4:
		print('usage: python %s train_filename test_filename number_of_trees' % sys.argv[0])
	else:
		main(sys.argv)
