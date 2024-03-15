import sys
import numpy as np


def lasso(X, y, s, reg, itercnt):
	"""
	Subgradient optimization for Lasso 

	Parameters
	----------
	X : an array of size (n, k)
		training input data for the regressor
	y : an array of size n
		training output for the regressor
	s : an array of size k
		initial weights
	reg : real
		weight for the regularization term
	itercnt : int
		number of iterations

	Returns
	-------
	w : an array of size k
		weights after itercnt iterations
	"""

	cnt, k = X.shape

	# make a copy of the initial vector
	# we can now change single elements of w without changing s
	w = s.copy()

	# place your code here

	return w


def main(argv):
	D = np.loadtxt(argv[1])
	y = D[:,0].copy() # copy is needed, otherwise next line will mess up the splice
	D[:,0] = 1 # replace the output column of D with constant, now the first feature gives us the bias term

	reg = float(argv[2])
	itercnt = int(argv[3])

	w = np.zeros(D.shape[1]) # starting point for weights

	print(lasso(D, y, w, reg, itercnt))



# This allows the script to be used as a module and as a standalone program
if __name__ == "__main__": 
	if len(sys.argv) != 4:
		print('usage: python %s filename' % sys.argv[0])
	else:
		main(sys.argv)
