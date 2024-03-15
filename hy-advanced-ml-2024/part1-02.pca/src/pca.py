import sys
from numpy import array, loadtxt, argsort
from numpy.linalg import eig


def covariance_matrix(X, bias=False):
	"""
	Computes covariance matrix.

	Parameters
	----------
	X : an array of size (n, k)
		input data
	bias: bool, optional
		If True, then the normalization should be n, otherwise it is n - 1.
		Default value is False.

	Returns
	-------
	C : an array of size (k, k)
		covariance matrix
	"""
	n, _ = X.shape

	dof = n if bias is True else (n - 1)
	diff = X - X.mean(axis=0)[None, :]
	# [n, k]
	cov = 1 / dof * (diff.T @ diff)
	return cov


def pca(X):
	"""
	Computes PCA with 2 components

	Parameters
	----------
	X : an array of size (n, k)
		input data

	Returns
	-------
	v1 : an array of size n
		ith element = first principal component of the ith data point
	v2 : an array of size n
		ith element = second principal component of the ith data point
	"""

	v1 = None
	v2 = None
	
	cov = covariance_matrix(X)
	# [k, k]

	eigenvalues, eigenvectors = eig(cov)
	idx_max = argsort(eigenvalues)[-2:]
	eig_vecs = eigenvectors[:, idx_max]
	# [k, 2]

	diff = X - X.mean(axis=0)[None, :]
	vecs_pca = diff @ eig_vecs
	# [n, 2]

	v1, v2 = vecs_pca[:, 1], vecs_pca[:, 0]

	# place your code here
	return v1, v2


def main(argv):
	X = loadtxt(argv[1])
	print(covariance_matrix(X))
	print(pca(X))



# This allows the script to be used as a module and as a standalone program
if __name__ == "__main__": 
	if len(sys.argv) != 2:
		print('usage: python %s filename' % sys.argv[0])
	else:
		main(sys.argv)
