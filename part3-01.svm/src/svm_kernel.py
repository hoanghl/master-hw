import sys
import numpy as np



class SVM:
	def __init__(self, C, kernel):
		""" 
		Initializes SVM with penalty C

		Parameters
		----------
		C : real
			penalty for violating the margin
		kernel : function
			kernel function, for example linear_kernel or rbf_kernel
		"""

		self.C = C
		self.kernel = kernel


	def step(self, i, j):
		""" 
		Optimizes alpha[i] and alpha[j] minimizing the dual program.
		Changes alpha[i], alpha[j]. Updates self.u

		Parameters
		----------
		i : int
			index referring to alpha[i]
		j : int
			index referring to alpha[j]

		Returns
		-------
		change : bool
			True if any alphas are updated
		"""

		# place your code here

		return True


	def optimize(self):
		"""
		Fits SVM weights to the training data.
		Uses self.X and self.y as the training data.
		"""
		cnt = self.X.shape[0]
		
		# find optimal alphas
		changes = True
		round = 0
		giveup = 1000 # This is just so that the exercise doesn't run forever, if incorrect

		while changes and round < giveup:
			changes = False
			round += 1
			for i in range(cnt):
				for j in range(cnt):
					if self.step(i, j):
						changes = True

		# find b
		for i in range(cnt):
			if self.alpha[i] > 0 and self.alpha[i] < self.C:
				self.b = self.y[i] - self.u[i]
				break


	
	def fit(self, X, y):
		"""
		Fits SVM weights to the training data

		Parameters
		----------
		X : an array of size (n, k)
			input matrix
		y : an array of size n
			labels
		"""

		cnt = X.shape[0]
		self.X = X
		self.y = y
		self.u = np.zeros(cnt) # cached svm output for the training data, equal to selft.score(self.X)
		self.alpha = np.zeros(cnt)
		self.b = 0
		self.optimize()
	

	def score(self, X):
		"""
		Computes the score for input X, negative scores indicate label -1,
		positive scores indicate label 1.

		Parameters
		----------
		X : an array of size (n, k)
			input matrix
		Returns
		-------
		p : an array of size n
			scores, p[i] = sum_j alpha[j] y[j] <x[j], z> + bias,
			where
			x[j], y[j] is the jth training sample,
			z is the ith sample in X,
			and < , > is the inner product.
		"""

		# place your code here
		return 0
		
		
	def predict(self, X):
		"""
		Predicts the labels for input X using self.score. Negative scores
		indicate label -1, positive scores indicate label 1.

		Parameters
		----------
		X : an array of size (n, k)
			input matrix

		Returns
		-------
		p : an array of size n
			predicted labels
		"""

		return np.sign(self.score(X))



def linear_kernel(x, y):
	return x @ y


def rbf_kernel(x, y):
	return np.exp(-0.5*np.sum((x - y)**2)) # RBF kernel with sigma = 1



def main(argv):
	D= np.loadtxt(argv[1])
	y = D[:,0]
	X = D[:,1:] 

	penalty = float(argv[2])

	svm = SVM(penalty, rbf_kernel)
	svm.fit(X, y)

	print('predictions for training data:')
	print(svm.predict(X))

	print('training error:')
	print(np.mean(svm.predict(X) != y))



# This allows the script to be used as a module and as a standalone program
if __name__ == "__main__": 
	if len(sys.argv) != 3:
		print('usage: python %s filename penalty' % sys.argv[0])
	else:
		main(sys.argv)
