#!/usr/bin/env python3

import sys
import unittest
from tmc import points
from tmc.utils import load, get_out

import numpy as np

from .utils import errormsg, ApproxTest



@points('4.1.1', '4.1.2', '4.1.3', '4.1.4', '4.1.5')
class EmTester(ApproxTest):

	def test_logsumrows(self):
		X = np.array([[-1000, -1001, -2000], [100, 103, 102]])
		logsumrows = load('src.em', 'logsumrows')
		S = logsumrows(X)
		S1 = np.array([-999.68673831, 103.34901222])
		self.assertApprox(S1, S, errormsg("Incorrect output for logsumrows (X = [[-1000, -1001, -2000], [100, 103, 102]])", S1, S))


	def test_responsibities(self):
		computeresponsibilities = load('src.em', 'computeresponsibilities')
		prior = np.array([0.9, 0.1])
		mu = np.array([[0], [1000]])
		C = np.array([[[1.0]], [[1.0]]])
		X = np.array([[500.001]])
		R = computeresponsibilities(X, prior, mu, C)
		R1 = np.array([[0.76803068, 0.23196932]])
		self.assertApprox(R1, R, errormsg("Incorrect output for computeresponsibilities\nprior = %s\nmu = %s\nC = %s\nX = %s\n" % (prior, mu, C, X), R1, R))

		prior = np.array([0.9, 0.1])
		mu = np.array([[0.0], [4.0]])
		C = np.array([[[1.0]], [[2.0]]])
		X = np.array([[2.0], [1], [3]])
		R = computeresponsibilities(X, prior, mu, C)
		R1 = np.array([[0.82401619, 0.17598381], [0.98653092, 0.01346908], [0.153657,   0.846343  ]])
		self.assertApprox(R1, R, errormsg("Incorrect output for computeresponsibilities\nprior = %s\nmu = %s\nC = %s\nX = %s\n" % (prior, mu, C, X), R1, R))



	def test_spherical(self):
		stats = load('src.em', 'computeparametersspherical')
		X = np.loadtxt('test/toy.txt')
		R = np.array([[0.01840738, 0.98159262], [0.69408228, 0.30591772], [0.5846214,  0.4153786 ], [0.58096681, 0.41903319], [0.55433175, 0.44566825]])
		prior, mu, C = stats(R, X)
		prior1 = np.array([0.48648192, 0.51351808])
		mu1 = np.array([[0.29458028, 0.5468195,  0.48396362], [0.24281914, 0.38164615, 0.64371726]])
		C1 = 0.08651849*np.array([np.eye(3), np.eye(3)])
		self.assertApprox(prior1, prior, errormsg("Incorrect prior for computeparametersspherical\nR = %s\nX = toy.txt\n" % (R), prior1, prior))
		self.assertApprox(mu1, mu, errormsg("Incorrect mu for computeparametersspherical\nR = %s\nX = toy.txt\n" % (R), mu1, mu))
		self.assertApprox(C1, C, errormsg("Incorrect C for computeparametersspherical\nR = %s\nX = toy.txt\n" % (R), C1, C))


	def test_diagonal(self):
		stats = load('src.em', 'computeparametersdiagonal')
		X = np.loadtxt('test/toy.txt')
		R = np.array([[0.01840738, 0.98159262], [0.69408228, 0.30591772], [0.5846214,  0.4153786 ], [0.58096681, 0.41903319], [0.55433175, 0.44566825]])
		prior, mu, C = stats(R, X)
		prior1 = np.array([0.48648192, 0.51351808])
		mu1 = np.array([[0.29458028, 0.5468195,  0.48396362], [0.24281914, 0.38164615, 0.64371726]])

		C1 = np.array([np.diag([0.10497885, 0.03838912, 0.11379078]), np.diag([0.05210405, 0.08631598, 0.12340598])])


		self.assertApprox(prior1, prior, errormsg("Incorrect prior for computeparametersdiagonal\nR = %s\nX = toy.txt\n" % (R), prior1, prior))
		self.assertApprox(mu1, mu, errormsg("Incorrect mu for computeparametersdiagonal\nR = %s\nX = toy.txt\n" % (R), mu1, mu))
		self.assertApprox(C1, C, errormsg("Incorrect C for computeparametersdiagonal\nR = %s\nX = toy.txt\n" % (R), C1, C))



	def test_same(self):
		stats = load('src.em', 'computeparameterssame')
		X = np.loadtxt('test/toy.txt')
		R = np.array([[0.01840738, 0.98159262], [0.69408228, 0.30591772], [0.5846214,  0.4153786 ], [0.58096681, 0.41903319], [0.55433175, 0.44566825]])
		prior, mu, C = stats(R, X)
		prior1 = np.array([0.48648192, 0.51351808])
		mu1 = np.array([[0.29458028, 0.5468195,  0.48396362], [0.24281914, 0.38164615, 0.64371726]])

		S = np.array([[0.07782669, -0.04597183,  0.02667775], [-0.04597183,  0.06300043, -0.04980006], [ 0.02667775, -0.04980006,  0.11872836]])
		C1 = np.array([S, S])

		self.assertApprox(prior1, prior, errormsg("Incorrect prior for computeparameterssame\nR = %s\nX = toy.txt\n" % (R), prior1, prior))
		self.assertApprox(mu1, mu, errormsg("Incorrect mu for computeparameterssame\nR = %s\nX = toy.txt\n" % (R), mu1, mu))
		self.assertApprox(C1, C, errormsg("Incorrect C for computeparameterssame\nR = %s\nX = toy.txt\n" % (R), C1, C))


	def test_normal(self):
		stats = load('src.em', 'computeparameters')
		X = np.loadtxt('test/toy.txt')
		R = np.array([[0.01840738, 0.98159262], [0.69408228, 0.30591772], [0.5846214,  0.4153786 ], [0.58096681, 0.41903319], [0.55433175, 0.44566825]])
		prior, mu, C = stats(R, X)
		prior1 = np.array([0.48648192, 0.51351808])
		mu1 = np.array([[0.29458028, 0.5468195,  0.48396362], [0.24281914, 0.38164615, 0.64371726]])

		C1 = np.array(
			[[[ 0.10497885, -0.05863552,  0.02955697],
			  [-0.05863552,  0.03838912, -0.02578816],
			  [ 0.02955697, -0.02578816,  0.11379078]],
			 [[ 0.05210405, -0.03397487,  0.02395011],
			  [-0.03397487,  0.08631598, -0.07254776],
			  [ 0.02395011, -0.07254776,  0.12340598]]])

		self.assertApprox(prior1, prior, errormsg("Incorrect prior for computeparameters\nR = %s\nX = toy.txt\n" % (R), prior1, prior))
		self.assertApprox(mu1, mu, errormsg("Incorrect mu for computeparameters\nR = %s\nX = toy.txt\n" % (R), mu1, mu))
		self.assertApprox(C1, C, errormsg("Incorrect C for computeparameters\nR = %s\nX = toy.txt\n" % (R), C1, C))


	def test_em_normal(self):
		X = np.loadtxt('test/toy.txt')
		em = load('src.em', 'em')
		S = np.array([[0.01840738, 0.98159262], [0.69408228, 0.30591772], [0.5846214,  0.4153786 ], [0.58096681, 0.41903319], [0.55433175, 0.44566825]])

		stats = load('src.em', 'computeparameters')
		R, prior, mu, C = em(X, S, 3, stats)

		prior1 = np.array([0.57368519, 0.42631481])
		mu1 = np.array([[0.36427589, 0.52459259, 0.54958373], [0.13844304, 0.37777013, 0.58809112]])

		C1 = np.array(
			[[[ 0.10339346, -0.06252964,  0.0180394 ],
			  [-0.06252964,  0.03981585, -0.01913735],
			  [ 0.0180394,  -0.01913735,  0.13048426]],
			 [[ 0.01573365, -0.03770214,  0.03844554],
			  [-0.03770214,  0.0978199,  -0.10328152],
			  [ 0.03844554, -0.10328152,  0.1170132 ]]])

		R1 = np.array(
			[[2.69991373e-39, 1.00000000e+00],
			 [1.00000000e+00, 2.14620464e-75],
			 [1.00000000e+00, 1.43240335e-27],
			 [4.14843864e-01, 5.85156136e-01],
			 [1.70913688e-03, 9.98290863e-01]])

		self.assertApprox(prior1, prior, errormsg("Incorrect prior for em\nR = %s\nX = toy.txt\nmode=normal\nitercnt=3\n" % (S), prior1, prior))
		self.assertApprox(mu1, mu, errormsg("Incorrect mu for em\nR = %s\nX = toy.txt\nmode=normal\nitercnt=3\n" % (S), mu1, mu))
		self.assertApprox(C1, C, errormsg("Incorrect C for em\nR = %s\nX = toy.txt\nmode=normal\nitercnt=3\n" % (S), C1, C))
		self.assertApprox(R1, R, errormsg("Incorrect output R for em\nR = %s\nX = toy.txt\nmode=normal\nitercnt=3\n" % (S), R1, R))


	def test_em_same(self):
		X = np.loadtxt('test/toy.txt')
		em = load('src.em', 'em')
		S = np.array([[0.01840738, 0.98159262], [0.69408228, 0.30591772], [0.5846214,  0.4153786 ], [0.58096681, 0.41903319], [0.55433175, 0.44566825]])

		stats = load('src.em', 'computeparameterssame')
		R, prior, mu, C = em(X, S, 3, stats)

		prior1 = np.array([0.53391506, 0.46608494])
		mu1 = np.array([[0.34084778, 0.53215143, 0.52775008], [0.15509532, 0.35327431, 0.62528245]])
		mu1 = np.array([[0.30107709, 0.53301858, 0.495302], [0.23010915, 0.38064597, 0.6469868 ]])
		C1 = np.array(
			[[[ 0.07724268, -0.04652695,  0.02729081],
			  [-0.04652695,  0.06403835, -0.05064043],
			  [ 0.02729081, -0.05064043,  0.1193784 ]],
			 [[ 0.07724268, -0.04652695,  0.02729081],
			  [-0.04652695,  0.06403835, -0.05064043],
			  [ 0.02729081, -0.05064043,  0.1193784 ]]])
		R1 = np.array(
			[[0.1136505,  0.8863495 ],
			 [0.78636137, 0.21363863],
			 [0.68284857, 0.31715143],
			 [0.71316785, 0.28683215],
			 [0.48574661, 0.51425339]])

		self.assertApprox(prior1, prior, errormsg("Incorrect prior for em\nR = %s\nX = toy.txt\nmode=same\nitercnt=3\n" % (S), prior1, prior))
		self.assertApprox(mu1, mu, errormsg("Incorrect mu for em\nR = %s\nX = toy.txt\nmode=same\nitercnt=3\n" % (S), mu1, mu))
		self.assertApprox(C1, C, errormsg("Incorrect C for em\nR = %s\nX = toy.txt\nmode=same\nitercnt=3\n" % (S), C1, C))
		self.assertApprox(R1, R, errormsg("Incorrect output R for em\nR = %s\nX = toy.txt\nmode=same\nitercnt=3\n" % (S), R1, R))
