#!/usr/bin/env python3

import sys
import unittest
from tmc import points
from tmc.utils import load, get_out

import numpy as np

from .utils import errormsg, ApproxTest


@points('3.1.1', '3.1.2', '3.1.3', '3.1.4', '3.1.5')
class SVMTester(ApproxTest):

	def test_svm(self):
		SVM = load('src.svm_kernel', 'SVM')
		rbf = load('src.svm_kernel', 'rbf_kernel')

		D = np.loadtxt('test/toy.txt')
		y = D[:,0]
		X = D[:,1:]

		svm = SVM(1, rbf)
		svm.fit(X, y)
		p = svm.score(X)
		palt = svm.score(X[0:3,:])

		alpha1 = np.array([2.0/3, 2.0/3, 2.0/3, 0., 1., 1.])
		u1 = np.array([2/3.0,  2/3.0,  2/3.0,  2/3.0, -1, -1])
		u2 = u1[0:3]

		self.assertApprox(alpha1, svm.alpha, errormsg("Incorrect self.alpha (input toy.txt, rbf, penalty=1)", alpha1, svm.alpha), tol=0.001)
		self.assertApprox(u1, svm.u, errormsg("Incorrect self.u (input toy.txt, rbf, penalty=1)", u1, svm.u), tol=0.001)
		self.assertApprox(u1 + 1/3.0, p, errormsg("Incorrect output of score() (input toy.txt, rbf, penalty=1, X = training data)", u1 + 1/3.0, p), tol=0.001)
		self.assertApprox(u2 + 1/3.0, palt, errormsg("Incorrect output of score() (input toy.txt, rbf, penalty=1, X = first 3 training datapoints)", u2 + 1/3.0, palt), tol=0.001)
