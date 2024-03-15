#!/usr/bin/env python3

import sys
import unittest
from tmc import points
from tmc.utils import load, get_out

from numpy import loadtxt, array
import numpy as np

from unittest.mock import MagicMock, patch, NonCallableMagicMock

from .utils import errormsg, ApproxTest



@points('1.2.1', '1.2.2', '1.2.3', '1.2.4', '1.2.5')
class CovarianceTester(ApproxTest):

	@patch('sklearn.decomposition.PCA', new=MagicMock)
	@patch('numpy.cov', new=MagicMock)
	def test_cov(self):
		covariance_matrix = load('src.pca', "covariance_matrix")
		X = loadtxt('test/toy.txt')
		C = covariance_matrix(X, bias=False)
		C1 = array([[0.09812, -0.054795, 0.030765], [-0.054795, 0.08727, -0.07049], [ 0.030765, -0.07049, 0.15638]])
		self.assertApprox(C, C1, errormsg("Incorrect covariance matrix (bias = False)", C1, C))

	@patch('sklearn.decomposition.PCA', new=MagicMock)
	@patch('numpy.cov', new=MagicMock)
	def test_cov_bias(self):
		covariance_matrix = load('src.pca', "covariance_matrix")
		X = loadtxt('test/toy.txt')
		C = covariance_matrix(X, bias=True)
		C1 = array([[0.078496, -0.043836, 0.024612], [-0.043836, 0.069816, -0.056392], [ 0.024612, -0.056392, 0.125104]])
		self.assertApprox(C, C1, errormsg("Incorrect covariance matrix (bias = True)", C1, C))


@points('1.2.1', '1.2.2', '1.2.3', '1.2.4', '1.2.5')
class PcaTester(ApproxTest):

	@patch('sklearn.decomposition.PCA', new=MagicMock)
	@patch('numpy.cov', new=MagicMock)
	def test_pca(self):
		pca = load('src.pca', "pca")
		X = loadtxt('test/toy.txt')
		u1 = array([-0.48326017, -0.27655539, -0.21992991,  0.66658711,  0.31315837])
		u2 = array([-0.1179138 ,  0.47968815, -0.30314492,  0.07743277, -0.1360622 ])
		o1, o2 = pca(X)
		v1 = o1 - o1.mean()
		v2 = o2 - o2.mean()
		if np.sign(u1[0]) != np.sign(v1[0]):
			v1 *= -1
		if np.sign(u2[0]) != np.sign(v2[0]):
			v2 *= -1
		self.assertApprox(u1, v1, errormsg("Incorrect first vector", u1, o1))
		self.assertApprox(u2, v2, errormsg("Incorrect second vector", u2, o2))
