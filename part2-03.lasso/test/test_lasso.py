#!/usr/bin/env python3

import sys
import unittest
from tmc import points
from tmc.utils import load, get_out

import numpy as np

from .utils import errormsg, ApproxTest


@points('2.3.1', '2.3.2', '2.3.3', '2.3.4', '2.3.5')
class LassoTester(ApproxTest):

	def test_lasso(self):
		lasso = load('src.lasso', 'lasso')

		D = np.loadtxt('test/toy.txt')
		y = D[:,0].copy() # copy is needed, otherwise next line will mess up the splice
		D[:,0] = 1 # replace the label column of D with constant, now the first feature gives us the bias term

		w = np.zeros(D.shape[1])

		w = lasso(D, y, w, 10, 100)
		
		w1 = np.zeros(11)
		w1[0] = 3.01421845
		w1[3] = 0.97828604
		w1[5] = 0.9834849

		self.assertApprox(w1, w, errormsg("Incorrect weights (input toy.txt, lam = 10, itercnt = 100)", w1, w))
