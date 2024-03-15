#!/usr/bin/env python3

import sys
import unittest
from tmc import points
from tmc.utils import load, get_out

import numpy as np

from .utils import errormsg, ApproxTest



@points('3.3.1', '3.3.2', '3.3.3', '3.3.4', '3.3.5')
class AdaboostTester(ApproxTest):

	def test_adaboost(self):
		adaboost = load('src.adaboost', 'adaboost')

		D = np.loadtxt('test/toy.txt')
		labels = D[:,0].copy() # copy is needed, otherwise next line will mess up the splice
		D[:,0] = 1 # replace the label column of D with constant, now the first feature gives us the bias term

		output, err_individual, err_ensemble, err_exponential = adaboost(D, labels, 100)

		output1 = np.loadtxt('test/output.txt')
		err_individual1 = np.loadtxt('test/err_individual.txt')
		err_ensemble1 = np.loadtxt('test/err_ensemble.txt')
		err_exponential1 = np.loadtxt('test/err_exponential.txt')


		self.assertApprox(output1, output,
			errormsg("Incorrect output (input toy.txt, itercnt = 100)", output1, output))
		self.assertApprox(err_individual1, err_individual,
			errormsg("Incorrect individual error (input toy.txt, itercnt = 100)", err_individual1, err_individual))
		self.assertApprox(err_ensemble1, err_ensemble,
			errormsg("Incorrect ensemble error (input toy.txt, itercnt = 100)", err_ensemble1, err_ensemble))
		self.assertApprox(err_exponential1, err_exponential,
			errormsg("Incorrect exponential error (input toy.txt, itercnt = 100)", err_exponential1, err_exponential))
