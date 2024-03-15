#!/usr/bin/env python3

import sys
import unittest
from tmc import points
from tmc.utils import load, get_out

import numpy as np

from .utils import errormsg, ApproxTest


@points('2.2.1', '2.2.2', '2.2.3', '2.2.4', '2.2.5')
class IrlsTester(ApproxTest):

	def test_irls(self):
		irls = load('src.irls', 'irls')

		D = np.loadtxt('test/toy.txt')
		labels = D[:,0].copy() # copy is needed, otherwise next line will mess up the splice
		D[:,0] = 1 # replace the label column of D with constant, now the first feature gives us the bias term

		w = np.zeros(D.shape[1])

		w, err, misclass = irls(D, labels, w, 10)

		w1 = np.array([
			14.2225575,  -3.83007457,  3.9215949,  -2.73087605,  0.9907367,  -1.9727234,
			-0.2329434,  -1.76999314, -5.60257676, -0.11626992, -5.22387472])
		err1 = np.array([
			138.62943611,  51.57260799, 31.2257094,   21.52724814,  16.25098254,
			13.29340328, 11.87038151, 11.4422655, 11.38955273,  11.38827445,
			11.38827308])
		misclass1 = np.array([0.5, 0.04, 0.035, 0.035, 0.03, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025])

		self.assertApprox(w1, w, errormsg("Incorrect weights (input toy.txt, w = 0, itercnt = 10)", w1, w))
		self.assertApprox(err1, err, errormsg("Incorrect error (input toy.txt, w = 0, itercnt = 10)", err1, err))
		self.assertApprox(misclass1, misclass, errormsg("Incorrect misclassication rate (input toy.txt, w = 0, itercnt = 10)", misclass1, misclass))
