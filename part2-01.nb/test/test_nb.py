#!/usr/bin/env python3

import sys
import unittest
from tmc import points
from tmc.utils import load, get_out

import numpy as np

from .utils import errormsg, ApproxTest


@points('2.1.1', '2.1.2', '2.1.3', '2.1.4', '2.1.5')
class NbTester(ApproxTest):

	def test_nb(self):
		D = np.loadtxt('test/toy.txt')
		X = D[:, 1:]
		labels = D[:,0]
		nb = load('src.nb', 'nb')
		w, b = nb(X, labels)
		w1 = np.array([-1.5040774 ,  1.38629436,  0.69314718])
		b1 = -0.040821994520255256
		self.assertApprox(w, w1, errormsg("Incorrect weight vector (input toy.txt)", w1, w))
		self.assertApprox(b, b1, errormsg("Incorrect bias (input toy.txt)", b1, b))
