#!/usr/bin/env python3

import sys
import unittest
from tmc import points
from tmc.utils import load, get_out

import numpy as np

from .utils import errormsg, ApproxTest


@points('4.3.1', '4.3.2', '4.3.3', '4.3.4', '4.3.5')
class NmfTester(ApproxTest):

	def test_nmf(self):
		D = np.loadtxt('test/toy.txt')
		nmf = load('src.nmf', 'nmf')
		W = np.array([
			[0.00935861, 0.49905781],
			[0.11338369, 0.04997402],
			[0.68540759, 0.48698807],
			[0.89765723, 0.64745207],
			[0.89696312, 0.72113493]])
		H = np.array([[0.83135342, 0.82756807, 0.83357958], [0.95704434, 0.36804444, 0.49483763]])

		W, H, err = nmf(D, W, H, 20)

		W1 = np.array([
			[1.30165413e-02, 1.51730811e+00],
			[2.71992608e-01, 1.26344551e+00],
			[7.43108401e-01, 1.45007967e+00],
			[1.01638264e+00, 4.35390265e-09],
			[8.17080346e-01, 4.41655312e-01]])
		H1 = np.array([[0.00814249, 0.75387463, 0.11264019], [0.27948174, 0.02955525, 0.55398798]])
		err1 = np.array([7.38628351,
			1.26212238, 0.83573853, 0.6878311,  0.58126599, 0.50481105, 0.45002629,
			0.41261534, 0.3884316,  0.37304945, 0.3630057,  0.35617026, 0.35135314,
			0.34787016, 0.34529876, 0.34336484, 0.34188605, 0.34073873, 0.33983732,
			0.33912138, 0.33854737])
		self.assertApprox(W1, W, errormsg("Incorrect W (input toy.txt, k = 2, itercnt = 20)", W1, W))
		self.assertApprox(H1, H, errormsg("Incorrect H (input toy.txt, k = 2, itercnt = 20)", H1, H))
		self.assertApprox(err1, err, errormsg("Incorrect err (input toy.txt, k = 2, itercnt = 20)", err1, err))
