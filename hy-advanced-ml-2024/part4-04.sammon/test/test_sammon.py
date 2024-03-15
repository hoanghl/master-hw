#!/usr/bin/env python3

import sys
import unittest
from tmc import points
from tmc.utils import load, get_out

import numpy as np
from sklearn.metrics import euclidean_distances

from .utils import errormsg, ApproxTest



@points('4.4.1', '4.4.2', '4.4.3', '4.4.4', '4.4.5')
class SammonTester(ApproxTest):

	def test_sammon(self):
		X = np.loadtxt('test/toy.txt')
		sammon = load('src.sammon', 'sammon')

		D = euclidean_distances(X)
		S = np.array([[-0.48326017, -0.1179138 ], [-0.27655539,  0.47968815], [-0.21992991, -0.30314492], [ 0.66658711,  0.07743277], [ 0.31315837, -0.1360622 ]])
		
		P = sammon(D, S, 1, 0.9, 1)
		P1 = np.array([[-0.58360472, -0.0814711 ], [-0.27228416,  0.50859955], [-0.16710242, -0.37715363], [ 0.67869649,  0.0807929 ], [ 0.34429482, -0.13076772]])
		self.assertApprox(P1, P, errormsg("Incorrect output (input toy.txt, P = PCA, eta = 1, tau = 0.9, itercnt = 1)", P1, P))
		
		S = np.array([[-0.48326017, -0.1179138 ], [-0.27655539,  0.47968815], [-0.21992991, -0.30314492], [ 0.66658711,  0.07743277], [ 0.31315837, -0.1360622 ]])
		P = sammon(D, S, 0.1, 0.99, 10)
		P1 = np.array([[-0.53095158, -0.09844152], [-0.27147441,  0.49127237], [-0.19519323, -0.34199692], [ 0.67225833,  0.08046224], [ 0.3253609,  -0.13129617]])
		self.assertApprox(P1, P, errormsg("Incorrect output (input toy.txt, P = PCA, eta = 1, tau = 0.9, itercnt = 10)", P1, P))
