#!/usr/bin/env python3

import sys
import unittest
from tmc import points
from tmc.utils import load, get_out

import numpy as np

from .utils import errormsg, ApproxTest



@points('1.1.1', '1.1.2', '1.1.3', '1.1.4', '1.1.5')
class GradTester(ApproxTest):

	def test_step(self):
		Gradient = load('src.gradient', 'Gradient')

		g = Gradient('determ', 0.1, 0.9, 3)

		s = g.step_size()
		s1 = 0.1
		self.assertApprox(s1, s, errormsg("Incorrect step size (mode = determ, round = 0, eta=0.1, tau=0.9, k=3)", s1, s))

		g.update_step(np.zeros(3))
		s = g.step_size()
		s1 = 0.09174311926605505
		self.assertApprox(s1, s, errormsg("Incorrect step size (mode = determ, round = 1, eta=0.1, tau=0.9, k=3)", s1, s))

		g.update_step(np.zeros(3))
		s = g.step_size()
		s1 = 0.08474576271186442
		self.assertApprox(s1, s, errormsg("Incorrect step size (mode = determ, round = 2, eta=0.1, tau=0.9, k=3)", s1, s))

		g = Gradient('adagrad', 0.1, 0.9, 3)

		s = g.step_size()
		s1 = 0.1 / 0.9 * np.ones(3)
		self.assertApprox(s1, s, errormsg("Incorrect step size (mode = adagrad, round = 0, eta=0.1, tau=0.9, k=3)", s1, s))

		g.update_step(np.array([0.1, 4, 1]))
		s = g.step_size()
		s1 = np.array([0.1, 0.02040816, 0.05263158])
		self.assertApprox(s1, s, errormsg("Incorrect step size (mode = adagrad, round = 1, eta=0.1, tau=0.9, k=3)", s1, s))

		g.update_step(np.array([-0.4, 2, 0.9]))
		s = g.step_size()
		s1 = np.array([0.07620147, 0.01861457, 0.04453624])
		self.assertApprox(s1, s, errormsg("Incorrect step size (mode = adagrad, round = 2, eta=0.1, tau=0.9, k=3)", s1, s))



	def test_grad(self):
		Gradient = load('src.gradient', 'Gradient')

		D = np.loadtxt('test/toy.txt')
		y = D[:,0]
		D = D[:,1:]

		g = Gradient('adagrad', 0.1, 0.9, 10)

		s = np.zeros(D.shape[1])
		w = g.gradient(D, y, s)
		w1 = np.array([-11.9610176, -0.3533132, -0.9462042, 1.38518016, 5.598737, 0.285209, -0.884304, -3.3476352, -1.426126, -1.000802])
		self.assertApprox(w1, w, errormsg("Incorrect weights (input toy.txt, w = %s)" % s, w1, w))

		s = w1
		w = g.gradient(D, y, s)
		w1 = np.array([-48.83610418,  -2.5810095, -5.81208385, 6.6755347, 16.77514095, 2.23376884,  -5.49589889, -15.15587321,  -7.81765637,  -6.42723113])
		self.assertApprox(w1, w, errormsg("Incorrect weights (input toy.txt, w = %s)" % s, w1, w))


	def test_epoch(self):
		Gradient = load('src.gradient', 'Gradient')

		D = np.loadtxt('test/toy.txt')
		y = D[:,0]
		D = D[:,1:]

		g = Gradient('adagrad', 0.01, 0.9, 10)
		
		s = np.zeros(D.shape[1])
		t = g.epoch(D, y, s, 100)
		t1 = np.array([[0.1329002, 0.0039257, 0.01051338, -0.01539089, -0.06220819, -0.00316899, 0.0098256, 0.03719595, 0.01584584, 0.01112002]])
		self.assertApprox(t1, t, errormsg("Incorrect trajectory (input toy.txt, w = %s, M=100, mode=adagrad)" % s, t1, t))
		s = t[-1,:]
		t = g.epoch(D, y, s, 75)
		t1 = np.array(
			[[ 0.14263503,  0.00168547,  0.01383354, -0.02312775, -0.07168872, -0.00234152,
			   0.01509373,  0.0447623,   0.01945453,  0.01279239],
			 [ 0.14736975,  0.01761528,  0.02227731, -0.02310112, -0.07537029, -0.01357497,
			   0.01737927,  0.0505323,   0.03073495,  0.02714526]])
		self.assertApprox(t1, t, errormsg("Incorrect trajectory (input toy.txt, w = %s, M=75, mode=adagrad)" % s, t1, t))

		g = Gradient('determ', 0.1, 0.9, 10)
		
		s = np.zeros(D.shape[1])
		t = g.epoch(D, y, s, 100)
		t1 = np.array([[1.19610176, 0.03533132,  0.09462042, -0.13851802, -0.5598737,  -0.0285209, 0.0884304, 0.33476352, 0.1426126, 0.1000802]])
		self.assertApprox(t1, t, errormsg("Incorrect trajectory (input toy.txt, w = %s, M=100, mode=determ)" % s, t1, t))
		s = t[-1,:]
		t = g.epoch(D, y, s, 75)
		t1 = np.array(
			[[ 2.02038883,  0.00960182,  0.12263136, -0.24354902, -1.02760878, -0.00740562,
			   0.13120268,  0.53091088,  0.18293202,  0.09916666],
			 [ 2.3845167,   0.08803739,  0.16472002, -0.23304163, -1.19795936, -0.08223335,
			   0.14441995,  0.66707909,  0.28253669,  0.20948627]])
		self.assertApprox(t1, t, errormsg("Incorrect trajectory (input toy.txt, w = %s, M=75, mode=determ)" % s, t1, t))
