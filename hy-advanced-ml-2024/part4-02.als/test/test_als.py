#!/usr/bin/env python3

import sys
import unittest
from tmc import points
from tmc.utils import load, get_out

import numpy as np

from .utils import errormsg, ApproxTest


@points('4.2.1', '4.2.2', '4.2.3', '4.2.4', '4.2.5')
class AlsTester(ApproxTest):

	def test_error(self):
		D = np.genfromtxt('test/toy.txt')
		error = load('src.als', 'error')
		W = np.array([
			[0.00935861, 0.49905781],
			[0.11338369, 0.04997402],
			[0.68540759, 0.48698807],
			[0.89765723, 0.64745207],
			[0.89696312, 0.72113493]])
		H = np.array(
			[[-0.23854025,  0.74591618, -0.47180022],
			 [ 0.48097106,  0.05636673,  1.57773169]])
		err = error(D, W, H, 0.1)
		err1 = 2.0650791165932234
		self.assertApprox(err, err1, errormsg("error() returns incorrect value\nX = toy.txt, reg = 0.1\nW =\n%s\nH =\n%s\n" % (str(W), str(H)), err1, err))


	def test_solve(self):
		D = np.genfromtxt('test/toy.txt')
		solve = load('src.als', 'solve')
		W = np.array([
			[0.00935861, 0.49905781],
			[0.11338369, 0.04997402],
			[0.68540759, 0.48698807],
			[0.89765723, 0.64745207],
			[0.89696312, 0.72113493]])
		H = solve(D, W, 0.1)

		H1 = np.array(
			[[-0.05901625,  0.59974042, -0.31824018],
			 [ 0.26972453,  0.16104422,  1.23531541]])
		self.assertApprox(H1, H, errormsg("solve() returns incorrect value\nX = toy.txt, reg = 0.1\nW =\n%s\n" % str(W), H1, H))


	def test_als(self):
		D = np.genfromtxt('test/toy.txt')
		als = load('src.als', 'als')
		W = np.array([
			[0.00935861, 0.49905781],
			[0.11338369, 0.04997402],
			[0.68540759, 0.48698807],
			[0.89765723, 0.64745207],
			[0.89696312, 0.72113493]])

		W, H, err = als(D, W, 0.001, 20)

		W1 = np.array(
			[[ 0.02345106,  0.5884415 ],
			 [ 0.2279559,   1.7483535 ],
			 [ 0.72591418,  0.83793902],
			 [-0.3162981,  -0.08118101],
			 [ 0.79808342,  0.41512675]])
		H1 = np.array(
			[[-0.23854025,  0.74591618, -0.47180022],
			 [ 0.48097106,  0.05636673,  1.57773169]])
		err1 = np.array(
			[0.0175272,  0.01583195, 0.01479687, 0.01409233, 0.01352739, 0.01285707,
			 0.01159863, 0.01076079, 0.0103009,  0.0099924,  0.00977729, 0.00962604,
			 0.0095192,  0.00944311, 0.00938812, 0.00934751, 0.00931662, 0.00929226,
			 0.00927228, 0.00925523])

		self.assertApprox(W1, W, errormsg("Incorrect W (input toy.txt, k = 2, reg = 0.001, itercnt = 20)", W1, W))
		self.assertApprox(H1, H, errormsg("Incorrect H (input toy.txt, k = 2, reg = 0.001, itercnt = 20)", W1, W))
		self.assertApprox(err1, err, errormsg("Incorrect err (input toy.txt, k = 2, reg = 0.001, itercnt = 20)", err1, err))
