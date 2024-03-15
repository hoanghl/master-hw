#!/usr/bin/env python3

import sys
import unittest
from tmc import points
from tmc.utils import load, get_out

import numpy as np

from .utils import errormsg, ApproxTest



@points('3.4.1', '3.4.2', '3.4.3', '3.4.4', '3.4.5')
class NnTester(ApproxTest):

	def test_activation(self):
		Neuron = load('src.nn', 'Neuron')
		n = Neuron()
		o = n.act(0)
		o1 = 0.5
		self.assertApprox(o, o1, errormsg('incorrect output for Neuron.act(0)', o1, o)) 
		o = n.act(1)
		o1 = 0.7310585786300049
		self.assertApprox(o, o1, errormsg('incorrect output for Neuron.act(1)', o1, o)) 
		o = n.act(-1)
		o1 = 0.2689414213699951
		self.assertApprox(o, o1, errormsg('incorrect output for Neuron.act(-1)', o1, o)) 

		o = n.der(0)
		o1 = 0.25
		self.assertApprox(o, o1, errormsg('incorrect output for Neuron.der(0)', o1, o)) 

		o = n.der(1)
		o1 = 0.19661193324148185
		self.assertApprox(o, o1, errormsg('incorrect output for Neuron.der(1)', o1, o)) 

		o = n.der(-1)
		o1 = 0.19661193324148185
		self.assertApprox(o, o1, errormsg('incorrect output for Neuron.der(-1)', o1, o)) 


	def test_output(self):
		NN = load('src.nn', 'NN')
		net = NN(3)
		net.join(0, 2)
		net.join(1, 2)
		net.edges[0].w = 1.3
		net.edges[0].w = 4.2
		net.neurons[2].bias = 2.3

		net.neurons[0].output = 0.3
		net.neurons[1].output = -0.2

		net.neurons[2].compute_output()
		o = net.neurons[2].input
		o1 = 3.56
		self.assertApprox(o, o1, errormsg('compute_output does not update Neuron.input correctly\ntwo input neurons with outputs 0.3 and -0.2\nedge weights = 1.3 and 4.2\nbias = 2.3', o1, o))

		o = net.neurons[2].output
		o1 = 0.9723475776771769
		self.assertApprox(o, o1, errormsg('compute_output does not update Neuron.output correctly\ntwo input neurons with outputs 0.3 and -0.2\nedge weights = 1.3 and 4.2\nbias = 2.3', o1, o))

	def build_nn(self):
		NN = load('src.nn', 'NN')
		net = NN(5)

		for i in range(2):
			for j in range(2):
				net.join(i, 2 + j)

		for i in range(2):
			net.join(2 + i, 4)

		np.random.seed(2022)
		net.randomize()
		return net

		
	def test_forward(self):
		net = self.build_nn()

		net.forward(np.array([0.1, 1.3]))
		o = np.zeros(5)
		for i in range(5):
			o[i] = net.neurons[i].output
		o1 = np.array([0.1, 1.3, 0.54794555, 0.94932965, 0.68710632])
		self.assertApprox(o, o1, errormsg('forward does not set correct output for the neurons', o1, o))


	def test_backward(self):
		net = self.build_nn()

		net.forward(np.array([0.1, 1.3]))
		net.backward(np.array(1))

		o = np.zeros(3)
		for i in range(3):
			o[i] = net.neurons[i + 2].delta
		o1 = np.array([-0.01502333, -0.00786805, -0.21499123])

		self.assertApprox(o, o1, errormsg('backward does not set correct delta for neurons', o1, o))


	def test_gradient(self):
		net = self.build_nn()

		X = np.array([[0.1, 1.3], [4, -1], [2, -1.4]])
		y = np.array([1, 0, 1])
		w, b, error = net.gradients(X, y)

		w1 = np.array([0.01062884, 0.00858572, -0.00442676, -0.00293692, -0.04038464, -0.06696546])
		self.assertApprox(w, w1, errormsg('gradient for edge weights is not correct', w1, w))

		b1 = np.array([0.00013137, -0.00077215, -0.00492857, -0.00190211, -0.07169663])
		self.assertApprox(b[2:], b1[2:], errormsg('gradient for biases is not correct', b1[2:], b[2:]))

		err1 = 0.4381352517958994
		self.assertApprox(error, err1, errormsg('error is not correct', err1, error))
