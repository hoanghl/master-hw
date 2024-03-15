import unittest
import numpy as np

def errormsg(msg, x1, x2):
	return "%s:\nexpected:\n%s\nreceived:\n%s" % (msg, str(x1), str(x2))


class ApproxTest(unittest.TestCase):
	def assertApprox(self, a, b, msg=None, tol=0.0001):
		f = None
		try:
			if np.isscalar(a) or np.isscalar(b):
				if not np.isscalar(a) or not np.isscalar(b):
					f = str(a) + ' and ' + str(b) + ' have different shapes'
			elif a.shape != b.shape:
				f = str(a) + ' and ' + str(b) + ' have different shapes'

			if not f and (np.abs(a - b)).max() > tol:
				f = 'Value ' + str(a) + ' is not close to value ' + str(b)
		except:
			f = 'Cannot compare ' + str(a) + ' and ' + str(b) 

		if f:
			if msg:
				f += ': ' + msg
			self.fail(f)
