import numpy as np
from collections import namedtuple
import math
prec = 53 #  double precision
eps = np.finfo(np.float64).eps # machine epsilon

thickening_constant = pow(2, -50) # depends on machine precision
error_constant = pow(2,-50) # compensation for rounding errors in computing f

p0 = 0.4
a = 3.6

class Interval(object):
	
	def __init__(self, left, right):
		self.left = left
		self.right = right

	def contains(self, other):
		return other.left >= self.left and other.right <= self.right

	def thicken(self, amount):
		return Interval(self.left-amount, self.right+amount)

	def __str__(self):
		return "["+str(self.left)+", "+str(self.right)+"]"


def f(x):
	if type(x) is Interval:
		return Interval(min(f(x.right), f(x.left)), max(f(x.right), f(x.left)))
	else:
		return a*x*(1-x)


def f_inv(x, side):
	if type(x) is Interval:
		if(side == 'left'):
			return Interval(0.5-math.sqrt(0.25-x.left/a), 0.5-math.sqrt(0.25-x.right/a))
		if(side == 'right'):
			return Interval(0.5+math.sqrt(0.25-x.right/a), 0.5+math.sqrt(0.25-x.left/a))


def prev_interval(I, side):
	I_thick = I.thicken(thickening_constant)
	I_prev = f_inv(I_thick, side) 
	while not f(I_prev).contains(I_thick):
		I_prev = I_prev.thicken(thickening_constant)
	I_prev.thicken(error_constant)
	return I_prev


def pseudo_orbit(N):
	ans = [p0]
	for i in range(N):
		ans.append(f(ans[-1]))
	return ans


def shadowing_bound(N):
	ps = pseudo_orbit(N)
	I = Interval(ps[-1], ps[-1]) # I_n = [p_n, p_n]
	beta = 0 # max distance between p_n and the shadow orbit
	for i in range(N-1,-1,-1):
		side = 'left' if ps[i] <= 0.5 else 'right'
		I = prev_interval(I, side)
		d = max(abs(ps[i]-I.left), abs(ps[i]-I.right)) # distance bound of p_n and shadow
		beta = max(beta, d)
	return beta

print shadowing_bound(1000000)







