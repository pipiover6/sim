# the theoretical background can be found in chapter 0 of the book "Probability with Martingales" by Williams

import numpy as np

class Distribution:
	def __init__(self, probabilities, values):
		self.pr = probabilities
		self.val = values

	def gen(self):
		return np.random.choice(self.val, p = self.pr)

	def genFun(self, theta):
		return sum([self.pr[i] * (theta**(self.val[i])) for i in range(len(self.pr))])
		

def branchingSim(distribution, generations):
	fathers = 1
	for _ in range(generations):
		children = 0
		for i in range(fathers):
			children += distribution.gen()
		fathers = children
	return fathers


def deathSim():
	dist = Distribution([1/6, 5/6], [0,2])

	gens = 7	# warning: increasing this value will have an exponential effect on runtime, and the result will barely change
	trials = 10**4
	deaths = 0
	for _ in range(trials):
		if branchingSim(dist,gens) == 0:
			deaths += 1
	theta = deaths/trials
	print(theta)
	print(dist.genFun(theta))
