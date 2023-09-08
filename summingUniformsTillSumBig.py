# The expected number of random uniform variables needed to sum to a value
# above x is exp(x) if x is in [0,1]

import numpy as np
import matplotlib.pyplot as plt

def simNumUniformsTillSumBig(bound):
	s = 0
	ret = 0
	while s <= bound:
		s += np.random.uniform()
		ret += 1
	return ret

def simExpectedNumUniformsTillSumBig(bound, rounds):
	ret = 0
	for _ in range(rounds):
		ret += simNumUniformsTillSumBig(bound)
	return ret/rounds

def test():
  numPts = 75
  numRounds = 5000
  xx = np.linspace(0,1,numPts)
  yy = [simExpectedNumUniformsTillSumBig(x, numRounds) for x in xx]
  plt.plot(xx,np.exp(xx),'b')
  plt.plot(xx,yy,'r')
  
  plt.show()
