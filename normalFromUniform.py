import numpy as np
import scipy as sc

def normalFromUniform(size=None):
  theta = 2 * np.pi * np.random.uniform(size=size)
  r = np.sqrt(-2 * np.log(np.random.uniform(size=size)))
  return np.multiply(r, np.cos(theta))
	
def test():
	M=8000
	a = np.sort(normalFromUniform(M))
	print(f"mean={np.mean(a)}, std={np.std(a)}")
	b = np.sort(np.random.normal(size=M))
	percentiles = [1,5,10,25,30,40,50]
	for i in percentiles:
		print(f"percentile = {i}, normal, np sim, cur sim = {np.sqrt(2)*(sc.special.erfinv(i/50 - 1))}, {b[M * i // 100]}, {a[M * i // 100]}")
