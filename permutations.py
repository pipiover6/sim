# we represent permutations as arrays of length n with elements [0,..,n-1] in some order

import numpy as np
import matplotlib.pyplot as plt

def getRandomPerm(n):
	arr = [i for i in range(n)]
	for j in range(n-1,0,-1):
		i = np.random.randint(j+1)
		arr[i], arr[j] = arr[j], arr[i]
	return arr

def getPermCycles(f):
	n = len(f)
	remaining = [i for i in range(n)]
	res = []
	while(remaining):
		a0 = remaining.pop()
		cyc = [a0]
		a = f[a0]
		while(a != a0):
			remaining.remove(a)
			cyc.append(a)
			a = f[a]
		res.append(cyc)
	return res

def insertRSK(decks, entry):
	index = 0
	while True:
		if (index == len(decks)):
			decks.append([entry])
			break

		if (entry > decks[index][-1]):
			decks[index].append(entry)
			break

		j = min([i for i in range(len(decks[index])) if decks[index][i] > entry])
		decks[index][j], entry = entry, decks[index][j]
		index += 1

	return decks

def plancherelShape(f):
	n = len(f)
	decks = []
	for i in range(n):
		decks = insertRSK(decks, f[i])
	return [len(d) for d in decks]

def plotRandomPlancherelShape(N=3600):
	y = np.asarray(plancherelShape(getRandomPerm(N)) / np.sqrt(N))
	x = [i/np.sqrt(N) for i in range(1,len(y)+1)]
	plt.plot(x, y, 'o')
	plt.show()

	
