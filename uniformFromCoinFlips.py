# Generating a uniform number in [1,...,n] using a finite (logarithmic) number of tosses of a fair coin and a coin of bias 1/n.

# Algorithm taken from the paper "On Dice and Coins: Models of Computation for Random Generation" by Feldman, Impagliazzo, Naor, Nisan, Rudich, and Shamir.

import numpy as np
import math

# returns True with probability 1/sides
def coinFlip(sides):
	r = np.random.uniform()
	return r * sides < 1

# return uniform integer in [0,2**nbits)
def uniformBits(nbits):
	rnd = 0
	for _ in range(nbits):
		rnd = rnd * 2 + int(coinFlip(2))
	return rnd	

def uniformFromCoins(n):
	assert(n>1)
	k = math.ceil(2 * np.log2(n))	# this is the number of times we will toss the fair coin
	K = 2**k
	Q = math.floor(K / (n-1))
	R = K % (n-1)
	assert(Q * (n-1) + R == K)
	assert(0 <= R < n)
	
	rnd = uniformBits(k)	# random number in [0,K)
	
	if coinFlip(n):		# we enter this branch with probability 1/n
		if rnd < R * (n-1):
			return 1 + (rnd % (n-1))
		return n

	else:			# we enter this branch with probability 1 - 1/n
		if rnd < Q * (n-1):
			return 1 + (rnd % (n-1))
		return n

def test(n, rounds=50000):
	pr = np.asarray([0 for _ in range(n)])
	s=0
	ss=0
	for _ in range(rounds):
		rnd = uniformFromCoins(n)
		s += rnd
		ss += rnd**2
		pr[rnd-1] += 1
	pr = pr/rounds
	s = s/rounds
	ss = ss/rounds
	

	print(f"{s} {ss} {pr}")
	print(f"expected {(n+1)/2} {(n+1)*(2*n+1)/6} {[1/n for _ in range(n)]}")
