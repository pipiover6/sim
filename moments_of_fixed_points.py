# Bell numbers as (lower) diagonal moments of the number of fixed points of a random permutation 

import itertools
import numpy as np


def fact(n):
    ret = 1
    for i in range(n):
        ret *= (i+1)
    return ret

def nfixed(pi):
    n = len(pi)
    ret = 0
    for i in range(n):
        if pi[i] == i:
            ret += 1

    return ret

def moment(n,k):
    perms = list(itertools.permutations([i for i in range(n)]))
    ret = 0
    for pi in perms:
        ret += nfixed(pi)**k
    return ret / fact(n)

A = np.zeros((10,10))
for i in range(10):
    for j in range(10):
        A[i][j] = moment(i,j)

np.set_printoptions(suppress=True)

print(A)
