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


# returns the number of partitions of [k] into at most n sets, using precomputed Bell values B
def get_bell(n,k,B):
    if k==0:
        ret = 1
    elif k > 0 and n == 0:
        ret = 0
    else:
        ret = 0
        for j in range(k):
            ret += (B[n-1][k-1-j] * fact(k-1) / (fact(j)*fact(k-1-j)))

    B[n][k] = ret
    return ret


B = np.zeros((10,10))

for n in range(10):
    for k in range(10):
        get_bell(n,k,B)

print(B)
print((A==B).all())
