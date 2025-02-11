import itertools
import torch
import snob as Snob2
import sys
import numpy as np
from irrep import Irrep
from itertools import permutations
from concurrent.futures import ProcessPoolExecutor
import time
import math
from fractions import Fraction

n = 8

permutationsTestSet = set()

def compose(p1,p2):
    return [p1[p2[i]-1] for i in range(0,len(p1))]

def inverse(p):
    return [p.index(i)+1 for i in range(1,len(p) + 1)]

def williamsCondition(p,n):
    i = p.index(n)
    r = p[(i % (n-1)) + 1]
    return r % (n-1) + 1 == p[0]

tau = [2,1] + [i for i in range(3,n+1)] # transposición (1,2)
sigma = [(i+1) for i in range(1,n)] + [1] # Ciclo completo
q = [(n-i+1) for i in range(1,n+1)] # El 1 va al n, el 2 al n-1, etc.
qtau = compose(q, tau)
qsigma = compose(q, sigma)
qsigmatau = compose(qsigma, tau)
invSigma = inverse(sigma)

p = compose(qsigma, tau)
permutationsTestSet.add(tuple(p))
print(p)
print(len(permutationsTestSet))

while(p != qtau):
    
    if(p != qsigmatau):
        if(williamsCondition(p,n) and p != qsigma):
            p = compose(p, tau)
        else:
            p = compose(p, invSigma)
    else:
        p = compose(p, invSigma)
    if(permutationsTestSet.__contains__(tuple(p))):
        print("Repetido")
        break
    permutationsTestSet.add(tuple(p))
    print(p)
    print(len(permutationsTestSet))


print("Elementos faltantes:", math.factorial(n) - len(permutationsTestSet))