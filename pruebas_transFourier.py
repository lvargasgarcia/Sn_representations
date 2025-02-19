import torch
import snob as Snob2
import cnine
import sys
import numpy as np
import math
from fractions import Fraction
from irrep import Irrep
from snob_base import SnFunction
# from fourierTransform import FourierTransform
import itertools
import random
import time
from fourierTransform import FourierTransform
from permutations_iterator import representar_decimal 
from gmpy2 import mpq

random.seed(44)

n = 5

i = 1
dict_gaussian = {}
for pi in itertools.permutations([i for i in range(1, n + 1)]):
    dict_gaussian[tuple(pi)] = round(random.normalvariate(10,20), 10)
    i += 1

def gaussian(pi):
    return dict_gaussian[tuple(pi)]

mode="YSR"

for pi in itertools.permutations([i for i in range(1, n + 1)]):
    print(pi,":", gaussian(pi))

t_0 = time.time()

ft = FourierTransform(n,gaussian,mode=mode)

t_f = time.time()
tiempo_ft = t_f - t_0
with open("tiempo.txt", "w") as f:
    f.write(str(time.time()) +  ":"+ str(t_f - t_0))
errores = []
i = 1
for pi in itertools.permutations([i for i in range(1, n + 1)]):
    gauss = gaussian(pi)
    print(gauss)
    gauss = mpq(representar_decimal(gauss))
    inv_ft = ft.inverseFourierTransform(pi)
    print(pi,":", gauss)
    print("---------------------------------------------")
    print(pi,":", inv_ft)
    print(inv_ft.numerator/inv_ft.denominator)
    print("----------------------------------------------")
    print("----------------------------------------------")
    errores.append(abs(gauss - inv_ft))
    i = i + 1
    if i > 20:
        break

print("error maximo:", max(errores).numerator/max(errores).denominator)
print("Tiempo:", tiempo_ft)
# print(np.iinfo(np.longlong).max)

# m = ft.images[(2,2)]

# print([[elem.numerator/elem.denominator for elem in row] for row in m])

# n = 4

# permutationsTestSet = set()

# def compose(p1,p2):
#     return [p1[p2[i]-1] for i in range(0,len(p1))]

# def inverse(p):
#     return [p.index(i)+1 for i in range(1,len(p) + 1)]

# def williamsCondition(p,n):
#     i = p.index(n)
#     r = p[(i % (n-1)) + 1]
#     return r % (n-1) + 1 == p[0]

# tau = [2,1] + [i for i in range(3,n+1)] # transposición (1,2)
# sigma = [(i+1) for i in range(1,n)] + [1] # Ciclo completo
# q = [(n-i+1) for i in range(1,n+1)] # El 1 va al n, el 2 al n-1, etc.
# qtau = compose(q, tau)
# qsigma = compose(q, sigma)
# qsigmatau = compose(qsigma, tau)
# invSigma = inverse(sigma)

# irrep = Irrep(Snob2.IntegerPartition([2,2]), mode="YKR")
# tau_matrix = np.array([[elem for elem in row] for row in irrep.evaluate(Snob2.SnElement(tau))], dtype=object)
# qsigmatau_matrix = np.array([[elem for elem in row] for row in irrep.evaluate(Snob2.SnElement(qsigmatau))], dtype=object)
# invSigma_matrix = np.array([[elem for elem in row] for row in irrep.evaluate(Snob2.SnElement(invSigma))], dtype=object)

# print("Tau")
# print(tau_matrix)
# print("QsigmaTau")
# print(qsigmatau_matrix)
# print("InvSigma")
# print(invSigma_matrix)

# p = compose(qsigma, tau)
# permutationsTestSet.add(tuple(p))
# print(p)
# print(len(permutationsTestSet))

# valores = "{{" + str(f"{Fraction(representar_decimal(dict_gaussian[tuple(p)])).numerator},") + str(f"{Fraction(representar_decimal(dict_gaussian[tuple(p)])).denominator}") + "}," 

# while(p != qtau):
    
#     if(p != qsigmatau):
#         if(williamsCondition(p,n) and p != qsigma):
#             p = compose(p, tau)
#         else:
#             p = compose(p, invSigma)
#     else:
#         p = compose(p, invSigma)
#     if(permutationsTestSet.__contains__(tuple(p))):
#         print("Repetido")
#         break
#     permutationsTestSet.add(tuple(p))
#     print(p)
#     valores += "{" + str(f"{Fraction(representar_decimal(dict_gaussian[tuple(p)])).numerator},") + str(f"{Fraction(representar_decimal(dict_gaussian[tuple(p)])).denominator}") + "},"
#     print(len(permutationsTestSet))


# print("Elementos faltantes:", math.factorial(n) - len(permutationsTestSet))
# print(valores)