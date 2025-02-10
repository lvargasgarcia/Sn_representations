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

n = 6
dict_gaussian = {tuple(pi): random.normalvariate(0, 1) for pi in itertools.permutations([i for i in range(1, n + 1)])}

def gaussian(pi):
    return dict_gaussian[tuple(pi)]

# def test_ft_creationtime():

#     n = 9
#     t_0 = time.time()
#     ft = FourierTransform(n, gaussian, mode="YOR")
#     t_f = time.time()
#     print("Tiempo:", t_f - t_0, "s")
#     # print(ft)

# test_ft_creationtime()

# rho = Irrep(Snob2.IntegerPartition([4,2,2]), mode="YOR")
# t_0 = time.time()
# # k = 0
# A = np.eye(rho.matrices[0].shape[0])
# for pi in itertools.permutations([i for i in range(1, n + 1)]):
#     A += gaussian(pi)*rho.evaluate(Snob2.SnElement(pi))
#     print(pi)
#     # k = k + 1
#     # if k > 10000:
#     #     break
# print("Tiempo:", time.time() - t_0, "s")

print("Función original")

for pi in itertools.permutations([i for i in range(1, n + 1)]):
    print(pi,":", gaussian(pi))

ft = FourierTransform(n,gaussian,mode="YSR")

errores = []
for pi in itertools.permutations([i for i in range(1, n + 1)]):
    gauss = gaussian(pi)
    inv_ft = ft.inverseFourierTransform(pi)
    inv_ft = inv_ft.numerator / inv_ft.denominator
    print(pi,":", gauss)
    print("---------------------------------------------")
    print(pi,":", inv_ft)
    print("----------------------------------------------")
    print("----------------------------------------------")
    errores.append(abs(gauss - inv_ft))
    break

print("Error máximo:", max(errores))