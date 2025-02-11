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

n = 5
dict_gaussian = {tuple(pi): random.normalvariate(0, 1) for pi in itertools.permutations([i for i in range(1, n + 1)])}

def gaussian(pi):
    return dict_gaussian[tuple(pi)]

mode="YKR"

for pi in itertools.permutations([i for i in range(1, n + 1)]):
    print(pi,":", gaussian(pi))

ft = FourierTransform(n,gaussian,mode=mode)

errores = []
for pi in itertools.permutations([i for i in range(1, n + 1)]):
    gauss = gaussian(pi)
    inv_ft = ft.inverseFourierTransform(pi)
    if mode != "YOR":
        inv_ft = inv_ft.numerator/inv_ft.denominator
    print(pi,":", gauss)
    print("---------------------------------------------")
    print(pi,":", inv_ft)
    print("----------------------------------------------")
    print("----------------------------------------------")
    errores.append(abs(gauss - inv_ft))

print("Error máximo:", max(errores))