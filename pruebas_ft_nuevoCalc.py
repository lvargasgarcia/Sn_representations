import torch
import snob as Snob2
import cnine
import sys
import numpy as np
import math
import argparse
from fractions import Fraction
from irrep import Irrep
from snob_base import SnFunction
import itertools
import random
import time
from fourierTransform import FourierTransform
from permutations_iterator import representar_decimal 
from gmpy2 import mpq
import matplotlib.pyplot as plt

n = 8

dict_gaussian = {pi: int(random.uniform(1000,50)) for pi in itertools.permutations([i for i in range(1, n + 1)])}

def gaussian(pi):
    return dict_gaussian[tuple(pi)]

mode="YSR"

t_0 = time.time()
ft = FourierTransform(n, gaussian, mode=mode)
t_f = time.time()
tiempo_ft = t_f - t_0
max_errores = []

for k in range(n):
    errores = []

    for pi in itertools.permutations([i for i in range(1, n + 1)]):
        gauss = gaussian(pi)
        gauss = mpq(representar_decimal(gauss))
        inv_ft = ft.inverseFourierTransform(pi, order=k)
        errores.append(abs(gauss - inv_ft))

    errormax = max(errores).numerator / max(errores).denominator
    max_errores.append(errormax)

max_errores = np.array(max_errores)/max(max_errores)

plt.plot(range(n), max_errores)
plt.grid()
plt.xlabel("Orden de truncamiento")
plt.ylabel("Error máximo")
plt.title("Error máximo vs. Orden de truncamiento")
plt.show()