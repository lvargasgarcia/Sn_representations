import snob as Snob2
import sys
import numpy as np
import math
from fractions import Fraction
from permutaciones import *
from permutaciones import get_list_from_SnElement
from snob_base import SnFunction

# f1(sigma) = sum_{i=1}^{n} i*sigma(i)
f = SnFunction.gaussian(4)
print(f)

fft = Snob2.ClausenFFT(4)
print(fft(f))