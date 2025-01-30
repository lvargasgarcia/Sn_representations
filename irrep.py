import torch
import snob as Snob2
import sys
import numpy as np
import math
from fractions import Fraction

# Aquí se implementará el objeto que representa una representación irreducible de una partición alpha

# El objeto consisirá en un array de n elementos, (alpha es una partición de n)

# Cada elemento será una matriz comprimida que contenga la info de rho^alpha(t_i), estos elementos pueden estar ya calculados o no

# Se irán calculando según se necesiten

# Además, se implementará una función evaluate que dada una permutación, la descomponga en ciclos adyacentes, tome las matrices que necesite y las multiplique

# Eso nos dará una matriz de la representación irreducible de la partición alpha en cualquier elemento de S_n