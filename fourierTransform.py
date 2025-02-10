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

def generate_partitions(n):
    if n == 1:
        return [[1]]
    else:
        resp = [[n]]
        for i in range(n-1,0, -1):
            for part in generate_partitions(n-i):
                resp.append([i]+part)
        resp = [sorted(part, reverse=True) for part in resp]
        appeared = set()
        i = 0
        while i < len(resp):
            if tuple(resp[i]) in appeared:
                resp.pop(i)
            else:
                appeared.add(tuple(resp[i]))
                i += 1
        return resp


class FourierTransform:
    
    def __init__(self, n , f, mode="YKR"):
        self.n = n
        self.nfact = math.factorial(n)
        self.f = f
        self.mode = mode
        partitions = generate_partitions(self.n)

        with ProcessPoolExecutor(max_workers=4) as executor:
            self.images = dict(executor.map(self._buildFT, partitions))
        
        self.irreps = {tuple(partition): Irrep(Snob2.IntegerPartition(partition), mode=self.mode) for partition in partitions}

    def _buildFT(self, partition):
        t_0 = time.time()
        irrep = Irrep(Snob2.IntegerPartition(partition), mode=self.mode)
        print("Irrep para particion", partition, "construido en", time.time() - t_0, "s")
        matrix = np.eye(irrep.matrices[0].shape[0], dtype=object) if self.mode != "YOR" else np.eye(irrep.matrices[0].shape[0])
        # k = 1
        for pi in permutations([i for i in range(1, self.n + 1)]):
            matrix += (self.f(pi) if self.mode == "YOR" else Fraction(self.f(pi)))*irrep.evaluate(Snob2.SnElement(pi))
            # k += 1
            # if k > 100:
            #     break

        print("Construida particion: ", partition, "en", time.time() - t_0, "s")
        
        return tuple(partition), matrix
    
    def _evaluate(self, partition):
        return self.images[tuple(partition)]
    
    def __str__(self):
        resp = "---- fourier transform ----\n"
        for key in self.images:
            resp += str(key) + "\n"
            resp += "-------------------\n"
            resp += str(self.images[key]) + "\n"
        return resp

    def inverseFourierTransform(self, pi):
        f = 0
        for partition in self.irreps.keys():
            irrep = self.irreps[partition]
            d_lambda = irrep.matrices[0].shape[0]
            f += d_lambda * (np.trace(self.images[partition] @ irrep.evaluate(Snob2.SnElement(pi).inv())))
        return f/self.nfact - (1 if pi == tuple([i for i in range(1, self.n + 1)]) else 0)
