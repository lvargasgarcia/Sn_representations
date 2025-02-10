import torch
import snob as Snob2
import sys
import numpy as np
from irrep import Irrep
from itertools import permutations
from concurrent.futures import ProcessPoolExecutor
import time

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
        self.f = f
        self.mode = mode
        #self.irreps = {tuple(partition): Irrep(Snob2.IntegerPartition(partition), mode=self.mode) for partition in generate_partitions(self.n)}
        self.irreps = {}

        # Usamos ProcessPoolExecutor para paralelizar la construcción de las matrices
        with ProcessPoolExecutor(max_workers=4) as executor:
            # Mapear cada partición a la función _buildFT
            self.images = dict(executor.map(self._buildFT, generate_partitions(self.n)))

    def _buildFT(self, partition):
        t_0 = time.time()
        irrep = Irrep(Snob2.IntegerPartition(partition), mode=self.mode)
        self.irreps[tuple(partition)] = irrep
        print("Irrep para particion", partition, "construido en", time.time() - t_0, "s")
        matrix = np.eye(irrep.matrices[0].shape[0])
        k = 1
        for pi in permutations([i for i in range(1, self.n + 1)]):
            matrix += self.f(pi)*irrep.evaluate(Snob2.SnElement(pi))
            k += 1
            if k > 100:
                break

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
        return None