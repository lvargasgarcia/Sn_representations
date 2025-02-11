import torch
import snob as Snob2
import sys
import numpy as np
from irrep import Irrep
from itertools import permutations
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor
import time
import math
from fractions import Fraction
from numba import jit, prange
from multiprocessing import Pool, cpu_count

def multiplicar_fila(i, A, B, cols_B):
    """Multiplica la fila i de A con toda la matriz B"""
    return [sum(A[i, k] * B[k, j] for k in range(A.shape[1])) for j in range(cols_B)]

def multiplicar_matrices_fraccionadas(A, B):
    filas_A, cols_A = A.shape
    filas_B, cols_B = B.shape

    if cols_A != filas_B:
        raise ValueError("Las matrices no son multiplicables")

    C = np.empty((filas_A, cols_B), dtype=object)

    with ThreadPoolExecutor() as executor:
        resultados = list(executor.map(lambda i: multiplicar_fila(i, A, B, cols_B), range(filas_A)))

    # Convertir resultados en matriz
    C[:] = np.array(resultados, dtype=object)

    return C


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

def compose(p1,p2):
    return [p1[p2[i]-1] for i in range(0,len(p1))]

def inverse(p):
    return [p.index(i)+1 for i in range(1,len(p) + 1)]

def williamsCondition(p,n):
    i = p.index(n)
    r = p[(i % (n-1)) + 1]
    return r % (n-1) + 1 == p[0]

class FourierTransform:
    
    def __init__(self, n , f, mode="YKR"):
        self.n = n
        self.nfact = math.factorial(n)
        self.f = f
        self.mode = mode
        partitions = generate_partitions(self.n)

        with ProcessPoolExecutor() as executor:
            self.images = dict(executor.map(self._buildFT, partitions))
        
        self.irreps = {tuple(partition): Irrep(Snob2.IntegerPartition(partition), mode=self.mode) for partition in partitions}

    def _buildFT_williams_YOR(self, irrep):

        n = self.n
        tau = [2,1] + [i for i in range(3,n+1)] # transposición (1,2)
        sigma = [(i+1) for i in range(1,n)] + [1] # Ciclo completo
        q = [(n-i+1) for i in range(1,n+1)] # El 1 va al n, el 2 al n-1, etc.
        qtau = compose(q, tau)
        qsigma = compose(q, sigma)
        qsigmatau = compose(qsigma, tau)
        invSigma = inverse(sigma)

        p = compose(qsigma, tau)

        tau_matrix = irrep.evaluate(Snob2.SnElement(tau))
        qsigmatau_matrix = irrep.evaluate(Snob2.SnElement(qsigmatau))
        invSigma_matrix = irrep.evaluate(Snob2.SnElement(invSigma))

        p_matrix = qsigmatau_matrix
        fourier_matrix = self.f(tuple(p))*p_matrix

        while(p != qtau):

            if(p != qsigmatau):
                if(williamsCondition(p,n) and p != qsigma):
                    p = compose(p, tau)
                    p_matrix = p_matrix @ tau_matrix
                    fourier_matrix += self.f(tuple(p))*p_matrix
                else:
                    p = compose(p, invSigma)
                    p_matrix = p_matrix @ invSigma_matrix
                    fourier_matrix += self.f(tuple(p))*p_matrix
            else:
                p = compose(p, invSigma)
                p_matrix = p_matrix @ invSigma_matrix
                fourier_matrix += self.f(tuple(p))*p_matrix
        
        return fourier_matrix

    def _buildFT_williams(self, irrep):

        n = self.n
        tau = [2,1] + [i for i in range(3,n+1)] # transposición (1,2)
        sigma = [(i+1) for i in range(1,n)] + [1] # Ciclo completo
        q = [(n-i+1) for i in range(1,n+1)] # El 1 va al n, el 2 al n-1, etc.
        qtau = compose(q, tau)
        qsigma = compose(q, sigma)
        qsigmatau = compose(qsigma, tau)
        invSigma = inverse(sigma)

        p = compose(qsigma, tau)

        tau_matrix = irrep.evaluate(Snob2.SnElement(tau))
        qsigmatau_matrix = irrep.evaluate(Snob2.SnElement(qsigmatau))
        invSigma_matrix = irrep.evaluate(Snob2.SnElement(invSigma))
        p_matrix = qsigmatau_matrix


        fourier_matrix = Fraction(self.f(tuple(p))) * p_matrix

        while(p != qtau):

            if(p != qsigmatau):
                if(williamsCondition(p,n) and p != qsigma):
                    p = compose(p, tau)
                    p_matrix = multiplicar_matrices_fraccionadas(p_matrix, tau_matrix)
                    f_p = Fraction(self.f(tuple(p)))
                    fourier_matrix += f_p * p_matrix
                else:
                    p = compose(p, invSigma)
                    p_matrix = multiplicar_matrices_fraccionadas(p_matrix, invSigma_matrix)
                    f_p = Fraction(self.f(tuple(p)))
                    fourier_matrix += f_p * p_matrix
            else:
                p = compose(p, invSigma)
                p_matrix = multiplicar_matrices_fraccionadas(p_matrix, invSigma_matrix)
                f_p = Fraction(self.f(tuple(p)))
                fourier_matrix += f_p * p_matrix
        
        return fourier_matrix
    
    def _buildFT(self, partition):
        irrep = Irrep(Snob2.IntegerPartition(partition), mode=self.mode)
        t_0 = time.time()
        matrix = self._buildFT_williams(irrep) if self.mode != "YOR" else self._buildFT_williams_YOR(irrep)
        # matrix = np.eye(irrep.matrices[0].shape[0], dtype=object) if self.mode != "YOR" else np.eye(irrep.matrices[0].shape[0])
        # # k = 1
        # for pi in permutations([i for i in range(1, self.n + 1)]):
        #     matrix += (self.f(pi) if self.mode == "YOR" else Fraction(self.f(pi)))*irrep.evaluate(Snob2.SnElement(pi))
        #     # k += 1
        #     # if k > 100:
        #     #     break

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
        return f/self.nfact
