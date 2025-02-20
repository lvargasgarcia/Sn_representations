from queue import PriorityQueue
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
from gmpy2 import mpq
from williams_interface import build_coefficient
from permutations_iterator import *

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
    
    def __init__(self, n , f, mode="YKR", truncating_order=None):

        self.truncating_order = truncating_order if truncating_order is not None else n
        self.n = n
        self.nfact = math.factorial(n)
        self.f = f
        self.mode = mode
        partitions = generate_partitions(self.n)
        self.images = {}
        inverseFT = np.array([mpq(0) for _ in range(self.nfact)], dtype=object)

        irreps = {tuple(partition): Irrep(Snob2.IntegerPartition(partition), mode=self.mode) for partition in partitions}
        d_lambda = {tuple(partition): irreps[tuple(partition)].matrices[0].shape[0] for partition in partitions}
        (self.f_nums, self.f_dens, self.williams_sequence) = generate_williams_list(f, n)
        print(self.williams_sequence)

        for partition, d in d_lambda.items():
            print(partition, ":", d)

        # Crear una PriorityQueue
        priority_queue = PriorityQueue()

        # Insertar las particiones con su prioridad en orden de d_lambda creciente
        for partition in partitions:
            # Prioridad es inversamente proporcional a d_lambda
            priority = d_lambda[tuple(partition)]
            priority_queue.put((priority, partition))
               
        tau = [2,1] + [i for i in range(3,n+1)] # transposición (1,2)
        sigma = [(i+1) for i in range(1,n)] + [1] # Ciclo completo
        q = [(n-i+1) for i in range(1,n+1)] # El 1 va al n, el 2 al n-1, etc.
        qtau = compose(q, tau)
        qsigma = compose(q, sigma)
        qsigmatau = compose(qsigma, tau)
        invSigma = inverse(sigma)
        
        # # Usar ProcessPoolExecutor con la cola de prioridades
        with ProcessPoolExecutor(max_workers=6) as executor:
            
            futures = []
            
            while not priority_queue.empty():
                _, partition = priority_queue.get()
                ord = self.n - partition[0]
                doInv = 1 if ord <= self.truncating_order else 0
                irrep = irreps[tuple(partition)]
                tau_matrix = np.array([[mpq(elem) for elem in row] for row in irrep.evaluate(Snob2.SnElement(tau))], dtype=object)
                qsigmatau_matrix = np.array([[mpq(elem) for elem in row] for row in irrep.evaluate(Snob2.SnElement(qsigmatau))], dtype=object)
                invSigma_matrix = np.array([[mpq(elem) for elem in row] for row in irrep.evaluate(Snob2.SnElement(invSigma))], dtype=object)
                future = executor.submit(self._buildFT, partition, tau_matrix, qsigmatau_matrix, invSigma_matrix, doInv)
                futures.append(future)
            
            for future in futures:
                partition, fourier_coef = future.result()
                ord = self.n - partition[0]
                doInv = 1 if ord <= self.truncating_order else 0
                if doInv == 0:
                    self.images[partition] = fourier_coef
                else:
                    fourier_coeff, inv_f = fourier_coef
                    self.images[partition] = fourier_coeff
                    inverseFT += inv_f

        # VERSIÓN SIN PARALELIZACIÓN (PARA DEBUGGING)

        # futures = []

        # while not priority_queue.empty():
        #     _, partition = priority_queue.get()
        #     ord = self.n - partition[0]
        #     doInv = 1 if ord <= self.truncating_order else 0
        #     irrep = irreps[tuple(partition)]
        #     tau_matrix = np.array([[mpq(elem) for elem in row] for row in irrep.evaluate(Snob2.SnElement(tau))], dtype=object)
        #     qsigmatau_matrix = np.array([[mpq(elem) for elem in row] for row in irrep.evaluate(Snob2.SnElement(qsigmatau))], dtype=object)
        #     invSigma_matrix = np.array([[mpq(elem) for elem in row] for row in irrep.evaluate(Snob2.SnElement(invSigma))], dtype=object)
        #     future = self._buildFT(partition, tau_matrix, qsigmatau_matrix, invSigma_matrix, doInv)
        #     futures.append(future)
        
        # for future in futures:
        #     partition, fourier_coef = future
        #     ord = self.n - partition[0]
        #     doInv = 1 if ord <= self.truncating_order else 0
        #     if doInv == 0:
        #         self.images[partition] = fourier_coef
        #     else:
        #         fourier_coeff, inv_f = fourier_coef
        #         self.images[partition] = fourier_coeff
        #         inverseFT += inv_f
        
        self.irreps = irreps
        self.d_lambda = d_lambda

        self.inverseFT = {}

        p = compose(qsigma, tau)

        k = 0

        self.inverseFT[tuple(inverse(p))] = inverseFT[k]
        k += 1

        while(p != qtau):
            
            if(p != qsigmatau):
                if(williamsCondition(p,n) and p != qsigma):
                    p = compose(p, tau)

                else:
                    p = compose(p, invSigma)

            else:
                p = compose(p, invSigma)

            self.inverseFT[tuple(inverse(p))] = inverseFT[k]
            k += 1


    def _buildFT_williams_YOR(self, tau_matrix, qsigmatau_matrix, invSigma_matrix):

        n = self.n
        tau = [2,1] + [i for i in range(3,n+1)] # transposición (1,2)
        sigma = [(i+1) for i in range(1,n)] + [1] # Ciclo completo
        q = [(n-i+1) for i in range(1,n+1)] # El 1 va al n, el 2 al n-1, etc.
        qsigma = compose(q, sigma)
        qtau = compose(q, tau)
        qsigmatau = compose(qsigma, tau)
        invSigma = inverse(sigma)

        p = compose(qsigma, tau)

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

    def _buildFT_williams(self, tau_matrix, qsigmatau_matrix, invSigma_matrix, doInv):

 
    
        # fourier_matrix = mpq(self.f(tuple(p))) * p_matrix

        return build_coefficient(self.n, tau_matrix.shape[0], tau_matrix, qsigmatau_matrix, invSigma_matrix, self.f_nums, self.f_dens, self.williams_sequence, doInv)

        # while(p != qtau):

        #     if(p != qsigmatau):
        #         if(williamsCondition(p,n) and p != qsigma):
        #             p = compose(p, tau)
        #             # p_matrix = multiplicar_matrices_fraccionadas(p_matrix, tau_matrix)
        #             # p_matrix = multiply_mpq_matrices(p_matrix, tau_matrix)
        #             p_matrix = p_matrix @ tau_matrix

        #             f_p = mpq(self.f(tuple(p)))
        #             fourier_matrix += f_p*p_matrix

        #         else:
        #             p = compose(p, invSigma)
        #             # p_matrix = multiplicar_matrices_fraccionadas(p_matrix, invSigma_matrix)
        #             # p_matrix = multiply_mpq_matrices(p_matrix, invSigma_matrix)
        #             p_matrix = p_matrix @ invSigma_matrix

        #             f_p = mpq(self.f(tuple(p)))
        #             fourier_matrix += f_p*p_matrix

        #     else:
        #         p = compose(p, invSigma)
        #         # p_matrix = multiplicar_matrices_fraccionadas(p_matrix, invSigma_matrix)
        #         # p_matrix = multiply_mpq_matrices(p_matrix, invSigma_matrix)
        #         p_matrix = p_matrix @ invSigma_matrix

        #         f_p = mpq(self.f(tuple(p)))
        #         fourier_matrix += f_p*p_matrix


    
    def _buildFT(self, partition, tau_matrix, qsigmatau_matrix, invSigma_matrix, doInv):
        
        t_0 = time.time()
        if self.mode != "YOR":
            matrix = self._buildFT_williams(tau_matrix, qsigmatau_matrix, invSigma_matrix, doInv)
        else:
            matrix = self._buildFT_williams_YOR(tau_matrix, qsigmatau_matrix, invSigma_matrix)

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
            matIrrep = irrep.evaluate(Snob2.SnElement(pi).inv())
            coefFourier = self.images[partition]
            d_lambda = mpq(irrep.matrices[0].shape[0])
            f += d_lambda * (np.trace(coefFourier @ matIrrep))
        return f/self.nfact
