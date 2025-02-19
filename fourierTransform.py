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
from permutations_iterator import generate_williams_list

# def multiplicar_fila(i, A, B, cols_B, mcm):
#     """Multiplica la fila i de A con toda la matriz B"""
#     C = np.empty(cols_B, dtype=object)
#     fila = A[i]
#     for j in range(cols_B):
#         sum = 0
#         for k in range(len(fila)):
#             sum += fila[k].numerator * (mcm // fila[k].denominator) * B[k][j].numerator * (mcm // B[k][j].denominator)
#         C[j] = Fraction(sum, mcm**2)
#     return C

def calcular_mcm_matriz(matrix):
    mcm = 1
    for row in matrix:
        for fraction in row:
            den = fraction.denominator
            mcm = mcm * (den // math.gcd(mcm, den))
    return mcm

def mat_mul(A,B):
    
    mcm_A = calcular_mcm_matriz(A)
    mcm_B = calcular_mcm_matriz(B)
    common_mcm = mcm_A * mcm_B // math.gcd(mcm_A, mcm_B)

    num_A = np.array([[elem.numerator for elem in row] for row in A], dtype=int)
    den_A = np.array([[elem.denominator for elem in row] for row in A], dtype=int)
    num_B = np.array([[elem.numerator for elem in row] for row in B], dtype=int)
    den_B = np.array([[elem.denominator for elem in row] for row in B], dtype=int)

    num_A = num_A * (common_mcm // den_A)
    num_B = num_B * (common_mcm // den_B)

    nums = num_A @ num_B

    return np.array([[mpq(nums[i][j], common_mcm**2) for j in range(B.shape[1])] for i in range(A.shape[0])], dtype=object)


# def multiplicar_matrices_fraccionadas(A, B):
#     filas_A, cols_A = A.shape
#     filas_B, cols_B = B.shape

#     if cols_A != filas_B:
#         raise ValueError("Las matrices no son multiplicables")

#     C = np.empty((filas_A, cols_B), dtype=object)

#     mcm_A= calcular_mcm_matriz(A)
#     mcm_B = calcular_mcm_matriz(B)
#     mcm = (mcm_A //math.gcd(mcm_A, mcm_B)) * mcm_B

#     with ThreadPoolExecutor() as executor:
#         resultados = list(executor.map(lambda i: multiplicar_fila(i, A, B, cols_B, mcm), range(filas_A)))

#     # Convertir resultados en matriz
#     C[:] = np.array(resultados, dtype=object)

#     return C

def multiply_mpq_fila(i, A, B, cols_B):
    C = np.empty(cols_B, dtype=object)
    fila = A[i]
    for j in range(cols_B):
        sum = 0
        for k in range(len(fila)):
            sum += fila[k] * B[k][j]
        C[j] = mpq(sum)
    return C

def _worker_multiply_mpq_fila(args):
    i, A, B, cols_B = args
    C = np.empty(cols_B, dtype=object)
    fila = A[i]
    for j in range(cols_B):
        sum = sum(fila[k] * B[k][j] for k in range(len(fila)))
        C[j] = mpq(sum)
    return C

def multiply_mpq_matrices(A, B, threshold=50000):
    filas_A, cols_A = A.shape
    filas_B, cols_B = B.shape

    if cols_A != filas_B:
        raise ValueError("Las matrices no son multiplicables")
    
    # Si la matriz es pequeña, usa un método secuencial más eficiente
    if filas_A * cols_B < threshold:
        return A@B
    
    # Si la matriz es grande, paraleliza
    C = np.empty((filas_A, cols_B), dtype=object)
    args_list = [(i, A, B, cols_B) for i in range(filas_A)]
    
    with ProcessPoolExecutor() as executor:
        resultados = list(executor.map(_worker_multiply_mpq_fila, args_list))
    
    C[:] = np.array(resultados, dtype=object)
    return C


# def multiply_mpq_matrix_scalar(matrix, scalar):
#     with ProcessPoolExecutor() as executor:
#         return np.array(list(executor.map(lambda row: [elem * scalar for elem in row], matrix)), dtype=object)

# def sum_mpq_matrices(A,B):
#     with ProcessPoolExecutor() as executor:
#         return np.array(list(executor.map(lambda i: [elem_A + elem_B for elem_A, elem_B in zip(A[i], B[i])], range(len(A))), dtype=object))


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

        irreps = {tuple(partition): Irrep(Snob2.IntegerPartition(partition), mode=self.mode) for partition in partitions}
        d_lambda = {tuple(partition): irreps[tuple(partition)].matrices[0].shape[0] for partition in partitions}
        (self.f_nums, self.f_dens) = generate_williams_list(f, n)
        print("f_nums:", self.f_nums)
        print("f_dens:", self.f_dens)

        for partition, d in d_lambda.items():
            print(partition, ":", d)

        # Crear una PriorityQueue
        priority_queue = PriorityQueue()

        # Insertar las particiones con su prioridad en orden de d_lambda creciente
        for partition in partitions:
            # Prioridad es inversamente proporcional a d_lambda
            priority = d_lambda[tuple(partition)]
            priority_queue.put((priority, partition))
        
        # Usar ProcessPoolExecutor con la cola de prioridades
        with ProcessPoolExecutor() as executor:
            futures = []
            while not priority_queue.empty():
                _, partition = priority_queue.get()
                future = executor.submit(self._buildFT, partition)
                futures.append(future)
            self.images = dict(future.result() for future in futures)
        
        # self.images = {tuple(partition): self._buildFT(partition) for partition in partitions}
        self.irreps = irreps
        self.d_lambda = d_lambda

        # with ProcessPoolExecutor() as executor:
        #     self.images = dict(executor.map(self._buildFT, partitions))
    
    # def __init__(self, n , f, mode="YKR"):
    #     self.n = n
    #     self.nfact = math.factorial(n)
    #     self.f = f
    #     self.mode = mode
    #     partitions = generate_partitions(self.n)

    #     self.d_lambda = {tuple(partition): Irrep(Snob2.IntegerPartition(partition), self.mode).matrices[0].shape[0] for partition in partitions}

    #     # Crear una PriorityQueue
    #     self.priority_queue = PriorityQueue()

    #     # Insertar las particiones con su prioridad en orden de d_lambda creciente
    #     for partition in partitions:
    #         # Prioridad es inversamente proporcional a d_lambda
    #         priority = self.d_lambda[tuple(partition)]
    #         self.priority_queue.put((priority, partition))
        
    #     # # Usar ProcessPoolExecutor con la cola de prioridades
    #     # with ProcessPoolExecutor() as executor:
    #     #     # # Lanzar los procesos en función de la prioridad
    #     #     # futures = []
    #     #     # while not self.priority_queue.empty():
    #     #     #     _, partition = self.priority_queue.get()
    #     #     #     future = executor.submit(self._buildFT, partition)
    #     #     #     futures.append(future)

    #     #     # Esperar a que todos los procesos terminen y almacenar los resultados
    #     #     self.images = dict(executor.map(self._buildFT, partitions))

    #     with ProcessPoolExecutor() as executor:
    #         self.images = dict(executor.map(self._buildFT, partitions))
    #     # with ProcessPoolExecutor() as executor:
    #     #     self.inverseFT = dict(executor.map(self._evaluate, partitions))
        
        
    #     self.irreps = {tuple(partition): Irrep(Snob2.IntegerPartition(partition), mode=self.mode) for partition in partitions}

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

        tau_matrix = np.array([[mpq(elem) for elem in row] for row in irrep.evaluate(Snob2.SnElement(tau))], dtype=object)
        qsigmatau_matrix = np.array([[mpq(elem) for elem in row] for row in irrep.evaluate(Snob2.SnElement(qsigmatau))], dtype=object)
        invSigma_matrix = np.array([[mpq(elem) for elem in row] for row in irrep.evaluate(Snob2.SnElement(invSigma))], dtype=object)
        p_matrix = qsigmatau_matrix
    
        # fourier_matrix = mpq(self.f(tuple(p))) * p_matrix

        fourier_matrix = build_coefficient(self.n, tau_matrix.shape[0], tau_matrix, qsigmatau_matrix, invSigma_matrix, self.f_nums, self.f_dens)

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

        
        return fourier_matrix
    
    def _buildFT(self, partition):
        
        # print("Construyendo particion: ", partition)
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
            matIrrep = irrep.evaluate(Snob2.SnElement(pi).inv())
            coefFourier = self.images[partition]
            d_lambda = mpq(irrep.matrices[0].shape[0])
            f += d_lambda * (np.trace(coefFourier @ matIrrep))
        return f/self.nfact
