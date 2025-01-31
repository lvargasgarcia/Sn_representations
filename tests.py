import torch
import snob as Snob2
import cnine
from irrep import Irrep
from permutaciones import *
from matrix_utils import *

def transponer_diagonal_secundaria(matriz):
    n = len(matriz)
    return [[matriz[n - 1 - j][n - 1 - i] for j in range(n)] for i in range(n)]

def compare_matrices(matrix1, matrix2):
    # Verificar si las matrices tienen el mismo tamaño
    if len(matrix1) != len(matrix2) or len(matrix1[0]) != len(matrix2[0]):
        return False

    # Comparar cada elemento de las matrices
    for i in range(len(matrix1)):
        for j in range(len(matrix1[i])):

            val1 = matrix1[i][j]
            val2 = matrix2[i][j]

            if abs(val1 - val2) > 1e-6:
                return False
            
    return True

def generate_transposition(i, n):
    initial = [i for i in range(1, n+1)]
    initial[i-2], initial[i-1] = initial[i-1], initial[i-2]
    return Snob2.SnElement(initial)

# Comparación de resultados del módulo que he creado con snob2, usaremos el modo "YOR" en mi módulo para comparar resultados

# Para S_n con n en [1..8], probaremos sigma^alpha (t_n) siendo alpha particiones con distinto orden

# rho = Snob2.SnIrrep([4,1,1,1])

# rep_snob = transponer_diagonal_secundaria(rho[Snob2.SnElement([1,2,3,5,4,6,7])].torch().tolist()) 
# mi_rep = build_irrep_of_transposition(Snob2.IntegerPartition([4,1,1,1]), 5, mode="YOR")
# mi_rep = decompress(mi_rep).tolist()

# print("------- Resultado de snob -------")
# print_matrix(rep_snob)
# print("------- Resultado de mi módulo -------")
# print_matrix(mi_rep)
# print("------- Comparación -------")
# print(compare_matrices(rep_snob, mi_rep))

# Pruebas con S_4 
partition = [2,1,1]
pi = Snob2.SnElement([2,1,4,3])
rho = Snob2.SnIrrep(partition)
mi_rho = Irrep(Snob2.IntegerPartition(partition), mode="YOR")
rep_snob = transponer_diagonal_secundaria(rho[pi].torch().tolist())
print("------ Matriz de representación de Snob ------")
print_matrix(rep_snob)
print("------ Descomposición de pi en producto de transposiciones adyacenytes ------")
print(express_into_adyacent_transpositions(pi))
print("------ Matrices de representación de Snob de t_2 y t_4")
ms_2 = np.array(transponer_diagonal_secundaria(rho[Snob2.SnElement([2,1,3,4])].torch().tolist()))
print("Matriz de t2:")
print_matrix(ms_2)
ms_4 = np.array(transponer_diagonal_secundaria(rho[Snob2.SnElement([1,2,4,3])].torch().tolist()))
print("Matriz de t4:")
print_matrix(ms_4)
print("------ Mis matrices de representación de t2 y t4 ------")
m_2 = mi_rho.evaluate(Snob2.SnElement([2,1,3,4]))
print("Mi matriz de t2:")
print_matrix(m_2)
print("Mi matriz de t4:")
m_4 = mi_rho.evaluate(Snob2.SnElement([1,2,4,3]))
print_matrix(m_4)
print("Matriz final de Snob:")
print_matrix(ms_2@ms_4)
print("Mi matriz final:")
print_matrix(mi_rho.evaluate(pi))





