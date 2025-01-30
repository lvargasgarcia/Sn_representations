import torch
import snob as Snob2
import cnine
from prueba import *

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

rho = Snob2.SnIrrep([5,1,1,1])

rep_snob = transponer_diagonal_secundaria(rho[Snob2.SnElement([2,1,3,4,5,6,7,8])].torch().tolist()) 
mi_rep = build_irrep_of_transposition(Snob2.IntegerPartition([5,1,1,1]), 2, mode="YOR")
mi_rep = decompress(mi_rep).tolist()

print("------- Resultado de snob -------")
print_matrix(rep_snob)
print("------- Resultado de mi módulo -------")
print_matrix(mi_rep)
print("------- Comparación -------")
print(compare_matrices(rep_snob, mi_rep))


