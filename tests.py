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

# Comparación de resultados del módulo que he creado con snob2, usaremos el modo "YOR" en mi módulo para comparar resultados

# Para S_n con n en [1..8], probaremos sigma^alpha (t_n) siendo alpha particiones con distinto orden

rho = Snob2.SnIrrep([4,2,1,1])

rep_snob = transponer_diagonal_secundaria(rho[Snob2.SnElement([1,2,3,4,5,6,8,7])].torch().tolist()) 
mi_rep = build_irrep(Snob2.IntegerPartition([4,2,1,1]), mode="YOR")
mi_rep = decompress(mi_rep).tolist()

print("------- Resultado de snob -------")
print(rep_snob)
print("------- Resultado de mi módulo -------")
print(mi_rep)
print("------- Comparación -------")
print(compare_matrices(rep_snob, mi_rep))

    
