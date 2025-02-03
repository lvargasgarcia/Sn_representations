import torch
import snob as Snob2
import cnine
from irrep import Irrep
from permutaciones import *
from matrix_utils import *



def compare_matrices(matrix1, matrix2):
    # Verificar si las matrices tienen el mismo tamaño
    if len(matrix1) != len(matrix2) or len(matrix1[0]) != len(matrix2[0]):
        return False
    maxerror = 0
    a = -1
    b = -1
    # Comparar cada elemento de las matrices
    for i in range(len(matrix1)):
        for j in range(len(matrix1[i])):

            val1 = matrix1[i][j]
            val2 = matrix2[i][j]

            if abs(val1 - val2) > maxerror:
                maxerror = abs(val1 - val2)
                a = i
                b = j
            
    print("Valores del error")
    print(matrix1[a][b], matrix2[a][b])
    return maxerror

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

# # Pruebas con S_4 
# partition = [2,1,1]
# pi = Snob2.SnElement([2,1,4,3])
# rho = Snob2.SnIrrep(partition)
# mi_rho = Irrep(Snob2.IntegerPartition(partition), mode="YOR")
# rep_snob = transponer_diagonal_secundaria(rho[pi].torch().tolist())
# print("------ Matriz de representación de Snob ------")
# print_matrix(rep_snob)
# print("------ Descomposición de pi en producto de transposiciones adyacenytes ------")
# print(express_into_adyacent_transpositions(pi))
# print("------ Matrices de representación de Snob de t_2 y t_4")
# ms_2 = np.array(transponer_diagonal_secundaria(rho[Snob2.SnElement([2,1,3,4])].torch().tolist()))
# print("Matriz de t2:")
# print_matrix(ms_2)
# ms_4 = np.array(transponer_diagonal_secundaria(rho[Snob2.SnElement([1,2,4,3])].torch().tolist()))
# print("Matriz de t4:")
# print_matrix(ms_4)
# print("------ Mis matrices de representación de t2 y t4 ------")
# m_2 = mi_rho.evaluate(Snob2.SnElement([2,1,3,4]))
# print("Mi matriz de t2:")
# print_matrix(m_2)
# print("Mi matriz de t4:")
# m_4 = mi_rho.evaluate(Snob2.SnElement([1,2,4,3]))
# print_matrix(m_4)
# print("Matriz final de Snob:")
# print_matrix(ms_2@ms_4)
# print("Mi matriz final:")
# print_matrix(mi_rho.evaluate(pi))

# #Prueba con S_8

# partition = [5,2,1]

# error = []

# # S_4

# partition = [2,1,1]
# G = Snob2.Sn(4)
# rho = Snob2.SnIrrep(partition)
# mi_rho = Irrep(Snob2.IntegerPartition(partition), mode="YOR")

# error = []

# for i in range(10):
    
#     pi = G[i]
    
#     if str(pi) == "[ 1 2 3 4 ]":
#         continue
    
#     print("Permutación:", pi)
    
#     print("Descomposicion en transposiciones adyacentes de pi:")
#     print(express_into_adyacent_transpositions(pi))

#     rep_snob = transponer_diagonal_secundaria(rho[pi].torch().tolist())
#     print("------ Matriz de representación de Snob ------")
#     print_matrix(rep_snob)

#     mi_rep = mi_rho.evaluate(pi)
#     print("------ Matriz de representación de mi módulo ------")
#     print_matrix(mi_rep)

#     error.append(compare_matrices(rep_snob, mi_rep))

#     print("------ Comparación ------")
#     print("Maxerror", error[len(error)-1])

# print("Errores")
# print(error)

# print("Error máximo")

# print(max(error))


# print("------ Comparación ------")
# print("Maxerror", compare_matrices(rep_snob, mi_rep))

# rep_snob_2 = transponer_diagonal_secundaria(rho[2,1,3,4,5,6,7,8].torch().tolist())
# mi_rep_2 = mi_rho.evaluate(Snob2.SnElement([2,1,3,4,5,6,7,8]))
# print("------ Comparación de t2 ------")
# print("Maxerror", compare_matrices(rep_snob_2, mi_rep_2))
# rep_snob_3 = transponer_diagonal_secundaria(rho[1,3,2,4,5,6,7,8].torch().tolist())
# mi_rep_3 = mi_rho.evaluate(Snob2.SnElement([1,3,2,4,5,6,7,8]))
# print("------ Comparación de t3 ------")
# print("Maxerror", compare_matrices(rep_snob_3, mi_rep_3))
# rep_snob_4 = transponer_diagonal_secundaria(rho[1,2,4,3,5,6,7,8].torch().tolist())
# mi_rep_4 = mi_rho.evaluate(Snob2.SnElement([1,2,4,3,5,6,7,8]))
# print("------ Comparación de t4 ------")
# print("Maxerror", compare_matrices(rep_snob_4, mi_rep_4))
# rep_snob_8 = transponer_diagonal_secundaria(rho[1,2,3,4,5,6,8,7].torch().tolist())
# mi_rep_8 = mi_rho.evaluate(Snob2.SnElement([1,2,3,4,5,6,8,7]))
# print("------ Comparación de t8 ------")
# print("Maxerror", compare_matrices(rep_snob_8, mi_rep_8))


# rep_snob_2 = np.array(rep_snob_2)
# rep_snob_3 = np.array(rep_snob_3)
# rep_snob_4 = np.array(rep_snob_4)
# rep_snob_8 = np.array(rep_snob_8)


# mat_final_snob = rep_snob_4@rep_snob_3@rep_snob_2@rep_snob_2@rep_snob_3@rep_snob_2@rep_snob_8
# mi_mat_final = mi_rep_4@mi_rep_3@mi_rep_2@mi_rep_2@mi_rep_3@mi_rep_2@mi_rep_8

# print("------ Comparación de la matriz final ------")
# print("Maxerror", compare_matrices(mat_final_snob, mi_mat_final))

# print("------ Comparación de la matriz final con la matriz de Snob ------")
# print("Maxerror", compare_matrices(rep_snob, mat_final_snob))

# print("------ Comparación de la matriz final con la matriz de mi módulo ------")
# print("Maxerror", compare_matrices(mi_rep, mi_mat_final))

partition = [2,1,1]
# pi = [2,3,4,1]
# pi = Snob2.SnElement(pi)

# rho = Snob2.SnIrrep(partition)
# mi_rho = Irrep(Snob2.IntegerPartition(partition), mode="YOR")

# print("Permutación:", pi)

# print("Descomposicion en transposiciones adyacentes de pi:")
# print(express_into_adyacent_transpositions(pi))

# rep_snob = (rho[pi].torch().tolist())
# print("------ Matriz de representación de Snob ------")
# print_matrix(rep_snob)

# mi_rep = mi_rho.evaluate(pi)

# print("------ Matriz de representación de mi módulo ------")
# print_matrix(apply_correction_transformation(mi_rep))

# print("------ Comparación ------")
# print("Maxerror", compare_matrices(rep_snob, mi_rep))

# mi_rep_2 = mi_rho.evaluate(Snob2.SnElement([2,1,3,4]))
# mi_rep_3 = mi_rho.evaluate(Snob2.SnElement([1,3,2,4]))
# mi_rep_4 = mi_rho.evaluate(Snob2.SnElement([1,2,4,3]))

# print_matrix(apply_correction_transformation(mi_rep_2@mi_rep_3@mi_rep_4))

# print("Maxerror", compare_matrices(rep_snob, mi_rep_4@mi_rep_3@mi_rep_2@mi_rep_3@mi_rep_2))

print("----- Pruebas en transposiciones -----")

error = []

rho = Snob2.SnIrrep(partition)
mi_rho = Irrep(Snob2.IntegerPartition(partition), mode="YOR")

for n in range(2,5):
    
    pi = [i for i in range(1, 5)]
    pi[n-2], pi[n-1] = pi[n-1], pi[n-2]
    print(pi)

    pi = Snob2.SnElement(pi)

    print("Descomposicion en transposiciones adyacentes de pi:")
    print(express_into_adyacent_transpositions(pi))

    rep_snob = rho[pi].torch().tolist()
    print("------ Matriz de representación de Snob ------")
    print_matrix(rep_snob)

    mi_rep = mi_rho.evaluate(pi)
    print("------ Matriz de representación de mi módulo ------")
    print_matrix(mi_rep)

    error.append(compare_matrices(rep_snob, mi_rep))

    print("------ Comparación ------")
    print("Maxerror", error[len(error)-1])

print("Errores")
print(error)

s_2 = np.array(rho[2,1,3,4].torch().tolist())
s_3 = np.array(rho[1,3,2,4].torch().tolist())
s_4 = np.array(rho[1,2,4,3].torch().tolist())

t_2 = mi_rho.evaluate(Snob2.SnElement([2,1,3,4]))
t_3 = mi_rho.evaluate(Snob2.SnElement([1,3,2,4]))
t_4 = mi_rho.evaluate(Snob2.SnElement([1,2,4,3]))

print(compare_matrices(s_2, t_2))
print(compare_matrices(s_3, t_3))
print(compare_matrices(s_4, t_4))

print("---------------------------------")
print_matrix(rho[2,3,4,1].torch().tolist())
print("-----------------------")
print_matrix(s_2@s_3@s_4)
print("-----------------------")
print_matrix(t_2@t_3@t_4)