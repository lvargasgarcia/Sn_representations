import torch
import snob as Snob2
import sys
import numpy as np
import math
from fractions import Fraction
from permutaciones import *

def print_matrix(matrix):
    for sublista in matrix:
        print(f"[{' '.join(map(str, sublista))}]")

def search_young_tableau(tableau):
    shape = tableau.shape()
    n = shape.getn()
    n1 = (-1,-1)
    n2 = (-1,-1)
    
    for i in range(n):
        for j in range(shape[i]):
            if tableau[i, j] == n:
                n1 = (i,j)
            if tableau[i, j] == n - 1:
                n2 = (i,j)
            if n1 != (-1,-1) and n2 != (-1,-1):
                return (n1, n2)  # Devolver tuplas

    return None  # Si no se encuentra, devolver None

def n_pos(tableau):
    shape = tableau.shape()
    height = shape.height()
    n = shape.getn()
    for i in range(height):
        for j in range(shape[i]):
            if tableau[i, j] == n:
                n1 = i
                return n1
    return None 

def get_subpartition(key, partition):
    resp = [] 
    for i in range(partition.height()):
          k = partition[i] 
          if key == i:
               k = k - 1
          resp.append(k)
    while resp and resp[-1] == 0:
        resp.pop()

    return Snob2.IntegerPartition(resp)

def branch(partition):
    tableaux_dict = {}
    T = Snob2.StandardYoungTableaux(partition)
    for i in range(len(T)):
        tableau = T[i]
        key = n_pos(tableau)
        if (key != None) and key not in tableaux_dict:
            tableaux_dict[key] = tableau
    
    ordered_1stlev = sorted(tableaux_dict.keys(), key=lambda t: t)
    partitions = []
    
    for key in ordered_1stlev:
        gamma = get_subpartition(key, partition)
        partitions.append(gamma)
    
    return partitions


def get_partition_from_gamma(y_symbol, partition):
    resp = [] 

    for i in range(partition.height()):
          k = partition[i]
          (a, b) = y_symbol # Desempaquetar la tupla
          if a == i:
               k = k - 1
          if b == i:
              k = k - 1
          resp.append(k)
    while resp and resp[-1] == 0:
        resp.pop()

    return Snob2.IntegerPartition(resp)

# ahora vamos a implementar el calculo de la fórmula de frame-robinson-thrall para una particion 

def cells_below_cell(partition, i, j):
    resp = 0
    for k in range(i + 1, partition.height()):
        if(partition[k] >= j + 1):
          resp = resp + 1
        else:
          break
    return resp

def hook_length_prod(partition):
     resp = 1
     for i in range(partition.height()):
          for j in range(partition[i]):
              cells_in_row = (partition[i] - j)
              cells_below = cells_below_cell(partition, i, j)
              hook_length = (cells_in_row + cells_below)
              resp = resp * hook_length
     return resp

def frame_robinson_thrall(partition):
     n = partition.getn()
     n_factorial = math.factorial(n)
     den = hook_length_prod(partition)
     return int(n_factorial / den)

# # Imprimir el diccionario resultante
# for key in ordered_2ndlev:
#     ((a,_),(b,_)) = key
#     print(f"Símbolo de Yamanouchi ({a},{b}):")
#     print("Tablero:")
#     tableau = tableaux_dict[key]
#     print(tableau)
#     print("Gamma asociada al tablero:")
#     gamma = get_partition_from_gamma((a,b), lamb)
#     print(gamma)
#     print("d_gamma")
#     print(frame_robinson_thrall(gamma))
#     print("--------------------")


# Ahora vamos a crear una estructura de datos para cada elemento en la matriz de representación, ya que esta es por bloques
# tiene nxn bloques y su dimension real es d_alpha, donde alpha es la partición dada 

# La matriz sera un np.array de dimension nxn y cada elemento será una dupla ((a,b) --> dimension del bloque, coef --> coeficiente de la identidad)
    
def generate_representation_matrix(partition, ordered_2ndlev, mode="YKR"):
    
    n = len(ordered_2ndlev)
    matrix = [[[-1,-1,None] for _ in range(n)] for _ in range(n)]
    # Primero hacemos la diagonal
    
    for i in range(n):

          ((a,b),(c,d)) = ordered_2ndlev[i]
          gamma = get_partition_from_gamma((a,c), partition)
          d_gamma = frame_robinson_thrall(gamma)
          coef = None
          # Cláusula 3 del teorema de young
          if a == c: # R-chain
               coef = 1
          if b == d: # C-chain
               coef = -1
          matrix[i][i] = [d_gamma, d_gamma, coef]

    for i in range(n):
        for j in range(n):
            
            if i < j:

               ((a,b),(c,d)) = ordered_2ndlev[i]
               ((a2,_),(c2,_)) = ordered_2ndlev[j]
               gamma = get_partition_from_gamma((a,c), partition)
               gamma2 = get_partition_from_gamma((a2,c2), partition)
               d_gamma = matrix[i][i][0]
               d_gamma2 = matrix[j][j][0]
               
               matrix[i][j][0] = d_gamma
               matrix[i][j][1] = d_gamma2

               matrix[j][i][0] = d_gamma2
               matrix[j][i][1] = d_gamma

               if str(gamma) == str(gamma2):
                   
                   chi = abs(a-c) + abs(b-d)
                   chi_sqr = chi**2
                   q_sqr_chi = Fraction((chi_sqr - 1), chi_sqr)
                   
                   if mode == "YKR":
                       if matrix[i][i][2] == None:
                           matrix[i][i][2] = Fraction(1,chi)
                       if matrix[j][j][2] == None:
                           matrix[j][j][2] = Fraction(-1,chi)
                       if matrix[i][j][2] == None:
                           matrix[i][j][2] = Fraction(1)
                       if matrix[j][i][2] == None:
                           matrix[j][i][2] = q_sqr_chi
                    
                   if mode == "YSR":
                         if matrix[i][i][2] == None:
                              matrix[i][i][2] = Fraction(1,chi)
                         if matrix[j][j][2] == None:
                              matrix[j][j][2] = Fraction(-1,chi)
                         if matrix[i][j][2] == None:
                              matrix[i][j][2] = q_sqr_chi
                         if matrix[j][i][2] == None:
                              matrix[j][i][2] = Fraction(1)
                   if mode == "YOR":
                         
                         sqrt = math.sqrt(q_sqr_chi.numerator/q_sqr_chi.denominator)

                         if matrix[i][i][2] == None:
                              matrix[i][i][2] = Fraction(1,chi)
                         if matrix[j][j][2] == None:
                              matrix[j][j][2] = Fraction(-1,chi)
                         if matrix[i][j][2] == None:
                              matrix[i][j][2] = sqrt
                         if matrix[j][i][2] == None:
                              matrix[j][i][2] = sqrt
               else:
                   matrix[i][j][2] = Fraction(0)
                   matrix[j][i][2] = Fraction(0)
                                         

    return matrix

def pretty_print(matrix):
    
    print("Matriz de representación:")
    
    print("---------Dimensiones---------")
    
    dimensions = [[str((cell[0], cell[1])) for cell in row] for row in matrix]
    col_widths = [max(len(row[i]) for row in dimensions) for i in range(len(dimensions[0]))]
    
    for row in dimensions:
        print("  ".join(cell.rjust(col_widths[i]) for i, cell in enumerate(row)))
    
    print("---------Coeficientes---------")

    coefficients = [[str(cell[2]) for cell in row] for row in matrix]
    col_widths = [max(len(row[i]) for row in coefficients) for i in range(len(coefficients[0]))]
    
    for row in coefficients:
        print("  ".join(cell.rjust(col_widths[i]) for i, cell in enumerate(row)))


#Función global

def build_irrep(partition, mode="YKR"):
     tableaux_dict = {}
     T = Snob2.StandardYoungTableaux(partition)
     for i in range(len(T)):
          tableau = T[i]
          key = search_young_tableau(tableau)
          if key and key not in tableaux_dict:
               tableaux_dict[key] = tableau
     ordered_2ndlev = sorted(tableaux_dict.keys(), key=lambda t: int(f"{t[0][0]}{t[1][0]}"))
     matrix = generate_representation_matrix(partition, ordered_2ndlev, mode)

     return matrix


# Método para convertir la matriz por bloques en una matriz de fracciones 


def decompress(matrix):
    # Inicializamos una lista vacía para los bloques
    blocks = []
    
    for row in matrix:
        block_row = []
        for block in row:
            # Obtenemos las dimensiones y el coeficiente
            block_rows, block_cols, coef = block

            # Convertimos el coeficiente a número flotante si es una fracción
            if isinstance(coef, Fraction):
                coef = coef.numerator / coef.denominator

            # Crear un bloque de ceros de tamaño block_rows x block_cols
            if block_rows == block_cols:
                # Si el bloque es cuadrado, inicializamos como la matriz identidad
                block_matrix = np.eye(block_rows) * coef
            else:
                # Si no es cuadrado, creamos una matriz de ceros
                block_matrix = np.zeros((block_rows, block_cols))

            # Agregar el bloque a la fila
            block_row.append(block_matrix)

        # Usamos hstack para juntar los bloques en la fila
        blocks.append(np.hstack(block_row))

    # Usamos vstack para juntar las filas completas
    result_matrix = np.vstack(blocks)
    
    return result_matrix

# lamb=Snob2.IntegerPartition([4,2,1,1])
# matrix = build_irrep(lamb, mode="YKR")
# pretty_print(matrix)

def direct_sum(matrices):
    n = len(matrices)
    
    total_rows = sum(len(mat) for mat in matrices)
    total_cols = sum(len(mat[0]) for mat in matrices)
    
    result = [[[0, 0, 0] for _ in range(total_cols)] for _ in range(total_rows)]
    
    row_start = 0
    col_start = 0
    
    for mat in matrices:
        mat_rows = len(mat)
        mat_cols = len(mat[0])
        
        for i in range(mat_rows):
            for j in range(mat_cols):
                result[row_start + i][col_start + j] = mat[i][j]
        
        row_start += mat_rows
        col_start += mat_cols
    
    for i in range(total_rows):
        for j in range(total_cols):
            if result[i][j] == [0, 0, 0]:
                result[i][j] = [result[i][i][0], result[j][j][0], 0]
    
    return result

def build_irrep_of_transposition(partition, t_n, mode="YKR"):
    
    # Caso base, la partición es de n elementos y queremos evaluar la representación en t_n
    
    n = partition.getn()
    if t_n == n:
        return build_irrep(partition, mode)
    else:
    
        # En otro caso, tenemos que expandir el árbol, es decir
        # Tenemos que calcular todas las transposiciones cuyos tableros de young salen de la actual
    
        partitions = branch(partition)
        submatrix_list = []
        for part in partitions:
            submatrix = build_irrep_of_transposition(part, t_n, mode)
            submatrix_list.append(submatrix)
        resp = direct_sum(submatrix_list)
        return resp


def build_irrep_of_permutation(partition, pi, mode="YKR"):
    transpositions = express_into_adyacent_transpositions(pi)
    currMatrix = np.array(decompress(build_irrep_of_transposition(partition, transpositions[0], mode)))
    for transposition in transpositions[1:]:
        currMatrix = currMatrix@np.array(decompress(build_irrep_of_transposition(partition, transposition, mode)).tolist())
    return currMatrix
