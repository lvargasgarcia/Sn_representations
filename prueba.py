import torch
import snob as Snob2
import sys
import numpy as np
import math
from fractions import Fraction

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

lamb=Snob2.IntegerPartition([4, 2, 1, 1])
T=Snob2.StandardYoungTableaux(lamb)

tableaux_dict = {}

for i in range(len(T)):
    tableau = T[i]
    key = search_young_tableau(tableau)
    
    if key and key not in tableaux_dict:  # Solo almacenar el primero encontrado
        tableaux_dict[key] = tableau


ordered_2ndlev = sorted(tableaux_dict.keys(), key=lambda t: int(f"{t[0][0]}{t[1][0]}"))
print(ordered_2ndlev)

# Aplicamos el teorema de Young del artículo de Clausen para construir la matriz de la representación 

# Antes, necesitamos una funcion que nos diga la particion asociada a una cadena gamma (usamos el símbolo de Yamanouchi)

def get_partition_from_gamma(y_symbol, partition):
    n = partition.getn()
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

# Imprimir el diccionario resultante
for key in ordered_2ndlev:
    ((a,_),(b,_)) = key
    print(f"Símbolo de Yamanouchi ({a},{b}):")
    print("Tablero:")
    tableau = tableaux_dict[key]
    print(tableau)
    print("Gamma asociada al tablero:")
    gamma = get_partition_from_gamma((a,b), lamb)
    print(gamma)
    print("d_gamma")
    print(frame_robinson_thrall(gamma))
    print("--------------------")


# Ahora vamos a crear una estructura de datos para cada elemento en la matriz de representación, ya que esta es por bloques
# tiene nxn bloques y su dimension real es d_alpha, donde alpha es la partición dada 

# La matriz sera un np.array de dimension nxn y cada elemento será una dupla ((a,b) --> dimension del bloque, coef --> coeficiente de la identidad)
    
def generate_representation_matrix(partition, mode="YKR"):
    
    n = len(ordered_2ndlev)
    matrix = [[[-1,-1,None] for _ in range(n)] for _ in range(n)]
    # Primero hacemos la diagonal
    
    for i in range(n):
          print(i)
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
               ((a2,b2),(c2,d2)) = ordered_2ndlev[j]
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
                   q_sqr_chi = math.sqrt(1 - (1/chi)**2)
                   print((i,j))
                   
                   if mode == "YKR":
                       if matrix[i][i][2] == None:
                           matrix[i][i][2] = Fraction(1/chi)
                       if matrix[j][j][2] == None:
                           matrix[j][j][2] = Fraction(-(1/chi))
                       if matrix[i][j][2] == None:
                           matrix[i][j][2] = Fraction(1)
                       if matrix[j][i][2] == None:
                           matrix[j][i][2] = Fraction(q_sqr_chi)
                    
                   if mode == "YSR":
                         if matrix[i][i][2] == None:
                              matrix[i][i][2] = Fraction(1/chi)
                         if matrix[j][j][2] == None:
                              matrix[j][j][2] = Fraction(-(1/chi))
                         if matrix[i][j][2] == None:
                              matrix[i][j][2] = Fraction(q_sqr_chi)
                         if matrix[j][i][2] == None:
                              matrix[j][i][2] = Fraction(1)
               else:
                   matrix[i][j][2] = Fraction(0)
                   matrix[j][i][2] = Fraction(0)
                                         

    return matrix

def pretty_print(matrix):
    print("Matriz de representación:")
    print("---------Dimensiones---------")
    for row in matrix:
        for elem in row:
            print(" ", (elem[0], elem[1]), " ", end="")
        print()
    print("---------Coeficientes---------")
    for row in matrix:
        for elem in row:
            print(" ", str(elem[2]), " ", end="")
        print()

matrix = generate_representation_matrix(lamb,mode="YSR")
print(pretty_print(matrix))
