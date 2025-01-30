import torch
import snob as Snob2


# En este fichero se implementa la descomposición de una permutación de S_n en producto de transposiciones adyacentes
# La idea es tener una funcion con esta signatura express_into_adyacent_transpositions :: SnElement -> [int], donde int
# Es la lista ordenada de enteros que representan los indices de las transposiciones
#                      (1 2 3 4 5)
# Por ejemplo, si pi = ----------- entonces express_into_adyacent_transpositions(pi) = [3,5] porque pi = (2 3) (4 5)
#                      (1 3 2 5 4)
