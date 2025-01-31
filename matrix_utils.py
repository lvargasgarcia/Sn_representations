import torch
import snob as Snob2
import sys
import numpy as np
from fractions import Fraction


def decompress(matrix):
    """
    Convierte una matriz de fracciones por bloques en una matriz de números reales.

    Args:
        matrix (list): La matriz de fracciones por bloques a convertir.

    Returns:
        np.ndarray: La matriz convertida a números reales.
    """
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

def pretty_print(matrix):
    """
    Imprime una matriz de representación por bloques de forma legible.

    Args:
        matrix (list): La matriz de representación por bloques a imprimir.
    """
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

def print_matrix(matrix):
    """
    Imprime una matriz de números reales de forma legible.

    Args:
        matrix (list): La matriz de números reales a imprimir.
    """
    for sublista in matrix:
        print(f"[{' '.join(map(str, sublista))}]")