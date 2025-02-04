import torch
import snob as Snob2
import cnine
from irrep import Irrep
from permutaciones import *
from matrix_utils import *
import unittest
import time
import datetime

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

    return maxerror <= 1e-6

def log_test_result(test_name, duration):
    with open("resultados.txt", "a") as file:
        file.write(f"{test_name}:{datetime.datetime.now()} passed in {duration:.2f} mins\n")

class TestIrrepRepresentation(unittest.TestCase):
    
    def _test_irrep_representation(self, partition):
        G = Snob2.Sn(partition.getn())
        rho = Snob2.SnIrrep(partition)
        mi_rho = Irrep(partition, mode="YOR")
        
        print("La partición es", partition)
        
        for i in range(len(G)):
            if i != 0:
                pi = G[i]
                # print("La permutación es", pi)
                # print("La descomposición de la permutación en transposiciones adyacentes es: ", express_into_adyacent_transpositions(pi))
                
                mi_matriz = mi_rho.evaluate(pi)
                matriz_snob = rho[pi].torch().tolist()
                
                self.assertTrue(np.allclose(mi_matriz, matriz_snob, atol=1e-6), "Las matrices no coinciden")

    
    def test_partitions_representations(self):
        
        t_0 = time.time()
        
        partitions = [[2], [1,1], [3], [2,1], [1,1,1], [4], [3,1], [2,2], [2,1,1], [5], [4,1], [3,2], [3,1,1], [6]
                      , [5,1], [4,2], [4,1,1], [7], [6,1], [5,2], [5,1,1], [8], [7,1], [6,2], [6,1,1]]
        partitions = [Snob2.IntegerPartition(partition) for partition in partitions]
        for partition in partitions:
            print("n =", partition.getn())
            print("Probando la partición", partition)
            self._test_irrep_representation(partition)
        
        duration = (time.time() - t_0) / 60
        print("Tiempo de ejecución:", duration, "minutos")
        log_test_result("test_partitions_representations", duration)
    
    def test_partitions_representations_fast(self):
        
        t_0 = time.time()

        partitions = [[2], [1,1], [3], [2,1], [1,1,1], [4], [3,1], [2,2], [2,1,1], [5], [4,1], [3,2], [3,1,1], [6]
                      , [5,1], [4,2], [4,1,1], [7], [6,1], [5,2], [5,1,1]]
        partitions = [Snob2.IntegerPartition(partition) for partition in partitions]
        for partition in partitions:
            print("n =", partition.getn())
            print("Probando la partición", partition)
            self._test_irrep_representation(partition)
        
        duration = (time.time() - t_0) / 60
        print("Tiempo de ejecución:", duration, "minutos")
        log_test_result("test_partitions_representations_fast", duration)

        
                
if __name__ == "__main__":
    unittest.main()

