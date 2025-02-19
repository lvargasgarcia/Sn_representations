#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>

// Función para inicializar una matriz de mpq_t
mpq_t **crear_matriz(int filas, int columnas) {
    mpq_t **matriz = malloc(filas * sizeof(mpq_t *));
    for (int i = 0; i < filas; i++) {
        matriz[i] = malloc(columnas * sizeof(mpq_t));
        for (int j = 0; j < columnas; j++) {
            mpq_init(matriz[i][j]);  // Inicializa cada elemento
        }
    }
    return matriz;
}

// Función para liberar la memoria de la matriz
void liberar_matriz(mpq_t **matriz, int filas, int columnas) {
    for (int i = 0; i < filas; i++) {
        for (int j = 0; j < columnas; j++) {
            mpq_clear(matriz[i][j]); // Limpia cada número racional
        }
        free(matriz[i]); // Libera la fila
    }
    free(matriz); // Libera la matriz
}

// Función para imprimir la matriz
void imprimir_matriz(mpq_t **matriz, int filas, int columnas) {
    printf("------------------------------------");
    printf("\n");
    for (int i = 0; i < filas; i++) {
        for (int j = 0; j < columnas; j++) {
            gmp_printf("%Qd ", matriz[i][j]); // Imprime como fracción
        }
        printf("\n");
    }
    printf("------------------------------------");
    printf("\n");
}


void sumar_matrices(mpq_t **A, mpq_t **B, mpq_t **C, int filas, int columnas) {
    for (int i = 0; i < filas; i++) {
        for (int j = 0; j < columnas; j++) {
            mpq_add(C[i][j], A[i][j], B[i][j]); 
        }
    }
}

void multiplicar_matrices(mpq_t **A, mpq_t **B, mpq_t **C, int filasA, int columnasA, int columnasB) {
    for (int i = 0; i < filasA; i++) {
        for (int j = 0; j < columnasB; j++) {
            mpq_set_ui(C[i][j], 0, 1); 
            for (int k = 0; k < columnasA; k++) {
                mpq_t temp;
                mpq_init(temp);
                mpq_mul(temp, A[i][k], B[k][j]); 
                mpq_add(C[i][j], C[i][j], temp); 
                mpq_clear(temp);
            }
        }
    }
}

void multiplicar_escalar_matriz(mpq_t escalar, mpq_t **A, mpq_t **C, int filas, int columnas) {
    for (int i = 0; i < filas; i++) {
        for (int j = 0; j < columnas; j++) {
            mpq_mul(C[i][j], escalar, A[i][j]); 
        }
    }
}

void compose(int *result, int *p1, int *p2, int n) {
    for (int i = 0; i < n; i++) {
        result[i] = p1[p2[i] - 1];
    }
}

void inverse(int *result, int *p, int n) {
    for (int i = 0; i < n; i++) {
        result[p[i] - 1] = i + 1;
    }
}

int williamsCondition(int *p, int n) {
    int i = 0;
    while (p[i] != n) i++;
    int r = p[(i % (n - 1)) + 1];
    return (r % (n - 1)) + 1 == p[0];
}

void print_permutation(int *p, int n) {
    for (int i = 0; i < n; i++) {
        printf("%d ", p[i]);
    }
    printf("\n");
}

int distintos(int *p1, int *p2, int n) {
    for (int i = 0; i < n; i++) {
        if (p1[i] != p2[i]) {
            return 1;
        }
    }
    return 0;
}

void williams_path(int n, mpq_t **tau_matrix, mpq_t **qsigmatau_matrix, mpq_t **invSigma_matrix, mpq_t *valores_f, mpq_t **fourier_matrix, int dim) {

    int *tau = malloc(n * sizeof(int));
    int *sigma = malloc(n * sizeof(int));
    int *q = malloc(n * sizeof(int));
    int *qtau = malloc(n * sizeof(int));
    int *qsigma = malloc(n * sizeof(int));
    int *qsigmatau = malloc(n * sizeof(int));
    int *invSigma = malloc(n * sizeof(int));
    int *p = malloc(n * sizeof(int));
    int *temp = malloc(n * sizeof(int));

    tau[0] = 2;
    tau[1] = 1;
    for (int i = 2; i < n; i++) {
        tau[i] = i + 1;
    }

    for (int i = 0; i < n - 1; i++) {
        sigma[i] = i + 2;
    }
    sigma[n - 1] = 1;

    for (int i = 0; i < n; i++) {
        q[i] = n - i;
    }

    compose(qtau, q, tau, n);
    compose(qsigma, q, sigma, n);
    compose(qsigmatau, qsigma, tau, n);
    inverse(invSigma, sigma, n);
    compose(p, qsigma, tau, n);

    print_permutation(p, n);

    mpq_t **aux_1 = crear_matriz(dim, dim);
    mpq_t **aux_2 = crear_matriz(dim, dim);
    mpq_t **p_matrix = crear_matriz(dim, dim);

    for(int i = 0; i < dim; i++) {
        for(int j = 0; j < dim; j++) {
            mpq_set(p_matrix[i][j], qsigmatau_matrix[i][j]);
        }
    }

    multiplicar_escalar_matriz(valores_f[0], p_matrix, fourier_matrix, dim, dim);

    int i = 1;

    while (1) {
        
        if (distintos(p, qtau, n) == 0) {
            break;
        }
        
        if (distintos(p, qsigmatau, n) == 1) {
            if (williamsCondition(p, n) && distintos(p, qsigma, n)) {
                
                compose(temp, p, tau, n);
                
                // multiplicar_matrices(p_matrix, tau_matrix, aux_1, dim, dim, dim);
                // multiplicar_escalar_matriz(valores_f[i], aux_1, aux_2, dim, dim);
                // sumar_matrices(fourier_matrix, aux_2, fourier_matrix, dim, dim);

                // liberar_matriz(p_matrix, dim, dim);
                p_matrix = aux_1;

                for (int i = 0; i < n; i++) {
                    p[i] = temp[i];
                }


            } else {
                
                compose(temp, p, invSigma, n);

                // multiplicar_matrices(p_matrix, invSigma_matrix, aux_1, dim, dim, dim);
                // multiplicar_escalar_matriz(valores_f[i], aux_1, aux_2, dim, dim);
                // sumar_matrices(fourier_matrix, aux_2, fourier_matrix, dim, dim);

                // liberar_matriz(p_matrix, dim, dim);
                p_matrix = aux_1;
                
                for (int i = 0; i < n; i++) {
                    p[i] = temp[i];
                }
            }
        } else {
            
            compose(temp, p, invSigma, n);
            
            for (int i = 0; i < n; i++) {
                p[i] = temp[i];
            }

            //multiplicar_matrices(p_matrix, invSigma_matrix, aux_1, dim, dim, dim);
            //multiplicar_escalar_matriz(valores_f[i], aux_1, aux_2, dim, dim);
            //sumar_matrices(fourier_matrix, aux_2, fourier_matrix, dim, dim);

            //liberar_matriz(p_matrix, dim, dim);
            //p_matrix = aux_1;
        }

        print_permutation(p, n);

        i++;

    }

    free(tau);
    free(sigma);
    free(q);
    free(qtau);
    free(qsigma);
    free(qsigmatau);
    free(invSigma);
    free(p);
    free(temp);
}

int leer_matriz(long** matrix, int dim){
    // Inicializar la matriz (dim x dim) y leer los valores desde stdin
    for (int i = 0; i < dim; i++) {
        matrix[i] = (long*)malloc(dim * sizeof(long));
        if (matrix[i] == NULL) {
            fprintf(stderr, "Error al reservar memoria para la fila %d.\n", i);
            return 1;
        }

        for (int j = 0; j < dim; j++) {
            if (fread(&matrix[i][j], sizeof(long), 1, stdin) != 1) {
                fprintf(stderr, "Error al leer los datos de la matriz.\n");
                return 1;
            }
        }
    }

    return 0;
}

int leer_f(mpq_t *valores_f, int length){
    
    long *f_nums = (long*)malloc(length * sizeof(long));
    if(f_nums == NULL){
        fprintf(stderr, "Error al reservar memoria para los numeradores de f");
        return 1;
    }

    for(int i = 0; i < length; i++){
        if (fread(&f_nums[i], sizeof(long), 1, stdin) != 1){
            fprintf(stderr, "Error al leer los numeradores de f");
            return 1;
        }
    }

    long *f_dens = (long*)malloc(length * sizeof(long));
    if(f_dens == NULL){
        fprintf(stderr, "Error al reservar memoria para los denominadores de f");
        return 1;
    }

    for(int i = 0; i < length; i++){
        if (fread(&f_dens[i], sizeof(long), 1, stdin) != 1){
            fprintf(stderr, "Error al leer los denominadores de f");
            return 1;
        }
    }


    for(int i = 0; i < length; i++){
        mpq_init(valores_f[i]);
        mpq_set_si(valores_f[i], f_nums[i], f_dens[i]);
        printf("f[%d] = ", i);
        printf("%ld/%ld\n", f_nums[i], f_dens[i]);
    }

}

void imprimir_matriz_entera(long** matrix, int dim){
    
    printf("Matriz (dim = %d):\n", dim);
    for (int i = 0; i < dim; i++) {
        for (int j = 0; j < dim; j++) {
            printf("%ld ", matrix[i][j]);
        }
        printf("\n");
    }

}

void limpiar_matriz_entera(long** matrix, int dim){
        
    for (int i = 0; i < dim; i++) {
            free(matrix[i]);
        }
    free(matrix);
}

mpq_t** cargar_matriz(int dim){
    
    mpq_t **matrix = crear_matriz(dim, dim);

    long **numeradores = (long**)malloc(dim * sizeof(long*));
    if (numeradores == NULL) {
        fprintf(stderr, "Error al reservar memoria para la matriz.\n");
        return NULL;
    }

    leer_matriz(numeradores, dim);
    
    long **denominadores = (long**)malloc(dim * sizeof(long*));
    if (denominadores == NULL) {
        fprintf(stderr, "Error al reservar memoria para la matriz.\n");
        return NULL;
    }

    leer_matriz(denominadores, dim);

    for (int i = 0; i < dim; i++) {
        for (int j = 0; j < dim; j++) {
            mpq_set_si(matrix[i][j], numeradores[i][j], denominadores[i][j]);
        }
    }

    limpiar_matriz_entera(numeradores, dim);
    limpiar_matriz_entera(denominadores, dim);
    
    return matrix;

}

void imprimir_vector(mpq_t *vector, int dim){
    printf("Vector (dim = %d):\n", dim);
    for (int i = 0; i < dim; i++) {
        gmp_printf("%Qd ", vector[i]);
    }
    printf("\n");
}


int main() {
    
    int dim;

    if (fread(&dim, sizeof(int), 1, stdin) != 1) {
        fprintf(stderr, "Error al leer la dimensión.\n");
        return 1;
    }

    mpq_t **tau_matrix = cargar_matriz(dim);

    mpq_t **qsigmatau_matrix = cargar_matriz(dim);

    mpq_t **invSigma_matrix = cargar_matriz(dim);

    int length;

    if (fread(&length, sizeof(int), 1, stdin) != 1) {
        fprintf(stderr, "Error al leer la long de f");
    }

    mpq_t *f = (mpq_t *)malloc(sizeof(long) * length);
    leer_f(f, length);

    int n;
    
    if (fread(&n, sizeof(int), 1, stdin) != 1) {
        fprintf(stderr, "Error al leer la long de f");
    }

    mpq_t **fourier_matrix = crear_matriz(dim, dim);

    williams_path(n, tau_matrix, qsigmatau_matrix, invSigma_matrix, f, fourier_matrix, dim);

    imprimir_matriz(fourier_matrix, dim, dim);

    liberar_matriz(fourier_matrix, dim, dim);
    liberar_matriz(tau_matrix, dim, dim);
    liberar_matriz(qsigmatau_matrix, dim, dim);
    liberar_matriz(invSigma_matrix, dim, dim);
    free(f);

    return 0;
}




