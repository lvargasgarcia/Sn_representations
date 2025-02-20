#include <gmp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>

// Función para inicializar una matriz de mpq_t
mpq_t **crear_matriz(int dim) {
    mpq_t **matriz = malloc(dim * sizeof(mpq_t *));
    for (int i = 0; i < dim; i++) {
        matriz[i] = malloc(dim * sizeof(mpq_t));
        for (int j = 0; j < dim; j++) {
            mpq_init(matriz[i][j]);  // Inicializa cada elemento
        }
    }
    return matriz;
}

mpq_t *crear_vector(int dim) {
    mpq_t *vector = malloc(dim * sizeof(mpq_t));
    for (int i = 0; i < dim; i++) {
        mpq_init(vector[i]);
    }
    return vector;
}

void copiar_matriz(mpq_t **A, mpq_t **B, int dim){
    for(int i = 0; i < dim; i++){
        for(int j = 0; j < dim; j++){
            mpq_set(B[i][j], A[i][j]);
        }
    }
}

// void multiplicar_matrices(mpq_t **A, mpq_t **B, mpq_t **C, int dim){
//     for(int i = 0; i < dim; i++){
//         for(int j = 0; j < dim; j++){
//             mpq_set_ui(C[i][j], 0, 1);
//             for(int k = 0; k < dim; k++){
//                 mpq_t temp;
//                 mpq_init(temp);
//                 mpq_mul(temp, A[i][k], B[k][j]);
//                 mpq_add(C[i][j], C[i][j], temp);
//                 //mpq_canonicalize(C[i][j]);
//                 mpq_clear(temp);
//             }
//         }
//     }
// }

void multiplicar_matrices(mpq_t **A, mpq_t **B, mpq_t **C, int dim){
    // Verificamos si la dimensión es mayor que 25 para paralelizar
    if (dim > 50) {
        #pragma omp parallel for collapse(2) num_threads(4)  // Paralelizamos los bucles i y j
        for (int i = 0; i < dim; i++) {
            for (int j = 0; j < dim; j++) {
                mpq_set_ui(C[i][j], 0, 1);
                for (int k = 0; k < dim; k++) {
                    mpq_t temp;
                    mpq_init(temp);
                    mpq_mul(temp, A[i][k], B[k][j]);
                    mpq_add(C[i][j], C[i][j], temp);
                    //mpq_canonicalize(C[i][j]);
                    mpq_clear(temp);
                }
            }
        }
    } else {
        // Si la dimensión es menor o igual a 25, hacemos la multiplicación de manera secuencial
        for (int i = 0; i < dim; i++) {
            for (int j = 0; j < dim; j++) {
                mpq_set_ui(C[i][j], 0, 1);
                for (int k = 0; k < dim; k++) {
                    mpq_t temp;
                    mpq_init(temp);
                    mpq_mul(temp, A[i][k], B[k][j]);
                    mpq_add(C[i][j], C[i][j], temp);
                    //mpq_canonicalize(C[i][j]);
                    mpq_clear(temp);
                }
            }
        }
    }
}


void multiplicar_matriz_escalar(mpq_t **A, mpq_t escalar, mpq_t **C, int dim){
    for(int i = 0; i < dim; i++){
        for(int j = 0; j < dim; j++){
            mpq_mul(C[i][j], A[i][j], escalar);
            //mpq_canonicalize(C[i][j]);
        }
    }
}

//Función para liberar una matriz de mpq_t
void liberar_matriz(mpq_t **matriz, int dim){
    for(int i = 0; i<dim; i++){
        for(int j = 0; j<dim; j++){
            mpq_clear(matriz[i][j]);
        }
        free(matriz[i]);
    }
    free(matriz);
}

void liberar_vector(mpq_t *vector, int dim){
    for(int i = 0; i<dim; i++){
        mpq_clear(vector[i]);
    }
    free(vector);
}

void imprimir_matriz(mpq_t **matriz, int dim){
    for(int i = 0; i<dim; i++){
        for(int j = 0; j<dim; j++){
            gmp_printf("%Qd ", matriz[i][j]);
        }
        printf("\n");
    }
}

void imprimir_matriz_double(mpq_t **matriz, int dim){
    for(int i = 0; i<dim; i++){
        for(int j = 0; j<dim; j++){
            printf("%lf ", mpq_get_d(matriz[i][j]));
        }
        printf("\n");
    }
}

void imprimir_vector(mpq_t *vector, int dim){
    for(int i = 0; i<dim; i++){
        gmp_printf("%Qd ", vector[i]);
    }
    printf("\n");
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

void williams_ft(int n, int dim, mpq_t **tau_matrix, mpq_t **qsigmatau_matrix, mpq_t **invSigma_matrix, mpq_t *f, 
    
    mpq_t **fourier_matrix, char *williams_sequence, int nfact) {
    
    mpq_t **aux_1 = crear_matriz(dim);
    mpq_t **aux_2 = crear_matriz(dim);

    mpq_t **p_matrix = crear_matriz(dim);

    copiar_matriz(qsigmatau_matrix, p_matrix, dim);

    multiplicar_matriz_escalar(p_matrix, f[0], fourier_matrix, dim);

    for(int k = 1; k < nfact; k++){
        // printf("%d\n", williams_sequence[k]);
        if(williams_sequence[k-1] == 't'){
            multiplicar_matrices(p_matrix, tau_matrix, aux_1, dim);
            copiar_matriz(aux_1, p_matrix, dim);
        }else{
            multiplicar_matrices(p_matrix, invSigma_matrix, aux_1, dim);
            copiar_matriz(aux_1, p_matrix, dim);
        }

        multiplicar_matriz_escalar(aux_1, f[k], aux_2, dim);

        //añadir aux_2 a fourier_matrix
        for(int i = 0; i < dim; i++){
            for(int j = 0; j < dim; j++){
                mpq_add(fourier_matrix[i][j], fourier_matrix[i][j], aux_2[i][j]);
                //mpq_canonicalize(fourier_matrix[i][j]);
            }
        }
    }

    liberar_matriz(aux_1, dim);
    liberar_matriz(aux_2, dim);
    liberar_matriz(p_matrix, dim);
}


// void traza(mpq_t **fourier_matrix, mpq_t **p_matrix, mpq_t *resp, int dim, mpq_t factor){

//     mpq_t traza;
//     mpq_init(traza);

//     mpq_set_si(traza, 0, 1);

//     for(int i = 0; i < dim; i++){
//         for(int j = 0; j < dim; j++){
//             mpq_t temp;
//             mpq_init(temp);
//             mpq_mul(temp, fourier_matrix[i][j], p_matrix[j][i]);
//             mpq_add(traza, traza, temp);
//             mpq_clear(temp);
//         }
//     }

//     mpq_mul(*resp, traza, factor);
//     mpq_clear(traza);
// }

mpq_t *williams_invft(int n, int dim, mpq_t **tau_matrix, mpq_t **qsigmatau_matrix, mpq_t **invSigma_matrix, 
    
    mpq_t **fourier_matrix, char *williams_sequence, int nfact) {
    
    mpq_t *invft_vector = crear_vector(nfact);
    mpq_t **aux_1 = crear_matriz(dim);
    mpq_t **aux_2 = crear_matriz(dim);

    //d_lambda/n!
    mpq_t factor;
    mpq_init(factor);
    mpq_set_si(factor, dim, nfact);
    mpq_canonicalize(factor);

    printf("Factor: ");
    gmp_printf("%Qd\n", factor);

    mpq_t **p_matrix = crear_matriz(dim);

    copiar_matriz(qsigmatau_matrix, p_matrix, dim);

    // // factor * traza(fourier_matrix * p_matrix)
    mpq_t traza;
    mpq_init(traza);
    mpq_set_si(traza, 0, 1);

    for(int i = 0; i < dim; i++){
        for(int j = 0; j < dim; j++){
            mpq_t temp;
            mpq_init(temp);
            mpq_mul(temp, fourier_matrix[i][j], p_matrix[j][i]);
            mpq_add(traza, traza, temp);
            mpq_clear(temp);
        }
    }

    mpq_mul(invft_vector[0], traza, factor);
    mpq_clear(traza);

    for(int k = 1; k < nfact; k++){
        // printf("%d\n", williams_sequence[k]);
        if(williams_sequence[k-1] == 't'){
            multiplicar_matrices(p_matrix, tau_matrix, aux_1, dim);
            copiar_matriz(aux_1, p_matrix, dim);
        }else{
            multiplicar_matrices(p_matrix, invSigma_matrix, aux_1, dim);
            copiar_matriz(aux_1, p_matrix, dim);
        }

        //factor * traza(fourier_matrix * p_matrix)
        mpq_t tr;
        mpq_init(tr);
        mpq_set_si(tr, 0, 1);

        for(int i = 0; i < dim; i++){
            for(int j = 0; j < dim; j++){
                mpq_t tempr;
                mpq_init(tempr);
                mpq_mul(tempr, fourier_matrix[i][j], p_matrix[j][i]);
                mpq_add(tr, tr, tempr);
                mpq_clear(tempr);
            }
        }

        mpq_mul(invft_vector[k], tr, factor);
        mpq_clear(tr);

    }

    mpq_clear(factor);
    liberar_matriz(aux_1, dim);
    liberar_matriz(aux_2, dim);
    liberar_matriz(p_matrix, dim);

    return invft_vector;

}

void inicializar_matriz(mpq_t **matriz, int dim, long **valores) {
    for (int i = 0; i < dim; i++) {
        for (int j = 0; j < dim; j++) {
            mpq_set_si(matriz[i][j], valores[i * dim + j][0], valores[i * dim + j][1]);
            mpq_canonicalize(matriz[i][j]);
        }
    }
}

void inicializar_vector(mpq_t *vector, int dim, long **valores) {
    for (int i = 0; i < dim; i++) {
        mpq_set_si(vector[i], valores[i][0], valores[i][1]);
        mpq_canonicalize(vector[i]);
    }
}

mpq_t *zip_vectores(long *num, long *den, int length){
    long **vectores = malloc(length * sizeof(long *));
    
    for(int i = 0; i < length; i++){
        vectores[i] = malloc(2 * sizeof(long));
        vectores[i][0] = num[i];
        vectores[i][1] = den[i];
    }
    
    mpq_t *resp = crear_vector(length);
    inicializar_vector(resp, length, vectores);
    
    for(int i = 0; i < length; i++){
        free(vectores[i]);
    }


    free(vectores);
    return resp;
}


mpq_t **zip_matrices(long *num, long *den, int dim){
    long **matrices = malloc(dim * dim * sizeof(long *));
    
    for(int i = 0; i < dim * dim; i++){
        matrices[i] = malloc(2 * sizeof(long));
        matrices[i][0] = num[i];
        matrices[i][1] = den[i];
    }
    
    mpq_t **resp = crear_matriz(dim);
    inicializar_matriz(resp, dim, matrices);
    
    for(int i = 0; i < dim * dim; i++){
        free(matrices[i]);
    }
    
    free(matrices);
    return resp;
}

mpq_t *matriz_a_vector(mpq_t **A, mpq_t *vector, int dim){
    
    for(int i = 0; i < dim; i++){
        for(int j = 0; j < dim; j++){
            mpq_set(vector[i * dim + j], A[i][j]);
        }
    }

    return vector;
}

char *matriz_a_vector_string(mpq_t **A, int dim) {
    char *result = NULL;  // Para almacenar la cadena final
    size_t offset = 0;    // Longitud actual de la cadena (sin incluir el '\0')

    // Recorremos la matriz
    for (int i = 0; i < dim; i++) {
        for (int j = 0; j < dim; j++) {
            char *number;
            // Convertimos el número a cadena
            gmp_asprintf(&number, "%Qd", A[i][j]);
            // printf("-------------\n");
            // printf("%s\n", number);
            // printf("-------------\n");

            size_t num_length = strlen(number);
            size_t extra = (offset == 0) ? num_length : (1 + num_length); // 1 para el espacio en iteraciones posteriores
            size_t new_length = offset + extra + 1; // +1 para el carácter nulo

            // Realocamos la memoria
            char *temp = realloc(result, new_length);
            if (!temp) {
                perror("Error reallocating memory");
                free(number);
                free(result);  // Liberamos lo que se haya asignado previamente
                return NULL;
            }
            result = temp;

            // Concatenamos el número (precedido de un espacio si no es el primero)
            if (offset == 0) {
                sprintf(result, "%s", number);
                offset = num_length;  // Actualizamos offset sin contar el '\0'
            } else {
                sprintf(result + offset, " %s", number);
                offset += 1 + num_length;
            }

            // printf("result: %s\n", result);

            free(number);
        }
    }

    return result;
}

char *vector_a_vector_string(mpq_t *vector, int length){
    
    char *result = NULL;  // Para almacenar la cadena final
    size_t offset = 0;    // Longitud actual de la cadena (sin incluir el '\0')
    
    for (int j = 0; j < length; j++) {
        char *number;
        // Convertimos el número a cadena
        gmp_asprintf(&number, "%Qd", vector[j]);
        // printf("-------------\n");
        // printf("%s\n", number);
        // printf("-------------\n");

        size_t num_length = strlen(number);
        size_t extra = (offset == 0) ? num_length : (1 + num_length); // 1 para el espacio en iteraciones posteriores
        size_t new_length = offset + extra + 1; // +1 para el carácter nulo

        // Realocamos la memoria
        char *temp = realloc(result, new_length);
        if (!temp) {
            perror("Error reallocating memory");
            free(number);
            free(result);  // Liberamos lo que se haya asignado previamente
            return NULL;
        }
        result = temp;

        // Concatenamos el número (precedido de un espacio si no es el primero)
        if (offset == 0) {
            sprintf(result, "%s", number);
            offset = num_length;  // Actualizamos offset sin contar el '\0'
        } else {
            sprintf(result + offset, " %s", number);
            offset += 1 + num_length;
        }

        // printf("result: %s\n", result);

        free(number);
    }

    return result;

}


char *williams_wrapper(int n, int nfact, int dim, long *tau_nums, long *tau_dens, long *qsigmatau_nums, long *qsigmatau_dens, long *invSigma_nums, long *invSigma_dens,
    
    long *f_nums, long *f_dens, char *williams_sequence, int doInv){

    mpq_t **tau_matrix = zip_matrices(tau_nums, tau_dens, dim);
    mpq_t **qsigmatau_matrix = zip_matrices(qsigmatau_nums, qsigmatau_dens, dim);
    mpq_t **invSigma_matrix = zip_matrices(invSigma_nums, invSigma_dens, dim);
    mpq_t *f = zip_vectores(f_nums, f_dens, nfact);

    mpq_t **fourier_matrix = crear_matriz(dim);

    williams_ft(n, dim, tau_matrix, qsigmatau_matrix, invSigma_matrix, f, fourier_matrix, williams_sequence, nfact);
    
    char *fourier_vector = matriz_a_vector_string(fourier_matrix, dim);

    char *invft_vector;
    char *result;
    
    if(doInv == 1){
        
        mpq_t *q_invft_vector = williams_invft(n, dim, tau_matrix, qsigmatau_matrix, invSigma_matrix, fourier_matrix, williams_sequence, nfact);
        invft_vector = vector_a_vector_string(q_invft_vector, nfact);
        liberar_vector(q_invft_vector, nfact);

        //hacer la concatenación fourier_vector + $ + invft_vector
        result = malloc(strlen(fourier_vector) + strlen(invft_vector) + 2);
        strcpy(result, fourier_vector);
        strcat(result, "$");
        strcat(result, invft_vector);
        free(invft_vector);
        free(fourier_vector);
    }

    // for(int i = 0; i < dim * dim; i++){
    //     mpz_t num;
    //     mpz_t den;
    //     mpz_init(num);
    //     mpz_init(den);
    //     mpq_get_num(num, fourier_vector[i]);
    //     mpq_get_den(den, fourier_vector[i]);
    //     fourier_vector_nums[i] = mpz_get_si(num);
    //     fourier_vector_dens[i] = mpz_get_si(den);
    //     mpz_clear(num);
    //     mpz_clear(den);
    //     mpq_clear(fourier_vector[i]);
    // }

    // free(fourier_vector);
    liberar_matriz(tau_matrix, dim);
    liberar_matriz(qsigmatau_matrix, dim);
    liberar_matriz(invSigma_matrix, dim);
    liberar_vector(f, nfact);
    liberar_matriz(fourier_matrix, dim);

    if(doInv == 1){
        return result;
        //return fourier_vector;
    }else{
        return fourier_vector;
    }

}

// int main(){
    
//     int dim = 2;
//     int n = 4;
//     int nfact = 24;
    
//     // Inicializar tau_matrix
//     long tau_valores[4][2] = {{-1, 1}, {0, 1}, {0, 1}, {1, 1}};
//     long *tau_valores_num = malloc(dim * dim * sizeof(long));
//     long *tau_valores_den = malloc(dim * dim * sizeof(long));
//     for(int i = 0; i < dim * dim; i++){
//         tau_valores_num[i] = tau_valores[i][0];
//         tau_valores_den[i] = tau_valores[i][1];
//     }
//     // mpq_t **tau_matrix = zip_matrices(tau_valores_num, tau_valores_den, dim);
//     // // mpq_t **tau_matrix = crear_matriz(dim);
//     // // inicializar_matriz(tau_matrix, dim, tau_valores);

//     // Inicializar qsigmatau_matrix
//     long qsigmatau_valores[4][2] = {{-1, 2}, {-3, 4}, {1, 1}, {-1, 2}};
//     long *qsigmatau_valores_num = malloc(dim * dim * sizeof(long));
//     long *qsigmatau_valores_den = malloc(dim * dim * sizeof(long));
//     for(int i = 0; i < dim * dim; i++){
//         qsigmatau_valores_num[i] = qsigmatau_valores[i][0];
//         qsigmatau_valores_den[i] = qsigmatau_valores[i][1];
//     }
//     // mpq_t **qsigmatau_matrix = zip_matrices(qsigmatau_valores_num, qsigmatau_valores_den, dim);
//     // mpq_t **qsigmatau_matrix = crear_matriz(dim);

//     // inicializar_matriz(qsigmatau_matrix, dim, qsigmatau_valores);

//     // Inicializar invSigma_matrix
//     long invSigma_valores[4][2] = {{1, 2}, {-3, 4}, {-1, 1}, {-1, 2}};
//     long *invSigma_valores_num = malloc(dim * dim * sizeof(long));
//     long *invSigma_valores_den = malloc(dim * dim * sizeof(long));
//     for(int i = 0; i < dim * dim; i++){
//         invSigma_valores_num[i] = invSigma_valores[i][0];
//         invSigma_valores_den[i] = invSigma_valores[i][1];
//     }
//     // mpq_t **invSigma_matrix = zip_matrices(invSigma_valores_num, invSigma_valores_den, dim);
//     // mpq_t **invSigma_matrix = crear_matriz(dim);

//     // inicializar_matriz(invSigma_matrix, dim, invSigma_valores);

//     long f_list[24][2] = {{1754535031782657,100000000000000},{5498636373984343,5000000000000},{613919104734999,125000000000000},{8617140607505623,100000000000000},{4465508510838359,1250000000000000},{1160356936430663,125000000000000},{3437005568598397,2500000000000000},{4232181797708361,10000000000000},{2798941837152681,12500000000000},{482427733963037,20000000000000},{8143547577358403,50000000000000},{6496219487737669,2500000000000000},{6267912893591607,100000000000000},{472471095927047,312500000000},{5698922332692467,125000000000000},{945040582286111,500000000000000},{3190522601659171,250000000000000},{5923432050678439,50000000000000},{2078573808259901,1000000000000},{8440173526572533,1250000000000000},{581841296241777,1000000000000},{3316213616154527,100000000000000},{1599833420178797,2000000000000},{7695982944381431,25000000000000}};
//     long *f_nums = malloc(nfact * sizeof(long));
//     long *f_dens = malloc(nfact * sizeof(long));
//     for(int i = 0; i < nfact; i++){
//         f_nums[i] = f_list[i][0];
//         f_dens[i] = f_list[i][1];
//     }
    
//     // mpq_t *f = zip_vectores(f_nums, f_dens, nfact);

//     // // Imprimir las matrices para verificar
//     // printf("Tau Matrix:\n");
//     // imprimir_matriz(tau_matrix, dim);
//     // printf("QsigmaTau Matrix:\n");
//     // imprimir_matriz(qsigmatau_matrix, dim);
//     // printf("InvSigma Matrix:\n");
//     // imprimir_matriz(invSigma_matrix, dim);
//     // imprimir_vector(f, nfact);





//     // mpq_t **fourier_matrix = crear_matriz(dim);
//     // williams(n, dim, tau_matrix, qsigmatau_matrix, invSigma_matrix, f, fourier_matrix);

//     // // multiplicar_matrices(tau_matrix, qsigmatau_matrix, fourier_matrix, dim);
//     // imprimir_matriz(fourier_matrix, dim);

//     // liberar_matriz(fourier_matrix, dim);

//     // liberar_matriz(tau_matrix, dim);
//     // liberar_matriz(qsigmatau_matrix, dim);
//     // liberar_matriz(invSigma_matrix, dim);
//     // liberar_vector(f, nfact);

//     char *seq = "iiititiiititiiitiiititi";
    
//     char *fourier_vector = williams_wrapper(n, nfact, dim, tau_valores_num, tau_valores_den, qsigmatau_valores_num, qsigmatau_valores_den, 
//         invSigma_valores_num, invSigma_valores_den, f_nums, f_dens, seq, 1);

//     printf("%s\n", fourier_vector);

//     free(fourier_vector);
//     free(tau_valores_num);
//     free(tau_valores_den);
//     free(qsigmatau_valores_num);
//     free(qsigmatau_valores_den);
//     free(invSigma_valores_num);
//     free(invSigma_valores_den);
//     free(f_nums);
//     free(f_dens);
    
//     // for(int i = 0; i < dim * dim; i++){
//     //     printf("%lf, ", (float)(fourier_vector_nums[i]/fourier_vector_dens[i]));
//     // }


    
//     return 0;
// }