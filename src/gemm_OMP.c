#include "mnblas.h"
#include "complexe2.h"
#include <omp.h>


void transpose_f_omp(float *matrice, int lignes, int colonnes, float *res) {

    float temp[lignes * colonnes];

    #pragma omp parallel for
    for (int i = 0; i < lignes; i++) {
        for (int j = 0; j < colonnes; j++) {
            temp[j * lignes + i] = matrice[i * colonnes + j];
        }
    }

    #pragma omp parallel for
    for (int i = 0; i < colonnes; i++) {
        for (int j = 0; j < lignes; j++) {
            res[i * lignes + j] = temp[i * lignes + j];
        }
    }
}

void transpose_f2_omp(complexe_float_t *matrice, int lignes, int colonnes, complexe_float_t *res) {

    complexe_float_t temp[lignes * colonnes];

    #pragma omp parallel for
    for (int i = 0; i < lignes; i++) {
        for (int j = 0; j < colonnes; j++) {
            temp[j * lignes + i] = matrice[i * colonnes + j];
        }
    }

    #pragma omp parallel for
    for (int i = 0; i < colonnes; i++) {
        for (int j = 0; j < lignes; j++) {
            res[i * lignes + j] = temp[i * lignes + j];
        }
    }
}

void transpose_d_omp(double *matrice, int lignes, int colonnes, double *res) {

    double temp[lignes * colonnes];

    #pragma omp parallel for
    for (int i = 0; i < lignes; i++) {
        for (int j = 0; j < colonnes; j++) {
            temp[j * lignes + i] = matrice[i * colonnes + j];
        }
    }

    #pragma omp parallel for
    for (int i = 0; i < colonnes; i++) {
        for (int j = 0; j < lignes; j++) {
            res[i * lignes + j] = temp[i * lignes + j];
        }
    }
}

void transpose_d2_omp(complexe_double_t *matrice, int lignes, int colonnes, complexe_double_t *res) {

    complexe_double_t temp[lignes * colonnes];

    #pragma omp parallel for
    for (int i = 0; i < lignes; i++) {
        for (int j = 0; j < colonnes; j++) {
            temp[j * lignes + i] = matrice[i * colonnes + j];
        }
    }

    #pragma omp parallel for
    for (int i = 0; i < colonnes; i++) {
        for (int j = 0; j < lignes; j++) {
            res[i * lignes + j] = temp[i * lignes + j];
        }
    }
}

void transpose_comp_f_omp(complexe_float_t *matrice, int lignes, int colonnes, complexe_float_t *res){
    complexe_float_t temp[lignes * colonnes];

    #pragma omp parallel for
    for (int i = 0; i < lignes; i++) {
        for (int j = 0; j < colonnes; j++) {
            temp[j * lignes + i] = matrice[i * colonnes + j];
        }
    }

    #pragma omp parallel for
    for (int i = 0; i < colonnes; i++) {
        for (int j = 0; j < lignes; j++) {
            res[i * lignes + j] = temp[i * lignes + j];
        }
    }

    #pragma omp parallel for
    for (int i = 0; i < lignes; i++) {
        for (int j = 0; j < colonnes; j++) {
            complexe_float_t c = {res[i * lignes + j].real, -res[i * lignes + j].imaginary};
            res[i * lignes + j] = c;
        }
    }
}

void transpose_comp_d_omp(complexe_double_t *matrice, int lignes, int colonnes, complexe_double_t *res){
    complexe_double_t temp[lignes * colonnes];

    #pragma omp parallel for
    for (int i = 0; i < lignes; i++) {
        for (int j = 0; j < colonnes; j++) {
            temp[j * lignes + i] = matrice[i * colonnes + j];
        }
    }

    #pragma omp parallel for
    for (int i = 0; i < colonnes; i++) {
        for (int j = 0; j < lignes; j++) {
            res[i * lignes + j] = temp[i * lignes + j];
        }
    }

    #pragma omp parallel for
    for (int i = 0; i < lignes; i++) {
        for (int j = 0; j < colonnes; j++) {
            complexe_double_t c = {res[i * lignes + j].real, -res[i * lignes + j].imaginary};
            res[i * lignes + j] = c;
        }
    }
}

void mncblas_sgemm_omp(MNCBLAS_LAYOUT layout, MNCBLAS_TRANSPOSE TransA,
                 MNCBLAS_TRANSPOSE TransB, const int M, const int N,
                 const int K, const float alpha, const float *A,
                 const int lda, const float *B, const int ldb,
                 const float beta, float *C, const int ldc){

 /*   if(TransA == 112){
        transpose_f(A,M,N,R);
    }
    if(TransA == 112){
        transpose_f(B,M,N,R);
    }*/
    register unsigned int i,j,k;

    #pragma omp parallel for
    for(i=0;i<M;i+=1){
        for(j=0;j<N;j++){
        C[i*M+j] = C[i*M+j]*beta;
        }
    }

    #pragma omp parallel for
    for(i=0;i<N;i+=1){
        for(k=0;k<N;k+=1){
            for(j=0;j<N;j+=1){
                C[i*M+j] += A[i*M+k] * B[k*N+j]*alpha; 
            }
        }
    }
    
}

void mncblas_dgemm_omp(MNCBLAS_LAYOUT layout, MNCBLAS_TRANSPOSE TransA,
                 MNCBLAS_TRANSPOSE TransB, const int M, const int N,
                 const int K, const double alpha, const double *A,
                 const int lda, const double *B, const int ldb,
                 const double beta, double *C, const int ldc){
    /*if(TransA == 112){
        transpose_d(A,M,N);
    }
    if(TransA == 112){
        transpose_d(B,M,N);
    }*/
    register unsigned int i,j,k;

    #pragma omp parallel for
    for(i=0;i<M;i+=1){
        for(j=0;j<N;j++){
        C[i*M+j] = C[i*M+j]*beta;
        }
    }

    #pragma omp parallel for
    for(i=0;i<M;i+=1){
        for(k=0;k<K;k+=1){
            for(j=0;j<N;j+=1){
                C[i*M+j] += A[i*M+k] * B[k*N+j]*alpha; 
            }
        }
    }
    
                 }

void mncblas_cgemm_omp(MNCBLAS_LAYOUT layout, MNCBLAS_TRANSPOSE TransA,
                 MNCBLAS_TRANSPOSE TransB, const int M, const int N,
                 const int K, const void *alpha, const void *A,
                 const int lda, const void *B, const int ldb,
                 const void *beta, void *C, const int ldc){
    register unsigned int i,j,k;
    complexe_float_t *a2 = (complexe_float_t*)A;
    complexe_float_t *b2 = (complexe_float_t*)B;
    complexe_float_t *c2 = (complexe_float_t*)C;
    register complexe_float_t alpha2 = *((complexe_float_t*)alpha);
    register complexe_float_t beta2 = *((complexe_float_t*)beta);

   /* if(TransA == 112){
        transpose_f2(A,M,N);
    }
    if(TransA == 112){
        transpose_f2(B,M,N);
    }
    if(TransA == 113){
        transpose_comp_f(A,M,N);
    }
    if(TransA == 113){
        transpose_comp_f(B,M,N);
    }*/

    #pragma omp parallel for
    for(i=0;i<M;i+=1){
        for(j=0;j<N;j++){
        c2[i*M+j] = mult_complexe_float(c2[i*M+j],beta2);
        }
    }

    #pragma omp parallel for
    for(i=0;i<N;i+=1){
        for(k=0;k<K;k+=1){
            for(j=0;j<N;j+=1){
                complexe_float_t ABa = mult_complexe_float(mult_complexe_float(a2[i*M+k], b2[k*N+j]), alpha2);
                c2[i*M+j] = add_complexe_float(c2[i*M+j], ABa); 
            }
        }
    }
    
                 }

void mncblas_zgemm_omp(MNCBLAS_LAYOUT layout, MNCBLAS_TRANSPOSE TransA,
                 MNCBLAS_TRANSPOSE TransB, const int M, const int N,
                 const int K, const void *alpha, const void *A,
                 const int lda, const void *B, const int ldb,
                 const void *beta, void *C, const int ldc){
    register unsigned int i,j,k;
    complexe_double_t *a2 = (complexe_double_t*)A;
    complexe_double_t *b2 = (complexe_double_t*)B;
    complexe_double_t *c2 = (complexe_double_t*)C;
    register complexe_double_t alpha2 = *((complexe_double_t*)alpha);
    register complexe_double_t beta2 = *((complexe_double_t*)beta);

/*    if(TransA == 112){
        transpose_d2(A,M,N);
    }
    if(TransA == 112){
        transpose_d2(B,M,N);
    }
    if(TransA == 113){
        transpose_comp_d(A,M,N);
    }
    if(TransA == 113){
        transpose_comp_d(B,M,N);
    }
*/
    #pragma omp parallel for
    for(i=0;i<M;i+=1){
        for(j=0;j<N;j++){
        c2[i*M+j] = mult_complexe_double(c2[i*M+j],beta2);
        }
    }

    #pragma omp parallel for
    for(i=0;i<N;i+=1){
        for(k=0;k<K;k+=1){
            for(j=0;j<N;j+=1){
                complexe_double_t ABa = mult_complexe_double(mult_complexe_double(a2[i*M+k], b2[k*N+j]), alpha2);
                c2[i*M+j] = add_complexe_double(c2[i*M+j], ABa); 
            }
        }
    }
    
}

