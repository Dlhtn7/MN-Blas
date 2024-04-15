#include "mnblas.h"
#include "complexe2.h"
#include <omp.h>

void mnblas_saxpy_omp(const int N, const float alpha, const float *X,const int incX, float *Y, const int incY){
    register unsigned int i = 0 ;
    #pragma omp parallel for
    for (i=0; i < N; i += incX){
      Y[i] += alpha*X[i] ;
    }
}

void mnblas_daxpy_omp(const int N, const double alpha, const double *X,
                 const int incX, double *Y, const int incY){
    register unsigned int i ;
    #pragma omp parallel for
    for (i=0; i < N; i += incX)
    {
      Y[i] += alpha*X[i] ;
    }
}

void mnblas_caxpy_omp(const int N, const void *alpha, const void *X,
                 const int incX, void *Y, const int incY){
    register unsigned int i;
    complexe_float_t *x2 = (complexe_float_t*)X;
    complexe_float_t *y2 = (complexe_float_t*)Y;
    register complexe_float_t a2 = *((complexe_float_t*)alpha);

    #pragma omp parallel for
    for (i=0; i < N; i += incX)
    {
      y2[i] = add_complexe_float(y2[i], mult_complexe_float(a2,x2[i]));
    }
}

void mnblas_zaxpy_omp(const int N, const void *alpha, const void *X,
                 const int incX, void *Y, const int incY){
    register unsigned int i = 0 ;
    complexe_double_t *x2 = (complexe_double_t*)X;
    complexe_double_t *y2 = (complexe_double_t*)Y;
    register complexe_double_t a2 = *((complexe_double_t*)alpha);

    #pragma omp parallel for
    for (i=0; i < N; i += incX)
    {
      y2[i] = add_complexe_double(y2[i], mult_complexe_double(a2,x2[i]));
    }
}
