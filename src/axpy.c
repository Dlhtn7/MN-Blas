#include "mnblas.h"
#include "complexe2.h"

void mnblas_saxpy(const int N, const float alpha, const float *X,const int incX, float *Y, const int incY){
    register unsigned int i = 0 ;
    register unsigned int j = 0 ;

    for (; ((i < N) && (j < N)) ; i += incX, j += incY){
      Y[j] += alpha*X[i] ;
    }
}

void mnblas_daxpy(const int N, const double alpha, const double *X,
                 const int incX, double *Y, const int incY){
    register unsigned int i = 0 ;
    register unsigned int j = 0 ;

    for (; ((i < N) && (j < N)) ; i += incX, j += incY){
      Y[j] += alpha*X[i] ;
    }
}

void mnblas_caxpy(const int N, const void *alpha, const void *X,
                 const int incX, void *Y, const int incY){
    register unsigned int i = 0 ;
    register unsigned int j = 0 ;
    complexe_float_t *x2 = (complexe_float_t*)X;
    complexe_float_t *y2 = (complexe_float_t*)Y;
    register complexe_float_t a2 = *((complexe_float_t*)alpha);

    for (; ((i < N) && (j < N)) ; i += incX, j += incY)
    {
      y2[j] = add_complexe_float(y2[j], mult_complexe_float(a2,x2[i]));
    }
}

void mnblas_zaxpy(const int N, const void *alpha, const void *X,
                 const int incX, void *Y, const int incY){
    register unsigned int i = 0 ;
    register unsigned int j = 0 ;
    complexe_double_t *x2 = (complexe_double_t*)X;
    complexe_double_t *y2 = (complexe_double_t*)Y;
    register complexe_double_t a2 = *((complexe_double_t*)alpha);

    for (; ((i < N) && (j < N)) ; i += incX, j += incY)
    {
      y2[j] = add_complexe_double(y2[j], mult_complexe_double(a2,x2[i]));
    }
}
