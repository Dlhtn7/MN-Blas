#include "mnblas.h"
#include "complexe2.h"
#include <omp.h>

void mncblas_scopy_omp(const int N, const float *X, const int incX, 
                 float *Y, const int incY)
{
  register unsigned int i;
  //register unsigned int j = 0 ;

  #pragma omp parallel for
  for (i=0; i < N; i += incX)    {
      Y [i] = X [i] ;
    }

  return ;
}

void mncblas_dcopy_omp(const int N, const double *X, const int incX, 
                 double *Y, const int incY)
{
  register unsigned int i ;
  #pragma omp parallel for
  for (i=0; i < N; i += incX)
    {
      Y [i] = X [i] ;
    }

  return ;

}

void mncblas_ccopy_omp(const int N, const void *X, const int incX, 
		                    void *Y, const int incY)
{
    register unsigned int i ;
    complexe_float_t *x2 = (complexe_float_t*)X;
    complexe_float_t *y2 = (complexe_float_t*)Y;

    #pragma omp parallel for
    for (i=0; i < N; i += incX)
    {
      y2[i] = x2[i] ;
    }

    return ;

}

void mncblas_zcopy_omp(const int N, const void *X, const int incX, 
		                    void *Y, const int incY)
{
    register unsigned int i ;
    complexe_double_t *x2 = (complexe_double_t*)X;
    complexe_double_t *y2 = (complexe_double_t*)Y;

    #pragma omp parallel for
    for (i=0; i < N; i += incX)    {
      y2[i] = x2[i] ;
    }

    return ;

}