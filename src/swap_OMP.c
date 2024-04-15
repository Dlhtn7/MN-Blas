#include "mnblas.h"
#include "complexe2.h"
#include <omp.h>


void mncblas_sswap_omp(const int N, float *X, const int incX, 
                 float *Y, const int incY)
{
  register unsigned int i ;
  register float save ;

  #pragma omp parallel for
  for (i=0; i < N; i += incX)

    {
      save = Y [i] ;
      Y [i] = X [i] ;
      X [i] = save ;
    }

  return ;
}

void mncblas_dswap_omp(const int N, double *X, const int incX, 
                 double *Y, const int incY)
{

  register unsigned int i;
  register double save ;
  
  #pragma omp parallel for
  for (i=0; i < N; i += incX)

    {
      save = Y [i] ;
      Y [i] = X [i] ;
      X [i] = save ;
    }

  return ;
}

void mncblas_cswap_omp(const int N, void *X, const int incX, 
		                    void *Y, const int incY)
{
    register unsigned int i ;
    register complexe_float_t *x2 = (complexe_float_t*)X ;
    register complexe_float_t *y2 = (complexe_float_t*)Y;
    register complexe_float_t save;


    #pragma omp parallel for
    for (i=0; i < N; i += incX)
    {
      save = y2[i] ;
      y2[i] = x2[i] ;
      x2[i] = save ;
    }

  return ;
}

void mncblas_zswap_omp(const int N, void *X, const int incX, 
		                    void *Y, const int incY)
{
    register unsigned int i ;
    register complexe_double_t *x2 = (complexe_double_t*)X ;
    register complexe_double_t *y2 = (complexe_double_t*)Y;
    register complexe_double_t save;


    #pragma omp parallel for
    for (i=0; i < N; i += incX)
    {
      save = y2[i] ;
      y2[i] = x2[i] ;
      x2[i] = save ;
    }

  return ;
}

