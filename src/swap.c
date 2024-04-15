#include "mnblas.h"
#include "complexe2.h"

void mncblas_sswap(const int N, float *X, const int incX, 
                 float *Y, const int incY)
{
  register unsigned int i = 0 ;
  register unsigned int j = 0 ;
  register float save ;
  
  for (; ((i < N) && (j < N)) ; i += incX, j+=incY)
    {
      save = Y [j] ;
      Y [j] = X [i] ;
      X [i] = save ;
    }

  return ;
}

void mncblas_dswap(const int N, double *X, const int incX, 
                 double *Y, const int incY)
{

  register unsigned int i = 0 ;
  register unsigned int j = 0 ;
  register double save ;
  
  for (; ((i < N) && (j < N)) ; i += incX, j+=incY)
    {
      save = Y [j] ;
      Y [j] = X [i] ;
      X [i] = save ;
    }

  return ;
}

void mncblas_cswap(const int N, void *X, const int incX, 
		                    void *Y, const int incY)
{
    register unsigned int i = 0 ;
    register unsigned int j = 0 ;
    register complexe_float_t *x2 = (complexe_float_t*)X ;
    register complexe_float_t *y2 = (complexe_float_t*)Y;
    register complexe_float_t save;

    for (; ((i < N) && (j < N)) ; i += incX, j+=incY)
    {
      save = y2[j] ;
      y2[j] = x2[i] ;
      x2[i] = save ;
    }

  return ;
}

void mncblas_zswap(const int N, void *X, const int incX, 
		                    void *Y, const int incY)
{
    register unsigned int i = 0 ;
    register unsigned int j = 0 ;
    register complexe_double_t *x2 = (complexe_double_t*)X ;
    register complexe_double_t *y2 = (complexe_double_t*)Y;
    register complexe_double_t save;

    for (; ((i < N) && (j < N)) ; i += incX, j+=incY)
    {
      save = y2[j] ;
      y2[j] = x2[i] ;
      x2[i] = save ;
    }

  return ;
}

