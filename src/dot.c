#include "mnblas.h"
#include <stdio.h>
#include "complexe2.h"


float mncblas_sdot(const int N, const float *X, const int incX, 
                 const float *Y, const int incY)
{
  register unsigned int i = 0 ;
  register unsigned int j = 0 ;
  register float dot = 0.0 ;
  
  for (; ((i < N) && (j < N)) ; i += incX, j+=incY)
    {
      dot = dot + X [i] * Y [j] ;
    }

  return dot ;
}

/*
float mncblas_sdot(const int N, const float *X, const int incX, 
                 const float *Y, const int incY)
{
  register unsigned int i = 0 ;
  register unsigned int j = 0 ;
  float dot = 0.0 ;

  
  for (i = 0 ; i < N ; i += incX)
    {
      dot += X [i] * Y [j] ;
      j+=incY ;
    }

  return dot ;
}*/

double mncblas_ddot(const int N, const double *X, const int incX, 
                 const double *Y, const int incY)
{
  register unsigned int i = 0 ;
  register unsigned int j = 0 ;
  register float dot = 0.0 ;
  
  for (; ((i < N) && (j < N)) ; i += incX, j+=incY)
    {
      dot = dot + X [i] * Y [j] ;
    }

  return dot ;
}

void   mncblas_cdotu_sub(const int N, const void *X, const int incX,
                       const void *Y, const int incY, void *dotu)
{
  register unsigned int i = 0 ;
  register unsigned int j = 0 ;
  complexe_float_t *x2 = (complexe_float_t*)X;
  complexe_float_t *y2 = (complexe_float_t*)Y;
  complexe_float_t *dot = (complexe_float_t*)dotu;

  
  for (; ((i < N) && (j < N)) ; i += incX, j += incY)
    {
      *dot= add_complexe_float(*dot, mult_complexe_float(x2[i], y2[j]));
      
    }

  return ;
}

void   mncblas_cdotc_sub(const int N, const void *X, const int incX,
                       const void *Y, const int incY, void *dotc)
{
  register unsigned int i = 0 ;
  register unsigned int j = 0 ;
  complexe_float_t *x2 = (complexe_float_t*)X;
  complexe_float_t *y2 = (complexe_float_t*)Y;
  complexe_float_t *dot = (complexe_float_t*)dotc;

  
  for (; ((i < N) && (j < N)) ; i += incX, j += incY)
    {
      complexe_float_t c = {x2[i].real, -x2[i].imaginary};
      *dot= add_complexe_float(*dot, mult_complexe_float(c, y2[j]));
      
    }

  return ;
}

void   mncblas_zdotu_sub(const int N, const void *X, const int incX,
                       const void *Y, const int incY, void *dotu)
{
  register unsigned int i = 0 ;
  register unsigned int j = 0 ;
  complexe_double_t *x2 = (complexe_double_t*)X;
  complexe_double_t *y2 = (complexe_double_t*)Y;
  complexe_double_t *dot = (complexe_double_t*)dotu;

  
  for (; ((i < N) && (j < N)) ; i += incX, j += incY)
    {
      *dot= add_complexe_double(*dot, mult_complexe_double(x2[i], y2[j]));
      
    }

  return ;
}
  
void   mncblas_zdotc_sub(const int N, const void *X, const int incX,
                       const void *Y, const int incY, void *dotc)
{
  register unsigned int i = 0 ;
  register unsigned int j = 0 ;
  complexe_double_t *x2 = (complexe_double_t*)X;
  complexe_double_t *y2 = (complexe_double_t*)Y;
  complexe_double_t *dot = (complexe_double_t*)dotc;

  for (; ((i < N) && (j < N)) ; i += incX, j += incY)
    {
      complexe_double_t c = {x2[i].real, -x2[i].imaginary};
      *dot= add_complexe_double(*dot, mult_complexe_double(c, y2[j]));
      
    }

  return ;
}
