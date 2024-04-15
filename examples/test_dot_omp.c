#include <stdio.h>

#include "mnblas.h"
#include "complexe2.h"

#include "flop.h"

#define VECSIZE    65536

#define NB_FOIS    10

typedef float vfloat [VECSIZE] ;
typedef double vdoublr [VECSIZE] ;
typedef complexe_float_t vfcomp [VECSIZE] ;
typedef complexe_double_t vdcomp [VECSIZE] ; 

vfloat vec1, vec2 ;
vdoublr d1, d2;
vfcomp v1, v2;
vdcomp v3, v4;

void vector_init (vfloat V, float x)
{
  register unsigned int i ;

  for (i = 0; i < VECSIZE; i++)
    V [i] = x ;

  return ;
}

void vector_print (vfloat V)
{
  register unsigned int i ;

  for (i = 0; i < VECSIZE; i++)
    printf ("%f ", V[i]) ;
  printf ("\n") ;
  
  return ;
}

void vector_init_d (vdoublr V, double x)
{
  register unsigned int i ;

  for (i = 0; i < VECSIZE; i++)
    V [i] = x ;

  return ;
}

void vector_print_d (vdoublr V)
{
  register unsigned int i ;

  for (i = 0; i < VECSIZE; i++)
    printf ("%f ", V[i]) ;
  printf ("\n") ;
  
  return ;
}

void vector_init_complexe_f (vfcomp V, complexe_float_t x)
{
  register unsigned int i ;

  for (i = 0; i < VECSIZE; i++)
    V [i] = x ;

  return ;
}

void vector_print_complexe_f (vfcomp V)
{
  register unsigned int i ;

  for (i = 0; i < VECSIZE; i++)
    printf ("%f, %f ", V[i].real, V[i].imaginary) ;
  printf ("\n") ;
  
  return ;
}

void vector_init_complexe_d (vdcomp V, complexe_double_t x)
{
  register unsigned int i ;

  for (i = 0; i < VECSIZE; i++)
    V [i] = x ;

  return ;
}

void vector_print_complexe_d (vdcomp V)
{
  register unsigned int i ;

  for (i = 0; i < VECSIZE; i++)
    printf ("%f, %f ", V[i].real, V[i].imaginary) ;
  printf ("\n") ;
  
  return ;
}

void s(){
  struct timespec start, end ;
  
  int i ;
  init_nano () ;
  double temp_moy = 0.0;
  for (i = 0 ; i < NB_FOIS; i++){
     vector_init (vec1, 1.0) ;
     vector_init (vec2, 2.0) ;
     
     TOP_NANO (start) ;
        mncblas_sdot_omp (VECSIZE, vec1, 1, vec2, 1) ;
     TOP_NANO (end);
     temp_moy = temp_moy+diff_nano(&start,&end);

     
     
  }
  printf ("sdot => nb flottants: %d en %e seconde\n", VECSIZE*2,temp_moy/NB_FOIS) ;


}

void d(){
  struct timespec start, end ;
  //double res2;
  int i ;
  init_nano () ;
  double temp_moy = 0.0;
  for (i = 0 ; i < NB_FOIS; i++){
     vector_init_d (d1, 1.0) ;
     vector_init_d (d2, 2.0) ;
     TOP_NANO (start) ;
        mncblas_ddot_omp (VECSIZE, d1, 1, d2, 1) ;
     TOP_NANO (end);
     temp_moy = temp_moy+diff_nano(&start,&end);
  }
  printf ("ddot => nb flottants: %d en %e seconde\n", VECSIZE*2,temp_moy/NB_FOIS) ;

}

void c(){
  struct timespec start, end ;
  complexe_float_t res3 = (complexe_float_t){0.0,0.0};
  int i ;
  init_nano () ;
  complexe_float_t c1 = {10.0,7.0};
  complexe_float_t c2 = {25.0,32.0}; 
  double temp_moy = 0.0;
  for (i = 0 ; i < NB_FOIS; i++){
     vector_init_complexe_f (v1, c1) ;
     vector_init_complexe_f (v2, c2) ;
     
     
     TOP_NANO (start) ;
        mncblas_cdotu_sub_omp (VECSIZE, v1, 1, v2, 1, &res3) ;
     TOP_NANO (end);
     temp_moy = temp_moy+diff_nano(&start,&end);
  }
  printf ("cdot => nb flottants: %d en %e seconde\n", VECSIZE*8,temp_moy/NB_FOIS) ;

}

void z(){
  struct timespec start, end ;
  complexe_double_t res4 = (complexe_double_t){0.0,0.0};
  int i ;
  init_nano () ;
  complexe_double_t c3 = {10.0,7.0};
  complexe_double_t c4 = {25.0,32.0}; 
  double temp_moy = 0.0;
  for (i = 0 ; i < NB_FOIS; i++){
     vector_init_complexe_d (v3, c3) ;
     vector_init_complexe_d (v4, c4) ;
     

     TOP_NANO (start) ;
        mncblas_zdotu_sub_omp (VECSIZE, v3, 1, v4, 1, &res4) ;
     TOP_NANO (end);
     temp_moy = temp_moy+diff_nano(&start,&end);
  }
  printf ("zdot => nb flottants: %d en %e seconde\n", VECSIZE*8,temp_moy/NB_FOIS) ;

}











int main (void)
{
  s();
  d();
  c();
  z();
  return 0;
}
