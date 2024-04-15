#include <stdio.h>

#include "mnblas.h"
#include "complexe2.h"

#include "flop.h"

#define VECSIZE    2048

#define NB_FOIS    10

typedef float vfloat [VECSIZE] ;
typedef double vdoublr [VECSIZE] ;
typedef complexe_float_t vfcomp [VECSIZE] ;
typedef complexe_double_t vdcomp [VECSIZE] ; 
typedef float mat_f [VECSIZE*VECSIZE];
typedef double mat_d [VECSIZE*VECSIZE];
typedef complexe_float_t mat_comp_f [VECSIZE*VECSIZE];
typedef complexe_double_t mat_comp_d [VECSIZE*VECSIZE];


vfloat vec1, vec2 ;
vdoublr d1, d2;
vfcomp v1, v2;
vdcomp v3, v4;
mat_f m1, m12, m13;
mat_d m2, m21,m22;
mat_comp_f m3,m31,m32;
mat_comp_d m4,m41,m42;


void m1_init (mat_f V, float x)
{
  register unsigned int i ;

  for (i = 0; i < VECSIZE*VECSIZE; i++)
    V [i] = x ;

  return ;
}

void m2_init (mat_d V, double x)
{
  register unsigned int i ;

  for (i = 0; i < VECSIZE*VECSIZE; i++)
    V [i] = x ;

  return ;
}

void m3_init (mat_comp_f V, complexe_float_t x)
{
  register unsigned int i ;

  for (i = 0; i < VECSIZE*VECSIZE; i++)
    V [i] = x ;

  return ;
}

void m4_init (mat_comp_d V, complexe_double_t x)
{
  register unsigned int i ;

  for (i = 0; i < VECSIZE*VECSIZE; i++)
    V [i] = x ;

  return ;
}

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
     m1_init (m1, 1.0) ;
     m1_init (m12, 2.0) ;
     m1_init (m13, 1.0);
     
     TOP_NANO (start) ;
        mncblas_sgemm_omp (101,111,111,VECSIZE,VECSIZE,VECSIZE,1.0,m1,1,m12,1,1.0,m13,1) ;
     TOP_NANO (end);
     temp_moy = temp_moy+diff_nano(&start,&end);

     
     
  }
  printf ("sgemm => %e seconde\n",temp_moy/NB_FOIS) ;


}

void d(){
  struct timespec start, end ;
  //double res2;
  int i ;
  init_nano () ;
  double temp_moy = 0.0;
  for (i = 0 ; i < NB_FOIS; i++){
     m2_init (m2, 1.0) ;
     m2_init (m21, 2.0) ;
     m2_init(m22, 3.0);
     TOP_NANO (start) ;
        mncblas_dgemm_omp (101,111,111,VECSIZE,VECSIZE,VECSIZE,1.0,m2,1,m21,1,1.0,m22,1) ;
     TOP_NANO (end);
     temp_moy = temp_moy+diff_nano(&start,&end);
  }
  printf ("dgemm => %e seconde\n", /*VECSIZE*VECSIZE*1+VECSIZE*VECSIZE*VECSIZE*3,*/temp_moy/NB_FOIS) ;



}

void c(){
  struct timespec start, end ;
  int i ;
  init_nano () ;
  complexe_float_t c1 = {10.0,7.0};
  complexe_float_t c2 = {25.0,32.0}; 
  complexe_float_t c3 = {1.0, 1.0};
  complexe_float_t a = (complexe_float_t){1.0,1.0};
  complexe_float_t b = (complexe_float_t){1.0,1.0};
  double temp_moy = 0.0;
  for (i = 0 ; i < NB_FOIS; i++){
     m3_init (m3, c1) ;
     m3_init (m31, c2) ;
     m3_init(m32, c3);
     
     
     TOP_NANO (start) ;
        mncblas_cgemm_omp (101,111,111,VECSIZE,VECSIZE,VECSIZE,&a,m3,1,m31,1,&b,m32,1) ;
     TOP_NANO (end);
     temp_moy = temp_moy+diff_nano(&start,&end);
  }
  printf ("cgemm => en %e seconde\n",temp_moy/NB_FOIS) ;


}

void z(){
  struct timespec start, end ;
  int i ;
  init_nano () ;
  complexe_double_t c3 = {10.0,7.0};
  complexe_double_t c4 = {25.0,32.0}; 
  complexe_double_t c5 = {1.0, 1.0};
  complexe_double_t a = (complexe_double_t){1.0,1.0};
  complexe_double_t b = (complexe_double_t){1.0,1.0};
  double temp_moy = 0.0;
  for (i = 0 ; i < NB_FOIS; i++){
     m4_init (m4, c3) ;
     m4_init (m41, c4) ;
     m4_init(m42, c5);
     

     TOP_NANO (start) ;
        mncblas_zgemm_omp (101,111,111,VECSIZE,VECSIZE,VECSIZE,&a,m4,1,m41,1,&b,m42,1) ;
     TOP_NANO (end);
     temp_moy = temp_moy+diff_nano(&start,&end);
  }
  printf ("zgemm ==> %e seconde\n", temp_moy/NB_FOIS) ;

}



int main(void){
    s();
    d();
   c();
   z();
    return 0;
}