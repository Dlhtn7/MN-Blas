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
mat_f m1;
mat_d m2;
mat_comp_f m3;
mat_comp_d m4;


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
     vector_init (vec1, 1.0) ;
     vector_init (vec2, 2.0) ;
     m1_init (m1, 1.0);
     
     TOP_NANO (start) ;
        mncblas_sgemv_omp (101,111,VECSIZE,VECSIZE,1.0,m1,1,vec1,1,1.0,vec2,1) ;
     TOP_NANO (end);
     temp_moy = temp_moy+diff_nano(&start,&end);

     
     
  }
  printf ("sgemv => nb flottants: %d en %e seconde\n", VECSIZE*3+VECSIZE*VECSIZE*2,temp_moy/NB_FOIS) ;


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
     m2_init(m2, 1.0);
     TOP_NANO (start) ;
        mncblas_dgemv_omp (101,111,VECSIZE,VECSIZE,1.0,m2,1,d1,1,1.0,d2,1) ;
     TOP_NANO (end);
     temp_moy = temp_moy+diff_nano(&start,&end);
  }
  printf ("dgemv => nb flottants: %d en %e seconde\n", VECSIZE*3+VECSIZE*VECSIZE*2,temp_moy/NB_FOIS) ;



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
     vector_init_complexe_f (v1, c1) ;
     vector_init_complexe_f (v2, c2) ;
     m3_init(m3, c3);
     
     
     TOP_NANO (start) ;
        mncblas_cgemv_omp (101,111,VECSIZE,VECSIZE,&a,m3,1,v1,1,&b,v2,1) ;
     TOP_NANO (end);
     temp_moy = temp_moy+diff_nano(&start,&end);
  }
  printf ("cgemv => nb flottants: %d en %e seconde\n", VECSIZE*14+VECSIZE*VECSIZE*8,temp_moy/NB_FOIS) ;


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
     vector_init_complexe_d (v3, c3) ;
     vector_init_complexe_d (v4, c4) ;
     m4_init(m4, c5);
     

     TOP_NANO (start) ;
        mncblas_zgemv_omp (101,111,VECSIZE,VECSIZE,&a,m4,1,v3,1,&b,v4,1) ;
     TOP_NANO (end);
     temp_moy = temp_moy+diff_nano(&start,&end);
  }
  printf ("zgemv => nb flottants: %d en %e seconde\n", VECSIZE*14+VECSIZE*VECSIZE*8,temp_moy/NB_FOIS) ;

}



int main(void){
    s();
    d();
    c();
    z();
    return 0;
}