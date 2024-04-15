#include <stdio.h>
#include <stdlib.h>
#include "mnblas.h"
#include "complexe2.h"
#include "flop.h"


void Voir_float(float* X, float* Y, int len){
    for(int i = 0; i<len; i++){
        printf("X[%d] = %f\n",i,X[i]);
    }
    for(int i = 0; i<len; i++){
        printf("Y[%d] = %f\n",i,Y[i]);
    }
}

void Voir_double(double* X, double* Y, int len){
    for(int i = 0; i<len; i++){
        printf("X[%d] = %f\n",i,X[i]);
    }
    for(int i = 0; i<len; i++){
        printf("Y[%d] = %f\n",i,Y[i]);
    }
}

void Voir_complexe_float(complexe_float_t* X, complexe_float_t* Y, int len){
    for(int i = 0; i<len; i++){
        printf("X[%d].imag = %f, X[%d].real= %f\n",i,X[i].imaginary,i,X[i].real);
    }
    for(int i = 0; i<len; i++){
        printf("Y[%d].imag = %f, Y[%d].real= %f\n",i,Y[i].imaginary,i,Y[i].real);
    }
}

void Voir_complexe_double(complexe_double_t* X, complexe_double_t* Y, int len){
    for(int i = 0; i<len; i++){
        printf("X[%d].imag = %f, X[%d].real= %f\n",i,X[i].imaginary,i,X[i].real);
    }
    for(int i = 0; i<len; i++){
        printf("Y[%d].imag = %f, Y[%d].real= %f\n",i,Y[i].imaginary,i,Y[i].real);
    }
}




void test_mncblas_scopy() {
    float X[] = {1.0, 2.0, 3.0, 4.0};
    float Y[] = {5.0, 6.0, 7.0, 8.0};
    int N = 4;

    Voir_float(X,Y,4);

    mncblas_scopy(N, X, 1, Y, 1);

    Voir_float(X,Y,4);


    printf("scopy réussi.\n");
}

void test_mncblas_dcopy() {
    double X[] = {1.0, 2.0, 3.0, 4.0};
    double Y[] = {5.0, 6.0, 7.0, 8.0};
    int N = 4;

    Voir_double(X,Y,4);

    mncblas_dcopy(N, X, 1, Y, 1);

    Voir_double(X,Y,4);


    printf("dcopy réussi.\n");
}

void test_mncblas_ccopy() {
    complexe_float_t X[3];
    complexe_float_t Y[3];
    X[0].real = 1, X[0].imaginary = 1, X[1].real = 2, X[1].imaginary = 2, X[2].real = 3, X[2].imaginary = 3;
    Y[0].real = 4, Y[0].imaginary = 4, Y[1].real = 5, Y[1].imaginary = 5, Y[2].real = 6, Y[2].imaginary = 6;

    int N = 3;
    Voir_complexe_float(X,Y,3);

    mncblas_ccopy(N, X, 1, Y, 1);

    Voir_complexe_float(X,Y,3);


    printf("mncblas_cswap réussi.\n");
}

void test_mncblas_zcopy() {
    complexe_double_t X[3];
    complexe_double_t Y[3];
    X[0].real = 1, X[0].imaginary = 1, X[1].real = 2, X[1].imaginary = 2, X[2].real = 3, X[2].imaginary = 3;
    Y[0].real = 4, Y[0].imaginary = 4, Y[1].real = 5, Y[1].imaginary = 5, Y[2].real = 6, Y[2].imaginary = 6;

    int N = 3;
    Voir_complexe_double(X,Y,3);

    mncblas_zcopy(N, X, 1, Y, 1);

    Voir_complexe_double(X,Y,3);


    printf("mncblas_cswap réussi.\n");
}



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
        mncblas_scopy (VECSIZE, vec1, 1, vec2, 1) ;
     TOP_NANO (end);
     temp_moy = temp_moy+diff_nano(&start,&end);

     
     
  }
  printf ("scopy => nb flottants: %f en %e seconde\n", 0.0,temp_moy/NB_FOIS) ;


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
        mncblas_dcopy (VECSIZE, d1, 1, d2, 1) ;
     TOP_NANO (end);
     temp_moy = temp_moy+diff_nano(&start,&end);
  }
  printf ("dcopy => nb flottants: %f en %e seconde\n", 0.0,temp_moy/NB_FOIS) ;

}

void c(){
  struct timespec start, end ;
  int i ;
  init_nano () ;
  complexe_float_t c1 = {10.0,7.0};
  complexe_float_t c2 = {25.0,32.0}; 
  double temp_moy = 0.0;
  for (i = 0 ; i < NB_FOIS; i++){
     vector_init_complexe_f (v1, c1) ;
     vector_init_complexe_f (v2, c2) ;
     
     
     TOP_NANO (start) ;
        mncblas_ccopy (VECSIZE, v1, 1, v2, 1) ;
     TOP_NANO (end);
     temp_moy = temp_moy+diff_nano(&start,&end);
  }
  printf ("ccopy => nb flottants: %f en %e seconde\n", 0.0,temp_moy/NB_FOIS) ;

}

void z(){
  struct timespec start, end ;
  int i ;
  init_nano () ;
  complexe_double_t c3 = {10.0,7.0};
  complexe_double_t c4 = {25.0,32.0}; 
  double temp_moy = 0.0;
  for (i = 0 ; i < NB_FOIS; i++){
     vector_init_complexe_d (v3, c3) ;
     vector_init_complexe_d (v4, c4) ;
     

     TOP_NANO (start) ;
        mncblas_zcopy (VECSIZE, v3, 1, v4, 1) ;
     TOP_NANO (end);
     temp_moy = temp_moy+diff_nano(&start,&end);
  }
  printf ("zcopy => nb flottants: %f en %e seconde\n", 0.0,temp_moy/NB_FOIS) ;

}

int main(){
    test_mncblas_scopy();
    test_mncblas_dcopy();
    test_mncblas_ccopy();
    test_mncblas_zcopy();

    s();
    d();
    c();
    z();


    return 0;
}