#include <stdio.h>

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

void test_mncblas_sswap() {
    float X[] = {1.0, 2.0, 3.0, 4.0};
    float Y[] = {5.0, 6.0, 7.0, 8.0};
    int N = 4;

    Voir_float(X,Y,4);

    mncblas_sswap(N, X, 1, Y, 1);

    Voir_float(X,Y,4);

    for (int i = 0; i < N; i++) {
        if (X[i] != 5.0 + i || Y[i] != 1.0 + i) {
            printf("Erreur dans mncblas_sswap à l'indice %d.\n", i);
            exit(1);
        }
    }

    printf("mncblas_sswap réussi.\n");
}




void test_mncblas_dswap() {
    double X[] = {1.0, 2.0, 3.0, 4.0};
    double Y[] = {5.0, 6.0, 7.0, 8.0};
    int N = 4;

    Voir_double(X,Y,4);

    mncblas_dswap(N, X, 1, Y, 1);

    Voir_double(X,Y,4);

    printf("mncblas_dswap réussi.\n");
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




void test_mncblas_cswap() {
    complexe_float_t X[3];
    complexe_float_t Y[3];
    X[0].real = 1, X[0].imaginary = 1, X[1].real = 2, X[1].imaginary = 2, X[2].real = 3, X[2].imaginary = 3;
    Y[0].real = 4, Y[0].imaginary = 4, Y[1].real = 5, Y[1].imaginary = 5, Y[2].real = 6, Y[2].imaginary = 6;

    int N = 3;
    Voir_complexe_float(X,Y,3);

    mncblas_cswap(N, X, 1, Y, 1);

    Voir_complexe_float(X,Y,3);


    printf("mncblas_cswap réussi.\n");
}

void test_mncblas_zswap() {
    complexe_double_t X[3];
    complexe_double_t Y[3];
    X[0].real = 1, X[0].imaginary = 1, X[1].real = 2, X[1].imaginary = 2, X[2].real = 3, X[2].imaginary = 3;
    Y[0].real = 4, Y[0].imaginary = 4, Y[1].real = 5, Y[1].imaginary = 5, Y[2].real = 6, Y[2].imaginary = 6;
    int N=3;

    Voir_complexe_double(X,Y,3);

    mncblas_zswap(N, X, 1, Y, 1);

    Voir_complexe_double(X,Y,3);

    printf("mncblas_zswap réussi.\n");
}









int main() {
    test_mncblas_sswap();
    test_mncblas_dswap();
    test_mncblas_cswap();
    test_mncblas_zswap();

    printf("Tous les tests ont réussi.\n");
    return 0;
}