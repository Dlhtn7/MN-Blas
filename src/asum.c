#include "mnblas.h"
#include "complexe2.h"


float abs_float(float x){
    if (x<0)
        return -x;
    else
        return x;
}

float abs_double(double x){
    if (x<0)
        return -x;
    else
        return x;
}

float  mnblas_sasum(const int N, const float *X, const int incX){
    register unsigned int i = 0 ;
    float res = 0.0;

    for (; ((i < N)) ; i += incX){
        res = res + abs_float(X[i]);
    }
    return res;
}

double mnblas_dasum(const int N, const double *X, const int incX){
    register unsigned int i = 0 ;
    double res = 0.0;

    for (; ((i < N)) ; i += incX){
        res = res + abs_double(X[i]);

    }
    return res;
}

float  mnblas_scasum(const int N, const void *X, const int incX){
    register unsigned int i = 0 ;
    complexe_float_t *x2 = (complexe_float_t*)X;
    float res = 0.0;

    for (; ((i < N)) ; i += incX){
        res = res + abs_float(x2[i].real) + abs_float(x2[i].imaginary);

    }
    return res;
}

double mnblas_dzasum(const int N, const void *X, const int incX){
    register unsigned int i = 0 ;
    complexe_double_t *x2 = (complexe_double_t*)X;
    double res = 0.0;

    for (; ((i < N)) ; i += incX){
        res = res + abs_double(x2[i].real) + abs_double(x2[i].imaginary);

    }
    return res;
}
