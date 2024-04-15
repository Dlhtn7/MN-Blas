#include "mnblas.h"
#include "complexe2.h"
#include <omp.h>

void mncblas_sgemv_omp(const MNCBLAS_LAYOUT layout,
                 const MNCBLAS_TRANSPOSE TransA, const int M, const int N,
                 const float alpha, const float *A, const int lda,
                 const float *X, const int incX, const float beta,
                 float *Y, const int incY){
                    int i,j;
                  //  if(TransA == 111){
                        #pragma omp parallel for
                        for(i=0;i<M;i+=incX){
                            float Ax = 0.0;
                            for(j=0;j<N;j+=incY){
                                Ax = Ax + A[i*N+j]*X[j];
                            }
                            Y[i] = Ax*alpha+ Y[i]*beta;
                        }
                 /*   }if (TransA == 112){
                        transpose_f(A,M,N);
                        #pragma omp parallel for
                        for(i=0;i<M;i+=incX){
                            float Ax = 0.0;
                            for(j=0;j<N;j+=incY){
                                Ax = Ax + A[i*N+j]*X[j];
                            }
                            Y[i] = Ax*alpha+ Y[i]*beta;
                        }
                    }*/
                 }

void mncblas_dgemv_omp(MNCBLAS_LAYOUT layout,
                 MNCBLAS_TRANSPOSE TransA, const int M, const int N,
                 const double alpha, const double *A, const int lda,
                 const double *X, const int incX, const double beta,
                 double *Y, const int incY){
                    int i,j;
               //     if(TransA == 111){
                        #pragma omp parallel for
                        for(i=0;i<M;i+=incX){
                            double Ax = 0.0;
                            for(j=0;j<N;j+=incY){
                                Ax = Ax + A[i*N+j]*X[j];
                            }
                            Y[i] = Ax*alpha+ Y[i]*beta;
                        }
                /*    }if (TransA == 112){
                        transpose_d(A,M,N);
                        #pragma omp parallel for
                        for(i=0;i<M;i+=incX){
                            double Ax = 0.0;
                            for(j=0;j<N;j+=incY){
                                Ax = Ax + A[i*N+j]*X[j];
                            }
                            Y[i] = Ax*alpha+ Y[i]*beta;
                        }
                    }*/
                 }

void mncblas_cgemv_omp(MNCBLAS_LAYOUT layout,
                 MNCBLAS_TRANSPOSE TransA, const int M, const int N,
                 const void *alpha, const void *A, const int lda,
                 const void *X, const int incX, const void *beta,
                 void *Y, const int incY){
                    int i,j;
                    register complexe_float_t * A2 = (complexe_float_t* ) A;
                    register complexe_float_t * X2 = (complexe_float_t*) X;
                    register complexe_float_t * Y2 = (complexe_float_t*) Y;
                    register complexe_float_t alpha2 = *((complexe_float_t*) alpha);
                    register complexe_float_t beta2 = *((complexe_float_t*) beta);

             //   if(TransA == 111){
                    #pragma omp parallel for
                    for (i=0;i<M;i+=incX){
                        complexe_float_t Ax = (complexe_float_t){0.0,0.0};
                        for(j=0;j<N;j+=incY){
                        Ax = add_complexe_float(Ax, mult_complexe_float(A2[i*N+j], X2[j]));
                        }
                        Y2[i] = add_complexe_float(mult_complexe_float(Y2[i], beta2) , mult_complexe_float(Ax, alpha2));
                    }
            /*    } if (TransA == 112){
                    transpose_f2(A,M,N);
                    #pragma omp parallel for
                    for (i=0;i<M;i+=incX){
                        complexe_float_t Ax = (complexe_float_t){0.0,0.0};
                        for(j=0;j<N;j+=incY){
                        Ax = add_complexe_float(Ax, mult_complexe_float(A2[i*N+j], X2[j]));
                        }
                        Y2[i] = add_complexe_float(mult_complexe_float(Y2[i], beta2) , mult_complexe_float(Ax, alpha2));
                    }

                 } if (TransA == 113){
                    transpose_comp_f(A,M,N);
                    #pragma omp parallel for
                    for (i=0;i<M;i+=incX){
                        complexe_float_t Ax = (complexe_float_t){0.0,0.0};
                        for(j=0;j<N;j+=incY){
                        Ax = add_complexe_float(Ax, mult_complexe_float(A2[i*N+j], X2[j]));
                        }
                        Y2[i] = add_complexe_float(mult_complexe_float(Y2[i], beta2) , mult_complexe_float(Ax, alpha2));
                    }
                 }*/
                }

void mncblas_zgemv_omp(MNCBLAS_LAYOUT layout,
                 MNCBLAS_TRANSPOSE TransA, const int M, const int N,
                 const void *alpha, const void *A, const int lda,
                 const void *X, const int incX, const void *beta,
                 void *Y, const int incY){int i,j;
                    register complexe_double_t * A2 = (complexe_double_t* ) A;
                    register complexe_double_t * X2 = (complexe_double_t*) X;
                    register complexe_double_t * Y2 = (complexe_double_t*) Y;
                    register complexe_double_t alpha2 = *((complexe_double_t*) alpha);
                    register complexe_double_t beta2 = *((complexe_double_t*) beta);

              //  if(TransA == 111){
                    #pragma omp parallel for
                    for (i=0;i<M;i+=incX){
                        complexe_double_t Ax = (complexe_double_t){0.0,0.0};
                        for(j=0;j<N;j+=incY){
                        Ax = add_complexe_double(Ax, mult_complexe_double(A2[i*N+j], X2[j]));
                        }
                        Y2[i] = add_complexe_double(mult_complexe_double(Y2[i], beta2) , mult_complexe_double(Ax, alpha2));
                    }
               /* } if (TransA == 112){
                    transpose_d2(A,M,N);
                    #pragma omp parallel for
                    for (i=0;i<M;i+=incX){
                        complexe_double_t Ax = (complexe_double_t){0.0,0.0};
                        for(j=0;j<N;j+=incY){
                        Ax = add_complexe_double(Ax, mult_complexe_double(A2[i*N+j], X2[j]));
                        }
                        Y2[i] = add_complexe_double(mult_complexe_double(Y2[i], beta2) , mult_complexe_double(Ax, alpha2));
                    }

                 } if (TransA == 113){
                    transpose_comp_d(A,M,N);
                    #pragma omp parallel for
                    for (i=0;i<M;i+=incX){
                        complexe_double_t Ax = (complexe_double_t){0.0,0.0};
                        for(j=0;j<N;j+=incY){
                        Ax = add_complexe_double(Ax, mult_complexe_double(A2[i*N+j], X2[j]));
                        }
                        Y2[i] = add_complexe_double(mult_complexe_double(Y2[i], beta2) , mult_complexe_double(Ax, alpha2));
                    }
                 }*/
                }
