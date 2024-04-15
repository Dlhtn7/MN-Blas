# une bibliotheque statique (libmnblas.a)
# une bibliotheque dynamique (libmnblasdyn.so)

all: libmnblas.a libmnblasdyn.so

LIB_DIR=../lib
INC_DIR=../include
OPTIM=-O2 -fopenmp

OPTIONS=-O2 -Wall -fPIC $(OPTIM) -I$(INC_DIR)
# Bibliothèque statique

libmnblas.a: swap_OMP.o dot_OMP.o copy_OMP.o axpy_OMP.o gemv_OMP.o gemm_OMP.o gemm.o gemv.o asum.o axpy.o copy.o swap.o dot.o complexe.o
	rm -f libmnblas.a $(LIB_DIR)/libmnblas.a
	ar -r libmnblas.a swap_OMP.o dot_OMP.o copy_OMP.o axpy_OMP.o gemv_OMP.o gemm_OMP.o gemm.o gemv.o asum.o axpy.o copy.o swap.o dot.o complexe.o
	cp libmnblas.a $(LIB_DIR)

# Bibliothèque dynamique

libmnblasdyn.so: swap_OMP.o dot_OMP.o copy_OMP.o axpy_OMP.o gemv_OMP.o gemm_OMP.o gemm.o gemv.o asum.o axpy.o copy.o swap.o dot.o complexe.o
	rm -f libmnblasdyn.so $(LIB_DIR)/libmnblasdyn.so
	gcc -shared -o libmnblasdyn.so swap_OMP.o dot_OMP.o copy_OMP.o axpy_OMP.o gemv_OMP.o gemm_OMP.o gemm.o gemv.o asum.o axpy.o copy.o swap.o dot.o complexe.o
	cp libmnblasdyn.so $(LIB_DIR)

complexe.o: $(INC_DIR)/complexe.h complexe.c
	gcc $(OPTIONS) -c complexe.c

# BLAS 1

copy.o : copy.c $(INC_DIR)/mnblas.h
	gcc $(OPTIONS) -c copy.c

swap.o : swap.c $(INC_DIR)/mnblas.h
	gcc  $(OPTIONS) -c swap.c

dot.o : dot.c $(INC_DIR)/mnblas.h
	gcc  $(OPTIONS) -c dot.c

axpy.o : axpy.c $(INC_DIR)/mnblas.h
	gcc  $(OPTIONS) -c axpy.c


asum.o : asum.c $(INC_DIR)/mnblas.h
	gcc  $(OPTIONS) -c asum.c

# BLAS 2
gemv.o : gemv.c $(INC_DIR)/mnblas.h
	gcc  $(OPTIONS) -c gemv.c
# BLAS 3

gemm.o : gemm.c $(INC_DIR)/mnblas.h
	gcc  $(OPTIONS) -c gemm.c

# OpenMP
copy_OMP.o : copy_OMP.c $(INC_DIR)/mnblas.h
	gcc $(OPTIONS) -c copy_OMP.c


dot_OMP.o : dot_OMP.c $(INC_DIR)/mnblas.h
	gcc  $(OPTIONS) -c dot_OMP.c

axpy_OMP.o : axpy_OMP.c $(INC_DIR)/mnblas.h
	gcc  $(OPTIONS) -c axpy_OMP.c

swap_OMP.o : swap_OMP.c $(INC_DIR)/mnblas.h
	gcc  $(OPTIONS) -c swap_OMP.c

gemv_OMP.o : gemv_OMP.c $(INC_DIR)/mnblas.h
	gcc  $(OPTIONS) -c gemv_OMP.c

gemm_OMP.o : gemm_OMP.c $(INC_DIR)/mnblas.h
	gcc  $(OPTIONS) -c gemm_OMP.c

clean:
	rm -f *.o *.a *.so *~  $(INC_DIR)/*~ $(LIB_DIR)/libmnblasdyn.so  $(LIB_DIR)/libmnblas.a ../lib/libmnblas.a ../lib/libmnblasdyn.so



