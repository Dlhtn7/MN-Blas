#include <stdio.h>
#include <stdlib.h>

#include "mnblas.h"
#include "complexe.h"

#include "flop.h"

#define    NB_FOIS        512

int main (int argc, char **argv)
{
 complexe_float_t c1 = {1.0, 2.0} ;
 complexe_float_t c2 = {3.0, 6.0} ;

 complexe_double_t cd1 = {1.0, 2.0};
 complexe_double_t cd2 = {3.0, 6.0};

 struct timespec start, end ;
 
 int i ;

 init_nano () ;
 
 c1 = add_complexe_float (c1, c2) ;

 printf ("AddF : c1.r %f c1.i %f\n", c1.real, c1.imaginary) ;

 cd1 = add_complexe_double (cd1, cd2) ;

 printf ("AddD : cd1.r %f cd1.i %f\n", cd1.real, cd1.imaginary) ;

 complexe_float_t c11 = mult_complexe_float (c1, c2) ;

 printf ("MultF : c1.r %f c1.i %f\n", c11.real, c11.imaginary) ;

 cd1 = mult_complexe_double (cd1, cd2) ;

 printf ("MultD : cd1.r %f cd1.i %f\n", cd1.real, cd1.imaginary) ;

 c1 = div_complexe_float (c1, c2) ;

 printf ("DivF : c1.r %f c1.i %f\n", c1.real, c1.imaginary) ;

 cd1 = div_complexe_double (cd1, cd2) ;

 printf ("DivD : cd1.r %f cd1.i %f\n", cd1.real, cd1.imaginary) ;



 cd1 = (complexe_double_t) {10.0, 7.0} ;
 cd2 = (complexe_double_t) {25.0, 32.0} ;


 TOP_NANO(start) ;
 
 for (i = 0 ; i < NB_FOIS; i++)
   {
     cd1 = add_complexe_double (cd1, cd2) ;
   }

 TOP_NANO(end) ;

  printf ("apres boucle cd1.real %f cd1.imaginary %fTemps %e seconde \n", cd1.real, cd1.imaginary, diff_nano (&start, &end)) ;

  TOP_NANO(start) ;
 
 for (i = 0 ; i < NB_FOIS; i++)
   {
     c1 = add_complexe_float (c1, c2) ;
   }

 TOP_NANO(end) ;

 printf ("apres boucle c1.real %f cd1.imaginary %fTemps %e seconde \n", c1.real, c1.imaginary, diff_nano (&start, &end)) ;

 exit (0) ;
}
