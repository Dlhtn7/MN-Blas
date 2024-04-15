#include <stdio.h>
#include <stdlib.h>
#include "mnblas.h"
#include "complexe2.h"


#define    NB_FOIS        512

#include "flop.h"

int main (int argc, char **argv)
{
 complexe_float_t c1= {1.0, 2.0} ;
 complexe_float_t c2= {3.0, 6.0} ;

 complexe_double_t cd1 ;
 complexe_double_t cd2 ;

 struct timespec start, end ;
 int i ;

 init_nano () ;
 
 c1 = add_complexe_float (c1, c2) ;

 printf ("c1.r %f c1.i %f\n", c1.real, c1.imaginary) ;

 cd1 = (complexe_double_t) {10.0, 7.0} ;
 cd2 = (complexe_double_t) {25.0, 32.0} ;


 TOP_NANO(start);
 
 for (i = 0 ; i < NB_FOIS; i++)
   {
     cd1 = add_complexe_double (cd1, cd2) ;
   }

 TOP_NANO(end) ;

 printf ("apres boucle cd1.real %f cd1.imaginary %f\nTemps %e seconde \n", cd1.real, cd1.imaginary, diff_nano(&start, &end)) ;


   TOP_NANO(start) ;

  for (i = 0 ; i < NB_FOIS; i++)
   {
     c1 = add_complexe_float (c1, c2) ;
   }

 TOP_NANO(end) ;

 printf ("apres boucle c1.real %f c1.imaginary %f\n Temps %e seconde \n", c1.real, c1.imaginary, diff_nano (&start, &end)) ;


complexe_float_t c11= {1.0, 2.0} ;
    TOP_NANO(start) ;

  for (i = 0 ; i < NB_FOIS; i++)
   {
     c11 = mult_complexe_float (c11, c2) ;
   }

 TOP_NANO(end) ;

 printf ("apres boucle c11.real %f c11.imaginary %f\n Temps %e seconde \n", c11.real, c11.imaginary, diff_nano (&start, &end)) ;


  complexe_double_t cd11 = (complexe_double_t) {10.0, 7.0} ;
 TOP_NANO(start) ;

  for (i = 0 ; i < NB_FOIS; i++)
   {
     cd11 = mult_complexe_double (cd11, cd2) ;
   }

 TOP_NANO(end) ;

 printf ("apres boucle cd11.real %f cd11.imaginary %f\n Temps %e seconde \n", cd11.real, cd11.imaginary, diff_nano (&start, &end)) ;


complexe_float_t c12= {1.0, 2.0} ;
    TOP_NANO(start) ;

  for (i = 0 ; i < NB_FOIS; i++)
   {
     c12 = div_complexe_float (c12, c2) ;
   }

 TOP_NANO(end) ;

 printf ("apres boucle c12.real %f c12.imaginary %f\n Temps %e seconde \n", c12.real, c12.imaginary, diff_nano (&start, &end)) ;


  complexe_double_t cd12 = (complexe_double_t) {10.0, 7.0} ;
 TOP_NANO(start) ;

  for (i = 0 ; i < NB_FOIS; i++)
   {
     cd12 = div_complexe_double (cd12, cd2) ;
   }

 TOP_NANO(end) ;

 printf ("apres boucle cd12.real %f cd12.imaginary %f\n Temps %e seconde \n", cd12.real, cd12.imaginary, diff_nano (&start, &end)) ;


 exit (0) ;
}





