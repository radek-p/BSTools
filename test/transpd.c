
/* //////////////////////////////////////////////////// */
/* This file is a part of the BSTools procedure package */
/* written by Przemyslaw Kiciak.                        */
/* //////////////////////////////////////////////////// */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "pkvaria.h"
#include "pknum.h"

double a[12] = {0.0, 1.0, 2.0,
                3.0, 4.0, 5.0,
                6.0, 7.0, 8.0,
                9.0,10.0,11.0};
double b[15], c[20];

int main ()
{
  int i, j;

  printf ( "a:\n" );
  for ( i = 0; i < 4; i++ ) {
    for ( j = 0; j < 3; j++ )
      printf ( "%6.1f", a[3*i+j] );
    printf ( "\n" );
  }
  printf ( "\n" );
  pkv_TransposeMatrixd ( 4, 3, 3, a, 5, b );
  printf ( "b:\n" );
  for ( i = 0; i < 3; i++ ) {
    for ( j = 0; j < 4; j++ )
      printf ( "%6.1f", b[5*i+j] );
    printf ( "\n" );
  }
  printf ( "\n" );
  pkn_MultMatrixd ( 4, 3, 3, a, 4, 5, b, 5, c );
  printf ( "c:\n" );
  for ( i = 0; i < 4; i++ ) {
    for ( j = 0; j < 4; j++ )
      printf ( "%6.1f", c[5*i+j] );
    printf ( "\n" );
  }
  printf ( "\n" );

  exit ( 0 );
} /*main*/

