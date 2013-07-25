
/* //////////////////////////////////////////////////// */
/* This file is a part of the BSTools procedure package */
/* written by Przemyslaw Kiciak.                        */
/* //////////////////////////////////////////////////// */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "pkvaria.h"
#include "pknum.h"

double a[8] =
  { 1,  0,
    1, -3,
    1, -3,
    1,  0 };
double b[8] = { -7, 14, 5, -10, 2, -4, 14, -28 };
double x[4];

int main ()
{
  int i;

  pkv_InitScratchMem ( 65536 );
  pkn_multiSolveRLSQd ( 4, 2, a, 2, 2, b, 2, x );
  for ( i = 0; i < 4; i++ )
    printf ( "%f\n", x[i] );
  pkv_DestroyScratchMem ();
  exit ( 0 );
} /*main*/

