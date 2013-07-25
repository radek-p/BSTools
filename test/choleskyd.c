
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2005                                  */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#include "pkvaria.h"
#include "pknum.h"

double a[6] = { 1.0,
                2.0, 5.0,
                3.0, 8.0, 14.0 };
double b[3] = { 1.0, -2.0, 3.0 };
double l[6], x[3];
double f[9];

int main ()
{
  pkv_InitScratchMem ( 65536 );
  memcpy ( l, a, 6*sizeof(double) );
  printf ( "%f\n%f %f\n%f %f %f\n\n", l[0], l[1], l[2], l[3], l[4], l[5] );
  pkn_CholeskyDecompd ( 3, l );
  printf ( "%f\n%f %f\n%f %f %f\n\n", l[0], l[1], l[2], l[3], l[4], l[5] );

  pkn_SymMatrixMultd ( 3, a, 1, 1, b, 1, x );
  printf ( "x: %f, %f, %f\n\n", x[0], x[1], x[2] );

  pkn_LowerTrMatrixMultd ( 3, l, 1, 1, b, 1, x );
  printf ( "x: %f, %f, %f\n\n", x[0], x[1], x[2] );

  pkn_UpperTrMatrixMultd ( 3, l, 1, 1, b, 1, x );
  printf ( "x: %f, %f, %f\n\n", x[0], x[1], x[2] );

  pkn_SymMatrixMultd ( 3, a, 1, 1, b, 1, x );
  pkn_LowerTrMatrixSolved ( 3, l, 1, 1, x, 1, x );
  pkn_UpperTrMatrixSolved ( 3, l, 1, 1, x, 1, x );
  printf ( "x: %f, %f, %f\n\n", x[0], x[1], x[2] );

  pkn_SymToFullMatrixd ( 3, a, 3, f );
  printf ( "%f %f %f\n", f[0], f[1], f[2] );
  printf ( "%f %f %f\n", f[3], f[4], f[5] );
  printf ( "%f %f %f\n\n", f[6], f[7], f[8] );

  pkn_FullToSymMatrixd ( 3, 3, f, l );
  printf ( "%f\n%f %f\n%f %f %f\n\n", l[0], l[1], l[2], l[3], l[4], l[5] );

  pkn_LTrToFullMatrixd ( 3, a, 3, f );
  printf ( "%f %f %f\n", f[0], f[1], f[2] );
  printf ( "%f %f %f\n", f[3], f[4], f[5] );
  printf ( "%f %f %f\n\n", f[6], f[7], f[8] );
/*
  pkn_UTrToFullMatrixd ( 3, a, 3, f );
  printf ( "%f %f %f\n", f[0], f[1], f[2] );
  printf ( "%f %f %f\n", f[3], f[4], f[5] );
  printf ( "%f %f %f\n\n", f[6], f[7], f[8] );
*/
  pkn_FullToUTrMatrixd ( 3, 3, f, l );
  printf ( "%f\n%f %f\n%f %f %f\n\n", l[0], l[1], l[2], l[3], l[4], l[5] );

  pkv_DestroyScratchMem ();
  exit ( 0 );
} /*main*/

