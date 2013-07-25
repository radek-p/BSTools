
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

float a[6] = { 1.0,
               2.0, 5.0,
               3.0, 8.0, 14.0 };
float b[3] = { 1.0, -2.0, 3.0 };
float l[6], x[3];
float f[9];

int main ()
{
  pkv_InitScratchMem ( 65536 );
  memcpy ( l, a, 6*sizeof(float) );
  printf ( "%f\n%f %f\n%f %f %f\n\n", l[0], l[1], l[2], l[3], l[4], l[5] );
  pkn_CholeskyDecompf ( 3, l );
  printf ( "%f\n%f %f\n%f %f %f\n\n", l[0], l[1], l[2], l[3], l[4], l[5] );

  pkn_SymMatrixMultf ( 3, a, 1, 1, b, 1, x );
  printf ( "x: %f, %f, %f\n\n", x[0], x[1], x[2] );

  pkn_LowerTrMatrixMultf ( 3, l, 1, 1, b, 1, x );
  printf ( "x: %f, %f, %f\n\n", x[0], x[1], x[2] );

  pkn_UpperTrMatrixMultf ( 3, l, 1, 1, b, 1, x );
  printf ( "x: %f, %f, %f\n\n", x[0], x[1], x[2] );

  pkn_SymMatrixMultf ( 3, a, 1, 1, b, 1, x );
  pkn_LowerTrMatrixSolvef ( 3, l, 1, 1, x, 1, x );
  pkn_UpperTrMatrixSolvef ( 3, l, 1, 1, x, 1, x );
  printf ( "x: %f, %f, %f\n\n", x[0], x[1], x[2] );

  pkn_SymToFullMatrixf ( 3, a, 3, f );
  printf ( "%f %f %f\n", f[0], f[1], f[2] );
  printf ( "%f %f %f\n", f[3], f[4], f[5] );
  printf ( "%f %f %f\n\n", f[6], f[7], f[8] );

  pkn_FullToSymMatrixf ( 3, 3, f, l );
  printf ( "%f\n%f %f\n%f %f %f\n\n", l[0], l[1], l[2], l[3], l[4], l[5] );

  pkn_LTrToFullMatrixf ( 3, a, 3, f );
  printf ( "%f %f %f\n", f[0], f[1], f[2] );
  printf ( "%f %f %f\n", f[3], f[4], f[5] );
  printf ( "%f %f %f\n\n", f[6], f[7], f[8] );
/*
  pkn_UTrToFullMatrixf ( 3, a, 3, f );
  printf ( "%f %f %f\n", f[0], f[1], f[2] );
  printf ( "%f %f %f\n", f[3], f[4], f[5] );
  printf ( "%f %f %f\n\n", f[6], f[7], f[8] );
*/
  pkn_FullToUTrMatrixf ( 3, 3, f, l );
  printf ( "%f\n%f %f\n%f %f %f\n\n", l[0], l[1], l[2], l[3], l[4], l[5] );

  pkv_DestroyScratchMem ();
  exit ( 0 );
} /*main*/

