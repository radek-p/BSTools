
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2005                                  */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

#include <math.h>
#include <string.h>

#include "pkvaria.h"

#undef CONST_
#define CONST_

#include "pknum.h"


/* /////////////////////////////////////////// */
/* rectangular matrix multiplication           */

void pkn_MultMatrixf ( int nrows_a, int rowlen_a, int pitch_a, CONST_ float *a,
                       int rowlen_b, int pitch_b, CONST_ float *b,
                       int pitch_c, float *c )
{
  int   i, k;
  float *aa, *bb, *cc;

  for ( i = 0, aa = a, cc = c;
        i < nrows_a;
        i++, aa += pitch_a, cc += pitch_c ) {
    bb = b;
    pkn_MultMatrixNumf ( 1, rowlen_b, 0, bb, aa[0], 0, cc );
    for ( k = 1, bb += pitch_b;
          k < rowlen_a;
          k++, bb += pitch_b )
      pkn_AddMatrixMf ( 1, rowlen_b, 0, cc, 0, bb, aa[k], 0, cc );
  }
} /*pkn_MultMatrixf*/

void pkn_MultMatrixAddf ( int nrows_a, int rowlen_a, int pitch_a, CONST_ float *a,
                          int rowlen_b, int pitch_b, CONST_ float *b,
                          int pitch_c, float *c )
{
  int   i, k;
  float *aa, *bb, *cc;

  for ( i = 0, aa = a, cc = c;
        i < nrows_a;
        i++, aa += pitch_a, cc += pitch_c ) {
    for ( k = 0, bb = b;
          k < rowlen_a;
          k++, bb += pitch_b )
      pkn_AddMatrixMf ( 1, rowlen_b, 0, cc, 0, bb, aa[k], 0, cc );
  }
} /*pkn_MultMatrixAddf*/

void pkn_MultMatrixSubf ( int nrows_a, int rowlen_a, int pitch_a, CONST_ float *a,
                          int rowlen_b, int pitch_b, CONST_ float *b,
                          int pitch_c, float *c )
{
  int   i, k;
  float *aa, *bb, *cc;

  for ( i = 0, aa = a, cc = c;
        i < nrows_a;
        i++, aa += pitch_a, cc += pitch_c ) {
    for ( k = 0, bb = b;
          k < rowlen_a;
          k++, bb += pitch_b )
      pkn_AddMatrixMf ( 1, rowlen_b, 0, cc, 0, bb, -aa[k], 0, cc );
  }
} /*pkn_MultMatrixSubf*/

