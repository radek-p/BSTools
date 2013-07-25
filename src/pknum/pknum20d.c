
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
/* multiplication of matrices, per coefficients */

void pkn_MultArrayd ( int nrows, int rowlen, int pitch_a, CONST_ double *a,
                      int pitch_b, CONST_ double *b,
                      int pitch_c, double *c )
{
  int    i, j;
  double *aa, *bb, *cc;

  for ( i = 0, aa = a, bb = b, cc = c;
        i < nrows;
        i++, aa += pitch_a, bb += pitch_b, cc += pitch_c )
    for ( j = 0; j < rowlen; j++ )
      cc[j] = aa[j]*bb[j];
} /*pkn_MultArrayd*/

