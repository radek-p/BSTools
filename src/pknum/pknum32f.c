
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2007                                  */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

#include <math.h>
#include <string.h>
#include <stdlib.h>

#include "pkvaria.h"

#undef CONST_
#define CONST_

#include "pknum.h"

/* ///////////////////////////////////////////////////////////// */
/* compute B = B - A*A^T, where B is a packed symmetric matrix,  */
/* and A is rectangular (full) */

void pkn_SymMatSubAATf ( int n, float *b, int m, int pitch_a, CONST_ float *a )
{
  int    i, j, k, l;
  float  *a1, *a2;
  double s;

  for ( i = l = 0, a1 = a;  i < n;  i++, a1 += pitch_a )
    for ( j = 0, a2 = a;  j <= i;  j++, l++, a2 += pitch_a ) {
      s = a1[0]*a2[0];
      for ( k = 1; k < m; k++ )
        s += a1[k]*a2[k];
      b[l] -= (float)s;
    }
} /*pkn_SymMatSubAATf*/

