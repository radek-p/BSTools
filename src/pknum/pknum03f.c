
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2005                                  */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

#include <stdio.h>
#include <string.h>
#include <math.h>

#include "pkvaria.h"
#include "pknum.h"


double pkn_detf ( int n, float *a )
{
  int   i, j, k, i0, j0;
  float m, s, d;
                              /* Gaussian elimination with pivoting */
  d = 1.0;
  for ( k = 0; k < n-1; k++ ) {
                              /* full pivoting */
    m = (float)fabs ( a[n*k+k] );  i0 = k;  j0 = k;
    for ( i = k; i < n; i++ )
      for ( j = k; j < n; j++ )
        if ( (s = (float)fabs ( a[n*i+j] )) > m ) { m = s;  i0 = i;  j0 = j; }
    if ( i0 != k ) {
      for ( j = k; j < n; j++ )
        { s = a[n*k+j];  a[n*k+j] = a[n*i0+j];  a[n*i0+j] = s; }
      d = -d;
    }
    if ( j0 != k ) {
      for ( i = k; i < n; i++ )
        { s = a[n*i+k];  a[n*i+k] = a[n*i+j0];  a[n*i+j0] = s; }
      d = -d;
    }
    if ( !a[n*k+k] )
      return 0.0;
                              /* elimination */
    s = a[n*k+k];
    for ( i = k+1; i < n; i++ ) {
      m = a[n*i+k]/s;
      for ( j = k+1; j < n; j++ )
        a[n*i+j] -= m*a[n*k+j];
    }
    d *= s;                   /* multiply by the diagonal element */
  }
  return d*a[n*n-1];       /* do not forget about the last one */
} /*pkn_detf*/

