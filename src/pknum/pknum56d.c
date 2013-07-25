
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2013                                  */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

#include "pkvaria.h"
#include "pknum.h"

void pkn_SymMatFindEigenvalueIntervald ( int n, double *a,
                                         double *lmin, double *lmax )
{
  int    i, j, k;
  double r, s, d, _lmin, _lmax;

  _lmin = _lmax = a[0];
  if ( n > 1 ) {
    r = fabs ( a[1] );
    for ( j = 2, k = 3;  j < n;  k += ++j )
      r += fabs ( a[k] );
    _lmin = a[0] - r;
    _lmax = a[0] + r;
    for ( i = 1; i < n; i++ ) {
      k = pkn_LowerTrMatIndex ( i, 0 );
      r = fabs ( a[k++] );
      for ( j = 1;  j < i;  j++ )
        r += fabs ( a[k++] );
      d = a[k];
      for ( j = i+1, k += j;  j < n;  k += ++j )
        r += fabs ( a[k] );
      s = d - r;  if ( s < _lmin ) _lmin = s;
      s = d + r;  if ( s > _lmax ) _lmax = s;
    }
  }
  *lmin = _lmin;
  *lmax = _lmax;
} /*pkn_SymMatFindEigenvalueIntervald*/

