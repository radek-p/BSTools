
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


void pkn_multiBandmMultInvUTMVectord ( int nrows,
                                       const bandm_profile *rprof,
                                       const double *r,
                                       int spdimen, const double *x, double *y )
{
  /* compute the product y = R^{-1}x, where R is an upper triangular band   */
  /* matrix and x, y are matrices nrows x spdimen, packed row by row.       */

  int    i, j, k, fr, ir, ix, jx;
  double a;

  if ( x != y )
    memcpy ( y, x, nrows*spdimen*sizeof(double) );
  for ( i = nrows-1, ix = i*spdimen;  i >= 0;  i--, ix -= spdimen ) {
    fr = rprof[i].firstnz;
    ir = rprof[i].ind-fr;
    a = r[ir+i];
    for ( k = 0; k < spdimen; k++ )
      y[ix+k] /= a;
    for ( j = fr, jx = j*spdimen;  j < i;  j++, jx += spdimen ) {
      a = r[ir+j];
      for ( k = 0; k < spdimen; k++ )
        y[jx+k] -= a*y[ix+k];
    }
  }
} /*pkn_multiBandmMultInvUTMVectord*/

