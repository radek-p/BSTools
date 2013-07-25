
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


void pkn_multiBandmMultInvTrUTMVectorf ( int nrows,
                                         const bandm_profile *rprof,
                                         const float *r,
                                         int spdimen, const float *x, float *y ) 
{
  /* compute the product y = R^{-T}x, where R is an upper triangular band */
  /* matrix and x, y are matrices nrows x spdimen, packed row by row.     */

  int   i, j, k, fr, ir, ix, jx;
  float a;

  if ( x != y )
    memcpy ( y, x, nrows*spdimen*sizeof(float) );
  for ( i = ix = 0;  i < nrows;  i++, ix += spdimen ) {
    fr = rprof[i].firstnz;
    ir = rprof[i].ind-fr;
    for ( j = fr, jx = j*spdimen;  j < i;  j++, jx += spdimen ) {
      a = r[ir+j];
      for ( k = 0; k < spdimen; k++ )
        y[ix+k] -= a*y[jx+k];
    }
    a = r[ir+i];
    for ( k = 0; k < spdimen; k++ )
      y[ix+k] /= a;
  }
} /*pkn_multiBandmMultInvTrUTMVectorf*/

