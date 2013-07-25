
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


void pkn_multiBandmInvReflectVectord ( int ncols,
                                       const bandm_profile *qprof,
                                       const double *q,
                                       int spdimen, double *b )
{
  int    i, j, k, k0, fr, lr;
  double wn;
  double *bb, *bbb;

  for ( i = ncols-1; i >= 0; i-- ) {
    fr = qprof[i].firstnz;
    k0 = qprof[i].ind;
    lr = qprof[i+1].ind-k0;
    for ( j = 0, bb = b;  j < spdimen;  j++, bb++ ) {
      wn = 0.0;
      for ( k = 0, bbb = bb + (fr*spdimen);  k < lr;  k++, bbb += spdimen )
        wn += q[k0+k]*(*bbb);
      wn *= q[i];
      for ( k = 0, bbb = bb + (fr*spdimen);  k < lr;  k++, bbb += spdimen )
        *bbb -= wn*q[k0+k];
    }
  }
} /*pkn_multiBandmInvReflectVectord*/

