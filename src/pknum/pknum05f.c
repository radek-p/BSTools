
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


void pkn_multiBandmReflectVectorf ( int ncols,
                                    const bandm_profile *qprof, const float *q,
                                    int spdimen, float *b )
{
  int    i, j, k, k0, fr, lr;
  double wn;
  float  *bb, *bbb;

  for ( i = 0; i < ncols; i++ ) {
    fr = qprof[i].firstnz;
    k0 = qprof[i].ind;
    lr = qprof[i+1].ind-k0;
    for ( j = 0, bb = b;  j < spdimen;  j++, bb++ ) {
      wn = 0.0;
      for ( k = 0, bbb = bb + (fr*spdimen);  k < lr;  k++, bbb += spdimen )
        wn += q[k0+k]*(*bbb);
      wn *= q[i];
      for ( k = 0, bbb = bb + (fr*spdimen);  k < lr;  k++, bbb += spdimen )
        *bbb -= (float)(wn*q[k0+k]);
    }
  }
} /*pkn_multiBandmReflectVectorf*/

