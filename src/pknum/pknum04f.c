
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


void pkn_BandmQRDecomposeMatrixf ( int nrows, int ncols,
                                   const bandm_profile *aprof,
                                   const float *a,
                                   bandm_profile *qprof, float *q,
                                   bandm_profile *rprof, float *r )
{
  int    i, j, k, fr, lr, k0, k1, k2, k3;
  float  *ac;
  double an, wn, b;
  int    size_ac;

  ac = (float*)pkv_GetScratchMem ( size_ac = nrows*sizeof(float) );

  qprof[0].ind = ncols;
  rprof[0].ind = 0;

  for ( i = 0; i < ncols; i++ ) {
           /* reconstruct the first column */
    fr = aprof[i].firstnz;                    /* first row with nonzero */
    lr = fr + (aprof[i+1].ind-aprof[i].ind);  /* last row with nonzero+1 */
    memset ( ac, 0, nrows*sizeof(float) );
    memcpy ( &ac[fr], &a[aprof[i].ind], (lr-fr)*sizeof(float) );

           /* apply the previously constructed Householder reflections */
    for ( j = 0; j < i; j++ ) {
      k0 = qprof[j].firstnz;
      k1 = k0 + qprof[j+1].ind-qprof[j].ind;
      if ( k0 < lr && k1 > fr) {
              /* do the reflection */
        k2 = max ( k0, fr );
        k3 = min ( k1, lr );
        wn = q[j]*pkn_ScalarProductf( k3-k2,
                      &q[qprof[j].ind+k2-k0], &ac[k2] );
        for ( k = 0; k < k1-k0; k++ )
          ac[k0+k] -= (float)(wn*q[qprof[j].ind+k]);

        fr = min ( fr, k0 );
      }
    }
    fr = min ( fr, i );

           /* construct the Householder reflection for the i-th column */
             /* copy data to destination arrays */
    memcpy ( &r[rprof[i].ind], &ac[fr], (i-fr)*sizeof(float) );
    memcpy ( &q[qprof[i].ind+1], &ac[i+1], (lr-i-1)*sizeof(float) );

             /* compute the norm of the column, in an */
    wn = pkn_ScalarProductf ( lr-i-1, &ac[i+1], &ac[i+1] );
    an = sqrt(wn+ac[i]*ac[i]);

             /* construct the reflection hyperplane normal vector */
    if ( ac[i] <= 0 )
         b = r[rprof[i].ind+i-fr] = (float)an;
    else b = r[rprof[i].ind+i-fr] = (float)(-an);
    b = q[qprof[i].ind] = (float)(ac[i]-b);
    wn += b*b;
    q[i] = (float)(2.0/wn);

             /* construct the profiles */
    qprof[i].firstnz = i;
    qprof[i+1].ind = qprof[i].ind + lr - i;
    rprof[i].firstnz = fr;
    rprof[i+1].ind = rprof[i].ind + i + 1 - fr;
  }
  qprof[ncols].firstnz = rprof[ncols].firstnz = 0;

  pkv_FreeScratchMem ( size_ac );
} /*pkn_BandmQRDecomposeMatrixf*/

