
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


void pkn_multiBandmMultTrVectorf ( int ncols,
                                   const bandm_profile *aprof, const float *a,
                                   int spdimen, const float *x, float *y )
{
  /* compute the product y = R^Tx, where R is a band matrix     */
  /* and x is the matrix nrows x spdimen, y is ncols x spdimen, */
  /* packed row by row.                                         */

  int i, j, k, fr, lr, jr, ix, kx;

  memset ( y, 0, ncols*spdimen*sizeof(float) );
  for ( i = ix = 0;  i < ncols;  i++, ix += spdimen ) {
    fr = aprof[i].firstnz;
    lr = fr + aprof[i+1].ind - aprof[i].ind;
    jr = aprof[i].ind-fr;
    for ( j = fr, kx = fr*spdimen;  j < lr;  j++, kx += spdimen )
      for ( k = 0; k < spdimen; k++ )
        y[ix+k] += a[jr+j]*x[kx+k];
  }
} /*pkn_multiBandmMultTrVectorf*/

