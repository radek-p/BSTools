
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2011                                  */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "pkvaria.h"
#include "pknum.h"

/* /////////////////////////////////////////////////////////////////////////// */
boolean pkn_MultSPMVectorf ( int nrows, int ncols, int nnz,
                             const index2 *ai, const float *ac,
                             int spdimen, const float *x,
                             float *y )
{
  int   k, i, j, l;
  float aij;

  memset ( y, 0, spdimen*nrows*sizeof(float) );
  for ( k = 0; k < nnz; k++ ) {
    i = ai[k].i;
    j = ai[k].j;
    if ( i >= 0 && i < nrows && j >= 0 && j < ncols ) {
      aij = ac[k];
      i *= spdimen;
      j *= spdimen;
      for ( l = 0; l < spdimen; l++ )
        y[i+l] += aij*x[j+l];
    }
    else
      return false;
  }
  return true;
} /*pkn_MultSPMVectorf*/

boolean pkn_MultSPMTVectorf ( int nrows, int ncols, int nnz,
                              const index2 *ai, const float *ac,
                              int spdimen, const float *x,
                              float *y )
{
  int   k, i, j, l;
  float aij;

  memset ( y, 0, spdimen*ncols*sizeof(float) );
  for ( k = 0; k < nnz; k++ ) {
    i = ai[k].i;
    j = ai[k].j;
    if ( i >= 0 && i < nrows && j >= 0 && j < ncols ) {
      aij = ac[k];
      i *= spdimen;
      j *= spdimen;
      for ( l = 0; l < spdimen; l++ )
        y[j+l] += aij*x[i+l];
    }
    else
      return false;
  }
  return true;
} /*pkn_MultSPMTVectorf*/

