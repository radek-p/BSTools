
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2005                                  */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "pkvaria.h"
#include "pknum.h"

#include "msgpool.h"


/* /////////////////////////////////////////// */
/* orthogonal-triangluar matrix decomposition  */

boolean pkn_QRDecomposeMatrixd ( int nrows, int ncols, double *a, double *aa )
{
  int    i, j, k, l, m, n, p;
  double b, c;

  for ( j = m = 0;  j < ncols;  j++, m += ncols+1 ) {
                                 /* construct the Householder reflection */
    b = 0.0;
    for ( i = j, l = m;  i < nrows;  i++, l += ncols )
      b += a[l]*a[l];
    if ( !b )
      return false;
    c = sqrt ( b );
    if ( a[m] >= 0 ) {
      aa[j] = a[m] + c;
      b += b + 2.0*a[m]*c;
      a[m] = -c;
    }
    else {
      aa[j] = a[m] - c;
      b += b - 2.0*a[m]*c;
      a[m] = c;
    }
    aa[j+ncols] = c = 2.0/b;

                                 /* reflect the matrix columns */
    for ( k = j+1, l = m+1;  k < ncols;  k++, l++ ) {
      b = aa[j]*a[l];
      for ( i = j+1, p = m+ncols, n = l+ncols;
            i < nrows;
            i++, p += ncols, n+= ncols )
        b += a[p]*a[n];
      b *= c;
      a[l] -= b*aa[j];
      for ( i = j+1, p = m+ncols, n = l+ncols;
            i < nrows;
            i++, p += ncols, n+= ncols )
        a[n] -= b*a[p];
    }
  }
  return true;
} /*pkn_QRDecomposeMatrixd*/

void pkn_multiReflectVectord ( int nrows, int ncols,
                               const double *a, const double *aa,
                               int spdimen, int pitch, double *b )
{
  int    i, j, l;
  double c, d;

  for ( j = 0; j < ncols; j++ ) {
    c = aa[j+ncols];
    for ( l = 0; l < spdimen; l++ ) {
      d = aa[j]*b[j*pitch+l];
      for ( i = j+1; i < nrows; i++ )
        d += a[i*ncols+j]*b[i*pitch+l];
      d *= c;
      b[j*pitch+l] -= d*aa[j];
      for ( i = j+1; i < nrows; i++ )
        b[i*pitch+l] -= d*a[i*ncols+j];
    }
  }
} /*pkn_multiReflectVectord*/

void pkn_multiInvReflectVectord ( int nrows, int ncols,
                                  const double *a, const double *aa,
                                  int spdimen, int pitch, double *b )
{
  int    i, j, l;
  double c, d;

  for ( j = ncols-1; j >= 0; j-- ) {
    c = aa[j+ncols];
    for ( l = 0; l < spdimen; l++ ) {
      d = aa[j]*b[j*pitch+l];
      for ( i = j+1; i < nrows; i++ )
        d += a[i*ncols+j]*b[i*pitch+l];
      d *= c;
      b[j*pitch+l] -= d*aa[j];
      for ( i = j+1; i < nrows; i++ )
        b[i*pitch+l] -= d*a[i*ncols+j];
    }
  }
} /*pkn_multiInvReflectVectord*/

void pkn_multiMultUTVectord ( int nrows, const double *a,
                              int spdimen, int bpitch, double *b,
                              int xpitch, double *x )
{
  int    i, j, k;
  double c;

  for ( i = 0; i < nrows; i++ )
    for ( j = 0; j < spdimen; j++ ) {
      c = 0.0;
      for ( k = i; k < nrows; k++ )
        c += a[i*nrows+k]*b[j+k*bpitch];
      x[i*xpitch+j] = c;
    }
} /*pkn_multiMultUTVectord*/

void pkn_multiMultInvUTVectord ( int nrows, const double *a,
                                 int spdimen, int bpitch, double *b,
                                 int xpitch, double *x )
{
  int    i, j, k;
  double c;

  for ( i = nrows-1; i >= 0; i-- )
    for ( k = 0; k < spdimen; k++ ) {
      c = b[i*bpitch+k];
      for ( j = i+1; j < nrows; j++ )
        c -= a[i*nrows+j]*x[j*xpitch+k];
      x[i*xpitch+k] = c/a[i*(nrows+1)];
    }
} /*pkn_multiMultInvUTVectord*/

void pkn_multiMultTrUTVectord ( int nrows, const double *a,
                                int spdimen, int bpitch, double *b,
                                int xpitch, double *x )
{
  int    i, j, k;
  double c;

  for ( i = 0; i < nrows; i++ )
    for ( j = 0; j < spdimen; j++ ) {
      c = 0.0;
      for ( k = 0; k <= i; k++ )
        c += a[k*nrows+i]*b[j+k*bpitch];
      x[i*xpitch+j] = c;
    }
} /*pkn_multiMultTrUTVectord*/

void pkn_multiMultInvTrUTVectord ( int nrows, const double *a,
                                   int spdimen, int bpitch, double *b,
                                   int xpitch, double *x )
{
  int    i, j, k;
  double c;

  for ( i = 0; i < nrows; i++ )
    for ( j = 0; j < spdimen; j++ ) {
      c = b[i*bpitch+j];
      for ( k = 0; k < i; k++ )
        c -= a[k*nrows+i]*x[j+k*xpitch];
      x[i*xpitch+j] = c/a[i*(nrows+1)];  
    }
} /*pkn_multiMultInvTrUTVectord*/

boolean pkn_multiSolveRLSQd ( int nrows, int ncols, double *a,
                              int spdimen, int bpitch, double *b,
                              int xpitch, double *x )
{
  void   *sp;
  double *aa;

  sp = pkv_GetScratchMemTop ();
  aa = pkv_GetScratchMemd ( 2*ncols );
  if ( !aa ) {
    pkv_SignalError ( LIB_PKNUM, 0, ERRMSG_0 );
    exit ( 1 );
  }

  if ( !pkn_QRDecomposeMatrixd ( nrows, ncols, a, aa ) ) {
    pkv_SetScratchMemTop ( sp );
    return false;
  }
  pkn_multiReflectVectord ( nrows, ncols, a, aa, spdimen, bpitch, b );
  pkn_multiMultInvUTVectord ( ncols, a, spdimen, bpitch, b, xpitch, x );

  pkv_SetScratchMemTop ( sp );
  return true;
} /*pkn_multiSolveRLSQd*/

void pkn_QRGetReflectiond ( int nrows, int ncols,
                            const double *a, const double *aa,
                            int nrefl, double *w, double *gamma )
{
  if ( gamma )
    *gamma = aa[ncols+nrefl];
  w[0] = aa[nrefl];
  pkv_Selectd ( nrows-nrefl-1, 1, ncols, 1,
                &a[(nrefl+1)*ncols+nrefl], &w[1] );
} /*pkn_QRGetReflectiond*/

