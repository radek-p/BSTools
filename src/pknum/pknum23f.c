
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2005, 2013                            */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

#include <math.h>
#include <string.h>

#include "pkvaria.h"
#include "pknum.h"

/* ///////////////////////////////////////////////////////////////////////// */
/* Gaussian elimination with full pivoting and solving systems of linear     */
/* equations. */

boolean pkn_GaussDecomposePLUQf ( int n, float *a, int *P, int *Q )
{
  int   i, j, k, i0, j0;
  float s, t;

  for ( k = 0; k < n-1; k++ ) {

        /* find the pivoting element */
    s = (float)fabs ( a[k*n+k] );  i0 = j0 = k;
    for ( i = k; i < n; i++ )
      for ( j = k; j < n; j++ ) {
        t = (float)fabs ( a[i*n+j] );
        if ( t > s )
          { s = t;  i0 = i;  j0 = j; }
      }
        /* exchange rows, if necessary */
    if ( i0 != k ) {
      pkv_Exchange ( &a[k*n], &a[i0*n], n*sizeof(float) );
      P[k] = i0;
    }
    else P[k] = k;
        /* exchange columns, if necessary */
    if ( j0 != k ) {
      for ( i = 0; i < n; i++ )
        { s = a[i*n+k];  a[i*n+k] = a[i*n+j0];  a[i*n+j0] = s; }
      Q[k] = j0;
    }
    else Q[k] = k;

        /* now the proper elimination */
    s = a[k*n+k];
    if ( !s )
      return false;
    for ( i = k+1; i < n; i++ ) {
      t = a[i*n+k] /= s;
      for ( j = k+1; j < n; j++ )
        a[i*n+j] -= t*a[k*n+j];
    }
  }
  return a[n*n-1] != 0.0;
} /*pkn_GaussDecomposePLUQf*/

void pkn_multiSolvePLUQf ( int n, const float *lu, const int *P, const int *Q,
                           int spdimen, int pitch, float *b )
{
  int   i, j, k;
  float t;

      /* permute the rows of b */
  for ( i = 0; i < n-1; i++ )
    if ( P[i] != i )
      pkv_Exchange ( &b[i*pitch], &b[P[i]*pitch], spdimen*sizeof(float) );

      /* solve the system with the lower triangular matrix */
  for ( i = 0; i < n-1; i++ )
    for ( j = i+1; j < n; j++ ) {
      t = lu[j*n+i];
      for ( k = 0; k < spdimen; k++ )
        b[j*pitch+k] -= t*b[i*pitch+k];
    }

      /* solve the system with the upper triangular matrix */
  for ( i = n-1; i >= 0; i-- ) {
    t = lu[i*n+i];
    for ( k = 0; k < spdimen; k++ )
      b[i*pitch+k] /= t;
    for ( j = i-1; j >= 0; j-- ) {
      t = lu[j*n+i];
      for ( k = 0; k < spdimen; k++ )
        b[j*pitch+k] -= t*b[i*pitch+k];
    }
  }

      /* permute the rows of the solution */
  for ( i = n-2; i >= 0; i-- )
    if ( Q[i] != i )
      pkv_Exchange ( &b[i*pitch], &b[Q[i]*pitch], spdimen*sizeof(float) );
} /*pkn_multiSolvePLUQf*/

boolean pkn_multiGaussSolveLinEqf ( int n, const float *a,
                                    int spdimen, int pitch, float *b )
{
  void  *sp;
  int   *P, *Q;
  float *lu;

  sp = pkv_GetScratchMemTop ();

  lu = pkv_GetScratchMemf ( n*n );
  P = pkv_GetScratchMem ( (n-1)*sizeof(int) );
  Q = pkv_GetScratchMem ( (n-1)*sizeof(int) );
  if ( !lu || !P || !Q )
    goto failure;

  memcpy ( lu, a, n*n*sizeof(float) );
  if ( !pkn_GaussDecomposePLUQf ( n, lu, P, Q ) )
    goto failure;
  pkn_multiSolvePLUQf ( n, lu, P, Q, spdimen, pitch, b );

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*pkn_multiGaussSolveLinEqf*/

boolean pkn_GaussInvertMatrixf ( int n, float *a )
{
  void  *sp;
  float *ai;
  int   i;

  sp = pkv_GetScratchMemTop ();
  if ( !(ai = pkv_GetScratchMemf ( n*n )) )
    goto failure;
        /* load the n x n identity matrix */
  memset ( ai, 0, n*n*sizeof(float) );
  for ( i = 0; i < n; i++ )
    ai[i*(n+1)] = 1.0;
  if ( !pkn_multiGaussSolveLinEqf ( n, a, n, n, ai ) )
    goto failure;
  memcpy ( a, ai, n*n*sizeof(float) );
  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*pkn_GaussInvertMatrixf*/

