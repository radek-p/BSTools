
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

boolean pkn_QRDecompUHessenbergf ( int n, float *ah )
{
  int   i, j;
  float c, s, xi;

  for ( i = 0; i < n-1; i++ ) {
    pkn_FindGivensRotXif ( ah[pkn_UHessenbergMatIndex(i,i)],
                           ah[pkn_UHessenbergMatIndex(i+1,i)], &c, &s, &xi );
    for ( j = i; j < n; j++ )
      pkn_ApplyGivensRotationf ( c, s, &ah[pkn_UHessenbergMatIndex(i,j)],
                                       &ah[pkn_UHessenbergMatIndex(i+1,j)] );
    if ( !ah[pkn_UHessenbergMatIndex(i,i)] )
      return false;
    ah[pkn_UHessenbergMatIndex(i+1,i)] = xi;
  }
  return true;
} /*pkn_QRDecompUHessenbergf*/

void pkn_multiSolveQRUHessenbergf ( int n, float *qrh,
                                    int spdimen, int pitch, float *b )
{
  int   i, j, k;
  float c, s, t, *b1, *b2;

        /* solve the system with the orthogonal matrix */
  b1 = b;
  b2 = b1 + pitch;
  for ( i = 0; i < n-1; i++ ) {
    pkn_FindXiGivensRotf ( qrh[pkn_UHessenbergMatIndex(i+1,i)], &c, &s );
    for ( k = 0; k < spdimen; k++ )
      pkn_ApplyGivensRotationf ( c, s, &b1[k], &b2[k] );
    b1 = b2;
    b2 = b1 + pitch;
  }
        /* solve the system with the upper triangular matrix */
  for ( i = n-1; i >= 0; i-- ) {
    t = qrh[pkn_UHessenbergMatIndex(i,i)];
    for ( k = 0; k < spdimen; k++ )
      b1[k] /= t;
    for ( j = i-1, b2 = b1-pitch;  j >= 0;  j--, b2 -= pitch ) {
      t = qrh[pkn_UHessenbergMatIndex(j,i)];
      for ( k = 0; k < spdimen; k++ )
        b2[k] -= t*b1[k];
    }    
    b1 -= pitch;
  }
} /*pkn_multiSolveQRUHessenbergf*/

boolean pkn_multiSolveUHessenbergLinEqf ( int n, float *ah,
                                          int spdimen, int pitch, float *b )
{
  void  *sp;
  float *qrh;

  sp = pkv_GetScratchMemTop ();
  qrh = pkv_GetScratchMemf ( ((n+1)*(n+2))/2-2 );
  if ( !qrh )
    goto failure;
  memcpy ( qrh, ah, (((n+1)*(n+2))/2-2)*sizeof(float) );
  if ( !pkn_QRDecompUHessenbergf ( n, qrh ) )
    goto failure;
  pkn_multiSolveQRUHessenbergf ( n, qrh, spdimen, pitch, b );
  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*pkn_multiSolveUHessenbergLinEqf*/

