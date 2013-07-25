
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

boolean pkn_QRDecompUHessenbergd ( int n, double *ah )
{
  int    i, j;
  double c, s, xi;

  for ( i = 0; i < n-1; i++ ) {
    pkn_FindGivensRotXid ( ah[pkn_UHessenbergMatIndex(i,i)],
                           ah[pkn_UHessenbergMatIndex(i+1,i)], &c, &s, &xi );
    for ( j = i; j < n; j++ )
      pkn_ApplyGivensRotationd ( c, s, &ah[pkn_UHessenbergMatIndex(i,j)],
                                       &ah[pkn_UHessenbergMatIndex(i+1,j)] );
    if ( !ah[pkn_UHessenbergMatIndex(i,i)] )
      return false;
    ah[pkn_UHessenbergMatIndex(i+1,i)] = xi;
  }
  return true;
} /*pkn_QRDecompUHessenbergd*/

void pkn_multiSolveQRUHessenbergd ( int n, double *qrh,
                                    int spdimen, int pitch, double *b )
{
  int    i, j, k;
  double c, s, t, *b1, *b2;

        /* solve the system with the orthogonal matrix */
  b1 = b;
  b2 = b1 + pitch;
  for ( i = 0; i < n-1; i++ ) {
    pkn_FindXiGivensRotd ( qrh[pkn_UHessenbergMatIndex(i+1,i)], &c, &s );
    for ( k = 0; k < spdimen; k++ )
      pkn_ApplyGivensRotationd ( c, s, &b1[k], &b2[k] );
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
} /*pkn_multiSolveQRUHessenbergd*/

boolean pkn_multiSolveUHessenbergLinEqd ( int n, double *ah,
                                          int spdimen, int pitch, double *b )
{
  void   *sp;
  double *qrh;

  sp = pkv_GetScratchMemTop ();
  qrh = pkv_GetScratchMemd ( ((n+1)*(n+2))/2-2 );
  if ( !qrh )
    goto failure;
  memcpy ( qrh, ah, (((n+1)*(n+2))/2-2)*sizeof(double) );
  if ( !pkn_QRDecompUHessenbergd ( n, qrh ) )
    goto failure;
  pkn_multiSolveQRUHessenbergd ( n, qrh, spdimen, pitch, b );
  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*pkn_multiSolveUHessenbergLinEqd*/

