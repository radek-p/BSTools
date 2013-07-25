
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2005                                  */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "pkvaria.h"
#include "pknum.h"
#include "msgpool.h"

void pkn_ComputeQTSQf ( int m, const float *s,
                        int n, const float *a, const float *aa,
                        float *b )
{
  void   *sp;
  float  *v, *p, gamma;
  double aux;
  int    i, j, k;

  sp = pkv_GetScratchMemTop ();
  v = pkv_GetScratchMemf ( m );
  p = pkv_GetScratchMemf ( m );
  if ( !v || !p )
    pkv_SignalError ( LIB_PKNUM, 15, ERRMSG_0 );

  if ( s != b )
    memcpy ( b, s, (m*(m+1)/2)*sizeof(float) );

        /* apply subsequent Householder reflections from both sides */
  for ( i = 0; i < n; i++ ) {
    if ( i )
      memset ( v, 0, i*sizeof(float) );
    pkn_QRGetReflectionf ( m, n, a, aa, i, &v[i], &gamma );
          /* compute w = S v gamma */
    for ( j = 0; j < m; j++ ) {
      aux = 0.0;
      for ( k = i; k < m; k++ )
        aux += b[pkn_SymMatIndex(j,k)]*v[k];
      p[j] = (float)(aux*gamma);
    }
          /* compute p = w - v v^T w gamma / 2 */
    aux = pkn_ScalarProductf ( m-i, &v[i], &p[i] )*gamma*0.5;
    for ( j = i; j < m; j++ )
      p[j] -= (float)(v[j]*aux);
          /* compute b = s - ( v p^T + p v^T ) */
    for ( j = 0; j < m; j++ )
      for ( k = 0; k <= j; k++ )
        b[pkn_SymMatIndex(j,k)] -= v[j]*p[k]+p[j]*v[k];
  }

  pkv_SetScratchMemTop ( sp );
} /*pkn_ComputeQTSQf*/

void pkn_ComputeQSQTf ( int m, const float *s,
                        int n, const float *a, const float *aa,
                        float *b )
{
  void   *sp;
  float  *v, *p, gamma;
  double aux;
  int    i, j, k;

  sp = pkv_GetScratchMemTop ();
  v = pkv_GetScratchMemf ( m );
  p = pkv_GetScratchMemf ( m );
  if ( !v || !p )
    pkv_SignalError ( LIB_PKNUM, 15, ERRMSG_0 );

  if ( s != b )
    memcpy ( b, s, (m*(m+1)/2)*sizeof(float) );

        /* apply subsequent Householder reflections from both sides */
  for ( i = n-1; i >= 0; i-- ) {
    if ( i )
      memset ( v, 0, i*sizeof(float) );
    pkn_QRGetReflectionf ( m, n, a, aa, i, &v[i], &gamma );
          /* compute w = S v gamma */
    for ( j = 0; j < m; j++ ) {
      aux = 0.0;
      for ( k = i; k < m; k++ )
        aux += b[pkn_SymMatIndex(j,k)]*v[k];
      p[j] = (float)(aux*gamma);
    }
          /* compute p = w - v v^T w gamma / 2 */
    aux = pkn_ScalarProductf ( m-i, &v[i], &p[i] )*gamma*0.5;
    for ( j = i; j < m; j++ )
      p[j] -= (float)(v[j]*aux);
          /* compute b = s - ( v p^T + p v^T ) */
    for ( j = 0; j < m; j++ )
      for ( k = 0; k <= j; k++ )
        b[pkn_SymMatIndex(j,k)] -= v[j]*p[k]+p[j]*v[k];
  }

  pkv_SetScratchMemTop ( sp );
} /*pkn_ComputeQSQTf*/

