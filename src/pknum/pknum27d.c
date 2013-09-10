
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2005, 2013                            */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "pkvaria.h"
#include "pknum.h"

boolean pkn_ComputeQTSQd ( int m, const double *s,
                           int n, const double *a, const double *aa,
                           double *b )
{
  void   *sp;
  double *v, *p, gamma;
  double aux;
  int    i, j, k;

  sp = pkv_GetScratchMemTop ();
  v = pkv_GetScratchMemd ( m );
  p = pkv_GetScratchMemd ( m );
  if ( !v || !p ) {
    PKV_SIGNALERROR ( LIB_PKNUM, 2, ERRMSG_2 );
    goto failure;
  }

  if ( s != b )
    memcpy ( b, s, (m*(m+1)/2)*sizeof(double) );

        /* apply subsequent Householder reflections from both sides */
  for ( i = 0; i < n; i++ ) {
    if ( i )
      memset ( v, 0, i*sizeof(double) );
    pkn_QRGetReflectiond ( m, n, a, aa, i, &v[i], &gamma );
          /* compute w = S v gamma */
    for ( j = 0; j < m; j++ ) {
      aux = 0.0;
      for ( k = i; k < m; k++ )
        aux += b[pkn_SymMatIndex(j,k)]*v[k];
      p[j] = aux*gamma;
    }
          /* compute p = w - v v^T w gamma / 2 */
    aux = pkn_ScalarProductd ( m-i, &v[i], &p[i] )*gamma*0.5;
    for ( j = i; j < m; j++ )
      p[j] -= v[j]*aux;
          /* compute b = s - ( v p^T + p v^T ) */
    for ( j = 0; j < m; j++ )
      for ( k = 0; k <= j; k++ )
        b[pkn_SymMatIndex(j,k)] -= v[j]*p[k]+p[j]*v[k];
  }

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*pkn_ComputeQTSQd*/

boolean pkn_ComputeQSQTd ( int m, const double *s,
                           int n, const double *a, const double *aa,
                           double *b )
{
  void   *sp;
  double  *v, *p, gamma;
  double aux;
  int    i, j, k;

  sp = pkv_GetScratchMemTop ();
  v = pkv_GetScratchMemd ( m );
  p = pkv_GetScratchMemd ( m );
  if ( !v || !p ) {
    PKV_SIGNALERROR ( LIB_PKNUM, 2, ERRMSG_2 );
    goto failure;
  }

  if ( s != b )
    memcpy ( b, s, (m*(m+1)/2)*sizeof(double) );

        /* apply subsequent Householder reflections from both sides */
  for ( i = n-1; i >= 0; i-- ) {
    if ( i )
      memset ( v, 0, i*sizeof(double) );
    pkn_QRGetReflectiond ( m, n, a, aa, i, &v[i], &gamma );
          /* compute w = S v gamma */
    for ( j = 0; j < m; j++ ) {
      aux = 0.0;
      for ( k = i; k < m; k++ )
        aux += b[pkn_SymMatIndex(j,k)]*v[k];
      p[j] = aux*gamma;
    }
          /* compute p = w - v v^T w gamma / 2 */
    aux = pkn_ScalarProductd ( m-i, &v[i], &p[i] )*gamma*0.5;
    for ( j = i; j < m; j++ )
      p[j] -= v[j]*aux;
          /* compute b = s - ( v p^T + p v^T ) */
    for ( j = 0; j < m; j++ )
      for ( k = 0; k <= j; k++ )
        b[pkn_SymMatIndex(j,k)] -= v[j]*p[k]+p[j]*v[k];
  }

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*pkn_ComputeQSQTd*/

