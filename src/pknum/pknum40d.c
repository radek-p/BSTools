
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2009, 2010                            */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "pkvaria.h"
#include "pknum.h"

/* ///////////////////////////////////////////////////////////////////////// */
boolean pkn_NRBComputeQTSQbld ( int n, int *prof, double *Amat, double **Arows,
                                int w, double *Bmat, double *bb,
                                int *qa11prof, double **QA11rows,
                                int *qa22prof, double **QA22rows,
                                double **QA21 )
{
  void   *sp;
  int    i, j, k, lv1, lv2, lp1, lp2;
  int    *prof11a, *prof11b, *prof22a, *prof22b;
  int    size11a, size11b, size22a, size22b, size21;
  double *v, *p, gamma, s;
  double **rows11a, **rows11b, **rows22a, **rows22b, *m21a, *m21b;
  void   *aa, *ab;

  if ( w <= 0 || w >= n )
    return false;
  sp = pkv_GetScratchMemTop ();
  aa = ab = NULL;

        /* allocate the memory */
  v = pkv_GetScratchMemd ( n );
  p = pkv_GetScratchMemd ( n );
  prof11a = pkv_GetScratchMemi ( n );
  prof11b = pkv_GetScratchMemi ( n );
  rows11a = pkv_GetScratchMem ( n*sizeof(double*) );
  rows11b = pkv_GetScratchMem ( n*sizeof(double*) );
  if ( !v || !p || !prof11a || !prof11b || !rows11a || !rows11b )
    goto failure;
  prof22a = &prof11a[w];
  prof22b = &prof11b[w];
  rows22a = &rows11a[w];
  rows22b = &rows11b[w];
          /* copy the profile of A11 */
  memcpy ( prof11a, prof, w*sizeof(int) );
          /* below it is the profile of A22 determined */
  for ( j = w; j < n; j++ ) {
    for ( k = max ( w, prof[j]); k < j && !Arows[j][k]; k++ )
      ;
    prof11a[j] = k-w;
  }
  size11a = pkn_NRBArraySize ( w, prof11a );
  size22a = pkn_NRBArraySize ( n-w, prof22a );
  size21 = w*(n-w);
  PKV_MALLOC ( aa, (size11a+size22a+size21)*sizeof(double) );
  if ( !aa )
    goto failure;
  rows11a[0] = aa;
  rows22a[0] = &rows11a[0][size11a];
  m21a = &rows22a[0][size22a];
  pkn_NRBFindRowsd ( w, prof11a, rows11a[0], rows11a );
  pkn_NRBFindRowsd ( n-w, prof22a, rows22a[0], rows22a );

        /* copy the matrix A, splitting it to blocks */
  memcpy ( rows11a[0], Amat, size11a*sizeof(double) );
  memset ( m21a, 0, size21*sizeof(double) );
  for ( j = w; j < n; j++ ) {
    k = prof[j];
    if ( k < w ) {
      memcpy ( &m21a[(j-w)*w+k], &Arows[j][k], (w-k)*sizeof(double) );
      k = prof11a[j];
      memcpy ( &rows22a[j-w][k], &Arows[j][w+k], (j-w-k+1)*sizeof(double) );
    }
    else
      memcpy ( &rows22a[j-w][k-w], &Arows[j][k], (j-k+1)*sizeof(double) );
  }

  for ( i = 0; i < w; i++ ) {
        /* use the Ortega-Householder algorithm to compute Q^TAQ */
          /* get the Householder reflection vector */
    pkn_QRGetReflectiond ( n, w, Bmat, bb, i, &v[i], &gamma );
          /* multiply A_{k-1} by v and by gamma */
    pkn_NRBSymMultd ( w, prof11a, rows11a[0], rows11a, 1, 1, v, 1, p );
    pkn_MultTMatrixAddd ( n-w, w, w, m21a, 1, 1, &v[w], 1, p );
    pkn_NRBSymMultd ( n-w, prof22a, rows22a[0], rows22a, 1, 1, &v[w], 1, &p[w] );
    pkn_MultMatrixAddd ( n-w, w, w, m21a, 1, 1, v, 1, &p[w] );
    pkn_MultMatrixNumd ( 1, n, 0, p, gamma, 0, p );
    s = pkn_ScalarProductd ( n-i, &v[i], &p[i] )*0.5*gamma;
    pkn_AddMatrixMd ( 1, n-i, 0, &p[i], 0, &v[i], -s, 0, &p[i] );
        /* now compute A_k = A_{k-1} - p*v^T - v*p^T */
          /* first it is necessary to find the profile of the result */
    memcpy ( prof11b, prof11a, n*sizeof(int) );
    for ( lv1 = i; lv1 < w; lv1++ )    
      if ( v[lv1] )
        break;
    for ( lv2 = w; lv2 < n; lv2++ )
      if ( v[lv2] )
        break;
    lv2 -= w;
    for ( lp1 = 0; lp1 < w; lp1++ )
      if ( p[lp1] )
        break;
    for ( lp2 = w; lp2 < n; lp2++ )
      if ( p[lp2] )
        break;
    lp2 -= w;
    for ( j = i; j < w; j++ )
      if ( v[j] && lp1 < prof11b[j] ) prof11b[j] = lp1;
    for ( j = 0; j < w; j++ )
      if ( p[j] && lv1 < prof11b[j] ) prof11b[j] = lv1;
    for ( j = 0; j < n-w; j++ ) {
      if ( v[j+w] && lp2 < prof22b[j] ) prof22b[j] = lp2;
      if ( p[j+w] && lv2 < prof22b[j] ) prof22b[j] = lv2;
    }
          /* allocate a suitable memory block */
    size11b = pkn_NRBArraySize ( w, prof11b );
    size22b = pkn_NRBArraySize ( n-w, prof22b );
    PKV_MALLOC ( ab, (size11b+size22b+size21)*sizeof(double) );
    if ( !ab )
      goto failure;
    rows11b[0] = ab;
    rows22b[0] = &rows11b[0][size11b];
    m21b = &rows22b[0][size22b];
    pkn_NRBFindRowsd ( w, prof11b, rows11b[0], rows11b );
    pkn_NRBFindRowsd ( n-w, prof22b, rows22b[0], rows22b );
          /* copy the matrix A_{k-1} */
    memset ( rows11b[0], 0, (size11b+size22b)*sizeof(double) );
    for ( j = 0; j < w; j++ )
      memcpy ( &rows11b[j][prof11a[j]], &rows11a[j][prof11a[j]],
               (j-prof11a[j]+1)*sizeof(double) );
    for ( j = 0; j < n-w; j++ )
      memcpy ( &rows22b[j][prof22a[j]], &rows22a[j][prof22a[j]],
               (j-prof22a[j]+1)*sizeof(double) );
    memcpy ( m21b, m21a, size21*sizeof(double) );
          /* add the two rank-1 matrices */
    for ( j = 0; j < w; j++ ) {
      if ( v[j] )
        for ( k = lp1; k <= j; k++ )
          rows11b[j][k] -= v[j]*p[k];
      if ( p[j] )
        for ( k = lv1; k <= j; k++ )
          rows11b[j][k] -= p[j]*v[k];
    }
    for ( j = 0; j < n-w; j++ ) {
      if ( v[j+w] ) {
        for ( k = lp1; k < w; k++ )
          m21b[j*w+k] -= v[j+w]*p[k];
        for ( k = lp2; k <= j; k++ )
          rows22b[j][k] -= v[j+w]*p[k+w];
      }
      if ( p[j+w] ) {
        for ( k = lv1; k < w; k++ )
          m21b[j*w+k] -= p[j+w]*v[k];
        for ( k = lv2; k <= j; k++ )
          rows22b[j][k] -= p[j+w]*v[k+w];
      }
    }
    v[i] = 0.0;
        /* now deallocate A_{k-1} and copy A_k */
    PKV_FREE ( rows11a[0] );
    memcpy ( prof11a, prof11b, n*sizeof(int) );
    memcpy ( rows11a, rows11b, n*sizeof(double*) );
    m21a = m21b;
    aa = rows11a[0];
    ab = NULL;
  }
        /* assign the result to output parameters */
  memcpy ( qa11prof, prof11a, w*sizeof(int) );
  memcpy ( QA11rows, rows11a, w*sizeof(double*) );
  memcpy ( qa22prof, prof22a, (n-w)*sizeof(int) );
  memcpy ( QA22rows, rows22a, (n-w)*sizeof(double*) );
  *QA21 = m21a;
  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  if ( aa ) PKV_FREE ( aa );
  if ( ab ) PKV_FREE ( ab );
  return false;
} /*pkn_NRBComputeQTSQbld*/

/* ///////////////////////////////////////////////////////////////////////// */
boolean pkn_NRBComputeQSQTbld ( int n, int *prof, double *Amat, double **Arows,
                                int w, double *Bmat, double *bb,
                                int *qa11prof, double **QA11rows,
                                int *qa22prof, double **QA22rows,
                                double **QA21 )
{
  void   *sp;
  int    i, j, k, lv1, lv2, lp1, lp2;
  int    *prof11a, *prof11b, *prof22a, *prof22b;
  int    size11a, size11b, size22a, size22b, size21;
  double *v, *p, gamma, s;
  double **rows11a, **rows11b, **rows22a, **rows22b, *m21a, *m21b;
  void   *aa, *ab;

  if ( w <= 0 || w >= n )
    return false;
  sp = pkv_GetScratchMemTop ();
  aa = ab = NULL;

        /* allocate the memory */
  v = pkv_GetScratchMemd ( n );
  p = pkv_GetScratchMemd ( n );
  prof11a = pkv_GetScratchMemi ( n );
  prof11b = pkv_GetScratchMemi ( n );
  rows11a = pkv_GetScratchMem ( n*sizeof(double*) );
  rows11b = pkv_GetScratchMem ( n*sizeof(double*) );
  if ( !v || !p || !prof11a || !prof11b || !rows11a || !rows11b )
    goto failure;
  prof22a = &prof11a[w];
  prof22b = &prof11b[w];
  rows22a = &rows11a[w];
  rows22b = &rows11b[w];
          /* copy the profile of A11 */
  memcpy ( prof11a, prof, w*sizeof(int) );
          /* below it is the profile of A22 determined */
  for ( j = w; j < n; j++ ) {
    for ( k = max ( w, prof[j]); k < j && !Arows[j][k]; k++ )
      ;
    prof11a[j] = k-w;
  }
  size11a = pkn_NRBArraySize ( w, prof11a );
  size22a = pkn_NRBArraySize ( n-w, prof22a );
  size21 = w*(n-w);
  PKV_MALLOC ( aa, (size11a+size22a+size21)*sizeof(double) );
  if ( !aa )
    goto failure;
  rows11a[0] = aa;
  rows22a[0] = &rows11a[0][size11a];
  m21a = &rows22a[0][size22a];
  pkn_NRBFindRowsd ( w, prof11a, rows11a[0], rows11a );
  pkn_NRBFindRowsd ( n-w, prof22a, rows22a[0], rows22a );

        /* copy the matrix A, splitting it to blocks */
  memcpy ( rows11a[0], Amat, size11a*sizeof(double) );
  memset ( m21a, 0, size21*sizeof(double) );
  for ( j = w; j < n; j++ ) {
    k = prof[j];
    if ( k < w ) {
      memcpy ( &m21a[(j-w)*w+k], &Arows[j][k], (w-k)*sizeof(double) );
      k = prof11a[j];
      memcpy ( &rows22a[j-w][k], &Arows[j][w+k], (j-w-k+1)*sizeof(double) );
    }
    else
      memcpy ( &rows22a[j-w][k-w], &Arows[j][k], (j-k+1)*sizeof(double) );
  }

  for ( i = w-1; i >= 0; i-- ) {
        /* use the Ortega-Householder algorithm to compute Q^TAQ */
          /* get the Householder reflection vector */
    pkn_QRGetReflectiond ( n, w, Bmat, bb, i, &v[i], &gamma );
          /* multiply A_{k-1} by v and by gamma */
    pkn_NRBSymMultd ( w, prof11a, rows11a[0], rows11a, 1, 1, v, 1, p );
    pkn_MultTMatrixAddd ( n-w, w, w, m21a, 1, 1, &v[w], 1, p );
    pkn_NRBSymMultd ( n-w, prof22a, rows22a[0], rows22a, 1, 1, &v[w], 1, &p[w] );
    pkn_MultMatrixAddd ( n-w, w, w, m21a, 1, 1, v, 1, &p[w] );
    pkn_MultMatrixNumd ( 1, n, 0, p, gamma, 0, p );
    s = pkn_ScalarProductd ( n-i, &v[i], &p[i] )*0.5*gamma;
    pkn_AddMatrixMd ( 1, n-i, 0, &p[i], 0, &v[i], -s, 0, &p[i] );
        /* now compute A_k = A_{k-1} - p*v^T - v*p^T */
          /* first it is necessary to find the profile of the result */
    memcpy ( prof11b, prof11a, n*sizeof(int) );
    for ( lv1 = i; lv1 < w; lv1++ )    
      if ( v[lv1] )
        break;
    for ( lv2 = w; lv2 < n; lv2++ )
      if ( v[lv2] )
        break;
    lv2 -= w;
    for ( lp1 = 0; lp1 < w; lp1++ )
      if ( p[lp1] )
        break;
    for ( lp2 = w; lp2 < n; lp2++ )
      if ( p[lp2] )
        break;
    lp2 -= w;
    for ( j = i; j < w; j++ )
      if ( v[j] && lp1 < prof11b[j] ) prof11b[j] = lp1;
    for ( j = 0; j < w; j++ )
      if ( p[j] && lv1 < prof11b[j] ) prof11b[j] = lv1;
    for ( j = 0; j < n-w; j++ ) {
      if ( v[j+w] && lp2 < prof22b[j] ) prof22b[j] = lp2;
      if ( p[j+w] && lv2 < prof22b[j] ) prof22b[j] = lv2;
    }
          /* allocate a suitable memory block */
    size11b = pkn_NRBArraySize ( w, prof11b );
    size22b = pkn_NRBArraySize ( n-w, prof22b );
    PKV_MALLOC ( ab, (size11b+size22b+size21)*sizeof(double) );
    if ( !ab )
      goto failure;
    rows11b[0] = ab;
    rows22b[0] = &rows11b[0][size11b];
    m21b = &rows22b[0][size22b];
    pkn_NRBFindRowsd ( w, prof11b, rows11b[0], rows11b );
    pkn_NRBFindRowsd ( n-w, prof22b, rows22b[0], rows22b );
          /* copy the matrix A_{k-1} */
    memset ( rows11b[0], 0, (size11b+size22b)*sizeof(double) );
    for ( j = 0; j < w; j++ )
      memcpy ( &rows11b[j][prof11a[j]], &rows11a[j][prof11a[j]],
               (j-prof11a[j]+1)*sizeof(double) );
    for ( j = 0; j < n-w; j++ )
      memcpy ( &rows22b[j][prof22a[j]], &rows22a[j][prof22a[j]],
               (j-prof22a[j]+1)*sizeof(double) );
    memcpy ( m21b, m21a, size21*sizeof(double) );
          /* add the two rank-1 matrices */
    for ( j = 0; j < w; j++ ) {
      if ( v[j] )
        for ( k = lp1; k <= j; k++ )
          rows11b[j][k] -= v[j]*p[k];
      if ( p[j] )
        for ( k = lv1; k <= j; k++ )
          rows11b[j][k] -= p[j]*v[k];
    }
    for ( j = 0; j < n-w; j++ ) {
      if ( v[j+w] ) {
        for ( k = lp1; k < w; k++ )
          m21b[j*w+k] -= v[j+w]*p[k];
        for ( k = lp2; k <= j; k++ )
          rows22b[j][k] -= v[j+w]*p[k+w];
      }
      if ( p[j+w] ) {
        for ( k = lv1; k < w; k++ )
          m21b[j*w+k] -= p[j+w]*v[k];
        for ( k = lv2; k <= j; k++ )
          rows22b[j][k] -= p[j+w]*v[k+w];
      }
    }
    v[i] = 0.0;
        /* now deallocate A_{k-1} and copy A_k */
    PKV_FREE ( rows11a[0] );
    memcpy ( prof11a, prof11b, n*sizeof(int) );
    memcpy ( rows11a, rows11b, n*sizeof(double*) );
    m21a = m21b;
    aa = rows11a[0];
    ab = NULL;
  }
        /* assign the result to output parameters */
  memcpy ( qa11prof, prof11a, w*sizeof(int) );
  memcpy ( QA11rows, rows11a, w*sizeof(double*) );
  memcpy ( qa22prof, prof22a, (n-w)*sizeof(int) );
  memcpy ( QA22rows, rows22a, (n-w)*sizeof(double*) );
  *QA21 = m21a;
  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  if ( aa ) PKV_FREE ( aa );
  if ( ab ) PKV_FREE ( ab );
  return false;
} /*pkn_NRBComputeQSQTbld*/

