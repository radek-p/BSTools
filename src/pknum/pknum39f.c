
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
boolean pkn_NRBComputeQTSQf ( int n, int *prof, float *Amat, float **Arows,
                              int w, float *Bmat, float *bb,
                              int *qaprof, float **QArows )
{
  void  *sp;
  int   i, j, k, lv, lp;
  int   *prof1, *prof2, size;
  float *v, *p, gamma, s;
  float **rows1, **rows2;

  sp = pkv_GetScratchMemTop ();
  rows1 = NULL;
  i = 0;
  if ( w <= 0 || w >= n )
    goto failure;

        /* allocate the memory */
  v = pkv_GetScratchMemf ( n );
  p = pkv_GetScratchMemf ( n );
  prof1 = pkv_GetScratchMemi ( n );
  prof2 = pkv_GetScratchMemi ( n );
  rows1 = pkv_GetScratchMem ( n*sizeof(float*) );
  rows2 = pkv_GetScratchMem ( n*sizeof(float*) );
  if ( !v || !p || !prof1 || !prof2 || !rows1 || !rows2 )
    goto failure;

  memcpy ( prof1, prof, n*sizeof(int) );
  memcpy ( rows1, Arows, n*sizeof(float*) );
  for ( i = 0; i < w; i++ ) {
        /* use the Ortega-Householder algorithm to compute Q^TAQ */
          /* get the Householder reflection vector */
    pkn_QRGetReflectionf ( n, w, Bmat, bb, i, &v[i], &gamma );
          /* multiply A_{k-1} by v and by gamma */
    pkn_NRBSymMultf ( n, prof1, rows1[0], rows1, 1, 1, v, 1, p );
    pkn_MultMatrixNumf ( 1, n, 0, p, gamma, 0, p );
    s = pkn_ScalarProductf ( n-i, &v[i], &p[i] )*0.5*gamma;
    pkn_AddMatrixMf ( 1, n-i, 0, &p[i], 0, &v[i], -s, 0, &p[i] );
        /* now compute A_k = A_{k-1} - p*v^T - v*p^T */
          /* first it is necessary to find the profile of the result */
    memcpy ( prof2, prof1, n*sizeof(int) );
    for ( lv = i; lv < n; lv++ )    
      if ( v[lv] )
        break;
    for ( lp = 0; lp < n; lp++ )
      if ( p[lp] )
        break;
    for ( j = i; j < n; j++ )
      if ( v[j] && lp < prof2[j] ) prof2[j] = lp;
    for ( j = 0; j < n; j++ )
      if ( p[j] && lv < prof2[j] ) prof2[j] = lv;
          /* allocate a suitable memory block */
    size = pkn_NRBArraySize ( n, prof2 );
    PKV_MALLOC ( rows2[0], size*sizeof(float) );
    if ( !rows2[0] )
      goto failure;
    if ( !pkn_NRBFindRowsf ( n, prof2, rows2[0], rows2 ) )
      goto failure;
          /* copy the matrix A_{k-1} */
    memset ( rows2[0], 0, size*sizeof(float) );
    for ( j = 0; j < n; j++ )
      memcpy ( &rows2[j][prof1[j]], &rows1[j][prof1[j]],
               (j-prof1[j]+1)*sizeof(float) );
          /* add the two rank-1 matrices */
    for ( j = 0; j < n; j++ ) {
      if ( v[j] )
        for ( k = lp; k <= j; k++ )
          rows2[j][k] -= v[j]*p[k];
      if ( p[j] )
        for ( k = lv; k <= j; k++ )
          rows2[j][k] -= p[j]*v[k];
    }
    v[i] = 0.0;
        /* now deallocate A_{k-1} and copy A_k */
    if ( i > 0 )
      PKV_FREE ( rows1[0] );
    memcpy ( prof1, prof2, n*sizeof(int) );
    memcpy ( rows1, rows2, n*sizeof(float*) );
  }
        /* assign the result to output parameters */
  memcpy ( qaprof, prof1, n*sizeof(int) );
  memcpy ( QArows, rows1, n*sizeof(float*) );

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  if ( i > 0 ) PKV_FREE ( rows1[0] );
  return false;
} /*pkn_NRBComputeQTSQf*/

boolean pkn_NRBComputeQSQTf ( int n, int *prof, float *Amat, float **Arows,
                              int w, float *Bmat, float *bb,
                              int *qaprof, float **QArows )
{
  void  *sp;
  int   i, j, k, lv, lp;
  int   *prof1, *prof2, size;
  float *v, *p, gamma, s;
  float **rows1, **rows2;

  sp = pkv_GetScratchMemTop ();
  rows1 = NULL;
  i = w-1;
  if ( w <= 0 || w >= n )
    goto failure;

        /* allocate the memory */
  v = pkv_GetScratchMemf ( n );
  p = pkv_GetScratchMemf ( n );
  prof1 = pkv_GetScratchMemi ( n );
  prof2 = pkv_GetScratchMemi ( n );
  rows1 = pkv_GetScratchMem ( n*sizeof(float*) );
  rows2 = pkv_GetScratchMem ( n*sizeof(float*) );
  if ( !v || !p || !prof1 || !prof2 || !rows1 || !rows2 )
    goto failure;

  memset ( v, 0, (w-1)*sizeof(float) );
  memcpy ( prof1, prof, n*sizeof(int) );
  memcpy ( rows1, Arows, n*sizeof(float*) );
  for ( i = w-1; i >= 0; i-- ) {
        /* use the Ortega-Householder algorithm to compute Q^TAQ */
          /* get the Householder reflection vector */
    pkn_QRGetReflectionf ( n, w, Bmat, bb, i, &v[i], &gamma );
          /* multiply A_{k-1} by v and by gamma */
    pkn_NRBSymMultf ( n, prof1, rows1[0], rows1, 1, 1, v, 1, p );
    pkn_MultMatrixNumf ( 1, n, 0, p, gamma, 0, p );
    s = pkn_ScalarProductf ( n-i, &v[i], &p[i] )*0.5*gamma;
    pkn_AddMatrixMf ( 1, n-i, 0, &p[i], 0, &v[i], -s, 0, &p[i] );
        /* now compute A_k = A_{k-1} - p*v^T - v*p^T */
          /* first it is necessary to find the profile of the result */
    memcpy ( prof2, prof1, n*sizeof(int) );
    for ( lv = i; lv < n; lv++ )    
      if ( v[lv] )
        break;
    for ( lp = 0; lp < n; lp++ )
      if ( p[lp] )
        break;
    for ( j = i; j < n; j++ )
      if ( v[j] && lp < prof2[j] ) prof2[j] = lp;
    for ( j = 0; j < n; j++ )
      if ( p[j] && lv < prof2[j] ) prof2[j] = lv;
          /* allocate a suitable memory block */
    size = pkn_NRBArraySize ( n, prof2 );
    PKV_MALLOC ( rows2[0], size*sizeof(float) );
    if ( !rows2[0] )
      goto failure;
    if ( !pkn_NRBFindRowsf ( n, prof2, rows2[0], rows2 ) )
      goto failure;
          /* copy the matrix A_{k-1} */
    memset ( rows2[0], 0, size*sizeof(float) );
    for ( j = 0; j < n; j++ )
      memcpy ( &rows2[j][prof1[j]], &rows1[j][prof1[j]],
               (j-prof1[j]+1)*sizeof(float) );
          /* add the two rank-1 matrices */
    for ( j = 0; j < n; j++ ) {
      if ( v[j] )
        for ( k = lp; k <= j; k++ )
          rows2[j][k] -= v[j]*p[k];
      if ( p[j] )
        for ( k = lv; k <= j; k++ )
          rows2[j][k] -= p[j]*v[k];
    }
        /* now deallocate A_{k-1} and copy A_k */
    if ( i < w-1 )
      PKV_FREE ( rows1[0] );
    memcpy ( prof1, prof2, n*sizeof(int) );
    memcpy ( rows1, rows2, n*sizeof(float*) );
  }
        /* assign the result to output parameters */
  memcpy ( qaprof, prof1, n*sizeof(int) );
  memcpy ( QArows, rows1, n*sizeof(float*) );

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  if ( i < w-1 ) PKV_FREE ( rows1[0] );
  return false;
} /*pkn_NRBComputeQSQTf*/

