
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2005                                  */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

#include <math.h>
#include <string.h>

#include "pkvaria.h"

#undef CONST_
#define CONST_

#include "pknum.h"

/* ////////////////////////////////////////////////////// */
/* Cholesky decomposition of a symmetric                  */
/* positive-definite matrix                               */

boolean pkn_CholeskyDecompf ( int n, float *a )
{
  int    i, j, k, ij;
  double ljj, ljk, aij;

  for ( j = 0; j < n; j++ ) {
    ljj = a[pkn_SymMatIndex(j,j)];
    for ( k = 0; k < j; k++ ) {
      ljk = a[pkn_SymMatIndex(j,k)];
      ljj -= ljk*ljk;
    }
    if ( ljj <= 0.0 )
      return false; /* this matrix is not positive-definite */

    a[pkn_SymMatIndex(j,j)] = (float)(ljj = sqrt(ljj));
    for ( i = j+1; i < n; i++ ) {
      ij = pkn_SymMatIndex(i,j);
      aij = a[ij];
      for ( k = 0; k < j; k++ )
        aij -= a[pkn_SymMatIndex(i,k)]*a[pkn_SymMatIndex(j,k)];
      a[ij] = (float)(aij/ljj);
    }
  }
  return true;
} /*pkn_CholeskyDecompf*/

void pkn_SymMatrixMultf ( int n, CONST_ float *a, int spdimen,
                          int bpitch, CONST_ float *b,
                          int xpitch, float *x )
{
  int   i, j, i0;
  float aa, *xx, *xxx, *bb, *bbb;

  for ( i = 0, xx = x;  i < n;  i++, xx += xpitch )
    memset ( xx, 0, spdimen*sizeof(float) );

  for ( i = 0, xx = x, bb = b;
        i < n;
        i++, xx += xpitch, bb += bpitch ) {
    i0 = pkn_SymMatIndex ( i, 0 );
    for ( j = 0, xxx = x, bbb = b;
          j < i;
          j++, xxx += xpitch, bbb += bpitch ) {
      aa = a[i0+j];
      pkn_AddMatrixMf ( 1, spdimen, 0, xx, 0, bbb, aa, 0, xx );
      pkn_AddMatrixMf ( 1, spdimen, 0, xxx, 0, bb, aa, 0, xxx );
    }
    aa = a[i0+i];
    pkn_AddMatrixMf ( 1, spdimen, 0, xx, 0, bb, aa, 0, xx );
  }
} /*pkn_SymMatrixMultf*/

void pkn_LowerTrMatrixMultf ( int n, CONST_ float *l, int spdimen,
                              int bpitch, CONST_ float *b,
                              int xpitch, float *x )
{
  int   i, j, i0;
  float ll, *xx, *bbb;

  for ( i = 0, xx = x;  i < n;  i++, xx += xpitch )
    memset ( xx, 0, spdimen*sizeof(float) );

  for ( i = 0, xx = x;  i < n;  i++, xx += xpitch ) {
    i0 = pkn_SymMatIndex ( i, 0 );
    for ( j = 0, bbb = b;  j <= i;  j++, bbb += bpitch ) {
      ll = l[i0+j];
      pkn_AddMatrixMf ( 1, spdimen, 0, xx, 0, bbb, ll, 0, xx );
    }
  }
} /*pkn_LowerTrMatrixMultf*/

void pkn_UpperTrMatrixMultf ( int n, CONST_ float *l, int spdimen,
                              int bpitch, CONST_ float *b,
                              int xpitch, float *x )
{
  int   i, j, i0;
  float ll, *xxx, *bb;

  for ( i = 0, xxx = x;  i < n;  i++, xxx += xpitch )
    memset ( xxx, 0, spdimen*sizeof(float) );

  for ( i = 0, bb = b;  i < n;  i++, bb += bpitch ) {
    i0 = pkn_SymMatIndex ( i, 0 );
    for ( j = 0, xxx = x;  j <= i;  j++, xxx += xpitch ) {
      ll = l[i0+j];
      pkn_AddMatrixMf ( 1, spdimen, 0, xxx, 0, bb, ll, 0, xxx );
    }
  }
} /*pkn_UpperTrMatrixMultf*/

void pkn_LowerTrMatrixSolvef ( int n, const float *l, int spdimen,
                               int bpitch, const float *b,
                               int xpitch, float *x )
{
  int   i, j, i0;
  float ll, *xx, *xxx;

  if ( x != b )
    pkv_Selectf ( n, spdimen, bpitch, xpitch, b, x );

  for ( i = 0, xx = x;  i < n;  i++, xx += xpitch ) {
    i0 = pkn_SymMatIndex ( i, 0 );
    for ( j = 0, xxx = x;  j < i;  j++, xxx += xpitch ) {
      ll = l[i0+j];
      pkn_AddMatrixMf ( 1, spdimen, 0, xx, 0, xxx, -ll, 0, xx );
    }
    ll = l[i0+i];
    pkn_MultMatrixNumf ( 1, spdimen, 0, xx, 1.0/ll, 0, xx );
  }
} /*pkn_LowerTrMatrixSolvef*/

void pkn_UpperTrMatrixSolvef ( int n, const float *l, int spdimen,
                               int bpitch, const float *b,
                               int xpitch, float *x )
{
  int   i, j, i0;
  float ll, *xx, *xxx;

  if ( x != b )
    pkv_Selectf ( n, spdimen, bpitch, xpitch, b, x );

  for ( i = n-1, xx = x+((n-1)*xpitch);  i >= 0;  i--, xx -= xpitch ) {
    i0 = pkn_SymMatIndex ( i, 0 );
    ll = l[i0+i];
    pkn_MultMatrixNumf ( 1, spdimen, 0, xx, 1.0/ll, 0, xx );
    for ( j = 0, xxx = x;  j < i;  j++, xxx += xpitch ) {
      ll = l[i0+j];
      pkn_AddMatrixMf ( 1, spdimen, 0, xxx, 0, xx, -ll, 0, xxx );
    }
  }
} /*pkn_UpperTrMatrixSolvef*/

/* ////////////////////////////////////////////////////// */
/* conversion between the packed and full representations */
/* of symmetric and triangular matrices                   */

void pkn_SymToFullMatrixf ( int n, const float *syma,
                            int pitch, float *fulla )
{
  int i, k, l;

  memcpy ( fulla, syma, sizeof(float) );
  for ( i = k = 1, l = pitch;  i < n;  i++, k += i, l += pitch ) {
    memcpy ( &fulla[l], &syma[k], (i+1)*sizeof(float) );
    pkv_Selectf ( i, 1, 1, pitch, &syma[k], &fulla[i] );
  }
} /*pkn_SymToFullMatrixf*/

void pkn_FullToSymMatrixf ( int n, int pitch, const float *fulla,
                            float *syma )
{
  int i, k, l;

  for ( i = k = l = 0;  i < n;  i++, k += i, l += pitch )
    memcpy ( &syma[k], &fulla[l], (i+1)*sizeof(float) );
} /*pkn_FullToSymMatrixf*/

void pkn_LTrToFullMatrixf ( int n, const float *ltra,
                            int pitch, float *fulla )
{
  int i, k, l;

  for ( i = k = l = 0;  i < n-1;  i++, k += i, l += pitch ) {
    memcpy ( &fulla[l], &ltra[k], (i+1)*sizeof(float) );
    memset ( &fulla[l+i+1], 0, (n-i-1)*sizeof(float) );
  }
  memcpy ( &fulla[l], &ltra[k], n*sizeof(float) );
} /*pkn_LTrToFullMatrixf*/

void pkn_UTrToFullMatrixf ( int n, const float *utra,
                            int pitch, float *fulla )
{
  int i, k, l;

  memcpy ( fulla, utra, sizeof(float) );
  for ( i = k = 1, l = pitch;  i < n;  i++, k += i, l += pitch ) {
    memset ( &fulla[l], 0, i*sizeof(float) );
    pkv_Selectf ( i+1, 1, 1, pitch, &utra[k], &fulla[i] );
  }
} /*pkn_UTrToFullMatrixf*/

void pkn_FullToUTrMatrixf ( int n, int pitch, const float *fulla,
                            float *utra )
{
  int i, k;

  for ( i = k = 0;  i < n;  i++, k += i )
    pkv_Selectf ( i+1, 1, pitch, 1, &fulla[i], &utra[k] );
} /*pkn_FullToUTrMatrixf*/

