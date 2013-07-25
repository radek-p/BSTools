
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2007                                  */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

#include <math.h>
#include <string.h>
#include <stdlib.h>

#include "pkvaria.h"

#undef CONST_
#define CONST_

#include "pknum.h"

/* ///////////////////////////////////////////////////////////// */
void pkn_MatrixLowerTrMultf ( int m, int n, int bpitch, CONST_ float *b,
                              CONST_ float *l, int xpitch, float *x )
{
  int    i, j, k;
  float  *xx, *bb;
  double s;

  for ( i = 0, xx = x, bb = b;  i < m;  i++, xx += xpitch, bb += bpitch )
    for ( j = 0; j < n; j++ ) {
      for ( s = 0.0, k = j;  k < n;  k++ )
        s += bb[k]*l[pkn_SymMatIndex(k,j)];
      xx[j] = (float)s;
    }
} /*pkn_MatrixLowerTrMultf*/

void pkn_MatrixUpperTrMultf ( int m, int n, int bpitch, CONST_ float *b,
                              CONST_ float *l, int xpitch, float *x )
{
  int    i, j, k;
  float  *xx, *bb;
  double s;

  for ( i = 0, xx = x, bb = b;  i < m;  i++, xx += xpitch, bb += bpitch )
    for ( j = 0; j < n; j++ ) {
      for ( s = 0.0, k = 0;  k <= j;  k++ )
        s += bb[k]*l[pkn_SymMatIndex(k,j)];
      xx[j] = (float)s;
    }
} /*pkn_MatrixUpperTrMultf*/

void pkn_MatrixLowerTrSolvef ( int m, int n, int bpitch, CONST_ float *b,
                               CONST_ float *l, int xpitch, float *x )
{
  int    i, j, k;
  float  *xx, *bb;
  double s;

  for ( k = 0, xx = x, bb = b;  k < m;  k++, xx += xpitch, bb += bpitch )
    for ( i = n-1; i >= 0; i-- ) {
      s = bb[i];
      for ( j = i+1; j < n; j++ )
        s -= xx[j]*l[pkn_SymMatIndex(j,i)];
      xx[i] = (float)(s/l[pkn_SymMatIndex(i,i)]);
    }
} /*pkn_MatrixLowerTrSolvef*/

void pkn_MatrixUpperTrSolvef ( int m, int n, int bpitch, CONST_ float *b,
                               CONST_ float *l, int xpitch, float *x )
{
  int    i, j, k;
  float  *xx, *bb;
  double s;
  
  for ( k = 0, xx = x, bb = b;  k < m;  k++, xx += xpitch, bb += bpitch )
    for ( i = 0; i < n; i++ ) {
      s = bb[i];
      for ( j = 0; j < i; j++ )
        s -= xx[j]*l[pkn_SymMatIndex(i,j)];
      xx[i] = (float)(s/l[pkn_SymMatIndex(i,i)]);
    }
} /*pkn_MatrixUpperTrSolvef*/

/* ///////////////////////////////////////////////////////////// */
void pkn_MatrixLowerTrMultAddf ( int m, int n, int bpitch, CONST_ float *b,
                                 CONST_ float *l, int xpitch, float *x )
{
  int    i, j, k;
  float  *xx, *bb;
  double s;

  for ( i = 0, xx = x, bb = b;  i < m;  i++, xx += xpitch, bb += bpitch )
    for ( j = 0; j < n; j++ ) {
      for ( s = 0.0, k = j;  k < n;  k++ )
        s += bb[k]*l[pkn_SymMatIndex(k,j)];
      xx[j] += (float)s;
    }
} /*pkn_MatrixLowerTrMultAddf*/

void pkn_MatrixUpperTrMultAddf ( int m, int n, int bpitch, CONST_ float *b,
                                 CONST_ float *l, int xpitch, float *x )
{
  int    i, j, k;
  float  *xx, *bb;
  double s;

  for ( i = 0, xx = x, bb = b;  i < m;  i++, xx += xpitch, bb += bpitch )
    for ( j = 0; j < n; j++ ) {
      for ( s = 0.0, k = 0;  k <= j;  k++ )
        s += bb[k]*l[pkn_SymMatIndex(k,j)];
      xx[j] += (float)s;
    }
} /*pkn_MatrixUpperTrMultAddf*/

boolean pkn_MatrixLowerTrSolveAddf ( int m, int n, int bpitch, CONST_ float *b,
                                     CONST_ float *l, int xpitch, float *x )
{
  void   *sp;
  int    i, j, k;
  float  *xx, *bb, *z;
  double s;

  sp = pkv_GetScratchMemTop ();
  if ( !(z = pkv_GetScratchMemf ( n )) )
    goto failure;
  for ( k = 0, xx = x, bb = b;  k < m;  k++, xx += xpitch, bb += bpitch ) {
    for ( i = n-1; i >= 0; i-- ) {
      s = bb[i];
      for ( j = i+1; j < n; j++ )
        s -= z[j]*l[pkn_SymMatIndex(j,i)];
      z[i] = (float)(s/l[pkn_SymMatIndex(i,i)]);
    }
    pkn_AddMatrixf ( 1, n, 0, xx, 0, z, 0, xx );
  }
  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*pkn_MatrixLowerTrSolveAddf*/

boolean pkn_MatrixUpperTrSolveAddf ( int m, int n, int bpitch, CONST_ float *b,
                                     CONST_ float *l, int xpitch, float *x )
{
  void   *sp;
  int    i, j, k;
  float  *xx, *bb, *z;
  double s;

  sp = pkv_GetScratchMemTop ();
  if ( !(z = pkv_GetScratchMemf ( n )) )
    goto failure;
  for ( k = 0, xx = x, bb = b;  k < m;  k++, xx += xpitch, bb += bpitch ) {
    for ( i = 0; i < n; i++ ) {
      s = bb[i];
      for ( j = 0; j < i; j++ )
        s -= z[j]*l[pkn_SymMatIndex(i,j)];
      z[i] = (float)(s/l[pkn_SymMatIndex(i,i)]);
    }
    pkn_AddMatrixf ( 1, n, 0, xx, 0, z, 0, xx );
  }
  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*pkn_MatrixUpperTrSolveAddf*/

/* ///////////////////////////////////////////////////////////// */
void pkn_MatrixLowerTrMultSubf ( int m, int n, int bpitch, CONST_ float *b,
                                 CONST_ float *l, int xpitch, float *x )
{
  int    i, j, k;
  float  *xx, *bb;
  double s;

  for ( i = 0, xx = x, bb = b;  i < m;  i++, xx += xpitch, bb += bpitch )
    for ( j = 0; j < n; j++ ) {
      for ( s = 0.0, k = j;  k < n;  k++ )
        s += bb[k]*l[pkn_SymMatIndex(k,j)];
      xx[j] -= (float)s;
    }
} /*pkn_MatrixLowerTrMultSubf*/

void pkn_MatrixUpperTrMultSubf ( int m, int n, int bpitch, CONST_ float *b,
                                 CONST_ float *l, int xpitch, float *x )
{
  int    i, j, k;
  float  *xx, *bb;
  double s;

  for ( i = 0, xx = x, bb = b;  i < m;  i++, xx += xpitch, bb += bpitch )
    for ( j = 0; j < n; j++ ) {
      for ( s = 0.0, k = 0;  k <= j;  k++ )
        s += bb[k]*l[pkn_SymMatIndex(k,j)];
      xx[j] -= (float)s;
    }
} /*pkn_MatrixUpperTrMultSubf*/

boolean pkn_MatrixLowerTrSolveSubf ( int m, int n, int bpitch, CONST_ float *b,
                                     CONST_ float *l, int xpitch, float *x )
{
  void   *sp;
  int    i, j, k;
  float  *xx, *bb, *z;
  double s;

  sp = pkv_GetScratchMemTop ();
  if ( !(z = pkv_GetScratchMemf ( n )) )
    goto failure;
  for ( k = 0, xx = x, bb = b;  k < m;  k++, xx += xpitch, bb += bpitch ) {
    for ( i = n-1; i >= 0; i-- ) {
      s = bb[i];
      for ( j = i+1; j < n; j++ )
        s -= z[j]*l[pkn_SymMatIndex(j,i)];
      z[i] = (float)(s/l[pkn_SymMatIndex(i,i)]);
    }
    pkn_SubtractMatrixf ( 1, n, 0, xx, 0, z, 0, xx );
  }
  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*pkn_MatrixLowerTrSolveSubf*/

boolean pkn_MatrixUpperTrSolveSubf ( int m, int n, int bpitch, CONST_ float *b,
                                     CONST_ float *l, int xpitch, float *x )
{
  void   *sp;
  int    i, j, k;
  float  *xx, *bb, *z;
  double s;

  sp = pkv_GetScratchMemTop ();
  if ( !(z = pkv_GetScratchMemf ( n )) )
    goto failure;
  for ( k = 0, xx = x, bb = b;  k < m;  k++, xx += xpitch, bb += bpitch ) {
    for ( i = 0; i < n; i++ ) {
      s = bb[i];
      for ( j = 0; j < i; j++ )
        s -= z[j]*l[pkn_SymMatIndex(i,j)];
      z[i] = (float)(s/l[pkn_SymMatIndex(i,i)]);
    }
    pkn_SubtractMatrixf ( 1, n, 0, xx, 0, z, 0, xx );
  }
  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*pkn_MatrixUpperTrSolveSubf*/

