
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
void pkn_MatrixLowerTrMultd ( int m, int n, int bpitch, CONST_ double *b,
                              CONST_ double *l, int xpitch, double *x )
{
  int    i, j, k;
  double *xx, *bb;
  long double s;

  for ( i = 0, xx = x, bb = b;  i < m;  i++, xx += xpitch, bb += bpitch )
    for ( j = 0; j < n; j++ ) {
      for ( s = 0.0, k = j;  k < n;  k++ )
        s += bb[k]*l[pkn_SymMatIndex(k,j)];
      xx[j] = (double)s;
    }
} /*pkn_MatrixLowerTrMultd*/

void pkn_MatrixUpperTrMultd ( int m, int n, int bpitch, CONST_ double *b,
                              CONST_ double *l, int xpitch, double *x )
{
  int    i, j, k;
  double *xx, *bb;
  long double s;

  for ( i = 0, xx = x, bb = b;  i < m;  i++, xx += xpitch, bb += bpitch )
    for ( j = 0; j < n; j++ ) {
      for ( s = 0.0, k = 0;  k <= j;  k++ )
        s += bb[k]*l[pkn_SymMatIndex(k,j)];
      xx[j] = (double)s;
    }
} /*pkn_MatrixUpperTrMultd*/

void pkn_MatrixLowerTrSolved ( int m, int n, int bpitch, CONST_ double *b,
                               CONST_ double *l, int xpitch, double *x )
{
  int    i, j, k;
  double *xx, *bb;
  long double s;

  for ( k = 0, xx = x, bb = b;  k < m;  k++, xx += xpitch, bb += bpitch )
    for ( i = n-1; i >= 0; i-- ) {
      s = bb[i];
      for ( j = i+1; j < n; j++ )
        s -= xx[j]*l[pkn_SymMatIndex(j,i)];
      xx[i] = (double)(s/l[pkn_SymMatIndex(i,i)]);
    }
} /*pkn_MatrixLowerTrSolved*/

void pkn_MatrixUpperTrSolved ( int m, int n, int bpitch, CONST_ double *b,
                               CONST_ double *l, int xpitch, double *x )
{
  int    i, j, k;
  double *xx, *bb;
  long double s;
  
  for ( k = 0, xx = x, bb = b;  k < m;  k++, xx += xpitch, bb += bpitch )
    for ( i = 0; i < n; i++ ) {
      s = bb[i];
      for ( j = 0; j < i; j++ )
        s -= xx[j]*l[pkn_SymMatIndex(i,j)];
      xx[i] = (double)(s/l[pkn_SymMatIndex(i,i)]);
    }
} /*pkn_MatrixUpperTrSolved*/

/* ///////////////////////////////////////////////////////////// */
void pkn_MatrixLowerTrMultAddd ( int m, int n, int bpitch, CONST_ double *b,
                                 CONST_ double *l, int xpitch, double *x )
{
  int    i, j, k;
  double *xx, *bb;
  long double s;

  for ( i = 0, xx = x, bb = b;  i < m;  i++, xx += xpitch, bb += bpitch )
    for ( j = 0; j < n; j++ ) {
      for ( s = 0.0, k = j;  k < n;  k++ )
        s += bb[k]*l[pkn_SymMatIndex(k,j)];
      xx[j] += (double)s;
    }
} /*pkn_MatrixLowerTrMultAddd*/

void pkn_MatrixUpperTrMultAddd ( int m, int n, int bpitch, CONST_ double *b,
                                 CONST_ double *l, int xpitch, double *x )
{
  int    i, j, k;
  double *xx, *bb;
  long double s;

  for ( i = 0, xx = x, bb = b;  i < m;  i++, xx += xpitch, bb += bpitch )
    for ( j = 0; j < n; j++ ) {
      for ( s = 0.0, k = 0;  k <= j;  k++ )
        s += bb[k]*l[pkn_SymMatIndex(k,j)];
      xx[j] += (double)s;
    }
} /*pkn_MatrixUpperTrMultAddd*/

boolean pkn_MatrixLowerTrSolveAddd ( int m, int n, int bpitch, CONST_ double *b,
                                     CONST_ double *l, int xpitch, double *x )
{
  void   *sp;
  int    i, j, k;
  double *xx, *bb, *z;
  long double s;

  sp = pkv_GetScratchMemTop ();
  if ( !(z = pkv_GetScratchMemd ( n )) )
    goto failure;
  for ( k = 0, xx = x, bb = b;  k < m;  k++, xx += xpitch, bb += bpitch ) {
    for ( i = n-1; i >= 0; i-- ) {
      s = bb[i];
      for ( j = i+1; j < n; j++ )
        s -= z[j]*l[pkn_SymMatIndex(j,i)];
      z[i] = (double)(s/l[pkn_SymMatIndex(i,i)]);
    }
    pkn_AddMatrixd ( 1, n, 0, xx, 0, z, 0, xx );
  }
  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*pkn_MatrixLowerTrSolveAddd*/

boolean pkn_MatrixUpperTrSolveAddd ( int m, int n, int bpitch, CONST_ double *b,
                                     CONST_ double *l, int xpitch, double *x )
{
  void   *sp;
  int    i, j, k;
  double *xx, *bb, *z;
  long double s;

  sp = pkv_GetScratchMemTop ();
  if ( !(z = pkv_GetScratchMemd ( n )) )
    goto failure;
  for ( k = 0, xx = x, bb = b;  k < m;  k++, xx += xpitch, bb += bpitch ) {
    for ( i = 0; i < n; i++ ) {
      s = bb[i];
      for ( j = 0; j < i; j++ )
        s -= z[j]*l[pkn_SymMatIndex(i,j)];
      z[i] = (double)(s/l[pkn_SymMatIndex(i,i)]);
    }
    pkn_AddMatrixd ( 1, n, 0, xx, 0, z, 0, xx );
  }
  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*pkn_MatrixUpperTrSolveAddd*/

/* ///////////////////////////////////////////////////////////// */
void pkn_MatrixLowerTrMultSubd ( int m, int n, int bpitch, CONST_ double *b,
                                 CONST_ double *l, int xpitch, double *x )
{
  int    i, j, k;
  double *xx, *bb;
  long double s;

  for ( i = 0, xx = x, bb = b;  i < m;  i++, xx += xpitch, bb += bpitch )
    for ( j = 0; j < n; j++ ) {
      for ( s = 0.0, k = j;  k < n;  k++ )
        s += bb[k]*l[pkn_SymMatIndex(k,j)];
      xx[j] -= (double)s;
    }
} /*pkn_MatrixLowerTrMultSubd*/

void pkn_MatrixUpperTrMultSubd ( int m, int n, int bpitch, CONST_ double *b,
                                 CONST_ double *l, int xpitch, double *x )
{
  int    i, j, k;
  double *xx, *bb;
  long double s;

  for ( i = 0, xx = x, bb = b;  i < m;  i++, xx += xpitch, bb += bpitch )
    for ( j = 0; j < n; j++ ) {
      for ( s = 0.0, k = 0;  k <= j;  k++ )
        s += bb[k]*l[pkn_SymMatIndex(k,j)];
      xx[j] -= (double)s;
    }
} /*pkn_MatrixUpperTrMultSubd*/

boolean pkn_MatrixLowerTrSolveSubd ( int m, int n, int bpitch, CONST_ double *b,
                                     CONST_ double *l, int xpitch, double *x )
{
  void   *sp;
  int    i, j, k;
  double *xx, *bb, *z;
  long double s;

  sp = pkv_GetScratchMemTop ();
  if ( !(z = pkv_GetScratchMemd ( n )) )
    goto failure;
  for ( k = 0, xx = x, bb = b;  k < m;  k++, xx += xpitch, bb += bpitch ) {
    for ( i = n-1; i >= 0; i-- ) {
      s = bb[i];
      for ( j = i+1; j < n; j++ )
        s -= z[j]*l[pkn_SymMatIndex(j,i)];
      z[i] = (double)(s/l[pkn_SymMatIndex(i,i)]);
    }
    pkn_SubtractMatrixd ( 1, n, 0, xx, 0, z, 0, xx );
  }
  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*pkn_MatrixLowerTrSolveSubd*/

boolean pkn_MatrixUpperTrSolveSubd ( int m, int n, int bpitch, CONST_ double *b,
                                     CONST_ double *l, int xpitch, double *x )
{
  void   *sp;
  int    i, j, k;
  double *xx, *bb, *z;
  long double s;

  sp = pkv_GetScratchMemTop ();
  if ( !(z = pkv_GetScratchMemd ( n )) )
    goto failure;
  for ( k = 0, xx = x, bb = b;  k < m;  k++, xx += xpitch, bb += bpitch ) {
    for ( i = 0; i < n; i++ ) {
      s = bb[i];
      for ( j = 0; j < i; j++ )
        s -= z[j]*l[pkn_SymMatIndex(i,j)];
      z[i] = (double)(s/l[pkn_SymMatIndex(i,i)]);
    }
    pkn_SubtractMatrixd ( 1, n, 0, xx, 0, z, 0, xx );
  }
  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*pkn_MatrixUpperTrSolveSubd*/

