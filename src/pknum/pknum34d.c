
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2008                                  */
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

/* ///////////////////////////////////////////////////////////////////////// */
boolean pkn_NRBFindRowsd ( int n, const int *prof, CONST_ double *a,
                           double **row )
{
  int    i;
  double *r;

  if ( prof[0] != 0 )
    return false;
  row[0] = r = a;
  for ( i = 1; i < n; i++ ) {
    if ( prof[i] < 0 || prof[i] > i )
      return false;
    r += i-prof[i]+1;
    row[i] = r-i;
  }
  return true;
} /*pkn_NRBFindRowsd*/

/* The parameter abort may be NULL, or it may point to a variable          */
/* initially set to 0 (false). If this variable is set to a nonzero value, */     
/* it will cause the almost immediate termination of the procedure, which  */
/* may be desirable in a multithread program.                              */         
boolean pkn_NRBSymCholeskyDecompd ( int n, const int *prof, double *a,
                                    double **row, boolean *abort )
{
  void    *sp;
  double  **r;
  int     i, j, k;
  double  ljj, ljk, aij;
  boolean _abort;

  sp = pkv_GetScratchMemTop ();
  if ( !row ) {
    r = pkv_GetScratchMem ( n*sizeof(double*) );
    if ( !r )
      goto failure;
    if ( !pkn_NRBFindRowsd ( n, prof, a, r ) )
      goto failure;
  }
  else
    r = row;
  if ( !abort ) {
    _abort = false;
    abort = &_abort;
  }
        /* it is assumed that if the parameter row passed to the procedure */
        /* is not NULL, then the profile is correct and the row array is */
        /* also correct */
        /* now the actual Cholesky's decomposition */
  for ( j = 0; j < n; j++ ) {
    ljj = r[j][j];
    for ( k = prof[j]; k < j; k++ ) {
      ljk = r[j][k];
      ljj -= ljk*ljk;
    }
    if ( ljj <= 0.0 )
      goto failure;  /* this matrix is not positive-definite */
    r[j][j] = (double)(ljj = sqrt(ljj));
    for ( i = j+1; i < n; i++ ) {
      if ( prof[i] <= j ) {
        aij = r[i][j];
        for ( k = max(prof[i],prof[j]); k < j; k++ )
          aij -= r[i][k]*r[j][k];
        r[i][j] = (double)(aij/ljj);
      }
    }
    if ( *abort )
      goto failure;
  }
  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*pkn_NRBSymCholeskyDecompd*/

boolean pkn_NRBSymMultd ( int n, const int *prof, CONST_ double *a,
                          double **row,
                          int spdimen, int xpitch, CONST_ double *x,
                          int ypitch, double *y )
{
  void   *sp;
  double **r;
  int    i, k;
  double *yy;

  sp = pkv_GetScratchMemTop ();
  if ( !row ) {
    r = pkv_GetScratchMem ( n*sizeof(double*) );
    if ( !r )
      goto failure;
    if ( !pkn_NRBFindRowsd ( n, prof, a, r ) )
      goto failure;
  }
  else
    r = row;
        /* it is assumed that if the parameter row passed to the procedure */
        /* is not NULL, then the profile is correct and the row array is */
        /* also correct */
        /* now the actual multiplication */
  for ( i = 0, yy = y;  i < n;  i++, yy += ypitch )
    memset ( yy, 0, spdimen*sizeof(double) );

  for ( i = 0; i < n; i++ ) {
    k = i-prof[i];
    pkn_MultMatrixAddd ( 1, k+1, 0, &r[i][prof[i]], spdimen,
                         xpitch, &x[prof[i]*xpitch], ypitch, &y[i*ypitch] );
    pkn_MultMatrixAddd ( k, 1, 1, &r[i][prof[i]], spdimen,
                         xpitch, &x[i*xpitch], ypitch, &y[prof[i]*ypitch] );
  }

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*pkn_NRBSymMultd*/

boolean pkn_NRBLowerTrMultd ( int n, const int *prof, CONST_ double *a,
                              double **row,
                              int spdimen, int xpitch, CONST_ double *x,
                              int ypitch, double *y )
{
  void   *sp;
  double **r;
  int    i;
  double *yy;

  sp = pkv_GetScratchMemTop ();
  if ( !row ) {
    r = pkv_GetScratchMem ( n*sizeof(double*) );
    if ( !r )
      goto failure;
    if ( !pkn_NRBFindRowsd ( n, prof, a, r ) )
      goto failure;
  }
  else
    r = row;
        /* it is assumed that if the parameter row passed to the procedure */
        /* is not NULL, then the profile is correct and the row array is */
        /* also correct */
        /* now the actual multiplication */
  for ( i = 0, yy = y;  i < n;  i++, yy += ypitch )
    memset ( yy, 0, spdimen*sizeof(double) );

  for ( i = 0; i < n; i++ )
    pkn_MultMatrixAddd ( 1, i-prof[i]+1, 0, &r[i][prof[i]], spdimen,
                         xpitch, &x[prof[i]*xpitch], ypitch, &y[i*ypitch] );

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*pkn_NRBLowerTrMultd*/

boolean pkn_NRBUpperTrMultd ( int n, const int *prof, CONST_ double *a,
                              double **row,
                              int spdimen, int xpitch, CONST_ double *x,
                              int ypitch, double *y )
{
  void   *sp;
  double **r;
  int    i;
  double *yy;

  sp = pkv_GetScratchMemTop ();
  if ( !row ) {
    r = pkv_GetScratchMem ( n*sizeof(double*) );
    if ( !r )
      goto failure;
    if ( !pkn_NRBFindRowsd ( n, prof, a, r ) )
      goto failure;
  }
  else
    r = row;
        /* it is assumed that if the parameter row passed to the procedure */
        /* is not NULL, then the profile is correct and the row array is */
        /* also correct */
        /* now the actual multiplication */
  for ( i = 0, yy = y;  i < n;  i++, yy += ypitch )
    memset ( yy, 0, spdimen*sizeof(double) );

  for ( i = 0; i < n; i++ )
    pkn_MultMatrixAddd ( i-prof[i]+1, 1, 1, &r[i][prof[i]], spdimen,
                         xpitch, &x[i*xpitch], ypitch, &y[prof[i]*ypitch] );

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*pkn_NRBUpperTrMultd*/

boolean pkn_NRBLowerTrSolved ( int n, const int *prof, CONST_ double *l,
                               double **row,
                               int spdimen, int bpitch, CONST_ double *b,
                               int xpitch, double *x )
{
  void   *sp;
  double **r;
  int    i;

  sp = pkv_GetScratchMemTop ();
  if ( !row ) {
    r = pkv_GetScratchMem ( n*sizeof(double*) );
    if ( !r )
      goto failure;
    if ( !pkn_NRBFindRowsd ( n, prof, l, r ) )
      goto failure;
  }
  else
    r = row;
        /* it is assumed that if the parameter row passed to the procedure */
        /* is not NULL, then the profile is correct and the row array is */
        /* also correct */
        /* now solving the system with the lower triangular matrix */
  if ( x != b )
    pkv_Selectd ( n, spdimen, bpitch, xpitch, b, x );

  for ( i = 0; i < n; i++ ) {
    if ( i > prof[i] )
      pkn_MultMatrixSubd ( 1, i-prof[i], 0, &r[i][prof[i]],
                           spdimen, xpitch, &x[prof[i]*xpitch],
                           xpitch, &x[i*xpitch] );
    pkn_MultMatrixNumd ( 1, spdimen, 0, &x[i*xpitch], 1.0/r[i][i],
                         0, &x[i*xpitch] );
  }

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*pkn_NRBLowerTrSolved*/

boolean pkn_NRBUpperTrSolved ( int n, const int *prof, CONST_ double *l,
                               double **row,
                               int spdimen, int bpitch, CONST_ double *b,
                               int xpitch, double *x )
{
  void   *sp;
  double **r;
  int    i;

  sp = pkv_GetScratchMemTop ();
  if ( !row ) {
    r = pkv_GetScratchMem ( n*sizeof(double*) );
    if ( !r )
      goto failure;
    if ( !pkn_NRBFindRowsd ( n, prof, l, r ) )
      goto failure;
  }
  else
    r = row;
        /* it is assumed that if the parameter row passed to the procedure */
        /* is not NULL, then the profile is correct and the row array is */
        /* also correct */
        /* now solving the system with the upper triangular matrix */
  if ( x != b )
    pkv_Selectd ( n, spdimen, bpitch, xpitch, b, x );

  for ( i = n-1; i >= 0; i-- ) {
    pkn_MultMatrixNumd ( 1, spdimen, 0, &x[i*xpitch], 1.0/r[i][i],
                         0, &x[i*xpitch] );
    if ( i > prof[i] )
      pkn_MultMatrixSubd ( i-prof[i], 1, 1, &r[i][prof[i]], spdimen,
                           0, &x[i*xpitch], xpitch, &x[prof[i]*xpitch] );
  }

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*pkn_NRBUpperTrSolved*/

