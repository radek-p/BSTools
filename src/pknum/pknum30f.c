
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

#include "msgpool.h"

/* ///////////////////////////////////////////////////////////// */
/* Symmetric and triangular matrices with Block2 block structure */

boolean pkn_Block2CholeskyDecompMf ( int k, int r, int s, int t, float *A )
{
#define BlPos(i,j) pkn_Block2FindBlockPos( k, r, s, t, i, j )
  int   i;
  float *Aii, *Aij, *Aik, *Ail;

       /* block columns from 0 to k-1 */
  for ( i = 0; i < k; i++ ) {
    Aii = &A[BlPos(i,i)];
    if ( !pkn_CholeskyDecompf ( r, Aii ) )
      return false;
    Aij = &A[BlPos(i+k,i)];
    pkn_MatrixUpperTrSolvef ( s, r, r, Aij, Aii, r, Aij );
    if ( i < k-1 )
      Aij = &A[BlPos(i+k+1,i)];
    else
      Aij = &A[BlPos(k,i)];
    pkn_MatrixUpperTrSolvef ( s, r, r, Aij, Aii, r, Aij );
    Aij = &A[BlPos(2*k,i)];
    pkn_MatrixUpperTrSolvef ( t, r, r, Aij, Aii, r, Aij );
  }

        /* block column k */
  Aii = &A[BlPos(k,k)];
  Aij = Aik = &A[BlPos(k,0)];
  pkn_SymMatSubAATf ( s, Aii, r, r, Aij );
  Aij = &A[BlPos(k,k-1)];
  pkn_SymMatSubAATf ( s, Aii, r, r, Aij );
  if ( !pkn_CholeskyDecompf ( s, Aii ) )
    return false;
  Aij = &A[BlPos(k+1,k)];
  Ail = &A[BlPos(k+1,0)];
  pkn_MultMatrixTSubf ( s, r, r, Ail, s, r, Aik, s, Aij );
  pkn_MatrixUpperTrSolvef ( s, s, s, Aij, Aii, s, Aij );
  Aij = &A[BlPos(2*k-1,k)];
  Aik = &A[BlPos(k,k-1)];
  Ail = &A[BlPos(2*k-1,k-1)];
  pkn_MultMatrixTSubf ( s, r, r, Ail, s, r, Aik, s, Aij );
  pkn_MatrixUpperTrSolvef ( s, s, s, Aij, Aii, s, Aij );
  Aij = &A[BlPos(2*k,k)];
  Aik = &A[BlPos(k,0)];
  Ail = &A[BlPos(2*k,0)];
  pkn_MultMatrixTSubf ( t, r, r, Ail, s, r, Aik, s, Aij );
  Aik = &A[BlPos(k,k-1)];
  Ail = &A[BlPos(2*k,k-1)];
  pkn_MultMatrixTSubf ( t, r, r, Ail, s, r, Aik, s, Aij );
  pkn_MatrixUpperTrSolvef ( t, s, s, Aij, Aii, s, Aij );

        /* block columns from k+1 to 2k-3 */
  for ( i = k+1; i < 2*k-2; i++ ) {
    Aii = &A[BlPos(i,i)];
    Aij = &A[BlPos(i,i-k-1)];
    pkn_SymMatSubAATf ( s, Aii, r, r, Aij );
    Aij = &A[BlPos(i,i-k)];
    pkn_SymMatSubAATf ( s, Aii, r, r, Aij );
    Aij = &A[BlPos(i,i-1)];
    pkn_SymMatSubAATf ( s, Aii, s, s, Aij );
    if ( !pkn_CholeskyDecompf ( s, Aii ) )
      return false;
    Aij = &A[BlPos(i+1,i)];
    Aik = &A[BlPos(i,i-k)];
    Ail = &A[BlPos(i+1,i-k)];
    pkn_MultMatrixTSubf ( s, r, r, Ail, s, r, Aik, s, Aij );
    pkn_MatrixUpperTrSolvef ( s, s, s, Aij, Aii, s, Aij );
    Aij = &A[BlPos(2*k-1,i)];
    Aik = &A[BlPos(i,i-1)];
    Ail = &A[BlPos(2*k-1,i-1)];
    pkn_MultMatrixTSubf ( s, s, s, Ail, s, s, Aik, s, Aij );
    pkn_MatrixUpperTrSolvef ( s, s, s, Aij, Aii, s, Aij );
    Aij = &A[BlPos(2*k,i)];
    Aik = &A[BlPos(i,i-k-1)];
    Ail = &A[BlPos(2*k,i-k-1)];
    pkn_MultMatrixTSubf ( t, r, r, Ail, s, r, Aik, s, Aij );
    Aik = &A[BlPos(i,i-k)];
    Ail = &A[BlPos(2*k,i-k)];
    pkn_MultMatrixTSubf ( t, r, r, Ail, s, r, Aik, s, Aij );
    Aik = &A[BlPos(i,i-1)];
    Ail = &A[BlPos(2*k,i-1)];
    pkn_MultMatrixTSubf ( t, s, s, Ail, s, s, Aik, s, Aij );
    pkn_MatrixUpperTrSolvef ( t, s, s, Aij, Aii, s, Aij );
  }

        /* block column 2k-2 */
  Aii = &A[BlPos(2*k-2,2*k-2)];
  Aij = &A[BlPos(2*k-2,k-3)];
  pkn_SymMatSubAATf ( s, Aii, r, r, Aij );
  Aij = &A[BlPos(2*k-2,k-2)];
  pkn_SymMatSubAATf ( s, Aii, r, r, Aij );
  Aij = &A[BlPos(2*k-2,2*k-3)];
  pkn_SymMatSubAATf ( s, Aii, s, s, Aij );
  if ( !pkn_CholeskyDecompf ( s, Aii ) )
    return false;
  Aij = &A[BlPos(2*k-1,2*k-2)];
  Aik = &A[BlPos(2*k-2,k-2)];
  Ail = &A[BlPos(2*k-1,k-2)];
  pkn_MultMatrixTSubf ( s, r, r, Ail, s, r, Aik, s, Aij );
  Aik = &A[BlPos(2*k-2,2*k-3)];
  Ail = &A[BlPos(2*k-1,2*k-3)];
  pkn_MultMatrixTSubf ( s, s, s, Ail, s, s, Aik, s, Aij );
  pkn_MatrixUpperTrSolvef ( s, s, s, Aij, Aii, s, Aij );
  Aij = &A[BlPos(2*k,2*k-2)];
  Aik = &A[BlPos(2*k-2,k-3)];
  Ail = &A[BlPos(2*k,k-3)];
  pkn_MultMatrixTSubf ( t, r, r, Ail, s, r, Aik, s, Aij );
  Aik = &A[BlPos(2*k-2,k-2)];
  Ail = &A[BlPos(2*k,k-2)];
  pkn_MultMatrixTSubf ( t, r, r, Ail, s, r, Aik, s, Aij );
  Aik = &A[BlPos(2*k-2,2*k-3)];
  Ail = &A[BlPos(2*k,2*k-3)];
  pkn_MultMatrixTSubf ( t, s, s, Ail, s, s, Aik, s, Aij );
  pkn_MatrixUpperTrSolvef ( t, s, s, Aij, Aii, s, Aij );

        /* block column 2k-1 */
  Aii = &A[BlPos(2*k-1,2*k-1)];
  Aij = &A[BlPos(2*k-1,k-2)];
  pkn_SymMatSubAATf ( s, Aii, r, r, Aij );
  Aij = &A[BlPos(2*k-1,k-1)];
  pkn_SymMatSubAATf ( s, Aii, r, r, Aij );
  for ( i = k; i < 2*k-1; i++ ) {
    Aij = &A[BlPos(2*k-1,i)];
    pkn_SymMatSubAATf ( s, Aii, s, s, Aij );
  }
  if ( !pkn_CholeskyDecompf ( s, Aii ) )
    return false;
  Aij = &A[BlPos(2*k,2*k-1)];
  Aik = &A[BlPos(2*k-1,k-2)];
  Ail = &A[BlPos(2*k,k-2)];
  pkn_MultMatrixTSubf ( t, r, r, Ail, s, r, Aik, s, Aij );
  Aik = &A[BlPos(2*k-1,k-1)];
  Ail = &A[BlPos(2*k,k-1)];
  pkn_MultMatrixTSubf ( t, r, r, Ail, s, r, Aik, s, Aij );
  for ( i = k; i < 2*k-1; i++ ) {
    Aik = &A[BlPos(2*k-1,i)];
    Ail = &A[BlPos(2*k,i)];
    pkn_MultMatrixTSubf ( t, s, s, Ail, s, s, Aik, s, Aij );
  }
  pkn_MatrixUpperTrSolvef ( t, s, s, Aij, Aii, s, Aij );

        /* block 2k,2k */
  Aii = &A[BlPos(2*k,2*k)];
  for ( i = 0; i < k; i++ ) {
    Aij = &A[BlPos(2*k,i)];
    pkn_SymMatSubAATf ( t, Aii, r, r, Aij );
  }
  for ( ; i < 2*k; i++ ) {
    Aij = &A[BlPos(2*k,i)];
    pkn_SymMatSubAATf ( t, Aii, s, s, Aij );
  }
  return pkn_CholeskyDecompf ( t, Aii );
#undef BlPos
} /*pkn_Block2CholeskyDecompMf*/

void pkn_Block2LowerTrMSolvef ( int k, int r, int s, int t, CONST_ float *L,
                                int spdimen, int xpitch, float *x )
{
#define BlPos(i,j) pkn_Block2FindBlockPos( k, r, s, t, i, j )
  int   i;
  float *Lij, *xx, *xxx, *xxxx;

  for ( i = 0, xx = x;  i < k;  i++, xx += r*xpitch ) {
    Lij = &L[BlPos(i,i)];
    pkn_LowerTrMatrixSolvef ( r, Lij, spdimen, xpitch, xx, xpitch, xx );
  }
  Lij = &L[BlPos(k,0)];
  pkn_MultMatrixSubf ( s, r, r, Lij, spdimen, xpitch, x, xpitch, xx );
  Lij = &L[BlPos(k,k-1)];
  pkn_MultMatrixSubf ( s, r, r, Lij, spdimen, xpitch, &x[(k-1)*r*xpitch],
                       xpitch, xx );
  Lij = &L[BlPos(k,k)];
  pkn_LowerTrMatrixSolvef ( s, Lij, spdimen, xpitch, xx, xpitch, xx );
  for ( i = k+1, xxxx = xx, xx += s*xpitch, xxx = x;
        i < 2*k-1;
        i++, xxxx = xx, xx += s*xpitch ) {
    Lij = &L[BlPos(i,i-k-1)];
    pkn_MultMatrixSubf ( s, r, r, Lij, spdimen, xpitch, xxx, xpitch, xx );
    xxx += r*xpitch;
    Lij = &L[BlPos(i,i-k)];
    pkn_MultMatrixSubf ( s, r, r, Lij, spdimen, xpitch, xxx, xpitch, xx );
    Lij = &L[BlPos(i,i-1)];
    pkn_MultMatrixSubf ( s, s, s, Lij, spdimen, xpitch, xxxx, xpitch, xx );
    Lij = &L[BlPos(i,i)];
    pkn_LowerTrMatrixSolvef ( s, Lij, spdimen, xpitch, xx, xpitch, xx );
  }
  Lij = &L[BlPos(2*k-1,k-2)];
  pkn_MultMatrixSubf ( s, r, r, Lij, spdimen, xpitch, xxx, xpitch, xx );
  xxx += r*xpitch;
  Lij = &L[BlPos(2*k-1,k-1)];
  pkn_MultMatrixSubf ( s, r, r, Lij, spdimen, xpitch, xxx, xpitch, xx );
  for ( i = k, xxx += r*xpitch;  i < 2*k-1;  i++, xxx += s*xpitch ) {
    Lij = &L[BlPos(2*k-1,i)];
    pkn_MultMatrixSubf ( s, s, s, Lij, spdimen, xpitch, xxx, xpitch, xx );
  }
  Lij = &L[BlPos(2*k-1,2*k-1)];
  pkn_LowerTrMatrixSolvef ( s, Lij, spdimen, xpitch, xx, xpitch, xx );
  for ( i = 0, xx += s*xpitch, xxx = x;  i < k;  i++, xxx += r*xpitch ) {
    Lij = &L[BlPos(2*k,i)];
    pkn_MultMatrixSubf ( t, r, r, Lij, spdimen, xpitch, xxx, xpitch, xx );
  }
  for ( ; i < 2*k; i++, xxx += s*xpitch ) {
    Lij = &L[BlPos(2*k,i)];
    pkn_MultMatrixSubf ( t, s, s, Lij, spdimen, xpitch, xxx, xpitch, xx );
  }
  Lij = &L[BlPos(2*k,2*k)];
  pkn_LowerTrMatrixSolvef ( t, Lij, spdimen, xpitch, xx, xpitch, xx );
#undef BlPos
} /*pkn_Block2LowerTrMSolvef*/

void pkn_Block2UpperTrMSolvef ( int k, int r, int s, int t, CONST_ float *L,
                                int spdimen, int xpitch, float *x )
{
#define BlPos(i,j) pkn_Block2FindBlockPos( k, r, s, t, i, j )
  int   i;
  float *Lij, *xx, *xxx, *xxxx, *xxxxx;

  Lij = &L[BlPos(2*k,2*k)];
  xx = &x[k*(r+s)*xpitch];
  pkn_UpperTrMatrixSolvef ( t, Lij, spdimen, xpitch, xx, xpitch, xx );
  xxxxx = xx;
  xxxx = xx -= s*xpitch;
  Lij = &L[BlPos(2*k,2*k-1)];
  pkn_MultTMatrixSubf ( t, s, s, Lij, spdimen, xpitch, xxxxx, xpitch, xx );
  Lij = &L[BlPos(2*k-1,2*k-1)];
  pkn_UpperTrMatrixSolvef ( s, Lij, spdimen, xpitch, xx, xpitch, xx );
  xxx = xx -= s*xpitch;
  Lij = &L[BlPos(2*k-1,2*k-2)];
  pkn_MultTMatrixSubf ( s, s, s, Lij, spdimen, xpitch, xxxx, xpitch, xx );
  Lij = &L[BlPos(2*k,2*k-2)];
  pkn_MultTMatrixSubf ( t, s, s, Lij, spdimen, xpitch, xxxxx, xpitch, xx );
  Lij = &L[BlPos(2*k-2,2*k-2)];
  pkn_UpperTrMatrixSolvef ( s, Lij, spdimen, xpitch, xx, xpitch, xx );
  for ( xx -= s*xpitch, i = 2*k-3; i >= k; i--, xxx = xx, xx -= s*xpitch ) {
    Lij = &L[BlPos(i+1,i)];
    pkn_MultTMatrixSubf ( s, s, s, Lij, spdimen, xpitch, xxx, xpitch, xx );
    Lij = &L[BlPos(2*k-1,i)];
    pkn_MultTMatrixSubf ( s, s, s, Lij, spdimen, xpitch, xxxx, xpitch, xx );
    Lij = &L[BlPos(2*k,i)];
    pkn_MultTMatrixSubf ( t, s, s, Lij, spdimen, xpitch, xxxxx, xpitch, xx );
    Lij = &L[BlPos(i,i)];
    pkn_UpperTrMatrixSolvef ( s, Lij, spdimen, xpitch, xx, xpitch, xx );
  }
  xx = &x[(k-1)*r*xpitch];
  Lij = &L[BlPos(k,k-1)];
  pkn_MultTMatrixSubf ( s, r, r, Lij, spdimen, xpitch, xxx, xpitch, xx );
  Lij = &L[BlPos(2*k-1,k-1)];
  pkn_MultTMatrixSubf ( s, r, r, Lij, spdimen, xpitch, xxxx, xpitch, xx );
  Lij = &L[BlPos(2*k,k-1)];
  pkn_MultTMatrixSubf ( t, r, r, Lij, spdimen, xpitch, xxxxx, xpitch, xx );
  Lij = &L[BlPos(k-1,k-1)];
  pkn_UpperTrMatrixSolvef ( r, Lij, spdimen, xpitch, xx, xpitch, xx );
  for ( xx -= r*xpitch, xxx = xxxx-s*xpitch, i = k-2;
        i >= 0;
        i--, xx -= r*xpitch, xxxx = xxx, xxx -= s*xpitch ) {
    Lij = &L[BlPos(i+k,i)];
    pkn_MultTMatrixSubf ( s, r, r, Lij, spdimen, xpitch, xxx, xpitch, xx );
    Lij = &L[BlPos(i+k+1,i)];
    pkn_MultTMatrixSubf ( s, r, r, Lij, spdimen, xpitch, xxxx, xpitch, xx );
    Lij = &L[BlPos(2*k,i)];
    pkn_MultTMatrixSubf ( t, r, r, Lij, spdimen, xpitch, xxxxx, xpitch, xx );
    Lij = &L[BlPos(i,i)];
    pkn_UpperTrMatrixSolvef ( r, Lij, spdimen, xpitch, xx, xpitch, xx );
  }
#undef BlPos
} /*pkn_Block2UpperTrMSolvef*/

void pkn_Block2SymMatrixMultf ( int k, int r, int s, int t, CONST_ float *A,
                                int spdimen, int xpitch, CONST_ float *x,
                                int ypitch, float *y )
{
#define BlPos(i,j) pkn_Block2FindBlockPos( k, r, s, t, i, j )
  int   i;
  float *Aij, *xx, *yy;

  for ( i = 0; i < k-1; i++ ) {
    yy = &y[i*r*ypitch];
    Aij = &A[BlPos(i,i)];
    xx = &x[i*r*xpitch];
    pkn_SymMatrixMultf ( r, Aij, spdimen, xpitch, xx, ypitch, yy );
    Aij = &A[BlPos(i+k,i)];
    xx = &x[(k*r+i*s)*xpitch];
    pkn_MultTMatrixAddf ( s, r, r, Aij, spdimen, xpitch, xx, ypitch, yy );
    Aij = &A[BlPos(i+k+1,i)];
    xx += s*xpitch;
    pkn_MultTMatrixAddf ( s, r, r, Aij, spdimen, xpitch, xx, ypitch, yy );
    Aij = &A[BlPos(2*k,i)];
    xx = &x[k*(r+s)*xpitch];
    pkn_MultTMatrixAddf ( t, r, r, Aij, spdimen, xpitch, xx, ypitch, yy );
  }
  yy = &y[(k-1)*r*ypitch];
  Aij = &A[BlPos(k-1,k-1)];
  xx = &x[(k-1)*r*xpitch];
  pkn_SymMatrixMultf ( r, Aij, spdimen, xpitch, xx, ypitch, yy );
  Aij = &A[BlPos(k,k-1)];
  xx += r*xpitch;
  pkn_MultTMatrixAddf ( s, r, r, Aij, spdimen, xpitch, xx, ypitch, yy );
  Aij = &A[BlPos(2*k-1,k-1)];
  xx = &x[(k*r+(k-1)*s)*xpitch];
  pkn_MultTMatrixAddf ( s, r, r, Aij, spdimen, xpitch, xx, ypitch, yy );
  Aij = &A[BlPos(2*k,k-1)];
  xx = &x[k*(r+s)*xpitch];
  pkn_MultTMatrixAddf ( t, r, r, Aij, spdimen, xpitch, xx, ypitch, yy );

  yy = &y[k*r*ypitch];
  Aij = &A[BlPos(k,k)];
  xx = &x[k*r*xpitch];
  pkn_SymMatrixMultf ( s, Aij, spdimen, xpitch, xx, ypitch, yy );
  Aij = &A[BlPos(k,0)];
  pkn_MultMatrixAddf ( s, r, r, Aij, spdimen, xpitch, x, ypitch, yy );
  Aij = &A[BlPos(k,k-1)];
  xx = &x[(k-1)*r*xpitch];
  pkn_MultMatrixAddf ( s, r, r, Aij, spdimen, xpitch, xx, ypitch, yy );
  Aij = &A[BlPos(k+1,k)];
  xx = &x[(k*r+s)*xpitch];
  pkn_MultTMatrixAddf ( s, s, s, Aij, spdimen, xpitch, xx, ypitch, yy );
  Aij = &A[BlPos(2*k-1,k)];
  xx = &x[(k*r+(k-1)*s)*xpitch];
  pkn_MultTMatrixAddf ( s, s, s, Aij, spdimen, xpitch, xx, ypitch, yy );
  Aij = &A[BlPos(2*k,k)];
  xx = &x[k*(r+s)*xpitch];
  pkn_MultTMatrixAddf ( t, s, s, Aij, spdimen, xpitch, xx, ypitch, yy );
  for ( i = k+1; i < 2*k-1; i++ ) {
    yy = &y[(k*r+(i-k)*s)*ypitch];
    Aij = &A[BlPos(i,i)];
    xx = &x[(k*r+(i-k)*s)*xpitch];
    pkn_SymMatrixMultf ( s, Aij, spdimen, xpitch, xx, ypitch, yy );
    Aij = &A[BlPos(i,i-k-1)];
    xx = &x[(i-k-1)*r*xpitch];
    pkn_MultMatrixAddf ( s, r, r, Aij, spdimen, xpitch, xx, ypitch, yy );
    Aij = &A[BlPos(i,i-k)];
    xx = &x[(i-k)*r*xpitch];
    pkn_MultMatrixAddf ( s, r, r, Aij, spdimen, xpitch, xx, ypitch, yy );
    Aij = &A[BlPos(i,i-1)];
    xx = &x[(k*r+(i-k-1)*s)*xpitch];
    pkn_MultMatrixAddf ( s, s, s, Aij, spdimen, xpitch, xx, ypitch, yy );
    if ( i < 2*k-2 ) {
      Aij = &A[BlPos(i+1,i)];
      xx = &x[(k*r+(i-k+1)*s)*xpitch];
      pkn_MultTMatrixAddf ( s, s, s, Aij, spdimen, xpitch, xx, ypitch, yy );
    }
    Aij = &A[BlPos(2*k-1,i)];
    xx = &x[(k*r+(k-1)*s)*xpitch];
    pkn_MultTMatrixAddf ( s, s, s, Aij, spdimen, xpitch, xx, ypitch, yy );
    Aij = &A[BlPos(2*k,i)];
    xx = &x[k*(r+s)*xpitch];
    pkn_MultTMatrixAddf ( t, s, s, Aij, spdimen, xpitch, xx, ypitch, yy );
  }
  yy = &y[(k*r+(k-1)*s)*ypitch];
  Aij = &A[BlPos(2*k-1,2*k-1)];
  xx = &x[(k*r+(k-1)*s)*xpitch];
  pkn_SymMatrixMultf ( s, Aij, spdimen, xpitch, xx, ypitch, yy );
  Aij = &A[BlPos(2*k-1,k-2)];
  xx = &x[(k-2)*r*xpitch];
  pkn_MultMatrixAddf ( s, r, r, Aij, spdimen, xpitch, xx, ypitch, yy );
  Aij = &A[BlPos(2*k-1,k-1)];
  xx = &x[(k-1)*r*xpitch];
  pkn_MultMatrixAddf ( s, r, r, Aij, spdimen, xpitch, xx, ypitch, yy );
  for ( i = k; i < 2*k-1; i++ ) {
    Aij = &A[BlPos(2*k-1,i)];
    xx = &x[(k*r+(i-k)*s)*xpitch];
    pkn_MultMatrixAddf ( s, s, s, Aij, spdimen, xpitch, xx, ypitch, yy );
  }
  Aij = &A[BlPos(2*k,2*k-1)];
  xx = &x[k*(r+s)*xpitch];
  pkn_MultTMatrixAddf ( t, s, s, Aij, spdimen, xpitch, xx, ypitch, yy );

  yy = &y[k*(r+s)*ypitch];
  Aij = &A[BlPos(2*k,2*k)];
  xx = &x[k*(r+s)*xpitch];
  pkn_SymMatrixMultf ( t, Aij, spdimen, xpitch, xx, ypitch, yy );
  for ( i = 0; i < k; i++ ) {
    Aij = &A[BlPos(2*k,i)];
    xx = &x[i*r*xpitch];
    pkn_MultMatrixAddf ( t, r, r, Aij, spdimen, xpitch, xx, ypitch, yy );
  }
  for ( ; i < 2*k; i++ ) {
    Aij = &A[BlPos(2*k,i)];
    xx = &x[(k*r+(i-k)*s)*xpitch];
    pkn_MultMatrixAddf ( t, s, s, Aij, spdimen, xpitch, xx, ypitch, yy );
  }
#undef BlPos
} /*pkn_Block2SymMatrixMultf*/

