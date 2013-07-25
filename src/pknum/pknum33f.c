
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

/* ///////////////////////////////////////////////////////////// */
/* processing matrices with the Block3 block structure           */

boolean pkn_Block3CholeskyDecompMf ( int k, int r, int s, float *A )
{
#define BlPos(i,j) pkn_Block3FindBlockPos( k, r, s, i, j )
  int   i;
  float *Aii, *Aij, *Aji, *Aki;

  for ( i = 0; i < k-1; i++ ) {
    Aii = &A[BlPos( i, i )];
    if ( i > 0 ) {
      Aij = &A[BlPos( i, i-1)];
      pkn_SymMatSubAATf ( r, Aii, r, r, Aij );
    }
    if ( !pkn_CholeskyDecompf ( r, Aii ) )
      return false;
    Aji = &A[BlPos( i+1, i )];
    pkn_MatrixUpperTrSolvef ( r, r, r, Aji, Aii, r, Aji );
    Aki = &A[BlPos( k, i )];
    if ( i > 0 ) {
      Aji = &A[BlPos( k, i-1 )];
      Aij = &A[BlPos( i, i-1 )];
      pkn_MultMatrixTSubf ( s, r, r, Aji, r, r, Aij, r, Aki );
    }
    pkn_MatrixUpperTrSolvef ( s, r, r, Aki, Aii, r, Aki );
  }
  Aii = &A[BlPos( k-1, k-1 )];
  Aij = &A[BlPos( k-1, k-2 )];
  pkn_SymMatSubAATf ( r, Aii, r, r, Aij );
  if ( !pkn_CholeskyDecompf ( r, Aii ) )
    return false;
  Aki = &A[BlPos( k, k-1 )];
  Aji = &A[BlPos( k, k-2 )];
  Aij = &A[BlPos( k-1, k-2 )];
  pkn_MultMatrixTSubf ( s, r, r, Aji, r, r, Aij, r, Aki );
  pkn_MatrixUpperTrSolvef ( s, r, r, Aki, Aii, r, Aki );
  Aii = &A[BlPos( k, k )];
  for ( i = 0; i < k; i++ ) {
    Aij = &A[BlPos( k, i )];
    pkn_SymMatSubAATf ( s, Aii, r, r, Aij );
  }
  return pkn_CholeskyDecompf ( s, Aii );
#undef BlPos
} /*pkn_Block3CholeskyDecompMf*/

void pkn_Block3LowerTrMSolvef ( int k, int r, int s, CONST_ float *L,
                                int spdimen, int xpitch, float *x )
{
#define BlPos(i,j) pkn_Block3FindBlockPos( k, r, s, i, j )
  int   i;
  float *Lij, *xx, *xxx;

  Lij = &L[BlPos(0,0)];
  pkn_LowerTrMatrixSolvef ( r, Lij, spdimen, xpitch, x, xpitch, x );
  for ( i = 1, xxx = x, xx = xxx + r*xpitch;
        i < k;
        i++, xxx = xx, xx += r*xpitch ) {
    Lij = &L[BlPos(i,i-1)];
    pkn_MultMatrixSubf ( r, r, r, Lij, spdimen, xpitch, xxx, xpitch, xx );
    Lij = &L[BlPos(i,i)];
    pkn_LowerTrMatrixSolvef ( r, Lij, spdimen, xpitch, xx, xpitch, xx );
  }
  for ( i = 0, xxx = x;  i < k;  i++, xxx += r*xpitch ) {
    Lij = &L[BlPos(k,i)];
    pkn_MultMatrixSubf ( s, r, r, Lij, spdimen, xpitch, xxx, xpitch, xx );
  }
  Lij = &L[BlPos(k,k)];
  pkn_LowerTrMatrixSolvef ( s, Lij, spdimen, xpitch, xx, xpitch, xx );
#undef BlPos
} /*pkn_Block3LowerTrMSolvef*/

void pkn_Block3UpperTrMSolvef ( int k, int r, int s, CONST_ float *L,
                                int spdimen, int xpitch, float *x )
{
#define BlPos(i,j) pkn_Block3FindBlockPos( k, r, s, i, j )
  int   i;
  float *Lij, *xx, *xxx;

  Lij = &L[BlPos(k,k)];
  x += k*r*xpitch;
  pkn_UpperTrMatrixSolvef ( s, Lij, spdimen, xpitch, x, xpitch, x );
  Lij = &L[BlPos(k,k-1)];
  xxx = x;  xx = xxx-r*xpitch;
  pkn_MultTMatrixSubf ( s, r, r, Lij, spdimen, xpitch, x, xpitch, xx );
  Lij = &L[BlPos(k-1,k-1)];
  pkn_UpperTrMatrixSolvef ( r, Lij, spdimen, xpitch, xx, xpitch, xx );
  for ( i = k-2, xxx = xx, xx -= r*xpitch;
        i >= 0;
        i--, xxx = xx, xx -= r*xpitch ) {
    Lij = &L[BlPos(k,i)];
    pkn_MultTMatrixSubf ( s, r, r, Lij, spdimen, xpitch, x, xpitch, xx );
    Lij = &L[BlPos(i+1,i)];
    pkn_MultTMatrixSubf ( r, r, r, Lij, spdimen, xpitch, xxx, xpitch, xx );
    Lij = &L[BlPos(i,i)];
    pkn_UpperTrMatrixSolvef ( r, Lij, spdimen, xpitch, xx, xpitch, xx );
  }
#undef BlPos
} /*pkn_Block3UpperTrMSolvef*/

void pkn_Block3SymMatrixMultf ( int k, int r, int s, CONST_ float *A,
                                int spdimen, int xpitch, CONST_ float *x,
                                int ypitch, float *y )
{
#define BlPos(i,j) pkn_Block3FindBlockPos( k, r, s, i, j )
  int   i;
  float *xx, *yy, *Aij;

  for ( i = 0; i < k; i++ ) {
    yy = y + i*r*ypitch;
    xx = x + i*r*xpitch;
    Aij = &A[BlPos(i,i)];
    pkn_SymMatrixMultf ( r, Aij, spdimen, xpitch, xx, ypitch, yy );
    if ( i > 0 ) {
      Aij = &A[BlPos(i,i-1)];
      pkn_MultMatrixAddf ( r, r, r, Aij,
                           spdimen, xpitch, xx-r*xpitch, ypitch, yy );
    }
    if ( i < k-1 ) {
      Aij = &A[BlPos(i+1,i)];
      pkn_MultTMatrixAddf ( r, r, r, Aij,
                            spdimen, xpitch, xx+r*xpitch, ypitch, yy );
    }
    Aij = &A[BlPos(k,i)];
    pkn_MultTMatrixAddf ( s, r, r, Aij,
                          spdimen, xpitch, x+k*r*xpitch, ypitch, yy );
  }
  Aij = &A[BlPos(k,k)];
  yy = y+k*r*ypitch;
  pkn_SymMatrixMultf ( s, Aij, spdimen, xpitch, x+k*r*xpitch, ypitch, yy );
  for ( i = 0, xx = x;  i < k;  i++, xx += r*xpitch ) {
    Aij = &A[BlPos(k,i)];
    pkn_MultMatrixAddf ( s, r, r, Aij, spdimen, xpitch, xx, ypitch, yy );
  }
#undef BlPos
} /*pkn_Block3SymMatrixMultf*/

