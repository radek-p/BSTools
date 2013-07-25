
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

boolean pkn_Block3CholeskyDecompMd ( int k, int r, int s, double *A )
{
#define BlPos(i,j) pkn_Block3FindBlockPos( k, r, s, i, j )
  int    i;
  double *Aii, *Aij, *Aji, *Aki;

  for ( i = 0; i < k-1; i++ ) {
    Aii = &A[BlPos( i, i )];
    if ( i > 0 ) {
      Aij = &A[BlPos( i, i-1)];
      pkn_SymMatSubAATd ( r, Aii, r, r, Aij );
    }
    if ( !pkn_CholeskyDecompd ( r, Aii ) )
      return false;
    Aji = &A[BlPos( i+1, i )];
    pkn_MatrixUpperTrSolved ( r, r, r, Aji, Aii, r, Aji );
    Aki = &A[BlPos( k, i )];
    if ( i > 0 ) {
      Aji = &A[BlPos( k, i-1 )];
      Aij = &A[BlPos( i, i-1 )];
      pkn_MultMatrixTSubd ( s, r, r, Aji, r, r, Aij, r, Aki );
    }
    pkn_MatrixUpperTrSolved ( s, r, r, Aki, Aii, r, Aki );
  }
  Aii = &A[BlPos( k-1, k-1 )];
  Aij = &A[BlPos( k-1, k-2 )];
  pkn_SymMatSubAATd ( r, Aii, r, r, Aij );
  if ( !pkn_CholeskyDecompd ( r, Aii ) )
    return false;
  Aki = &A[BlPos( k, k-1 )];
  Aji = &A[BlPos( k, k-2 )];
  Aij = &A[BlPos( k-1, k-2 )];
  pkn_MultMatrixTSubd ( s, r, r, Aji, r, r, Aij, r, Aki );
  pkn_MatrixUpperTrSolved ( s, r, r, Aki, Aii, r, Aki );
  Aii = &A[BlPos( k, k )];
  for ( i = 0; i < k; i++ ) {
    Aij = &A[BlPos( k, i )];
    pkn_SymMatSubAATd ( s, Aii, r, r, Aij );
  }
  return pkn_CholeskyDecompd ( s, Aii );
#undef BlPos
} /*pkn_Block3CholeskyDecompMd*/

void pkn_Block3LowerTrMSolved ( int k, int r, int s, CONST_ double *L,
                                int spdimen, int xpitch, double *x )
{
#define BlPos(i,j) pkn_Block3FindBlockPos( k, r, s, i, j )
  int    i;
  double *Lij, *xx, *xxx;

  Lij = &L[BlPos(0,0)];
  pkn_LowerTrMatrixSolved ( r, Lij, spdimen, xpitch, x, xpitch, x );
  for ( i = 1, xxx = x, xx = xxx + r*xpitch;
        i < k;
        i++, xxx = xx, xx += r*xpitch ) {
    Lij = &L[BlPos(i,i-1)];
    pkn_MultMatrixSubd ( r, r, r, Lij, spdimen, xpitch, xxx, xpitch, xx );
    Lij = &L[BlPos(i,i)];
    pkn_LowerTrMatrixSolved ( r, Lij, spdimen, xpitch, xx, xpitch, xx );
  }
  for ( i = 0, xxx = x;  i < k;  i++, xxx += r*xpitch ) {
    Lij = &L[BlPos(k,i)];
    pkn_MultMatrixSubd ( s, r, r, Lij, spdimen, xpitch, xxx, xpitch, xx );
  }
  Lij = &L[BlPos(k,k)];
  pkn_LowerTrMatrixSolved ( s, Lij, spdimen, xpitch, xx, xpitch, xx );
#undef BlPos
} /*pkn_Block3LowerTrMSolved*/

void pkn_Block3UpperTrMSolved ( int k, int r, int s, CONST_ double *L,
                                int spdimen, int xpitch, double *x )
{
#define BlPos(i,j) pkn_Block3FindBlockPos( k, r, s, i, j )
  int    i;
  double *Lij, *xx, *xxx;

  Lij = &L[BlPos(k,k)];
  x += k*r*xpitch;
  pkn_UpperTrMatrixSolved ( s, Lij, spdimen, xpitch, x, xpitch, x );
  Lij = &L[BlPos(k,k-1)];
  xxx = x;  xx = xxx-r*xpitch;
  pkn_MultTMatrixSubd ( s, r, r, Lij, spdimen, xpitch, x, xpitch, xx );
  Lij = &L[BlPos(k-1,k-1)];
  pkn_UpperTrMatrixSolved ( r, Lij, spdimen, xpitch, xx, xpitch, xx );
  for ( i = k-2, xxx = xx, xx -= r*xpitch;
        i >= 0;
        i--, xxx = xx, xx -= r*xpitch ) {
    Lij = &L[BlPos(k,i)];
    pkn_MultTMatrixSubd ( s, r, r, Lij, spdimen, xpitch, x, xpitch, xx );
    Lij = &L[BlPos(i+1,i)];
    pkn_MultTMatrixSubd ( r, r, r, Lij, spdimen, xpitch, xxx, xpitch, xx );
    Lij = &L[BlPos(i,i)];
    pkn_UpperTrMatrixSolved ( r, Lij, spdimen, xpitch, xx, xpitch, xx );
  }
#undef BlPos
} /*pkn_Block3UpperTrMSolved*/

void pkn_Block3SymMatrixMultd ( int k, int r, int s, CONST_ double *A,
                                int spdimen, int xpitch, CONST_ double *x,
                                int ypitch, double *y )
{
#define BlPos(i,j) pkn_Block3FindBlockPos( k, r, s, i, j )
  int    i;
  double *xx, *yy, *Aij;

  for ( i = 0; i < k; i++ ) {
    yy = y + i*r*ypitch;
    xx = x + i*r*xpitch;
    Aij = &A[BlPos(i,i)];
    pkn_SymMatrixMultd ( r, Aij, spdimen, xpitch, xx, ypitch, yy );
    if ( i > 0 ) {
      Aij = &A[BlPos(i,i-1)];
      pkn_MultMatrixAddd ( r, r, r, Aij,
                           spdimen, xpitch, xx-r*xpitch, ypitch, yy );
    }
    if ( i < k-1 ) {
      Aij = &A[BlPos(i+1,i)];
      pkn_MultTMatrixAddd ( r, r, r, Aij,
                            spdimen, xpitch, xx+r*xpitch, ypitch, yy );
    }
    Aij = &A[BlPos(k,i)];
    pkn_MultTMatrixAddd ( s, r, r, Aij,
                          spdimen, xpitch, x+k*r*xpitch, ypitch, yy );
  }
  Aij = &A[BlPos(k,k)];
  yy = y+k*r*ypitch;
  pkn_SymMatrixMultd ( s, Aij, spdimen, xpitch, x+k*r*xpitch, ypitch, yy );
  for ( i = 0, xx = x;  i < k;  i++, xx += r*xpitch ) {
    Aij = &A[BlPos(k,i)];
    pkn_MultMatrixAddd ( s, r, r, Aij, spdimen, xpitch, xx, ypitch, yy );
  }
#undef BlPos
} /*pkn_Block3SymMatrixMultd*/

