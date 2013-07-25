
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2005, 2008                            */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "pkvaria.h"

#ifdef CONST_
#undef CONST_
#endif
#define CONST_

#include "pknum.h"

/* ///////////////////////////////////////////////////////////////////////// */
/* The procedure pkn_Block1CholeskyDecompMd computes the Cholesky            */
/* decomposition of a symmetric positive-definite matrix A, having the       */
/* following block structure:                                                */
/* A00                Ak0^T */
/*     A11            Ak1^T */
/*         A22        Ak2^T */
/*             ...    ...   */
/* Ak0 Ak1 Ak2 ...    Akk,  */
/* with the symmetric blocks Aii (i=0...k-1) of dimensions r x r, the        */
/* symmetric block Akk of dimensions s x s, and full blocks Aki (i=0...k-1)  */
/* of dimensions s x r. Only the lower block of each symmetric block is      */
/* stored in memory, row by row in r(r+1)/2 subsequent entries of the array  */
/* Aii (for i < k) or in the array Akk of length s(s+1)/2. Each full block   */
/* Aki is stored rowwise in s x r subsequent entries of the array Aki        */
/* The result of decomposition is the lower triangular matrix L, which under */
/* the diagonal has the same block structure that the matrix A.              */
/* The coefficients of L are stored in the arrays Aii, Aki and Akk in the    */
/* same way that the coefficients of the matrix A.                           */
/* ///////////////////////////////////////////////////////////////////////// */

boolean pkn_Block1CholeskyDecompMd ( int k, int r, int s, double *A )
{
  int    b, i, j, t;
  int    br, lbb, lkb, li0, lj0;
  double *Akk, *Aki;
  long double ajj, aij;

  Akk = &A[pkn_Block1FindBlockPos ( k, r, s, k, k )];
  Aki = &A[pkn_Block1FindBlockPos ( k, r, s, k, 0 )];
        /* the decomposition of the diagonal blocks */
  for ( b = br = 0;  b < k;  b++, br += r ) {
    lbb = b*r*(r+1)/2;  /* index of the first element of the block Abb */
                        /* in the array A */
    lkb = b*s*r;        /* index of the first element of the block Akb */
                        /* in the array Akb */
    for ( j = 0; j < r; j++ ) {      /* this concerns the column b*r+j */
      lj0 = lbb + pkn_SymMatIndex ( j, 0 );
      ajj = A[lj0+j];
      for ( t = 0; t < j; t++ )
        ajj -= A[lj0+t]*A[lj0+t];
      if ( ajj <= 0.0 )
        return false;   /* the matrix is not positive-definite */
      A[lj0+j] = (double)(ajj = sqrt ( (double)ajj ));
      for ( i = j+1; i < r; i++ ) {      /* elimination in Abb */
        li0 = lbb + pkn_SymMatIndex ( i, 0 );
        aij = A[li0+j];
        for ( t = 0; t < j; t++ )
          aij -= A[li0+t]*A[lj0+t];
        A[li0+j] = (double)(aij/ajj);
      }
      for ( i = 0; i < s; i++ ) {        /* elimination in Akb */
        li0 = lkb + i*r;
        aij = Aki[li0+j];
        for ( t = 0; t < j; t++ )
          aij -= Aki[li0+t]*A[lj0+t];
        Aki[li0+j] = (double)(aij/ajj);
      }
    }
  }
        /* the decomposition of the last diagonal block, Akk */
  for ( j = 0; j < s; j++ ) {
    lj0 = pkn_SymMatIndex ( j, 0 );
    ajj = Akk[lj0+j];
    for ( b = 0; b < k; b++ ) {
      lkb = b*s*r;
      li0 = lkb + r*j;
      for ( t = 0; t < r; t++ )
        ajj -= Aki[li0+t]*Aki[li0+t]; 
    }
    for ( t = 0; t < j; t++ )
      ajj -= Akk[lj0+t]*Akk[lj0+t];
    if ( ajj <= 0.0 )
      return false; /* the matrix is not positive-definite */
    Akk[lj0+j] = (double)(ajj = sqrt ( (double)ajj ));
    for ( i = j+1; i < s; i++ ) {
      li0 = pkn_SymMatIndex ( i, 0 );
      aij = Akk[li0+j];
      for ( b = 0; b < k; b++ ) {
        lkb = b*s*r;
        lj0 = lkb + r*j;
        li0 = lkb + r*i;
        for ( t = 0; t < r; t++ )
          aij -= Aki[li0+t]*Aki[lj0+t];
      }
      lj0 = pkn_SymMatIndex ( j, 0 );
      li0 = pkn_SymMatIndex ( i, 0 );
      for ( t = 0; t < j; t++ )
        aij -= Akk[li0+t]*Akk[lj0+t];
      Akk[li0+j] = (double)(aij/ajj);
    }
  }
  return true;
} /*pkn_Block1CholeskyDecompMd*/

/* ///////////////////////////////////////////////////////////////////////// */
void pkn_Block1LowerTrMSolved ( int k, int r, int s, CONST_ double *A,
                                int spdimen, int xpitch, double *x )
{
  int    b, i, j, lbb, lkb, li0, lj0;
  double *Akk, *Aki;

  Akk = &A[pkn_Block1FindBlockPos ( k, r, s, k, k )];
  Aki = &A[pkn_Block1FindBlockPos ( k, r, s, k, 0 )];
  for ( b = 0; b < k; b++ ) {
    lbb = b*r*(r+1)/2;  /* index of the first element of the block Abb */
                        /* in the array A */
    lkb = b*s*r;        /* index of the first element of the block Akb */
                        /* in the array Akb */
    for ( j = 0; j < r; j++ ) {      /* this concerns the column b*r+j */
      lj0 = lbb + pkn_SymMatIndex ( j, 0 );
      pkn_MultMatrixNumd ( 1, spdimen, 0, &x[(b*r+j)*xpitch],
                       1.0/A[lj0+j], 0, &x[(b*r+j)*xpitch] );
      for ( i = j+1; i < r; i++ ) {
        li0 = lbb + pkn_SymMatIndex ( i, 0 );
        pkn_AddMatrixMd ( 1, spdimen, 0, &x[(b*r+i)*xpitch],
                          0, &x[(b*r+j)*xpitch], -A[li0+j],
                          0, &x[(b*r+i)*xpitch] );
      }
      for ( i = 0; i < s; i++ ) {
        li0 = lkb + r*i;
        pkn_AddMatrixMd ( 1, spdimen, 0, &x[(k*r+i)*xpitch],
                          0, &x[(b*r+j)*xpitch], -Aki[li0+j],
                          0, &x[(k*r+i)*xpitch] );
      }
    }
  }
  for ( j = 0; j < s; j++ ) {
    lj0 = pkn_SymMatIndex ( j, 0 );
    pkn_MultMatrixNumd ( 1, spdimen, 0, &x[(k*r+j)*xpitch], 1.0/Akk[lj0+j],
                         0, &x[(k*r+j)*xpitch] );
    for ( i = j+1; i < s; i++ ) {
      li0 = pkn_SymMatIndex ( i, 0 );
      pkn_AddMatrixMd ( 1, spdimen, 0, &x[(k*r+i)*xpitch],
                        0, &x[(k*r+j)*xpitch], -Akk[li0+j],
                        0, &x[(k*r+i)*xpitch] );
    }
  }
} /*pkn_Block1LowerTrMSolved*/

void pkn_Block1UpperTrMSolved ( int k, int r, int s, CONST_ double *A,
                                int spdimen, int xpitch, double *x )
{
  int    b, i, j, lbb, lkb, lj0;
  double *Akk, *Aki;

  Akk = &A[pkn_Block1FindBlockPos ( k, r, s, k, k )];
  Aki = &A[pkn_Block1FindBlockPos ( k, r, s, k, 0 )];
  for ( j = s-1; j >= 0; j-- ) {
    lj0 = pkn_SymMatIndex ( j, 0 );
    pkn_MultMatrixNumd ( 1, spdimen, 0, &x[(k*r+j)*xpitch],
                         1.0/Akk[lj0+j], 0, &x[(k*r+j)*xpitch] );
    for ( i = j-1; i >= 0; i-- )
      pkn_AddMatrixMd ( 1, spdimen, 0, &x[(k*r+i)*xpitch],
                        0, &x[(k*r+j)*xpitch], -Akk[lj0+i],
                        0, &x[(k*r+i)*xpitch] );
    for ( b = k-1; b >= 0; b-- ) {
      lkb = b*s*r;
      for ( i = r-1; i >= 0; i-- )
        pkn_AddMatrixMd ( 1, spdimen, 0, &x[(b*r+i)*xpitch],
                          0, &x[(k*r+j)*xpitch], -Aki[lkb+j*r+i], /** ???? **/
                          0, &x[(b*r+i)*xpitch] );
    }
  }
  for ( b = k-1; b >= 0; b-- ) {
    lbb = b*r*(r+1)/2;  /* index of the first element of the block Abb */
                        /* in the array Aii */
    for ( j = r-1; j >= 0; j-- ) {
      lj0 = lbb + pkn_SymMatIndex ( j, 0 );
      pkn_MultMatrixNumd ( 1, spdimen, 0, &x[(b*r+j)*xpitch], 1.0/A[lj0+j],
                           0, &x[(b*r+j)*xpitch] );
      for ( i = j-1; i >= 0; i-- ) {
        pkn_AddMatrixMd ( 1, spdimen, 0, &x[(b*r+i)*xpitch],
                          0, &x[(b*r+j)*xpitch], -A[lj0+i],
                          0, &x[(b*r+i)*xpitch] );
      }
    }
  }
} /*pkn_Block1UpperTrMSolved*/

void pkn_Block1SymMatrixMultd ( int k, int r, int s, CONST_ double *A,
                                int spdimen, int xpitch, double *x,
                                int ypitch, double *y )
{
  void   *sp; 
  int    i, t;
  double *z, *Aki, *Akk;

  sp = pkv_GetScratchMemTop ();
  z = pkv_GetScratchMemd ( max(r,s)*spdimen );
  if ( !z )
    exit ( 1 ); 
  Akk = &A[pkn_Block1FindBlockPos ( k, r, s, k, k )];
  Aki = &A[pkn_Block1FindBlockPos ( k, r, s, k, 0 )];
  t = r*(r+1)/2;
  for ( i = 0; i < k; i++ )
    pkn_SymMatrixMultd ( r, &A[i*t], spdimen,
                         xpitch, &x[i*r*xpitch], ypitch, &y[i*r*ypitch] );
  pkn_SymMatrixMultd ( s, Akk, spdimen,
                       xpitch, &x[k*r*xpitch], ypitch, &y[k*r*ypitch] );
  for ( i = 0; i < k; i++ ) {
    pkn_MultMatrixd ( s, r, r, &Aki[i*r*s], spdimen,
                      xpitch, &x[i*r*xpitch], spdimen, z );
    pkn_AddMatrixd ( s, spdimen, spdimen, z,
                     ypitch, &y[k*r*ypitch], ypitch, &y[k*r*ypitch] );
    pkn_MultTMatrixd ( s, r, r, &Aki[i*r*s], spdimen,
                       xpitch, &x[k*r*xpitch], spdimen, z );
    pkn_AddMatrixd ( r, spdimen, spdimen, z,
                     ypitch, &y[i*r*ypitch], ypitch, &y[i*r*ypitch] );
  }

  pkv_SetScratchMemTop ( sp );
} /*pkn_Block1SymMatrixMultd*/

