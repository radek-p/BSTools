
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
#include "pknum.h"

#include "msgpool.h"

/* ///////////////////////////////////////////////////////////// */
/* Symmetric and triangular matrices with Block1 block structure */
/* computing array size and indexes */

int pkn_Block1ArraySize ( int k, int r, int s )
{
  if ( k < 2 || r < 1 || s < 1 )
    return -1;
  else
    return k*((r*(r+1))/2+r*s)+(s*(s+1))/2;
} /*pkn_Block1ArraySize*/

int pkn_Block1FindBlockPos ( int k, int r, int s, int i, int j )
{
  if ( i < 0 || j < 0 || i > k || j > k || i < j )
    return -1;
  if ( i == j )       /* a diagonal block */
    return i*((r*(r+1))/2);
  else if ( i == k )  /* a block from the lowest row */
    return k*((r*(r+1))/2) + (s*(s+1))/2 + j*r*s;
  else
    return -1;
} /*pkn_Block1FindBlockPos*/

int pkn_Block1FindElemPos ( int k, int r, int s, int i, int j )
{
  int bi, bj, n;

      /* are the indices within the range? */
  if ( i < 0 || i >= k*r+s || j < 0 || j >= k*r+s )
    return -1;
      /* a symmetric matrix, with the lower triangle stored, */
      /* swap the indexes if i < j */
  if ( i < j ) { bi = i;  i = j;  j = bi; }
      /* compute the block indices */
  bi = i / r;  if ( bi > k ) bi = k;
  bj = j / r;  if ( bj > k ) bj = k;
  n = pkn_Block1FindBlockPos ( k, r, s, bi, bj );
  if ( n >= 0 ) {
    if ( bi == bj )
      return n + pkn_SymMatIndex ( i-bi*r, j-bj*r );
    else  /* there must be i >= k*r */
      return n + (i-k*r)*r + j-bj*r;
  }
  else
    return -1;
} /*pkn_Block1FindElemPos*/

