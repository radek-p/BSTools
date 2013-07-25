
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
/* Symmetric and triangular matrices with Block2 block structure */
/* computing array size and indexes */

int pkn_Block2ArraySize ( int k, int r, int s, int t )
{
  if ( k < 3 || r < 1 || s < 1 || t < 1 )
    return -1;
  else
    return k*(((r+1)*r+(s+1)*s)/2+(r+s)*(t+2*s))+((t+1)*t)/2-3*s*s;
} /*pkn_Block2ArraySize*/

int pkn_Block2FindBlockPos ( int k, int r, int s, int t, int i, int j )
{
  int n;

  if ( i < 0 || j < 0 || i > 2*k || j > 2*k || i < j )
    return -1;
  if ( i == j ) {  /* a diagonal block */
    if ( i < k )
      return (i*(r+1)*r)/2;
    else if ( i < 2*k )
      return (k*(r+1)*r+(i-k)*(s+1)*s)/2;
    else return (k*((r+1)*r+(s+1)*s))/2;
  }
        /* a subdiagonal block */
  n = (k*((r+1)*r+(s+1)*s)+(t+1)*t)/2;
  if ( i < k )
    return -1;
  else if ( i < 2*k ) {
    if ( j < k ) {
      if ( i == k ) {
        if ( j == k-1 ) return n;
        else if ( j == 0 ) return n+r*s;
        else return -1;
      }
      else {
        if ( j == i-k-1 ) return n + 2*(i-k)*r*s;
        else if ( j == i-k ) return n + (2*(i-k)+1)*r*s;
        else return -1;
      }
    }
    else {
      n += 2*k*r*s;
      if ( i < 2*k-1 ) {
        if ( j == i-1 )
          return n + (i-k-1)*s*s;
          else return -1;
      }
        /* here i == 2*k-1 */
      return n + (j-2)*s*s;
    }
  }
        /* here i == 2*k */
  n += (2*k*r+(2*k-3)*s)*s;
  if ( j < k )
    return n + j*r*t;
  else
    return n + (k*r+(j-k)*s)*t;
} /*pkn_Block2FindBlockPos*/

int pkn_Block2FindElemPos ( int k, int r, int s, int t, int i, int j )
{
  int br, bc, bb, n;

      /* are the indices within the range? */
  if ( i < 0 || i >= k*(r+s)+t || j < 0 || j >= k*(r+s)+t )
    return -1;
      /* a symmetric matrix, with the lower triangle stored, */
      /* swap the indexes if i < j */
  if ( i < j ) { br = i;  i = j;  j = br; }
      /* compute the block indices */
  if ( i < k*r ) { br = i/r;  i -= br*r; }
  else {
    i -= k*r;
    if ( i < k*s ) { br = i/s;  i -= br*s;  br += k; }
    else { br = 2*k;  i -= k*s; }
  }
  if ( j < k*r ) { bc = j/r;  j -= bc*r; }
  else {
    j -= k*r;
    if ( j < k*s ) { bc = j/s;  j -= bc*s;  bc += k; }
    else { bc = 2*k;  j -= k*s; }
  }
  bb = pkn_Block2FindBlockPos ( k, r, s, t, br, bc );
  if ( br == bc )
    return bb + pkn_SymMatIndex ( i, j );
  else if ( bb >= 0 ) {
    if ( bc < k ) n = r;
      else n = s;
    return bb + i*n+j;
  }
  else
    return -1;
} /*pkn_Block2FindElemPos*/

