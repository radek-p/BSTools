
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
#include "pknum.h"

#include "msgpool.h"

/* ///////////////////////////////////////////////////////////// */
/* Symmetric and triangular matrices with Block3 block structure */
/* computing array size and indexes */

int pkn_Block3ArraySize ( int k, int r, int s )
{
  return k*((r*(r+1))/2 + r*s) + (k-1)*r*r + (s*(s+1))/2;
} /*pkn_Block3ArraySize*/

int pkn_Block3FindBlockPos ( int k, int r, int s, int i, int j )
{
  int s1;

  if ( i >= 0 && i <= k && j >= 0 && j <= i ) {
    s1 = (r*(r+1))/2;
    if ( i == j )
      return i*s1;
    s1 = k*s1 + (s*(s+1))/2;
    if ( i < k && i == j+1 )
      return s1 + j*r*r;
    else if ( i == k )
      return s1 + ((k-1)*r + j*s)*r;
  }
  return -1;
} /*pkn_Block3FindBlockPos*/

int pkn_Block3FindElemPos ( int k, int r, int s, int i, int j )
{
  int ib, jb, ie, je, bp;

  if ( i < j ) { ib = i;  i = j;  j = ib; }
  if ( i >= 0 && i < k*r+s && j >= 0 && j <= i ) {
    ib = i / r;  if ( ib > k ) ib = k;
    jb = j / r;  if ( jb > k ) jb = k;
    bp = pkn_Block3FindBlockPos ( k, r, s, ib, jb );
    if ( bp >= 0 ) {
      ie = i-r*ib;
      je = j-r*jb;
      if ( ib == jb )
        return bp + pkn_SymMatIndex ( ie, je );
      else
        return bp + ie*r+je;
    }
  }
  return -1;
} /*pkn_Block3FindElemPos*/

