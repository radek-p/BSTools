
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Mateusz Markowski                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2010                                  */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "pkvaria.h"
#include "pknum.h"
#include "pkgeom.h"
#include "multibs.h"
#include "g1blendingd.h"

#include "g1blprivated.h"

int _g1bl_Asym3MatIndex ( int i, int j, int k, boolean *neg )
{
  int ii, jj, kk;

        /* sort and find the permutation sign */
  if ( i < j ) {
    if ( j < k )
      { ii = i;  jj = j;  kk = k;  *neg = false; }
    else {
      if ( i < k )
        { ii = i;  jj = k;  kk = j;  *neg = true; }
      else
        { ii = k;  jj = i;  kk = j;  *neg = false; }
    }
  }
  else {
    if ( j < k ) {
      if ( i < k )
        { ii = j;  jj = i;  kk = k;  *neg = true; }
      else
        { ii = j;  jj = k;  kk = i;  *neg = false; }
    }
    else
      { ii = k;  jj = j;  kk = i;  *neg = true; }
  }
        /* compute the index */
  jj -= 1;
  kk -= 2;
  return (((kk+3)*kk+2)*kk)/6 + ((jj+1)*jj)/2 + ii;
} /*_g1bl_Asym3MatIndex*/

int g1bl_NiSize ( int nkn )
{
  return 9*5*nkn*nkn;
} /*g2bl_NiSize*/

int g1bl_NijSize ( int nkn )
{
  return 45*3*nkn*nkn; /* 9 * 10/2 - liczba par dwuwskaźników t. że: (i0,i1) =< (j0,j1) */
} /*g2bl_NijSize*/

int g1bl_MijSize ( int nkn ) /* 9 * 8/2 - liczba par dwuwskaźników t. że: (i0,i1) < (j0,j1) */
{
  return 36*7*nkn*nkn;
} /*g2bl_MijSize*/

