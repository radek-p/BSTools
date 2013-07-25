
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2010                                  */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

#include <math.h>
#include <string.h>

#include "pkvaria.h"


int pkv_Binom ( int n, int k )
{
  int  b, i;

  if ( k > n-k ) k = n-k;
  b = 1;
  for ( i = 1; i <= k; i++ )
    b = (b*(n+1-i))/i;
  return b;
} /*pkv_Binom*/

