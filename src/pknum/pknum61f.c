
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2015                                  */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

#include "pkvaria.h"
#include "pknum.h"

boolean pkn_VarSignf ( int n, const float *a )
{
  int  i;
  char s, z;

  s = pkv_signf ( a[0] );
  for ( i = 1; i < n; i++ ) {
    z = pkv_signf ( a[i] );
    if ( z != s )
      return true;
  }
  return false;
} /*pkn_VarSignf*/

boolean pkn_OneSignChangef ( int n, const float *a, boolean *nonzero )
{
  int  i, k;
  char sa, sb;

  k = 0;
  sa = pkv_signf ( a[0] );
  *nonzero = (boolean)(sa != 0);
  for ( i = 1; i < n; i++ ) {
    sb = pkv_signf ( a[i] );
    *nonzero = (boolean)(*nonzero || sb != 0);
    if ( sb != sa ) {
      sa = sb;
      k++;
      if ( k > 1 )
        return false;
    }
  }
  return true;
} /*pkn_OneSignChangef*/

