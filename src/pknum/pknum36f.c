
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2009                                  */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdarg.h>

#include "pkvaria.h"

#undef CONST_
#define CONST_

#include "pknum.h"

/* ///////////////////////////////////////////////////////////////////////// */
boolean pkn_QuadRectanglesf ( float a, float b, int n,
                              float *qknots, float *qcoeff )
{
  int   i;
  float c, d;

  if ( n < 1 )
    return false;
  d = b-a;
  if ( qknots ) {
    for ( i = 0; i < n; i++ )
      qknots[i] = a + d*(float)(i+i+1)/(float)(n+n);
  }
  if ( qcoeff ) {
    c = d/(float)n;
    for ( i = 0; i < n; i++ )
      qcoeff[i] = c;
  }
  return true;
} /*pkn_QuadRectanglesf*/

boolean pkn_QuadSimpsonf ( float a, float b, int n,    
                           float *qknots, float *qcoeff )
{
  int   i;
  float c, d;

  if ( n < 3 || !(n & 0x01) )
    return false;
  d = b-a;
  if ( qknots ) {
    for ( i = 0; i < n; i++ )
      qknots[i] = a + d*(float)i/(float)(n-1);
  }
  if ( qcoeff ) {
    c = (float)(d/(3.0*(float)(n-1)));
    qcoeff[0] = qcoeff[n-1] = c;
    c += c;
    d = c+c;
    for ( i = 1; i < n-2; i += 2 ) {
      qcoeff[i] = d;
      qcoeff[i+1] = c;
    }
    qcoeff[n-2] = d;
  }
  return true;
} /*pkn_QuadSimpsonf*/

