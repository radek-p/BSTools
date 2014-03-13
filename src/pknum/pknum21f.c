
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2005, 2014                            */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

#include <math.h>
#include <string.h>

#include "pkvaria.h"
#include "pknum.h"


/* ///////////////////////////////////////////////// */
/* finding the minimum of a function of one variable */
/* Golden ratio division algorithm                   */

float pkn_GoldenRatf ( float (*f)(void*,float), void *usrptr,
                       float a, float b, float eps, boolean *error )
{
#define GP 0.6180339887498949
  float c, d, fc, fd;

  if ( a > b )
    { c = b;  b = a;  a = c; }

  c = (float)(GP*a + (1.0-GP)*b);  fc = f ( usrptr, c );
  d = (float)((1.0-GP)*a + GP*b);  fd = f ( usrptr, d );
  while ( b-a > eps ) {
    if ( fc > fd ) {
      a = c;  c = d;  fc = fd;
      d = (float)((1.0-GP)*a + GP*b);  fd = f ( usrptr, d );
    }
    else {
      b = d;  d = c;  fd = fc;
      c = (float)(GP*a + (1.0-GP)*b);  fc = f ( usrptr, c );
    }
  }
  *error = false;
  return (float)0.5*(a+b);
#undef GP
} /*pkn_GoldenRatf*/

