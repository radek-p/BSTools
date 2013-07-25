
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2005                                  */
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

double pkn_GoldenRatd ( double (*f) (double), double a, double b, double eps,
                        boolean *error )
{
#define GP 0.6180339887498949
  double c, d, fc, fd;

  if ( a > b )
    { c = b;  b = a;  a = c; }

  c = GP*a + (1.0-GP)*b;  fc = f(c);
  d = (1.0-GP)*a + GP*b;  fd = f(d);
  while ( b-a > eps ) {
    if ( fc > fd ) {
      a = c;  c = d;  fc = fd;
      d = (1.0-GP)*a + GP*b;  fd = f(d);
    }
    else {
      b = d;  d = c;  fd = fc;
      c = GP*a + (1.0-GP)*b;  fc = f(c);
    }
  }
  *error = false;
  return 0.5*(a+b);
#undef GP
} /*pkn_GoldenRatd*/

