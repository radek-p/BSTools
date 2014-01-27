
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2013, 2014                            */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

#include "pkvaria.h"
#include "pknum.h"

#define GRDIV(a,b) (exp((1.0-TAU)*log(a)+TAU*log(b)))

/* ///////////////////////////////////////////////////////////////////////// */
boolean _pkn_DivideIntervald ( double *ga, double *gc, double *gd, double *gb,
                               double *fa, double *fc, double *fd, double *fb )
{
  double  x[4], y[4], a, b, c, x1, x2;
  int     i, j;
  boolean ok;
  
        /* construct the cubic polynomial of interpolation */
  x[0] = log(*ga);  x[1] = log(*gc);  x[2] = log(*gd);  x[3] = log(*gb);
  y[0] = *fa;       y[1] = *fc;       y[2] = *fd;       y[3] = *fb;
          /* compute the divided differences */
  for ( i = 1; i <= 3; i++ )
    for ( j = 3; j >= i; j-- )
      y[j] = (y[j]-y[j-1])/(x[j]-x[j-i]);
          /* find the derivative coefficients */
  a = 6.0*y[3];
  b = 2.0*y[2] - 4.0*y[3]*(x[0]+x[1]+x[2]);
  c = y[1] - y[2]*(x[0]+x[1]) + 2.0*y[3]*(x[0]*(x[1]+x[2])+x[1]*x[2]);
  ok = false;
  if ( a ) {
    if ( pkn_SolveSqEqd ( b/(2.0*a), c/a, &x1, &x2 ) ) {
      if ( a > 0.0 )
        x1 = x2;
      ok = true;
    }
  }
  else {
    x1 = -c/b;
    ok = true;
  }
  if ( ok ) {
    x2 = y[3];
    for ( i = 2; i >= 0; i-- )
      x2 = x2*(x1-x[i]) + y[i];
    if ( x2 <= 0.25*(*fa+*fb) )
      ok = false;
  }
  if ( *fc < *fd ) {
    *gb = *gd;  *fb = *fd;
    if ( ok && x1 > x[0] && x1 < x[2] ) {
      if ( x1 > x[1] ) {
        *gd = exp ( x1 );
        return false;
      }
      *gd = *gc;  *fd = *fc;
      if ( x1 < x[1] )
        *gc = exp ( x1 );
      else
        *gc = GRDIV ( *ga, *gd );
    }
    else {
      *gd = *gc;  *fd = *fc;
      *gc = GRDIV ( *ga, *gd );
    }
    return true;
  }
  else {
    *ga = *gc;  *fa = *fc;
    if ( ok && x1 > x[1] && x1 < x[3] ) {
      if ( x1 < x[2] ) {
        *gc = exp ( x1 );
        return true;
      }
      *gc = *gd;  *fc = *fd;
      if ( x1 > x[2] )
        *gd = exp ( x1 );
      else
        *gd = GRDIV ( *gb, *gc );
    }
    else {
      *gc = *gd;  *fc = *fd;
      *gd = GRDIV ( *gb, *gc );
    }
    return false;
  }
} /*_pkn_DivideIntervald*/

