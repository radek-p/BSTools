
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


/* /////////////////////////////////////////// */
/* compute zeros of the polynomial x^2+2px+q   */

boolean pkn_SolveSqEqd ( double p, double q, double *x1, double *x2 )
{
  double delta;

  delta = p*p-q;
  if ( delta >= 0.0 ) {  /* real zeros */
    delta = sqrt ( delta );
    if ( p > 0.0 ) {
      *x1 = -p - delta;
      *x2 = q / *x1;
    }
    else if ( p < 0.0 ) {
      *x2 = -p + delta;
      *x1 = q / *x2;
    }
    else {
      *x2 = delta;
      *x1 = -*x2;
    }
    return true;
  }
  else {                 /* complex zeros */
    *x1 = -p;
    *x2 = sqrt ( -delta );
    return false;
  }
} /*pkn_SolveSqEqd*/

