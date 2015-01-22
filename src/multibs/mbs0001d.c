
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2005, 2015                            */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

#include <math.h>
#include <string.h>
#include <stdio.h>  /* for testing */

#include "pkvaria.h"
#include "pkgeom.h"
#include "multibs.h"

/* /////////////////////////////////////////// */
/* Horner scheme for Bezier curves and patches */
boolean mbs_multiBCHornerd ( int degree, int ncurves, int spdimen, int pitch,
                             const double *ctlpoints, double t, double *cpoints )
{
  int         i, j, k;
  long double s, e, be;
  double      *ci, *cci, *di;

  s = 1.0-t;
  pkv_Selectd ( ncurves, spdimen, pitch, spdimen, ctlpoints, cpoints );
  e = t;
  if ( degree <= 29 ) {
    int binom;

    binom = degree;
    for ( i = 1, ci = (double*)&ctlpoints[spdimen];
          i <= degree;
          i++, ci += spdimen ) {
      be = (double)binom*e;
      for ( j = 0, cci = ci, di = cpoints;
            j < ncurves;
            j++, cci += pitch, di += spdimen )
        for ( k = 0;  k < spdimen;  k++ )
          di[k] = s*di[k] + be*cci[k];
      e *= t;
      binom = (binom*(degree-i))/(i+1);
    }
    return true;
  }
  else if ( degree <= 61 ) {
          /* for high degrees 64 bit integers are needed to avoid */
          /* overflow in computing binomial coefficients */
#if __WORDSIZE == 64
    long int binom;
#else
    long long int binom;
#endif

    binom = degree;
    for ( i = 1, ci = (double*)&ctlpoints[spdimen];
          i <= degree;
          i++, ci += spdimen ) {
      be = (long double)binom*e;
      for ( j = 0, cci = ci, di = cpoints;
            j < ncurves;
            j++, cci += pitch, di += spdimen )
        for ( k = 0;  k < spdimen;  k++ )
          di[k] = s*di[k] + be*cci[k];
      e *= t;
      binom = (binom*(degree-i))/(i+1);
    }
    return true;
  }
  else {
    PKV_SIGNALERROR ( LIB_MULTIBS, ERRCODE_8, ERRMSG_8 );
    return false;
  }
} /*mbs_multiBCHornerd*/

