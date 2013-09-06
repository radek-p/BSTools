
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2005, 2013                            */
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
void mbs_multiBCHornerd ( int degree, int ncurves, int spdimen, int pitch,
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
  }
  else
    PKV_SIGNALERROR ( LIB_MULTIBS, ERRCODE_8, ERRMSG_8 );
} /*mbs_multiBCHornerd*/

void mbs_BCHornerC2Rd ( int degree, const point3d *ctlpoints, double t,
                        point2d *cpoint )
{
  point3d hcpoint;

  mbs_multiBCHornerd ( degree, 1, 3, 0, (double*)ctlpoints, t,
                       (double*)&hcpoint );
  Point3to2d ( &hcpoint, cpoint );
} /*mbs_BCHornerC2Rd*/

void mbs_BCHornerC3Rd ( int degree, const point4d *ctlpoints, double t,
                        point3d *cpoint )
{
  point4d hcpoint;

  mbs_multiBCHornerd ( degree, 1, 4, 0, (double*)ctlpoints, t,
                       (double*)&hcpoint );
  Point4to3d ( &hcpoint, cpoint );
} /*mbs_BCHornerC3Rd*/

void mbs_BCHornerPd ( int degreeu, int degreev, int spdimen,
                      const double *ctlpoints,
                      double u, double v, double *ppoint )
{
  double *auxc;
  int   auxc_size;

  auxc_size = (degreev+1)*spdimen*sizeof(double);
  if ( (auxc = pkv_GetScratchMem ( auxc_size )) ) {
    mbs_multiBCHornerd ( degreeu, 1, spdimen*(degreev+1), 0, ctlpoints,
                         u, auxc );
    mbs_multiBCHornerd ( degreev, 1, spdimen, 0, auxc, v, ppoint );
    pkv_FreeScratchMem ( auxc_size );
  }
} /*mbs_BCHornerPd*/

void mbs_BCHornerP3Rd ( int degreeu, int degreev, const point4d *ctlpoints,
                        double u, double v, point3d *ppoint )
{
  point4d auxp;

  mbs_BCHornerP4d ( degreeu, degreev, ctlpoints, u, v, &auxp );
  Point4to3d ( &auxp, ppoint );
} /*mbs_BCHornerP3Rd*/

