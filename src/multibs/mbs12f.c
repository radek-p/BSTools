
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

#include "msgpool.h"

/* /////////////////////////////////////////// */
/* Horner scheme for Bezier curves and patches */

void mbs_multiBCHornerf ( int degree, int ncurves, int spdimen, int pitch,
                          const float *ctlpoints, float t, float *cpoints )
{
  int    i, j, k;
  double s, e, be;
  float  *ci, *cci, *di;

  s = 1.0-t;
  pkv_Selectf ( ncurves, spdimen, pitch, spdimen, ctlpoints, cpoints );
  e = t;
  if ( degree <= 29 ) {
    int binom;

    binom = degree;
    for ( i = 1, ci = (float*)&ctlpoints[spdimen];
          i <= degree;
          i++, ci += spdimen ) {
      be = (double)binom*e;
      for ( j = 0, cci = ci, di = cpoints;
            j < ncurves;
            j++, cci += pitch, di += spdimen )
        for ( k = 0;  k < spdimen;  k++ )
          di[k] = (float)(s*di[k] + be*cci[k]);
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
    for ( i = 1, ci = (float*)&ctlpoints[spdimen];
          i <= degree;
          i++, ci += spdimen ) {
      be = (double)binom*e;
      for ( j = 0, cci = ci, di = cpoints;
            j < ncurves;
            j++, cci += pitch, di += spdimen )
        for ( k = 0;  k < spdimen;  k++ )
          di[k] = (float)(s*di[k] + be*cci[k]);
      e *= t;
      binom = (binom*(degree-i))/(i+1);
    }
  }
  else
    pkv_SignalError ( LIB_MULTIBS, 65, ERRMSG_4 );
} /*mbs_multiBCHornerf*/

void mbs_BCHornerC2Rf ( int degree, const point3f *ctlpoints, float t,
                        point2f *cpoint )
{
  point3f hcpoint;

  mbs_multiBCHornerf ( degree, 1, 3, 0, (float*)ctlpoints, t,
                       (float*)&hcpoint );
  Point3to2f ( &hcpoint, cpoint );
} /*mbs_BCHornerC2Rf*/

void mbs_BCHornerC3Rf ( int degree, const point4f *ctlpoints, float t,
                        point3f *cpoint )
{
  point4f hcpoint;

  mbs_multiBCHornerf ( degree, 1, 4, 0, (float*)ctlpoints, t,
                       (float*)&hcpoint );
  Point4to3f ( &hcpoint, cpoint );
} /*mbs_BCHornerC3Rf*/

void mbs_BCHornerPf ( int degreeu, int degreev, int spdimen,
                      const float *ctlpoints,
                      float u, float v, float *ppoint )
{
  float *auxc;
  int   auxc_size;

  auxc_size = (degreev+1)*spdimen*sizeof(float);
  if ( (auxc = pkv_GetScratchMem ( auxc_size )) ) {
    mbs_multiBCHornerf ( degreeu, 1, spdimen*(degreev+1), 0, ctlpoints,
                         u, auxc );
    mbs_multiBCHornerf ( degreev, 1, spdimen, 0, auxc, v, ppoint );
    pkv_FreeScratchMem ( auxc_size );
  }
} /*mbs_BCHornerPf*/

void mbs_BCHornerP3Rf ( int degreeu, int degreev, const point4f *ctlpoints,
                        float u, float v, point3f *ppoint )
{
  point4f auxp;

  mbs_BCHornerP4f ( degreeu, degreev, ctlpoints, u, v, &auxp );
  Point4to3f ( &auxp, ppoint );
} /*mbs_BCHornerP3Rf*/

