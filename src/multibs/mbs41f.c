
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2005, 2009                            */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdio.h>  /* for testing */

#include "pkvaria.h"
#include "pkgeom.h"
#include "multibs.h"

#include "msgpool.h"

/* /////////////////////////////////////////////////// */
/* computing values and derivatives of B-spline curves */

int mbs_multideBoorDer2f ( int degree, int lastknot, const float *knots,
                           int ncurves, int spdimen,
                           int pitch, const float *ctlpoints,
                           float t, float *p, float *d1, float *d2 )
{
  void   *sp;
  int    i, k, r;
  int    dpitch;
  float  *d, kn1[6], kn2[3];
  float  u, s;

  if ( degree < 2 ) {
    memset ( d2, 0, ncurves*spdimen*sizeof(float) );
    return mbs_multideBoorDerf ( degree, lastknot, knots, ncurves, spdimen,
                                 pitch, ctlpoints, t, p, d1 );
  }

  sp = pkv_GetScratchMemTop ();

  /* Find the proper interval between the knots and the multiplicity of */
  /* the left-side knot. Knots are numbered from 0 to lastknot */

  k = mbs_FindKnotIntervalf ( degree, lastknot, knots, t, &r );

  if ( r <= degree-2 ) {
    dpitch = (degree+1-r)*spdimen;
    d = pkv_GetScratchMemf ( ncurves*dpitch );
    if ( !d ) {
      pkv_SignalError ( LIB_MULTIBS, 41, ERRMSG_0 );
      exit ( 1 );
    }
    _mbs_multideBoorKernelf ( degree, knots, ncurves, spdimen,
                               pitch, ctlpoints, t, k, r, 2, dpitch, d );
    memcpy ( &kn1[0], &knots[k-r-2], 3*sizeof(float) );
  }
  else {
    dpitch = 3*spdimen;
    d = pkv_GetScratchMemf ( ncurves*dpitch );
    if ( !d ) {
      pkv_SignalError ( LIB_MULTIBS, 42, ERRMSG_0 );
      exit ( 1 );
    }
    pkv_Selectf ( ncurves, dpitch, pitch, dpitch,
                  &ctlpoints[(k-degree)*spdimen], &d[0] );
    memcpy ( &kn1[0], &knots[k-degree], 3*sizeof(float) );
  }
  memcpy ( &kn1[3], &knots[k+1], 3*sizeof(float) );
  kn2[0] = kn2[1] = kn2[2] = kn1[2];
  mbs_multiBSChangeLeftKnotsf ( ncurves, spdimen, 2, kn1,
                                dpitch, d, kn2 );
  kn2[0] = kn2[1] = kn2[2] = kn1[3];
  mbs_multiBSChangeRightKnotsf ( ncurves, spdimen, 2, 5, kn1,
                                 dpitch, d, kn2 );

      /* compute the derivatives as differences, with the */
      /* de Casteljau algorithm */
  u = kn1[3]-kn1[2];
  t = (t-kn1[2])/u;  s = (float)(1.0-t);

  pkn_AddMatrixf ( ncurves, spdimen,
                   dpitch, &d[spdimen], dpitch, &d[spdimen],
                   spdimen, d2 );
  pkn_SubtractMatrixf ( ncurves, spdimen,
                        dpitch, d, spdimen, d2, spdimen, d2 );
  pkn_AddMatrixf ( ncurves, spdimen, spdimen, d2, dpitch, &d[2*spdimen],
                   spdimen, d2 );
  pkn_MultMatrixNumf ( 1, ncurves*spdimen, 0, d2,
                       (float)(degree*(degree-1))/(u*u), 0, d2 );

  for ( i = 0; i < 2; i++ )
    pkn_MatrixLinCombf ( ncurves, spdimen, dpitch, &d[i*spdimen], s,
                         dpitch, &d[(i+1)*spdimen], t, dpitch, &d[i*spdimen] );

  pkn_SubtractMatrixf ( ncurves, spdimen, dpitch, &d[spdimen],
                        dpitch, d, spdimen, d1 );
  pkn_MultMatrixNumf ( 1, ncurves*spdimen, 0, d1,
                       (float)(degree)/u, 0, d1 );

  pkn_MatrixLinCombf ( ncurves, spdimen,
                       dpitch, d, s, dpitch, &d[spdimen], t, spdimen, p );

  pkv_SetScratchMemTop ( sp );
  return degree-r;
} /*mbs_multideBoorDer2f*/

int mbs_multideBoorDer3f ( int degree, int lastknot, const float *knots,
                           int ncurves, int spdimen,
                           int pitch, const float *ctlpoints,
                           float t, float *p, float *d1, float *d2, float *d3 )
{
  void   *sp;
  int    i, k, r;
  int    dpitch;
  float  *d, kn1[8], kn2[4];
  float  u, s;

  if ( degree < 3 ) {
    memset ( d3, 0, ncurves*spdimen*sizeof(float) );
    return mbs_multideBoorDer2f ( degree, lastknot, knots, ncurves, spdimen,
                                  pitch, ctlpoints, t, p, d1, d2 );
  }

  sp = pkv_GetScratchMemTop ();

  /* Find the proper interval between the knots and the multiplicity of */
  /* the left-side knot. Knots are numbered from 0 to lastknot */

  k = mbs_FindKnotIntervalf ( degree, lastknot, knots, t, &r );

  if ( r <= degree-3 ) {
    dpitch = (degree+1-r)*spdimen;
    d = pkv_GetScratchMemf ( ncurves*dpitch );
    if ( !d ) {
      pkv_SignalError ( LIB_MULTIBS, 43, ERRMSG_0 );
      exit ( 1 );
    }
    _mbs_multideBoorKernelf ( degree, knots, ncurves, spdimen,
                               pitch, ctlpoints, t, k, r, 3, dpitch, d );
    memcpy ( &kn1[0], &knots[k-r-3], 4*sizeof(float) );
  }
  else {
    dpitch = 4*spdimen;
    d = pkv_GetScratchMemf ( ncurves*dpitch );
    if ( !d ) {
      pkv_SignalError ( LIB_MULTIBS, 44, ERRMSG_0 );
      exit ( 1 );
    }
    pkv_Selectf ( ncurves, dpitch, pitch, dpitch,
                  &ctlpoints[(k-degree)*spdimen], &d[0] );
    memcpy ( &kn1[0], &knots[k-degree], 4*sizeof(float) );
  }
  memcpy ( &kn1[4], &knots[k+1], 4*sizeof(float) );
  kn2[0] = kn2[1] = kn2[2] = kn2[3] = kn1[3];
  mbs_multiBSChangeLeftKnotsf ( ncurves, spdimen, 3, kn1,
                                dpitch, d, kn2 );
  kn2[0] = kn2[1] = kn2[2] = kn2[3] = kn1[4];
  mbs_multiBSChangeRightKnotsf ( ncurves, spdimen, 3, 7, kn1,
                                 dpitch, d, kn2 );

      /* compute the derivatives as differences, with the */
      /* de Casteljau algorithm */
  u = kn1[4]-kn1[3];
  t = (t-kn1[3])/u;  s = (float)(1.0-t);

  pkn_SubtractMatrixf ( ncurves, spdimen,
                        dpitch, &d[spdimen], dpitch, &d[2*spdimen],
                        spdimen, d2 );
  pkn_SubtractMatrixf ( ncurves, spdimen,
                        dpitch, &d[3*spdimen], dpitch, d,
                        spdimen, d3 );
  pkn_AddMatrixMf ( 1, ncurves*spdimen, 0, d3, 0, d2, 3.0, 0, d3 );
  pkn_MultMatrixNumf ( 1, ncurves*spdimen, 0, d3,
                       (float)(degree*(degree-1)*(degree-2))/(u*u*u), 0, d3 );

  for ( i = 0; i < 3; i++ )
    pkn_MatrixLinCombf ( ncurves, spdimen, dpitch, &d[i*spdimen], s,
                         dpitch, &d[(i+1)*spdimen], t, dpitch, &d[i*spdimen] );
  pkn_AddMatrixf ( ncurves, spdimen,
                   dpitch, &d[spdimen], dpitch, &d[spdimen],
                   spdimen, d2 );
  pkn_SubtractMatrixf ( ncurves, spdimen,
                        dpitch, d, spdimen, d2, spdimen, d2 );
  pkn_AddMatrixf ( ncurves, spdimen, spdimen, d2, dpitch, &d[2*spdimen],
                   spdimen, d2 );
  pkn_MultMatrixNumf ( 1, ncurves*spdimen, 0, d2,
                       (float)(degree*(degree-1))/(u*u), 0, d2 );

  for ( i = 0; i < 2; i++ )
    pkn_MatrixLinCombf ( ncurves, spdimen, dpitch, &d[i*spdimen], s,
                         dpitch, &d[(i+1)*spdimen], t, dpitch, &d[i*spdimen] );
  pkn_SubtractMatrixf ( ncurves, spdimen, dpitch, &d[spdimen],
                        dpitch, d, spdimen, d1 );
  pkn_MultMatrixNumf ( 1, ncurves*spdimen, 0, d1,
                       (float)(degree)/u, 0, d1 );

  pkn_MatrixLinCombf ( ncurves, spdimen,
                       dpitch, d, s, dpitch, &d[spdimen], t, spdimen, p );

  pkv_SetScratchMemTop ( sp );
  return degree-r;
} /*mbs_multideBoorDer3f*/

