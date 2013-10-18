
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2005, 2013                            */
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

/* /////////////////////////////////////////////////// */
/* computing values and derivatives of B-spline curves */

int mbs_multideBoorDer2d ( int degree, int lastknot, const double *knots,
                           int ncurves, int spdimen,
                           int pitch, const double *ctlpoints,
                           double t, double *p, double *d1, double *d2 )
{
  void   *sp;
  int    i, k, r;
  int    dpitch;
  double  *d, kn1[6], kn2[3];
  double  u, s;

  if ( degree < 2 ) {
    memset ( d2, 0, ncurves*spdimen*sizeof(double) );
    return mbs_multideBoorDerd ( degree, lastknot, knots, ncurves, spdimen,
                                 pitch, ctlpoints, t, p, d1 );
  }

  sp = pkv_GetScratchMemTop ();

  /* Find the proper interval between the knots and the multiplicity of */
  /* the left-side knot. Knots are numbered from 0 to lastknot */

  k = mbs_FindKnotIntervald ( degree, lastknot, knots, t, &r );

  if ( r <= degree-2 ) {
    dpitch = (degree+1-r)*spdimen;
    d = pkv_GetScratchMemd ( ncurves*dpitch );
    if ( !d ) {
      PKV_SIGNALERROR ( LIB_MULTIBS, ERRCODE_2, ERRMSG_2 );
      goto failure;
    }
    _mbs_multideBoorKerneld ( degree, knots, ncurves, spdimen,
                               pitch, ctlpoints, t, k, r, 2, dpitch, d );
    memcpy ( &kn1[0], &knots[k-r-2], 3*sizeof(double) );
  }
  else {
    dpitch = 3*spdimen;
    d = pkv_GetScratchMemd ( ncurves*dpitch );
    if ( !d ) {
      PKV_SIGNALERROR ( LIB_MULTIBS, ERRCODE_2, ERRMSG_2 );
      goto failure;
    }
    pkv_Selectd ( ncurves, dpitch, pitch, dpitch,
                  &ctlpoints[(k-degree)*spdimen], &d[0] );
    memcpy ( &kn1[0], &knots[k-degree], 3*sizeof(double) );
  }
  memcpy ( &kn1[3], &knots[k+1], 3*sizeof(double) );
  kn2[0] = kn2[1] = kn2[2] = kn1[2];
  if ( !mbs_multiBSChangeLeftKnotsd ( ncurves, spdimen, 2, kn1,
                                      dpitch, d, kn2 ) )
    goto failure;
  kn2[0] = kn2[1] = kn2[2] = kn1[3];
  if ( !mbs_multiBSChangeRightKnotsd ( ncurves, spdimen, 2, 5, kn1,
                                       dpitch, d, kn2 ) )
    goto failure;

      /* compute the derivatives as differences, with the */
      /* de Casteljau algorithm */
  u = kn1[3]-kn1[2];
  t = (t-kn1[2])/u;  s = 1.0-t;

  pkn_AddMatrixd ( ncurves, spdimen,
                   dpitch, &d[spdimen], dpitch, &d[spdimen],
                   spdimen, d2 );
  pkn_SubtractMatrixd ( ncurves, spdimen,
                        dpitch, d, spdimen, d2, spdimen, d2 );
  pkn_AddMatrixd ( ncurves, spdimen, spdimen, d2, dpitch, &d[2*spdimen],
                   spdimen, d2 );
  pkn_MultMatrixNumd ( 1, ncurves*spdimen, 0, d2,
                       (double)(degree*(degree-1))/(u*u), 0, d2 );

  for ( i = 0; i < 2; i++ )
    pkn_MatrixLinCombd ( ncurves, spdimen, dpitch, &d[i*spdimen], s,
                         dpitch, &d[(i+1)*spdimen], t, dpitch, &d[i*spdimen] );

  pkn_SubtractMatrixd ( ncurves, spdimen, dpitch, &d[spdimen],
                        dpitch, d, spdimen, d1 );
  pkn_MultMatrixNumd ( 1, ncurves*spdimen, 0, d1,
                       (double)(degree)/u, 0, d1 );

  pkn_MatrixLinCombd ( ncurves, spdimen,
                       dpitch, d, s, dpitch, &d[spdimen], t, spdimen, p );

  pkv_SetScratchMemTop ( sp );
  return degree-r;

failure:
  pkv_SetScratchMemTop ( sp );
  return -1;
} /*mbs_multideBoorDer2d*/

int mbs_multideBoorDer3d ( int degree, int lastknot, const double *knots,
                           int ncurves, int spdimen,
                           int pitch, const double *ctlpoints,
                           double t, double *p, double *d1, double *d2, double *d3 )
{
  void   *sp;
  int    i, k, r;
  int    dpitch;
  double  *d, kn1[8], kn2[4];
  double  u, s;

  if ( degree < 3 ) {
    memset ( d3, 0, ncurves*spdimen*sizeof(double) );
    return mbs_multideBoorDer2d ( degree, lastknot, knots, ncurves, spdimen,
                                  pitch, ctlpoints, t, p, d1, d2 );
  }

  sp = pkv_GetScratchMemTop ();

  /* Find the proper interval between the knots and the multiplicity of */
  /* the left-side knot. Knots are numbered from 0 to lastknot */

  k = mbs_FindKnotIntervald ( degree, lastknot, knots, t, &r );

  if ( r <= degree-3 ) {
    dpitch = (degree+1-r)*spdimen;
    d = pkv_GetScratchMemd ( ncurves*dpitch );
    if ( !d ) {
      PKV_SIGNALERROR ( LIB_MULTIBS, ERRCODE_2, ERRMSG_2 );
      goto failure;
    }
    _mbs_multideBoorKerneld ( degree, knots, ncurves, spdimen,
                               pitch, ctlpoints, t, k, r, 3, dpitch, d );
    memcpy ( &kn1[0], &knots[k-r-3], 4*sizeof(double) );
  }
  else {
    dpitch = 4*spdimen;
    d = pkv_GetScratchMemd ( ncurves*dpitch );
    if ( !d ) {
      PKV_SIGNALERROR ( LIB_MULTIBS, ERRCODE_2, ERRMSG_2 );
      goto failure;
    }
    pkv_Selectd ( ncurves, dpitch, pitch, dpitch,
                  &ctlpoints[(k-degree)*spdimen], &d[0] );
    memcpy ( &kn1[0], &knots[k-degree], 4*sizeof(double) );
  }
  memcpy ( &kn1[4], &knots[k+1], 4*sizeof(double) );
  kn2[0] = kn2[1] = kn2[2] = kn2[3] = kn1[3];
  if ( !mbs_multiBSChangeLeftKnotsd ( ncurves, spdimen, 3, kn1,
                                      dpitch, d, kn2 ) )
    goto failure;
  kn2[0] = kn2[1] = kn2[2] = kn2[3] = kn1[4];
  if ( !mbs_multiBSChangeRightKnotsd ( ncurves, spdimen, 3, 7, kn1,
                                       dpitch, d, kn2 ) )
    goto failure;

      /* compute the derivatives as differences, with the */
      /* de Casteljau algorithm */
  u = kn1[4]-kn1[3];
  t = (t-kn1[3])/u;  s = 1.0-t;

  pkn_SubtractMatrixd ( ncurves, spdimen,
                        dpitch, &d[spdimen], dpitch, &d[2*spdimen],
                        spdimen, d2 );
  pkn_SubtractMatrixd ( ncurves, spdimen,
                        dpitch, &d[3*spdimen], dpitch, d,
                        spdimen, d3 );
  pkn_AddMatrixMd ( 1, ncurves*spdimen, 0, d3, 0, d2, 3.0, 0, d3 );
  pkn_MultMatrixNumd ( 1, ncurves*spdimen, 0, d3,
                       (double)(degree*(degree-1)*(degree-2))/(u*u*u), 0, d3 );

  for ( i = 0; i < 3; i++ )
    pkn_MatrixLinCombd ( ncurves, spdimen, dpitch, &d[i*spdimen], s,
                         dpitch, &d[(i+1)*spdimen], t, dpitch, &d[i*spdimen] );
  pkn_AddMatrixd ( ncurves, spdimen,
                   dpitch, &d[spdimen], dpitch, &d[spdimen],
                   spdimen, d2 );
  pkn_SubtractMatrixd ( ncurves, spdimen,
                        dpitch, d, spdimen, d2, spdimen, d2 );
  pkn_AddMatrixd ( ncurves, spdimen, spdimen, d2, dpitch, &d[2*spdimen],
                   spdimen, d2 );
  pkn_MultMatrixNumd ( 1, ncurves*spdimen, 0, d2,
                       (double)(degree*(degree-1))/(u*u), 0, d2 );

  for ( i = 0; i < 2; i++ )
    pkn_MatrixLinCombd ( ncurves, spdimen, dpitch, &d[i*spdimen], s,
                         dpitch, &d[(i+1)*spdimen], t, dpitch, &d[i*spdimen] );
  pkn_SubtractMatrixd ( ncurves, spdimen, dpitch, &d[spdimen],
                        dpitch, d, spdimen, d1 );
  pkn_MultMatrixNumd ( 1, ncurves*spdimen, 0, d1,
                       (double)(degree)/u, 0, d1 );

  pkn_MatrixLinCombd ( ncurves, spdimen,
                       dpitch, d, s, dpitch, &d[spdimen], t, spdimen, p );

  pkv_SetScratchMemTop ( sp );
  return degree-r;

failure:
  pkv_SetScratchMemTop ( sp );
  return -1;
} /*mbs_multideBoorDer3d*/

