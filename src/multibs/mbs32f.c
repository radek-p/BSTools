
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

/* /////////////////////////////////////////// */
/* Horner scheme for Bezier curves and patches */

boolean mbs_multiBCHornerDerf ( int degree, int ncurves, int spdimen, int pitch,
                                const float *ctlpoints,
                                float t, float *p, float *d )
{
  int   i;
  float *a;

  if ( degree == 0 ) {
    pkv_Selectf ( ncurves, spdimen, pitch, spdimen, ctlpoints, p );
    memset ( d, 0, ncurves*spdimen*sizeof(float) );
    return true;
  }
  else if ( (a = pkv_GetScratchMemf ( 2*spdimen )) ) {
    for ( i = 0;
          i < ncurves;
          i++, ctlpoints += pitch, p += spdimen, d += spdimen ) {
      mbs_multiBCHornerf ( degree-1, 2, spdimen, spdimen, ctlpoints, t, a );
      pkn_MatrixMDifferencef ( 1, spdimen, 0, &a[spdimen], 0, &a[0],
                               (float)degree, 0, d );
      pkn_MatrixLinCombf ( 1, spdimen, 0, &a[0], 1.0-t, 0, &a[spdimen], t,
                           0, p );
    }
    pkv_FreeScratchMemf ( 2*spdimen );
    return true;
  }
  else
    return false;
} /*mbs_multiBCHornerDerf*/

boolean mbs_BCHornerDerC2Rf ( int degree, const point3f *ctlpoints, float t,
                              point2f *p, vector2f *d )
{
  point3f hp, hd;

  if ( !mbs_BCHornerDerC3f ( degree, ctlpoints, t, &hp, &hd ) )
    return false;
  Point3to2f ( &hp, p );
  memcpy ( d, &hd, sizeof(vector2f) );
  AddVector2Mf ( d, p, -hd.z, d );
  MultVector2f ( 1.0/hp.z, d, d );
  return true;
} /*mbs_BCHornerDerC2Rf*/

boolean mbs_BCHornerDerC3Rf ( int degree, const point4f *ctlpoints, float t,
                              point3f *p, vector3f *d )
{
  point4f hp, hd;

  if ( !mbs_BCHornerDerC4f ( degree, ctlpoints, t, &hp, &hd ) )
    return false;
  Point4to3f ( &hp, p );
  memcpy ( d, &hd, sizeof(vector3f) );
  AddVector3Mf ( d, p, -hd.w, d );
  MultVector3f ( 1.0/hp.w, d, d );
  return true;
} /*mbs_BCHornerDerC3Rf*/

boolean mbs_BCHornerDerPf ( int degreeu, int degreev, int spdimen,
                            const float *ctlpoints,
                            float u, float v,
                            float *p, float *du, float *dv )
{
  void *sp;
  float *scr, *q;
  int   scr_size;

  sp = pkv_GetScratchMemTop ();
  if ( degreeu == 0 ) {
    if ( degreev == 0 ) {
      memcpy ( p, ctlpoints, spdimen*sizeof(float) );
      memset ( dv, 0, spdimen*sizeof(float) );
    }
    else {
      if ( !mbs_multiBCHornerDerf ( degreev, 1, spdimen, 0, ctlpoints, v, p, dv ) )
        goto failure;
    }
    memset ( du, 0, spdimen*sizeof(float) );
  }
  else if ( degreev == 0 ) {
    if ( !mbs_multiBCHornerDerf ( degreeu, 1, spdimen, 0, ctlpoints, u, p, du ) )
      goto failure;
    memset ( dv, 0, spdimen*sizeof(float) );
  }
  else {
    q = pkv_GetScratchMemf ( scr_size = (6+2*(degreev+1))*spdimen );
    if ( !q ) {
      PKV_SIGNALERROR ( LIB_MULTIBS, ERRCODE_2, ERRMSG_2 );
      goto failure;
    }
    scr = &q[6*spdimen];
    if ( !mbs_multiBCHornerf ( degreeu-1, 2, spdimen*(degreev+1), spdimen*(degreev+1),
                               ctlpoints, u, scr ) )
      goto failure; 
    if ( !mbs_multiBCHornerf ( degreev-1, 2, spdimen, spdimen,
                               &scr[0], v, &q[0] ) )
      goto failure;
    if ( !mbs_multiBCHornerf ( degreev-1, 2, spdimen, spdimen,
                               &scr[spdimen*(degreev+1)], v, &q[2*spdimen] ) )
      goto failure;
    pkn_SubtractMatrixf ( 2, spdimen, spdimen, &q[2*spdimen], spdimen, &q[0],
                          spdimen, &q[4*spdimen] );
    pkn_MatrixLinCombf ( 1, spdimen, 0, &q[4*spdimen], (1.0-v)*(float)degreeu,
                         0, &q[5*spdimen], v*(float)degreeu, 0, du );
    pkn_MatrixLinCombf ( 2, spdimen, spdimen, &q[0], 1.0-u,
                         spdimen, &q[2*spdimen], u, spdimen, &q[0] );
    pkn_MatrixMDifferencef ( 1, spdimen, 0, &q[spdimen], 0, &q[0],
                             (float)degreev, 0, dv );
    pkn_MatrixLinCombf ( 1, spdimen, 0, &q[0], 1.0-v, 0, &q[spdimen], v,
                         0, p );
    pkv_FreeScratchMemf ( scr_size );
  }
  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*mbs_BCHornerDerPf*/

boolean mbs_BCHornerDerP3Rf ( int degreeu, int degreev, const point4f *ctlpoints,
                              float u, float v,
                              point3f *p, vector3f *du, vector3f *dv )
{
  point4f hp, hdu, hdv;

  if ( mbs_BCHornerDerP4f ( degreeu, degreev, ctlpoints, u, v, &hp, &hdu, &hdv ) ) {
    Point4to3f ( &hp, p );
    memcpy ( du, &hdu, sizeof(vector3f) );
    AddVector3Mf ( du, p, -hdu.w, du );
    MultVector3f ( 1.0/hp.w, du, du );
    memcpy ( dv, &hdv, sizeof(vector3f) );
    AddVector3Mf ( dv, p, -hdv.w, dv );
    MultVector3f ( 1.0/hp.w, dv, dv );
    return true;
  }
  else
    return false;
} /*mbs_BCHornerDerP3Rf*/

boolean mbs_BCHornerNvP3f ( int degreeu, int degreev, const point3f *ctlpoints,
                            float u, float v,
                            point3f *p, vector3f *nv )
{
  vector3f du, dv;

  if ( mbs_BCHornerDerP3f ( degreeu, degreev, ctlpoints, u, v, p, &du, &dv ) ) {
    CrossProduct3f ( &du, &dv, nv );
    NormalizeVector3f ( nv );
    return true;
  }
  else
    return false;
} /*mbs_BCHornerNvP3f*/

boolean mbs_BCHornerNvP3Rf ( int degreeu, int degreev, const point4f *ctlpoints,
                             float u, float v,
                             point3f *p, vector3f *nv )
{
  point4f  P;
  vector4f Du, Dv;

  if ( mbs_BCHornerDerP4f ( degreeu, degreev, ctlpoints, u, v, &P, &Du, &Dv ) ) {
    Point4to3f ( &P, p );
    CrossProduct4P3f ( &P, &Du, &Dv, nv );
    NormalizeVector3f ( nv );
    return true;
  }
  else
    return false;
} /*mbs_BCHornerNvP3Rf*/

