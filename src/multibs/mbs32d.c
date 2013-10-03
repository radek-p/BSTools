
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

boolean mbs_multiBCHornerDerd ( int degree, int ncurves, int spdimen, int pitch,
                                const double *ctlpoints,
                                double t, double *p, double *d )
{
  int   i;
  double *a;

  if ( degree == 0 ) {
    pkv_Selectd ( ncurves, spdimen, pitch, spdimen, ctlpoints, p );
    memset ( d, 0, ncurves*spdimen*sizeof(double) );
    return true;
  }
  else if ( (a = pkv_GetScratchMemd ( 2*spdimen )) ) {
    for ( i = 0;
          i < ncurves;
          i++, ctlpoints += pitch, p += spdimen, d += spdimen ) {
      mbs_multiBCHornerd ( degree-1, 2, spdimen, spdimen, ctlpoints, t, a );
      pkn_MatrixMDifferenced ( 1, spdimen, 0, &a[spdimen], 0, &a[0],
                               (double)degree, 0, d );
      pkn_MatrixLinCombd ( 1, spdimen, 0, &a[0], 1.0-t, 0, &a[spdimen], t,
                           0, p );
    }
    pkv_FreeScratchMemd ( 2*spdimen );
    return true;
  }
  else
    return false;
} /*mbs_multiBCHornerDerd*/

void mbs_BCHornerDerC2Rd ( int degree, const point3d *ctlpoints, double t,
                           point2d *p, vector2d *d )
{
  point3d hp, hd;

  mbs_BCHornerDerC3d ( degree, ctlpoints, t, &hp, &hd );
  Point3to2d ( &hp, p );
  memcpy ( d, &hd, sizeof(vector2d) );
  AddVector2Md ( d, p, -hd.z, d );
  MultVector2d ( 1.0/hp.z, d, d );
} /*mbs_BCHornerDerC2Rd*/

void mbs_BCHornerDerC3Rd ( int degree, const point4d *ctlpoints, double t,
                           point3d *p, vector3d *d )
{
  point4d hp, hd;

  mbs_BCHornerDerC4d ( degree, ctlpoints, t, &hp, &hd );
  Point4to3d ( &hp, p );
  memcpy ( d, &hd, sizeof(vector3d) );
  AddVector3Md ( d, p, -hd.w, d );
  MultVector3d ( 1.0/hp.w, d, d );
} /*mbs_BCHornerDerC3Rd*/

boolean mbs_BCHornerDerPd ( int degreeu, int degreev, int spdimen,
                            const double *ctlpoints,
                            double u, double v,
                            double *p, double *du, double *dv )
{
  void   *sp;
  double *scr, *q;
  int    scr_size;

  sp = pkv_GetScratchMemTop ();
  if ( degreeu == 0 ) {
    if ( degreev == 0 ) {
      memcpy ( p, ctlpoints, spdimen*sizeof(double) );
      memset ( dv, 0, spdimen*sizeof(double) );
    }
    else {
      mbs_multiBCHornerDerd ( degreev, 1, spdimen, 0, ctlpoints, v, p, dv );
    }
    memset ( du, 0, spdimen*sizeof(double) );
  }
  else if ( degreev == 0 ) {
    mbs_multiBCHornerDerd ( degreeu, 1, spdimen, 0, ctlpoints, u, p, du );
    memset ( dv, 0, spdimen*sizeof(double) );
  }
  else {
    q = pkv_GetScratchMemd ( scr_size = (6+2*(degreev+1))*spdimen );
    if ( !q ) {
      PKV_SIGNALERROR ( LIB_MULTIBS, ERRCODE_2, ERRMSG_2 );
      goto failure;
    }
    scr = &q[6*spdimen];
    mbs_multiBCHornerd ( degreeu-1, 2, spdimen*(degreev+1), spdimen*(degreev+1),
                         ctlpoints, u, scr );
    mbs_multiBCHornerd ( degreev-1, 2, spdimen, spdimen,
                         &scr[0], v, &q[0] );
    mbs_multiBCHornerd ( degreev-1, 2, spdimen, spdimen,
                         &scr[spdimen*(degreev+1)], v, &q[2*spdimen] );
    pkn_SubtractMatrixd ( 2, spdimen, spdimen, &q[2*spdimen], spdimen, &q[0],
                          spdimen, &q[4*spdimen] );
    pkn_MatrixLinCombd ( 1, spdimen, 0, &q[4*spdimen], (1.0-v)*(double)degreeu,
                         0, &q[5*spdimen], v*(double)degreeu, 0, du );
    pkn_MatrixLinCombd ( 2, spdimen, spdimen, &q[0], 1.0-u,
                         spdimen, &q[2*spdimen], u, spdimen, &q[0] );
    pkn_MatrixMDifferenced ( 1, spdimen, 0, &q[spdimen], 0, &q[0],
                             (double)degreev, 0, dv );
    pkn_MatrixLinCombd ( 1, spdimen, 0, &q[0], 1.0-v, 0, &q[spdimen], v,
                         0, p );
  }
  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*mbs_BCHornerDerPd*/

boolean mbs_BCHornerDerP3Rd ( int degreeu, int degreev, const point4d *ctlpoints,
                              double u, double v,
                              point3d *p, vector3d *du, vector3d *dv )
{
  point4d hp, hdu, hdv;

  if ( mbs_BCHornerDerP4d ( degreeu, degreev, ctlpoints, u, v, &hp, &hdu, &hdv ) ) {
    Point4to3d ( &hp, p );
    memcpy ( du, &hdu, sizeof(vector3d) );
    AddVector3Md ( du, p, -hdu.w, du );
    MultVector3d ( 1.0/hp.w, du, du );
    memcpy ( dv, &hdv, sizeof(vector3d) );
    AddVector3Md ( dv, p, -hdv.w, dv );
    MultVector3d ( 1.0/hp.w, dv, dv );
    return true;
  }
  else
    return false;
} /*mbs_BCHornerDerP3Rd*/

boolean mbs_BCHornerNvP3d ( int degreeu, int degreev, const point3d *ctlpoints,
                            double u, double v,
                            point3d *p, vector3d *nv )
{
  vector3d du, dv;

  if ( mbs_BCHornerDerP3d ( degreeu, degreev, ctlpoints, u, v, p, &du, &dv ) ) {
    CrossProduct3d ( &du, &dv, nv );
    NormalizeVector3d ( nv );
    return true;
  }
  else
    return false;
} /*mbs_BCHornerNvP3d*/

boolean mbs_BCHornerNvP3Rd ( int degreeu, int degreev, const point4d *ctlpoints,
                             double u, double v,
                             point3d *p, vector3d *nv )
{
  point4d  P;
  vector4d Du, Dv;

  if ( mbs_BCHornerDerP4d ( degreeu, degreev, ctlpoints, u, v, &P, &Du, &Dv ) ) {
    Point4to3d ( &P, p );
    CrossProduct4P3d ( &P, &Du, &Dv, nv );
    NormalizeVector3d ( nv );
    return true;
  }
  else
    return false;
} /*mbs_BCHornerNvP3Rd*/

