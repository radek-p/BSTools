
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

boolean mbs_multiBCHornerDer2d ( int degree, int ncurves, int spdimen, int pitch,
                                 const double *ctlpoints,
                                 double t, double *p, double *d1, double *d2 )
{
  void  *sp;
  int   i;
  double *a, s;

  sp = pkv_GetScratchMemTop ();
  if ( degree <= 1 ) {
    if ( !mbs_multiBCHornerDerd ( degree, ncurves, spdimen, pitch, ctlpoints,
                                  t, p, d1 ) )
      goto failure;
    memset ( d2, 0, ncurves*spdimen*sizeof(double) );
  }
  else if ( (a = pkv_GetScratchMemd ( 3*spdimen )) ) {
    s = 1.0-t;
    for ( i = 0;
          i < ncurves;
          i++, ctlpoints += pitch,
          p += spdimen, d1 += spdimen, d2 += spdimen ) {
      mbs_multiBCHornerd ( degree-2, 3, spdimen, spdimen, ctlpoints, t, a );
      pkn_AddMatrixd ( 1, spdimen, 0, a, 0, &a[2*spdimen], 0, d2 );
      pkn_AddMatrixMd ( 1, spdimen, 0, d2, 0, &a[spdimen], -2.0, 0, d2 );
      pkn_MultMatrixNumd ( 1, spdimen, 0, d2,
                          (double)(degree*(degree-1)), 0, d2 );
      pkn_MatrixLinCombd ( 1, spdimen, 0, a, s, 0, &a[spdimen], t, 0, a );
      pkn_MatrixLinCombd ( 1, spdimen, 0, &a[spdimen], s,
                           0, &a[2*spdimen], t, 0, &a[spdimen] );
      pkn_MatrixMDifferenced ( 1, spdimen, 0, &a[spdimen], 0, a,
                               (double)degree, 0, d1 );
      pkn_MatrixLinCombd ( 1, spdimen, 0, a, s, 0, &a[spdimen], t, 0, p );
    }
  }
  else
    goto failure;

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*mbs_multiBCHornerDer2d*/

boolean mbs_BCHornerDer2C2Rd ( int degree, const point3d *ctlpoints, double t,
                               point2d *p, vector2d *d1, vector2d *d2 )
{
  point3d hp, hd1, hd2;

  if ( !mbs_BCHornerDer2C3d ( degree, ctlpoints, t, &hp, &hd1, &hd2 ) )
    return false;
  Point3to2d ( &hp, p );
  memcpy ( d1, &hd1, sizeof(vector2d) );
  AddVector2Md ( d1, p, -hd1.z, d1 );
  MultVector2d ( 1.0/hp.z, d1, d1 );
  memcpy ( d2, &hd2, sizeof(vector2d) );
  AddVector2Md ( d2, d1, -2.0*hd1.z, d2 );
  AddVector2Md ( d2, p, -hd2.z, d2 );
  MultVector2d ( 1.0/hp.z, d2, d2 );
  return true;
} /*mbs_BCHornerDer2C2Rd*/

boolean mbs_BCHornerDer2C3Rd ( int degree, const point4d *ctlpoints, double t,
                               point3d *p, vector3d *d1, vector3d *d2 )
{
  point4d hp, hd1, hd2;

  if ( !mbs_BCHornerDer2C4d ( degree, ctlpoints, t, &hp, &hd1, &hd2 ) )
    return false;
  Point4to3d ( &hp, p );
  memcpy ( d1, &hd1, sizeof(vector3d) );
  AddVector3Md ( d1, p, -hd1.w, d1 );
  MultVector3d ( 1.0/hp.w, d1, d1 );
  memcpy ( d2, &hd2, sizeof(vector3d) );
  AddVector3Md ( d2, d1, -2.0*hd1.w, d2 );
  AddVector3Md ( d2, p, -hd2.w, d2 );
  MultVector3d ( 1.0/hp.w, d2, d2 );
  return true;
} /*mbs_BCHornerDer2C3Rd*/

boolean mbs_BCHornerDer2Pd ( int degreeu, int degreev, int spdimen,   
                             const double *ctlpoints,   
                             double u, double v,   
                             double *p, double *du, double *dv,
                             double *duu, double *duv, double *dvv )
{
  void   *sp;
  double *q, *r;
  int    n, m, i, k, pitch;
  int    size;

  sp = pkv_GetScratchMemTop ();
  pitch = (degreev+1)*spdimen;
  if ( (r = pkv_GetScratchMemd ( size = 36*spdimen+3*pitch )) ) {
    q = &r[36*spdimen];
    if ( degreeu <= 2 ) {
      n = degreeu+1;
      memcpy ( q, ctlpoints, n*pitch*sizeof(double) );
    }
    else {
      n = 3;
      if ( !mbs_multiBCHornerd ( degreeu-2, 3, pitch, pitch, ctlpoints, u, q ) )
        goto failure;
    }
    if ( degreev <= 2 ) {
      m = degreev+1;
      pkv_Selectd ( n, m*spdimen, pitch, 3*spdimen, q, r );
    }
    else {
      m = 3;
      for ( i = 0; i < n; i++ )
        if ( !mbs_multiBCHornerd ( degreev-2, 3, spdimen, spdimen,
                                   &q[i*pitch], v, &r[i*3*spdimen] ) )
          goto failure;
    }

                    /* now the last steps of the de Casteljau algorithm */
    switch ( n ) {
case 3:
      pkn_MatrixLinCombd ( 2, m*spdimen, 3*spdimen, r, 1.0-u,
          3*spdimen, &r[3*spdimen], u, 3*spdimen, &r[9*spdimen] );
      pkn_MatrixLinCombd ( 1, m*spdimen, 0, &r[9*spdimen], 1.0-u,
          0, &r[12*spdimen], u, 0, &r[15*spdimen] );
      k = 6;
      break;
case 2:
      memcpy ( &r[9*spdimen], r, 6*spdimen*sizeof(double) );
      pkn_MatrixLinCombd ( 1, m*spdimen, 0, &r[9*spdimen], 1.0-u,
          0, &r[12*spdimen], u, 0, &r[15*spdimen] );
      k = 3;
      break;
case 1:
      memcpy ( &r[15*spdimen], r, 3*spdimen*sizeof(double) );
      k = 1;
      break;
default: k = 0;  /* to suppress a warning */
    }
    switch ( m ) {
case 3:
      pkn_MatrixLinCombd ( k, spdimen, 3*spdimen,
          &r[3*(6-k)*spdimen], 1.0-v,
          3*spdimen, &r[(3*(6-k)+1)*spdimen], v,
          2*spdimen, &r[(2*(6-k)+18)*spdimen] );
      pkn_MatrixLinCombd ( k, spdimen, 3*spdimen,
          &r[(3*(6-k)+1)*spdimen], 1.0-v,
          3*spdimen, &r[(3*(6-k)+2)*spdimen], v,
          2*spdimen, &r[(2*(6-k)+19)*spdimen] );
      pkn_MatrixLinCombd ( k, spdimen, 2*spdimen,
          &r[(2*(6-k)+18)*spdimen], 1.0-v,
          2*spdimen, &r[(2*(6-k)+19)*spdimen], v,
          spdimen, &r[(6-k+30)*spdimen] );
      break;
case 2:
      pkv_Selectd ( k, 2*spdimen, 3*spdimen, 2*spdimen,
                    &r[(3*(6-k))*spdimen], &r[(2*(6-k)+18)*spdimen] );
      pkn_MatrixLinCombd ( k, spdimen, 2*spdimen,
          &r[(2*(6-k)+18)*spdimen], 1.0-v,
          2*spdimen, &r[(2*(6-k)+19)*spdimen], v,
          spdimen, &r[(6-k+30)*spdimen] );
      break;
case 1:
      pkv_Selectd ( k, spdimen, 3*spdimen, spdimen,
                    &r[(3*(6-k))*spdimen], &r[(6-k+30)*spdimen] );
      break;
    }
                    /* compute the patch point and derivatives */
    memcpy ( p, &r[35*spdimen], spdimen*sizeof(double) );
    if ( degreeu > 0 ) {
      pkn_MatrixMDifferenced ( 1, spdimen, 0, &r[34*spdimen],
                               0, &r[33*spdimen], (double)degreeu, 0, du );
      if ( degreeu > 1 ) {
        pkn_AddMatrixd ( 1, spdimen, 0, &r[30*spdimen], 0, &r[32*spdimen],
                         0, duu );
        pkn_AddMatrixMd ( 1, spdimen, 0, duu, 0, &r[31*spdimen], -2.0,
                          0, duu );
        pkn_MultMatrixNumd ( 1, spdimen, 0, duu, (double)(degreeu*(degreeu-1)),
                             0, duu );
      }
      else
        memset ( duu, 0, spdimen*sizeof(double) );
      if ( degreev > 0 ) {
        pkn_AddMatrixd ( 1, spdimen, 0, &r[24*spdimen], 0, &r[27*spdimen],
                         0, duv );
        pkn_SubtractMatrixd ( 1, spdimen, 0, duv, 0, &r[25*spdimen], 0, duv );
        pkn_SubtractMatrixd ( 1, spdimen, 0, duv, 0, &r[26*spdimen], 0, duv );
        pkn_MultMatrixNumd ( 1, spdimen, 0, duv, (double)(degreeu*degreev),
                             0, duv );
      }
      else
        memset ( duv, 0, spdimen*sizeof(double) );
    }
    else {
      memset ( du, 0, spdimen*sizeof(double) );
      memset ( duu, 0, spdimen*sizeof(double) );
      memset ( duv, 0, spdimen*sizeof(double) );
    }
    if ( degreev > 0 ) {
      pkn_MatrixMDifferenced ( 1, spdimen, 0, &r[29*spdimen],
                               0, &r[28*spdimen], (double)degreev, 0, dv );
      if ( degreev > 1 ) {
        pkn_AddMatrixd ( 1, spdimen, 0, &r[15*spdimen], 0, &r[17*spdimen],
                         0, dvv );
        pkn_AddMatrixMd ( 1, spdimen, 0, dvv, 0, &r[16*spdimen], -2.0,
                          0, dvv );
        pkn_MultMatrixNumd ( 1, spdimen, 0, dvv, (double)(degreev*(degreev-1)),
                             0, dvv );
      }
      else
        memset ( dvv, 0, spdimen*sizeof(double) );
    }
    else {
      memset ( dv, 0, spdimen*sizeof(double) );
      memset ( dvv, 0, spdimen*sizeof(double) );
    }
  }
  else
    goto failure;
  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*mbs_BCHornerDer2Pd*/

boolean mbs_BCHornerDer2P3Rd ( int degreeu, int degreev, const point4d *ctlpoints,
                               double u, double v,
                               point3d *p, vector3d *du, vector3d *dv,
                               vector3d *duu, vector3d *duv, vector3d *dvv )
{
  point4d hp, hdu, hdv, hduu, hduv, hdvv;
  double   iw;

  if ( !mbs_BCHornerDer2P4d ( degreeu, degreev, ctlpoints, u, v,
                              &hp, &hdu, &hdv, &hduu, &hduv, &hdvv ) )
    return false;
  Point4to3d ( &hp, p );

  iw = 1.0/hp.w;
  memcpy ( du, &hdu, sizeof(vector3d) );
  AddVector3Md ( du, p, -hdu.w, du );
  MultVector3d ( iw, du, du );
  memcpy ( dv, &hdv, sizeof(vector3d) );
  AddVector3Md ( dv, p, -hdv.w, dv );
  MultVector3d ( iw, dv, dv );

  memcpy ( duu, &hduu, sizeof(vector3d) );
  AddVector3Md ( duu, du, -2.0*hdu.w, duu );
  AddVector3Md ( duu, p, -hduu.w, duu );
  MultVector3d ( iw, duu, duu );
  memcpy ( duv, &hduv, sizeof(vector3d) );
  AddVector3Md ( duv, du, -hdu.w, duv );
  AddVector3Md ( duv, dv, -hdv.w, duv );
  AddVector3Md ( duv, p, -hduv.w, duv );
  MultVector3d ( iw, duv, duv );
  memcpy ( dvv, &hdvv, sizeof(vector3d) );
  AddVector3Md ( dvv, dv, -2.0*hdv.w, dvv );
  AddVector3Md ( dvv, p, -hdvv.w, dvv );
  MultVector3d ( iw, dvv, dvv );
  return true;
} /*mbs_BCHornerDer2P3Rd*/

