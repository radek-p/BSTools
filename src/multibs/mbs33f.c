
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

boolean mbs_multiBCHornerDer2f ( int degree, int ncurves, int spdimen, int pitch,
                                 const float *ctlpoints,
                                 float t, float *p, float *d1, float *d2 )
{
  int   i;
  float *a, s;

  if ( degree <= 1 ) {
    mbs_multiBCHornerDerf ( degree, ncurves, spdimen, pitch, ctlpoints,
                            t, p, d1 );
    memset ( d2, 0, ncurves*spdimen*sizeof(float) );
    return true;
  }
  else if ( (a = pkv_GetScratchMemf ( 3*spdimen )) ) {
    s = (float)(1.0-t);
    for ( i = 0;
          i < ncurves;
          i++, ctlpoints += pitch,
          p += spdimen, d1 += spdimen, d2 += spdimen ) {
      mbs_multiBCHornerf ( degree-2, 3, spdimen, spdimen, ctlpoints, t, a );
      pkn_AddMatrixf ( 1, spdimen, 0, a, 0, &a[2*spdimen], 0, d2 );
      pkn_AddMatrixMf ( 1, spdimen, 0, d2, 0, &a[spdimen], -2.0, 0, d2 );
      pkn_MultMatrixNumf ( 1, spdimen, 0, d2,
                          (float)(degree*(degree-1)), 0, d2 );
      pkn_MatrixLinCombf ( 1, spdimen, 0, a, s, 0, &a[spdimen], t, 0, a );
      pkn_MatrixLinCombf ( 1, spdimen, 0, &a[spdimen], s,
                           0, &a[2*spdimen], t, 0, &a[spdimen] );
      pkn_MatrixMDifferencef ( 1, spdimen, 0, &a[spdimen], 0, a,
                               (float)degree, 0, d1 );
      pkn_MatrixLinCombf ( 1, spdimen, 0, a, s, 0, &a[spdimen], t, 0, p );
    }
    pkv_FreeScratchMemf ( 3*spdimen );
    return true;
  }
  else
    return false;
} /*mbs_multiBCHornerDer2f*/

boolean mbs_BCHornerDer2C2Rf ( int degree, const point3f *ctlpoints, float t,
                               point2f *p, vector2f *d1, vector2f *d2 )
{
  point3f hp, hd1, hd2;

  if ( !mbs_BCHornerDer2C3f ( degree, ctlpoints, t, &hp, &hd1, &hd2 ) )
    return false;
  Point3to2f ( &hp, p );
  memcpy ( d1, &hd1, sizeof(vector2f) );
  AddVector2Mf ( d1, p, -hd1.z, d1 );
  MultVector2f ( 1.0/hp.z, d1, d1 );
  memcpy ( d2, &hd2, sizeof(vector2f) );
  AddVector2Mf ( d2, d1, -2.0*hd1.z, d2 );
  AddVector2Mf ( d2, p, -hd2.z, d2 );
  MultVector2f ( 1.0/hp.z, d2, d2 );
  return true;
} /*mbs_BCHornerDer2C2Rf*/

boolean mbs_BCHornerDer2C3Rf ( int degree, const point4f *ctlpoints, float t,
                               point3f *p, vector3f *d1, vector3f *d2 )
{
  point4f hp, hd1, hd2;

  if ( !mbs_BCHornerDer2C4f ( degree, ctlpoints, t, &hp, &hd1, &hd2 ) )
    return false;
  Point4to3f ( &hp, p );
  memcpy ( d1, &hd1, sizeof(vector3f) );
  AddVector3Mf ( d1, p, -hd1.w, d1 );
  MultVector3f ( 1.0/hp.w, d1, d1 );
  memcpy ( d2, &hd2, sizeof(vector3f) );
  AddVector3Mf ( d2, d1, -2.0*hd1.w, d2 );
  AddVector3Mf ( d2, p, -hd2.w, d2 );
  MultVector3f ( 1.0/hp.w, d2, d2 );
  return true;
} /*mbs_BCHornerDer2C3Rf*/

boolean mbs_BCHornerDer2Pf ( int degreeu, int degreev, int spdimen,   
                             const float *ctlpoints,   
                             float u, float v,   
                             float *p, float *du, float *dv,
                             float *duu, float *duv, float *dvv )
{
  void  *sp;
  float *q, *r;
  int   n, m, i, k, pitch;
  int   size;

  sp = pkv_GetScratchMemTop ();
  pitch = (degreev+1)*spdimen;
  if ( (r = pkv_GetScratchMemf ( size = 36*spdimen+3*pitch )) ) {
    q = &r[36*spdimen];
    if ( degreeu <= 2 ) {
      n = degreeu+1;
      memcpy ( q, ctlpoints, n*pitch*sizeof(float) );
    }
    else {
      n = 3;
      if ( !mbs_multiBCHornerf ( degreeu-2, 3, pitch, pitch, ctlpoints, u, q ) )
        goto failure;
    }
    if ( degreev <= 2 ) {
      m = degreev+1;
      pkv_Selectf ( n, m*spdimen, pitch, 3*spdimen, q, r );
    }
    else {
      m = 3;
      for ( i = 0; i < n; i++ )
        if ( !mbs_multiBCHornerf ( degreev-2, 3, spdimen, spdimen,
                                   &q[i*pitch], v, &r[i*3*spdimen] ) )
          goto failure;
    }

                    /* now the last steps of the de Casteljau algorithm */
    switch ( n ) {
case 3:
      pkn_MatrixLinCombf ( 2, m*spdimen, 3*spdimen, r, 1.0-u,
          3*spdimen, &r[3*spdimen], u, 3*spdimen, &r[9*spdimen] );
      pkn_MatrixLinCombf ( 1, m*spdimen, 0, &r[9*spdimen], 1.0-u,
          0, &r[12*spdimen], u, 0, &r[15*spdimen] );
      k = 6;
      break;
case 2:
      memcpy ( &r[9*spdimen], r, 6*spdimen*sizeof(float) );
      pkn_MatrixLinCombf ( 1, m*spdimen, 0, &r[9*spdimen], 1.0-u,
          0, &r[12*spdimen], u, 0, &r[15*spdimen] );
      k = 3;
      break;
case 1:
      memcpy ( &r[15*spdimen], r, 3*spdimen*sizeof(float) );
      k = 1;
      break;
default: k = 0;  /* to suppress a warning */
    }
    switch ( m ) {
case 3:
      pkn_MatrixLinCombf ( k, spdimen, 3*spdimen,
          &r[3*(6-k)*spdimen], 1.0-v,
          3*spdimen, &r[(3*(6-k)+1)*spdimen], v,
          2*spdimen, &r[(2*(6-k)+18)*spdimen] );
      pkn_MatrixLinCombf ( k, spdimen, 3*spdimen,
          &r[(3*(6-k)+1)*spdimen], 1.0-v,
          3*spdimen, &r[(3*(6-k)+2)*spdimen], v,
          2*spdimen, &r[(2*(6-k)+19)*spdimen] );
      pkn_MatrixLinCombf ( k, spdimen, 2*spdimen,
          &r[(2*(6-k)+18)*spdimen], 1.0-v,
          2*spdimen, &r[(2*(6-k)+19)*spdimen], v,
          spdimen, &r[(6-k+30)*spdimen] );
      break;
case 2:
      pkv_Selectf ( k, 2*spdimen, 3*spdimen, 2*spdimen,
                    &r[(3*(6-k))*spdimen], &r[(2*(6-k)+18)*spdimen] );
      pkn_MatrixLinCombf ( k, spdimen, 2*spdimen,
          &r[(2*(6-k)+18)*spdimen], 1.0-v,
          2*spdimen, &r[(2*(6-k)+19)*spdimen], v,
          spdimen, &r[(6-k+30)*spdimen] );
      break;
case 1:
      pkv_Selectf ( k, spdimen, 3*spdimen, spdimen,
                    &r[(3*(6-k))*spdimen], &r[(6-k+30)*spdimen] );
      break;
    }
                    /* compute the patch point and derivatives */
    memcpy ( p, &r[35*spdimen], spdimen*sizeof(float) );
    if ( degreeu > 0 ) {
      pkn_MatrixMDifferencef ( 1, spdimen, 0, &r[34*spdimen],
                               0, &r[33*spdimen], (float)degreeu, 0, du );
      if ( degreeu > 1 ) {
        pkn_AddMatrixf ( 1, spdimen, 0, &r[30*spdimen], 0, &r[32*spdimen],
                         0, duu );
        pkn_AddMatrixMf ( 1, spdimen, 0, duu, 0, &r[31*spdimen], -2.0,
                          0, duu );
        pkn_MultMatrixNumf ( 1, spdimen, 0, duu, (float)(degreeu*(degreeu-1)),
                             0, duu );
      }
      else
        memset ( duu, 0, spdimen*sizeof(float) );
      if ( degreev > 0 ) {
        pkn_AddMatrixf ( 1, spdimen, 0, &r[24*spdimen], 0, &r[27*spdimen],
                         0, duv );
        pkn_SubtractMatrixf ( 1, spdimen, 0, duv, 0, &r[25*spdimen], 0, duv );
        pkn_SubtractMatrixf ( 1, spdimen, 0, duv, 0, &r[26*spdimen], 0, duv );
        pkn_MultMatrixNumf ( 1, spdimen, 0, duv, (float)(degreeu*degreev),
                             0, duv );
      }
      else
        memset ( duv, 0, spdimen*sizeof(float) );
    }
    else {
      memset ( du, 0, spdimen*sizeof(float) );
      memset ( duu, 0, spdimen*sizeof(float) );
      memset ( duv, 0, spdimen*sizeof(float) );
    }
    if ( degreev > 0 ) {
      pkn_MatrixMDifferencef ( 1, spdimen, 0, &r[29*spdimen],
                               0, &r[28*spdimen], (float)degreev, 0, dv );
      if ( degreev > 1 ) {
        pkn_AddMatrixf ( 1, spdimen, 0, &r[15*spdimen], 0, &r[17*spdimen],
                         0, dvv );
        pkn_AddMatrixMf ( 1, spdimen, 0, dvv, 0, &r[16*spdimen], -2.0,
                          0, dvv );
        pkn_MultMatrixNumf ( 1, spdimen, 0, dvv, (float)(degreev*(degreev-1)),
                             0, dvv );
      }
      else
        memset ( dvv, 0, spdimen*sizeof(float) );
    }
    else {
      memset ( dv, 0, spdimen*sizeof(float) );
      memset ( dvv, 0, spdimen*sizeof(float) );
    }
  }
  else
    goto failure;

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*mbs_BCHornerDer2Pf*/

boolean mbs_BCHornerDer2P3Rf ( int degreeu, int degreev, const point4f *ctlpoints,
                               float u, float v,
                               point3f *p, vector3f *du, vector3f *dv,
                               vector3f *duu, vector3f *duv, vector3f *dvv )
{
  point4f hp, hdu, hdv, hduu, hduv, hdvv;
  float   iw;

  if ( !mbs_BCHornerDer2P4f ( degreeu, degreev, ctlpoints, u, v,
                              &hp, &hdu, &hdv, &hduu, &hduv, &hdvv ) )
    return false;
  Point4to3f ( &hp, p );

  iw = (float)(1.0/hp.w);
  memcpy ( du, &hdu, sizeof(vector3f) );
  AddVector3Mf ( du, p, -hdu.w, du );
  MultVector3f ( iw, du, du );
  memcpy ( dv, &hdv, sizeof(vector3f) );
  AddVector3Mf ( dv, p, -hdv.w, dv );
  MultVector3f ( iw, dv, dv );

  memcpy ( duu, &hduu, sizeof(vector3f) );
  AddVector3Mf ( duu, du, -2.0*hdu.w, duu );
  AddVector3Mf ( duu, p, -hduu.w, duu );
  MultVector3f ( iw, duu, duu );
  memcpy ( duv, &hduv, sizeof(vector3f) );
  AddVector3Mf ( duv, du, -hdu.w, duv );
  AddVector3Mf ( duv, dv, -hdv.w, duv );
  AddVector3Mf ( duv, p, -hduv.w, duv );
  MultVector3f ( iw, duv, duv );
  memcpy ( dvv, &hdvv, sizeof(vector3f) );
  AddVector3Mf ( dvv, dv, -2.0*hdv.w, dvv );
  AddVector3Mf ( dvv, p, -hdvv.w, dvv );
  MultVector3f ( iw, dvv, dvv );
  return true;
} /*mbs_BCHornerDer2P3Rf*/

