
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

boolean _mbs_BCHornerDer2P3Rf ( int degreeu, int degreev, const point4f *ctlpoints,
                                float u, float v,
                                point3f *p, vector3f *du, vector3f *dv,
                                vector3f *duu, vector3f *duv, vector3f *dvv,
                                float *workspace )
{
  point4f hp, hdu, hdv, hduu, hduv, hdvv;
  float   iw;

  if ( !_mbs_BCHornerDer2P4f ( degreeu, degreev, ctlpoints, u, v,
                               &hp, &hdu, &hdv, &hduu, &hduv, &hdvv,
                               workspace ) )
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
} /*_mbs_BCHornerDer2P3Rf*/

