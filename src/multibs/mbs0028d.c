
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

boolean _mbs_BCHornerDer2P3Rd ( int degreeu, int degreev, const point4d *ctlpoints,
                                double u, double v,
                                point3d *p, vector3d *du, vector3d *dv,
                                vector3d *duu, vector3d *duv, vector3d *dvv,
                                double *workspace )
{
  point4d hp, hdu, hdv, hduu, hduv, hdvv;
  double  iw;

  if ( !_mbs_BCHornerDer2P4d ( degreeu, degreev, ctlpoints, u, v,
                               &hp, &hdu, &hdv, &hduu, &hduv, &hdvv,
                               workspace ) )
    return false;
  Point4to3d ( &hp, p );

  iw = (double)(1.0/hp.w);
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
} /*_mbs_BCHornerDer2P3Rd*/

