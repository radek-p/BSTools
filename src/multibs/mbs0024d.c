
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

boolean _mbs_BCHornerDer2C3Rd ( int degree, const point4d *ctlpoints, double t,
                                point3d *p, vector3d *d1, vector3d *d2,
                                double *workspace )
{
  point4d hp, hd1, hd2;

  if ( !_mbs_BCHornerDer2C4d ( degree, ctlpoints, t, &hp, &hd1, &hd2,
                               workspace ) )
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
} /*_mbs_BCHornerDer2C3Rd*/

