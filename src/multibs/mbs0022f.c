
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

boolean _mbs_BCHornerDer2C2Rf ( int degree, const point3f *ctlpoints, float t,
                                point2f *p, vector2f *d1, vector2f *d2,
                                float *workspace )
{
  point3f hp, hd1, hd2;

  if ( !_mbs_BCHornerDer2C3f ( degree, ctlpoints, t, &hp, &hd1, &hd2,
                               workspace ) )
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
} /*_mbs_BCHornerDer2C2Rf*/

