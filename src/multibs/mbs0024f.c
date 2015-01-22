
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

boolean _mbs_BCHornerDer2C3Rf ( int degree, const point4f *ctlpoints, float t,
                                point3f *p, vector3f *d1, vector3f *d2,
                                float *workspace )
{
  point4f hp, hd1, hd2;

  if ( !_mbs_BCHornerDer2C4f ( degree, ctlpoints, t, &hp, &hd1, &hd2,
                               workspace ) )
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
} /*_mbs_BCHornerDer2C3Rf*/

