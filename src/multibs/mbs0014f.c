
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2005, 2015                            */
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

boolean _mbs_BCHornerDerP3Rf ( int degreeu, int degreev, const point4f *ctlpoints,
                               float u, float v,
                               point3f *p, vector3f *du, vector3f *dv,
                               float *workspace )
{
  point4f hp, hdu, hdv;

  if ( _mbs_BCHornerDerP4f ( degreeu, degreev, ctlpoints, u, v,
                             &hp, &hdu, &hdv, workspace ) ) {
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
} /*_mbs_BCHornerDerP3Rf*/

