
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

boolean _mbs_BCHornerDerP3Rd ( int degreeu, int degreev, const point4d *ctlpoints,
                               double u, double v,
                               point3d *p, vector3d *du, vector3d *dv,
                               double *workspace )
{
  point4d hp, hdu, hdv;

  if ( _mbs_BCHornerDerP4d ( degreeu, degreev, ctlpoints, u, v,
                             &hp, &hdu, &hdv, workspace ) ) {
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
} /*_mbs_BCHornerDerP3Rd*/

