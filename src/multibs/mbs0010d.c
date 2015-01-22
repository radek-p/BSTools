
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

boolean mbs_BCHornerDerC2Rd ( int degree, const point3d *ctlpoints, double t,
                              point2d *p, vector2d *d )
{
  point3d hp, hd;

  if ( !mbs_BCHornerDerC3d ( degree, ctlpoints, t, &hp, &hd ) )
    return false;
  Point3to2d ( &hp, p );
  memcpy ( d, &hd, sizeof(vector2d) );
  AddVector2Md ( d, p, -hd.z, d );
  MultVector2d ( 1.0/hp.z, d, d );
  return true;
} /*mbs_BCHornerDerC2Rd*/

