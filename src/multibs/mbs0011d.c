
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

boolean mbs_BCHornerDerC3Rd ( int degree, const point4d *ctlpoints, double t,
                              point3d *p, vector3d *d )
{
  point4d hp, hd;

  if ( !mbs_BCHornerDerC4d ( degree, ctlpoints, t, &hp, &hd ) )
    return false;
  Point4to3d ( &hp, p );
  memcpy ( d, &hd, sizeof(vector3d) );
  AddVector3Md ( d, p, -hd.w, d );
  MultVector3d ( 1.0/hp.w, d, d );
  return true;
} /*mbs_BCHornerDerC3Rd*/

