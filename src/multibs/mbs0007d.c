
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

boolean _mbs_BCHornerP3Rd ( int degreeu, int degreev, const point4d *ctlpoints,
                            double u, double v, point3d *ppoint,
                            double *workspace )
{
  point4d auxp;

  if ( _mbs_BCHornerP4d ( degreeu, degreev, ctlpoints, u, v, &auxp,
                          workspace ) ) {
    Point4to3d ( &auxp, ppoint );
    return true;
  }
  else
    return false;
} /*_mbs_BCHornerP3Rd*/
