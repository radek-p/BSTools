
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2015                                  */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

#include <math.h>
#include <string.h>
#include <stdio.h>  /* for testing */

#include "pkvaria.h"
#include "pkgeom.h"
#include "multibs.h"

boolean _mbs_BCHornerP3Rf ( int degreeu, int degreev, const point4f *ctlpoints,
                            float u, float v, point3f *ppoint,
                            float *workspace )
{
  point4f auxp;

  if ( _mbs_BCHornerP4f ( degreeu, degreev, ctlpoints, u, v, &auxp,
                          workspace ) ) {
    Point4to3f ( &auxp, ppoint );
    return true;
  }
  else
    return false;
} /*_mbs_BCHornerP3Rf*/

