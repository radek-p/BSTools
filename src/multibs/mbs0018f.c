
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

boolean _mbs_BCHornerNvP3Rf ( int degreeu, int degreev, const point4f *ctlpoints,
                              float u, float v,
                              point3f *p, vector3f *nv, float *workspace )
{
  point4f  P;
  vector4f Du, Dv;

  if ( _mbs_BCHornerDerP4f ( degreeu, degreev, ctlpoints,
                             u, v, &P, &Du, &Dv, workspace ) ) {
    Point4to3f ( &P, p );
    CrossProduct4P3f ( &P, &Du, &Dv, nv );
    NormalizeVector3f ( nv );
    return true;
  }
  else
    return false;
} /*_mbs_BCHornerNvP3Rf*/

