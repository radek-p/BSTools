
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

boolean _mbs_BCHornerNvP3f ( int degreeu, int degreev, const point3f *ctlpoints,
                             float u, float v,
                             point3f *p, vector3f *nv, float *workspace )
{
  vector3f du, dv;

  if ( _mbs_BCHornerDerP3f ( degreeu, degreev, ctlpoints, u, v,
                             p, &du, &dv, workspace ) ) {
    CrossProduct3f ( &du, &dv, nv );
    NormalizeVector3f ( nv );
    return true;
  }
  else
    return false;
} /*mbs_BCHornerNvP3f*/

