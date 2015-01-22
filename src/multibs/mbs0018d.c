
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

boolean _mbs_BCHornerNvP3Rd ( int degreeu, int degreev, const point4d *ctlpoints,
                              double u, double v,
                              point3d *p, vector3d *nv, double *workspace )
{
  point4d  P;
  vector4d Du, Dv;

  if ( _mbs_BCHornerDerP4d ( degreeu, degreev, ctlpoints,
                             u, v, &P, &Du, &Dv, workspace ) ) {
    Point4to3d ( &P, p );
    CrossProduct4P3d ( &P, &Du, &Dv, nv );
    NormalizeVector3d ( nv );
    return true;
  }
  else
    return false;
} /*_mbs_BCHornerNvP3Rd*/

