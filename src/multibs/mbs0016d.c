
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

boolean _mbs_BCHornerNvP3d ( int degreeu, int degreev, const point3d *ctlpoints,
                             double u, double v,
                             point3d *p, vector3d *nv, double *workspace )
{
  vector3d du, dv;

  if ( _mbs_BCHornerDerP3d ( degreeu, degreev, ctlpoints, u, v,
                             p, &du, &dv, workspace ) ) {
    CrossProduct3d ( &du, &dv, nv );
    NormalizeVector3d ( nv );
    return true;
  }
  else
    return false;
} /*mbs_BCHornerNvP3d*/

