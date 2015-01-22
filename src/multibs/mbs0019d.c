
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

boolean mbs_BCHornerNvP3Rd ( int degreeu, int degreev, const point4d *ctlpoints,
                             double u, double v,
                             point3d *p, vector3d *nv )
{
  double  *workspace;
  boolean result;

  workspace = pkv_GetScratchMemd ( (6+2*(degreev+1))*4 );
  if ( workspace ) {
    result = _mbs_BCHornerNvP3Rd ( degreeu, degreev, ctlpoints, u, v,
                                   p, nv, workspace );
    pkv_SetScratchMemTop ( workspace );
    return result;
  }
  else
    return false;
} /*_mbs_BCHornerNvP3Rd*/

