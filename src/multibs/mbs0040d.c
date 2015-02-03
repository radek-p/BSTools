
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

boolean _mbs_FundFormsBP3d ( int degreeu, int degreev, const point3d *ctlpoints,
                             double u, double v,
                             double *firstform, double *secondform,
                             double *workspace )
{
  point3d  p;
  vector3d pu, pv, puu, puv, pvv, normal;

  if ( !_mbs_BCHornerDer2P3d ( degreeu, degreev, ctlpoints, u, v,
                               &p, &pu, &pv, &puu, &puv, &pvv, workspace ) )
    return false;

  firstform[0] = DotProduct3d ( &pu, &pu );
  firstform[1] = DotProduct3d ( &pu, &pv );
  firstform[2] = DotProduct3d ( &pv, &pv );

  CrossProduct3d ( &pu, &pv, &normal );
  NormalizeVector3d ( &normal );
  secondform[0] = DotProduct3d ( &normal, &puu );
  secondform[1] = DotProduct3d ( &normal, &puv );
  secondform[2] = DotProduct3d ( &normal, &pvv );
  return true;
} /*_mbs_FundFormsBP3d*/

