
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

boolean _mbs_FundFormsBP3f ( int degreeu, int degreev, const point3f *ctlpoints,
                             float u, float v,
                             float *firstform, float *secondform,
                             float *workspace )
{
  point3f  p;
  vector3f pu, pv, puu, puv, pvv, normal;

  if ( !_mbs_BCHornerDer2P3f ( degreeu, degreev, ctlpoints, u, v,
                               &p, &pu, &pv, &puu, &puv, &pvv, workspace ) )
    return false;

  firstform[0] = DotProduct3f ( &pu, &pu );
  firstform[1] = DotProduct3f ( &pu, &pv );
  firstform[2] = DotProduct3f ( &pv, &pv );

  CrossProduct3f ( &pu, &pv, &normal );
  NormalizeVector3f ( &normal );
  secondform[0] = DotProduct3f ( &normal, &puu );
  secondform[1] = DotProduct3f ( &normal, &puv );
  secondform[2] = DotProduct3f ( &normal, &pvv );
  return true;
} /*_mbs_FundFormsBP3f*/

