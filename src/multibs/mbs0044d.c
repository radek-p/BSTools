
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

boolean _mbs_PrincipalDirectionsBP3d ( int degreeu, int degreev,
                                       const point3d *ctlpoints, 
                                       double u, double v,
                                       double *k1, vector2d *v1,
                                       double *k2, vector2d *v2,
                                       double *workspace )
{
  double first[3], second[3];
  double a, b, c;

  if ( !_mbs_FundFormsBP3d ( degreeu, degreev, ctlpoints, u, v,
                             first, second, workspace ) )
    return false;
  a = first[0]*first[2] - first[1]*first[1];
  b = 2.0*second[1]*first[1] - second[0]*first[2] - second[2]*first[0];
  c = second[0]*second[2] - second[1]*second[1];

  pkn_SolveSqEqd ( 0.5*b/a, c/a, k1, k2 );

  a = second[0] - *k1*first[0];
  b = second[1] - *k1*first[1];
  c = second[2] - *k1*first[2];
  if ( fabs(a) > fabs(c) )
    SetVector2d ( v1, -b, a );
  else
    SetVector2d ( v1, c, -b );
  NormalizeVector2d ( v1 );

  a = second[0] - *k2*first[0];
  b = second[1] - *k2*first[1];
  c = second[2] - *k2*first[2];
  if ( fabs(a) > fabs(c) )
    SetVector2d ( v2, -b, a );
  else
    SetVector2d ( v2, c, -b );
  NormalizeVector2d ( v2 );
  return true;
} /*_mbs_PrincipalDirectionsBP3d*/

