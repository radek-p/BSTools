
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

boolean _mbs_PrincipalDirectionsBP3Rf ( int degreeu, int degreev,
                                        const point4f *ctlpoints, 
                                        float u, float v,
                                        float *k1, vector2f *v1,
                                        float *k2, vector2f *v2,
                                        float *workspace )
{
  float first[3], second[3];
  float a, b, c;

  if ( !_mbs_FundFormsBP3Rf ( degreeu, degreev, ctlpoints, u, v,
                              first, second, workspace ) )
    return false;
  a = first[0]*first[2] - first[1]*first[1];
  b = 2.0*second[1]*first[1] - second[0]*first[2] - second[2]*first[0];
  c = second[0]*second[2] - second[1]*second[1];

  pkn_SolveSqEqf ( 0.5*b/a, c/a, k1, k2 );

  a = second[0] - *k1*first[0];
  b = second[1] - *k1*first[1];
  c = second[2] - *k1*first[2];
  if ( fabs(a) > fabs(c) )
    SetVector2f ( v1, -b, a );
  else
    SetVector2f ( v1, c, -b );
  NormalizeVector2f ( v1 );

  a = second[0] - *k2*first[0];
  b = second[1] - *k2*first[1];
  c = second[2] - *k2*first[2];
  if ( fabs(a) > fabs(c) )
    SetVector2f ( v2, -b, a );
  else
    SetVector2f ( v2, c, -b );
  NormalizeVector2f ( v2 );
  return true;
} /*_mbs_PrincipalDirectionsBP3Rf*/
