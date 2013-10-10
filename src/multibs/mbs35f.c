
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2005, 2013                            */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

#include <math.h>
#include <string.h>
#include <stdio.h>  /* for testing */

#include "pkvaria.h"
#include "pkgeom.h"
#include "multibs.h"


/* /////////////////////////////////////////// */
/* computing the matrices of the fundamental */
/* forms of B\'{e}zier patches */

boolean mbs_FundFormsBP3Rf ( int degreeu, int degreev, const point4f *ctlpoints,
                             float u, float v,
                             float *firstform, float *secondform )
{
  point3f  p;
  vector3f pu, pv, puu, puv, pvv, normal;

  if ( !mbs_BCHornerDer2P3Rf ( degreeu, degreev, ctlpoints, u, v,
                               &p, &pu, &pv, &puu, &puv, &pvv ) )
    return false;

  firstform[0] = (float)DotProduct3f ( &pu, &pu );
  firstform[1] = (float)DotProduct3f ( &pu, &pv );
  firstform[2] = (float)DotProduct3f ( &pv, &pv );

  CrossProduct3f ( &pu, &pv, &normal );
  NormalizeVector3f ( &normal );
  secondform[0] = (float)DotProduct3f ( &normal, &puu );
  secondform[1] = (float)DotProduct3f ( &normal, &puv );
  secondform[2] = (float)DotProduct3f ( &normal, &pvv );
  return true;
} /*mbs_FundFormsBP3Rf*/


/* computing Gaussian and mean curvature */

boolean mbs_GMCurvaturesBP3Rf ( int degreeu, int degreev, const point4f *ctlpoints, 
                                float u, float v,
                                float *gaussian, float *mean )
{
  float first[3], second[3];
  float a, b, c;

  if ( !mbs_FundFormsBP3Rf ( degreeu, degreev, ctlpoints, u, v, first, second ) )
    return false;
  a = first[0]*first[2] - first[1]*first[1];
  b = (float)2.0*second[1]*first[1] - second[0]*first[2] - second[2]*first[0];
  c = second[0]*second[2] - second[1]*second[1];
  *gaussian = c/a;
  *mean = (float)0.5*b/a;
  return true;
} /*mbs_GMCurvaturesBP3Rf*/


/* computing principal curvatures and directions */

boolean mbs_PrincipalDirectionsBP3Rf ( int degreeu, int degreev,
                                       const point4f *ctlpoints, 
                                       float u, float v,
                                       float *k1, vector2f *v1,
                                       float *k2, vector2f *v2 )
{
  float first[3], second[3];
  float a, b, c;

  if ( !mbs_FundFormsBP3Rf ( degreeu, degreev, ctlpoints, u, v, first, second ) )
    return false;
  a = first[0]*first[2] - first[1]*first[1];
  b = (float)2.0*second[1]*first[1] - second[0]*first[2] - second[2]*first[0];
  c = second[0]*second[2] - second[1]*second[1];

  pkn_SolveSqEqf ( (float)0.5*b/a, c/a, k1, k2 );

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
} /*mbs_PrincipalDirectionsBP3Rf*/

