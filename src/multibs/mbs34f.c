
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2005, 2009                            */
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

void mbs_FundFormsBP3f ( int degreeu, int degreev, const point3f *ctlpoints,
                         float u, float v,
                         float *firstform, float *secondform )
{
  point3f  p;
  vector3f pu, pv, puu, puv, pvv, normal;

  mbs_BCHornerDer2P3f ( degreeu, degreev, ctlpoints, u, v,
                        &p, &pu, &pv, &puu, &puv, &pvv );

  firstform[0] = (float)DotProduct3f ( &pu, &pu );
  firstform[1] = (float)DotProduct3f ( &pu, &pv );
  firstform[2] = (float)DotProduct3f ( &pv, &pv );

  CrossProduct3f ( &pu, &pv, &normal );
  NormalizeVector3f ( &normal );
  secondform[0] = (float)DotProduct3f ( &normal, &puu );
  secondform[1] = (float)DotProduct3f ( &normal, &puv );
  secondform[2] = (float)DotProduct3f ( &normal, &pvv );
} /*mbs_FundFormsBP3f*/


/* computing Gaussian and mean curvature */

void mbs_GMCurvaturesBP3f ( int degreeu, int degreev, const point3f *ctlpoints, 
                            float u, float v,
                            float *gaussian, float *mean )
{
  float first[3], second[3];
  float a, b, c;

  mbs_FundFormsBP3f ( degreeu, degreev, ctlpoints, u, v, first, second );
  a = first[0]*first[2] - first[1]*first[1];
  b = (float)2.0*second[1]*first[1] - second[0]*first[2] - second[2]*first[0];
  c = second[0]*second[2] - second[1]*second[1];
  *gaussian = c/a;
  *mean = (float)0.5*b/a;
} /*mbs_GMCurvaturesBP3f*/


/* computing principal curvatures and directions */

void mbs_PrincipalDirectionsBP3f ( int degreeu, int degreev,
                                   const point3f *ctlpoints, 
                                   float u, float v,
                                   float *k1, vector2f *v1,
                                   float *k2, vector2f *v2 )
{
  float first[3], second[3];
  float a, b, c;

  mbs_FundFormsBP3f ( degreeu, degreev, ctlpoints, u, v, first, second );
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
} /*mbs_PrincipalDirectionsBP3f*/

