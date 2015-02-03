
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

boolean _mbs_GMCurvaturesBP3Rd ( int degreeu, int degreev, const point4d *ctlpoints, 
                             double u, double v,
                             double *gaussian, double *mean, double *workspace )
{
  double first[3], second[3];
  double a, b, c;

  if ( !_mbs_FundFormsBP3Rd ( degreeu, degreev, ctlpoints, u, v,
                              first, second, workspace ) )
    return false;
  a = first[0]*first[2] - first[1]*first[1];
  b = 2.0*second[1]*first[1] - second[0]*first[2] - second[2]*first[0];
  c = second[0]*second[2] - second[1]*second[1];
  *gaussian = c/a;
  *mean = 0.5*b/a;
  return true;
} /*_mbs_GMCurvaturesBP3Rd*/

