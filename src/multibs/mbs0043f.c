
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

boolean mbs_GMCurvaturesBP3f ( int degreeu, int degreev, const point3f *ctlpoints, 
                               float u, float v,
                               float *gaussian, float *mean )
{
  void    *sp;
  float  *workspace;
  boolean result;

  sp = workspace = pkv_GetScratchMemf ( (36+3*(degreev+1))*3 );
  if ( !workspace )
    return false;
  result = _mbs_GMCurvaturesBP3f ( degreeu, degreev, ctlpoints, u, v,
                                   gaussian, mean, workspace );
  pkv_SetScratchMemTop ( sp );
  return result;
} /*mbs_GMCurvaturesBP3f*/

