
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2005, 2015                            */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include "pkvaria.h"
#include "pknum.h"
#include "pkgeom.h"
#include "multibs.h"

boolean mbs_TabBezCurveDer3f ( int spdimen, int degree, const float *cp,
                               int nkn, const float *kn,
                               int ppitch,
                               float *p, float *dp, float *ddp, float *dddp )
{
  void    *sp;
  float   *workspace;
  boolean result;

  sp = workspace = pkv_GetScratchMemf ( 4*spdimen );
  if ( !workspace )
    return false;
  result = _mbs_TabBezCurveDer3f ( spdimen, degree, cp, nkn, kn,
                                   ppitch, p, dp, ddp, dddp, workspace );
  pkv_SetScratchMemTop ( sp );
  return result;
} /*mbs_TabBezCurveDer3f*/

