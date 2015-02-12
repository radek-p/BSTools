
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

boolean mbs_TabBezCurveDer2d ( int spdimen, int degree, const double *cp,
                               int nkn, const double *kn,
                               int ppitch,
                               double *p, double *dp, double *ddp )
{
  void    *sp;
  double  *workspace;
  boolean result;

  sp = workspace = pkv_GetScratchMem ( 3*spdimen );
  if ( !workspace )
    return false;
  result = _mbs_TabBezCurveDer2d ( spdimen, degree, cp, nkn, kn,
                                   ppitch, p, dp, ddp, workspace );
  pkv_SetScratchMemTop ( sp );
  return result;
} /*mbs_TabBezCurveDer2d*/

