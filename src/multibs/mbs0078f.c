
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

boolean _mbs_TabBezCurveDer3f ( int spdimen, int degree, const float *cp,
                                int nkn, const float *kn,
                                int ppitch,
                                float *p, float *dp, float *ddp, float *dddp,
                                float *workspace )
{
  int i, j;

  for ( i = j = 0;  i < nkn;  i++, j += ppitch )
    if ( !_mbs_multiBCHornerDer3f ( degree, 1, spdimen, 0, cp, kn[i],
                    &p[j], &dp[j], &ddp[j], &dddp[j], workspace ) )
      return false;
  return true;
} /*_mbs_TabBezCurveDer3f*/

