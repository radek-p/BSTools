
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

boolean mbs_TabBezC1Coons0Der2f ( int spdimen,
      int nknu, const float *knu, const float *hfuncu,
      const float *dhfuncu, const float *ddhfuncu,
      int nknv, const float *knv, const float *hfuncv,
      const float *dhfuncv, const float *ddhfuncv,
      int degc00, const float *c00,
      int degc01, const float *c01,
      int degd00, const float *d00,
      int degd01, const float *d01,
      float *p, float *pu, float *pv, float *puu, float *puv, float *pvv )
{
  void    *sp;
  float   *workspace;
  boolean result;

  sp = workspace = pkv_GetScratchMemf ( (10*(nknu+nknv)+4)*spdimen );
  if ( !workspace )
    return false;
  result = _mbs_TabBezC1Coons0Der2f ( spdimen,
                nknu, knu, hfuncu, dhfuncu, ddhfuncu,
                nknv, knv, hfuncv, dhfuncv, ddhfuncv,
                degc00, c00, degc01, c01, degd00, d00, degd01, d01,
                p, pu, pv, puu, puv, pvv, workspace );
  pkv_SetScratchMemTop ( sp );
  return result;
} /*mbs_TabBezC1Coons0Der2f*/

