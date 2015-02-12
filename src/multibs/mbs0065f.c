
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

boolean mbs_TabBezC1CoonsDer2f ( int spdimen,
      int nknu, const float *knu, const float *hfuncu,
      const float *dhfuncu, const float *ddhfuncu,
      int nknv, const float *knv, const float *hfuncv,
      const float *dhfuncv, const float *ddhfuncv,
      int degc00, const float *c00,
      int degc01, const float *c01,
      int degc10, const float *c10,
      int degc11, const float *c11,
      int degd00, const float *d00,
      int degd01, const float *d01,
      int degd10, const float *d10,
      int degd11, const float *d11,
      float *p, float *pu, float *pv, float *puu, float *puv, float *pvv )
{
  void   *sp;
  float  *workspace;
  boolean result;

  sp = workspace = pkv_GetScratchMemf ( (18*(nknu+nknv) + 16 +
                                         4*max(nknu,nknv))*spdimen );
  if ( !workspace )
    return false;
  result = _mbs_TabBezC1CoonsDer2f ( spdimen,
                nknu, knu, hfuncu, dhfuncu, ddhfuncu,
                nknv, knv, hfuncv, dhfuncv, ddhfuncv,
                degc00, c00, degc01, c01, degc10, c10, degc11, c11,
                degd00, d00, degd01, d01, degd10, d10, degd11, d11,
                p, pu, pv, puu, puv, pvv, workspace );
  pkv_SetScratchMemTop ( sp );
  return result;
} /*mbs_TabBezC1CoonsDer2f*/

