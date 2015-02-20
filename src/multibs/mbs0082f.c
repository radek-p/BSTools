
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

boolean mbs_TabBezC2CoonsDer3f ( int spdimen,
      int nknu, const float *knu, const float *hfuncu,
      const float *dhfuncu, const float *ddhfuncu, const float *dddhfuncu,
      int nknv, const float *knv, const float *hfuncv,
      const float *dhfuncv, const float *ddhfuncv, const float *dddhfuncv,
      int degc00, const float *c00,
      int degc01, const float *c01,
      int degc02, const float *c02,
      int degc10, const float *c10,
      int degc11, const float *c11,
      int degc12, const float *c12,
      int degd00, const float *d00,
      int degd01, const float *d01,
      int degd02, const float *d02,
      int degd10, const float *d10,
      int degd11, const float *d11,
      int degd12, const float *d12,
      float *p, float *pu, float *pv, float *puu, float *puv, float *pvv,
      float *puuu, float *puuv, float *puvv, float *pvvv )
{
  void *sp;
  float *workspace;
  boolean result;

  sp = workspace = pkv_GetScratchMemf ( (24*(nknu+nknv) +
                                        6*max(nknu,nknv) + 36)*spdimen );
  if ( !workspace )
    return false;
  result = _mbs_TabBezC2CoonsDer3f ( spdimen,
                nknu, knu, hfuncu, dhfuncu, ddhfuncu, dddhfuncu,
                nknv, knv, hfuncv, dhfuncv, ddhfuncv, dddhfuncv,
                degc00, c00, degc01, c01, degc02, c02,
                degc10, c10, degc11, c11, degc12, c12,
                degd00, d00, degd01, d01, degd02, d02,
                degd10, d10, degd11, d11, degd12, d12,
                p, pu, pv, puu, puv, pvv, puuu, puuv, puvv, pvvv, workspace );
  pkv_SetScratchMemTop ( sp );
  return result;
} /*mbs_TabBezC2CoonsDer3f*/
