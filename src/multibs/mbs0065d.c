
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

boolean mbs_TabBezC1CoonsDer2d ( int spdimen,
      int nknu, const double *knu, const double *hfuncu,
      const double *dhfuncu, const double *ddhfuncu,
      int nknv, const double *knv, const double *hfuncv,
      const double *dhfuncv, const double *ddhfuncv,
      int degc00, const double *c00,
      int degc01, const double *c01,
      int degc10, const double *c10,
      int degc11, const double *c11,
      int degd00, const double *d00,
      int degd01, const double *d01,
      int degd10, const double *d10,
      int degd11, const double *d11,
      double *p, double *pu, double *pv, double *puu, double *puv, double *pvv )
{
  void   *sp;
  double *workspace;
  boolean result;

  sp = workspace = pkv_GetScratchMemd ( (18*(nknu+nknv) + 16 +
                                         4*max(nknu,nknv))*spdimen );
  if ( !workspace )
    return false;
  result = _mbs_TabBezC1CoonsDer2d ( spdimen,
                nknu, knu, hfuncu, dhfuncu, ddhfuncu,
                nknv, knv, hfuncv, dhfuncv, ddhfuncv,
                degc00, c00, degc01, c01, degc10, c10, degc11, c11,
                degd00, d00, degd01, d01, degd10, d10, degd11, d11,
                p, pu, pv, puu, puv, pvv, workspace );
  pkv_SetScratchMemTop ( sp );
  return result;
} /*mbs_TabBezC1CoonsDer2d*/

