
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

boolean mbs_TabBezC1Coons0Der3d ( int spdimen,
      int nknu, const double *knu, const double *hfuncu,
      const double *dhfuncu, const double *ddhfuncu, const double *dddhfuncu,
      int nknv, const double *knv, const double *hfuncv,
      const double *dhfuncv, const double *ddhfuncv, const double *dddhfuncv,
      int degc00, const double *c00,
      int degc01, const double *c01,
      int degd00, const double *d00,
      int degd01, const double *d01,
      double *p, double *pu, double *pv, double *puu, double *puv, double *pvv,
      double *puuu, double *puuv, double *puvv, double *pvvv )
{
  void   *sp;
  double *workspace;
  boolean result;

  sp = workspace = pkv_GetScratchMemd ( (24*(nknu+nknv)+20)*spdimen );
  if ( !workspace )
    return false;
  result = _mbs_TabBezC1Coons0Der3d ( spdimen,
                nknu, knu, hfuncu, dhfuncu, ddhfuncu, dddhfuncu,
                nknv, knv, hfuncv, dhfuncv, ddhfuncv, dddhfuncv,
                degc00, c00, degc01, c01, degd00, d00, degd01, d01,
                p, pu, pv, puu, puv, pvv, puuu, puuv, puvv, pvvv, workspace );
  pkv_SetScratchMemTop ( sp );
  return result;
} /*mbs_TabBezC1Coons0Der3d*/

