
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

boolean _mbs_TabBezC1Coons0Der3d ( int spdimen,
      int nknu, const double *knu, const double *hfuncu,
      const double *dhfuncu, const double *ddhfuncu, const double *dddhfuncu,
      int nknv, const double *knv, const double *hfuncv,
      const double *dhfuncv, const double *ddhfuncv, const double *dddhfuncv,
      int degc00, const double *c00,
      int degc01, const double *c01,
      int degd00, const double *d00,
      int degd01, const double *d01,
      double *p, double *pu, double *pv, double *puu, double *puv, double *pvv,
      double *puuu, double *puuv, double *puvv, double *pvvv,
      double *workspace )
{
  int    ku, kv;
  double *c, *dc, *ddc, *dddc, *d, *dd, *ddd, *dddd, *pcorners, *aux;

  ku = 6*spdimen*nknu;
  kv = 6*spdimen*nknv;
  c = workspace;  dc = &c[ku];  ddc = &dc[ku];  dddc = &ddc[ku];
  d = &dddc[ku];  dd = &d[kv];  ddd = &dd[kv];  dddd = &ddd[kv];
  pcorners = &dddd[kv];
  aux = &pcorners[16*spdimen];

  if ( !_mbs_multiBCHornerDerd ( degc00, 1, spdimen, 0, c00, 0.0,
      &pcorners[0], &pcorners[spdimen*2], aux ) )
    return false;
  if ( !_mbs_multiBCHornerDerd ( degc01, 1, spdimen, 0, c01, 0.0,
      &pcorners[spdimen*1], &pcorners[spdimen*3], aux ) )
    return false;

  if ( !_mbs_TabBezCurveDer3d ( spdimen, degc00, c00, nknu, knu, 2*spdimen,
            &c[0], &dc[0], &ddc[0], &dddc[0], aux ) )
    return false;
  if ( !_mbs_TabBezCurveDer3d ( spdimen, degc01, c01, nknu, knu, 2*spdimen,
            &c[spdimen], &dc[spdimen], &ddc[spdimen], &dddc[spdimen], aux ) )
    return false;
  if ( !_mbs_TabBezCurveDer3d ( spdimen, degd00, d00, nknv, knv, 2*spdimen,
            &d[0], &dd[0], &ddd[0], &dddd[0], aux ) )
    return false;
  if ( !_mbs_TabBezCurveDer3d ( spdimen, degd01, d01, nknv, knv, 2*spdimen,
            &d[spdimen], &dd[spdimen], &ddd[spdimen], &dddd[spdimen], aux ) )
    return false;

  if ( p )
    _mbs_TabBezC1Coons0d ( spdimen, nknu, nknv, c, d, pcorners,
                                 hfuncu, hfuncv, p, aux );
  if ( pu )
    _mbs_TabBezC1Coons0d ( spdimen, nknu, nknv, dc, d, pcorners,
                                 dhfuncu, hfuncv, pu, aux );
  if ( pv )
    _mbs_TabBezC1Coons0d ( spdimen, nknu, nknv, c, dd, pcorners,
                                 hfuncu, dhfuncv, pv, aux );
  if ( puu )
    _mbs_TabBezC1Coons0d ( spdimen, nknu, nknv, ddc, d, pcorners,
                                 ddhfuncu, hfuncv, puu, aux );
  if ( puv )
    _mbs_TabBezC1Coons0d ( spdimen, nknu, nknv, dc, dd, pcorners,
                                 dhfuncu, dhfuncv, puv, aux );
  if ( pvv )
    _mbs_TabBezC1Coons0d ( spdimen, nknu, nknv, c, ddd, pcorners,
                                 hfuncu, ddhfuncv, pvv, aux );
  if ( puuu )
    _mbs_TabBezC1Coons0d ( spdimen, nknu, nknv, dddc, d, pcorners,
                                 dddhfuncu, hfuncv, puuu, aux );
  if ( puuv )
    _mbs_TabBezC1Coons0d ( spdimen, nknu, nknv, ddc, dd, pcorners,
                                 ddhfuncu, dhfuncv, puuv, aux );
  if ( puvv )
    _mbs_TabBezC1Coons0d ( spdimen, nknu, nknv, dc, ddd, pcorners,
                                 dhfuncu, ddhfuncv, puvv, aux );
  if ( pvvv )
    _mbs_TabBezC1Coons0d ( spdimen, nknu, nknv, c, dddd, pcorners,
                                 hfuncu, dddhfuncv, pvvv, aux );
  return true;
} /*_mbs_TabBezC1Coons0Der3d*/

