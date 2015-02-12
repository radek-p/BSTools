
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

boolean _mbs_TabBezC2Coons0Der3f ( int spdimen,
      int nknu, const float *knu, const float *hfuncu,
      const float *dhfuncu, const float *ddhfuncu, const float *dddhfuncu,
      int nknv, const float *knv, const float *hfuncv,
      const float *dhfuncv, const float *ddhfuncv, const float *dddhfuncv,
      int degc00, const float *c00,
      int degc01, const float *c01,
      int degc02, const float *c02,
      int degd00, const float *d00,
      int degd01, const float *d01,
      int degd02, const float *d02,
      float *p, float *pu, float *pv, float *puu, float *puv, float *pvv,
      float *puuu, float *puuv, float *puvv, float *pvvv,
      float *workspace )
{
  int    ku, kv;
  float *c, *dc, *ddc, *dddc, *d, *dd, *ddd, *dddd, *pcorners, *aux;

  ku = 3*spdimen*nknu;
  kv = 3*spdimen*nknv;
  c = workspace;  dc = &c[ku];  ddc = &dc[ku];  dddc = &ddc[ku];
  d = &dddc[ku];  dd = &d[kv];  ddd = &dd[kv];  dddd = &ddd[kv];
  pcorners = &dddd[kv];
  aux = &pcorners[36*spdimen];

  if ( !_mbs_multiBCHornerDer2f ( degc00, 1, spdimen, 0, c00, 0.0,
          &pcorners[0], &pcorners[spdimen*3], &pcorners[spdimen*6], aux ) )
    return false;
  if ( !_mbs_multiBCHornerDer2f ( degc01, 1, spdimen, 0, c01, 0.0,
          &pcorners[spdimen*1], &pcorners[spdimen*4], &pcorners[spdimen*7], aux ) )
    return false;
  if ( !_mbs_multiBCHornerDer2f ( degc02, 1, spdimen, 0, c02, 0.0,
          &pcorners[spdimen*2], &pcorners[spdimen*5], &pcorners[spdimen*8], aux ) )
    return false;

  if ( !_mbs_TabBezCurveDer3f ( spdimen, degc00, c00, nknu, knu, 3*spdimen,
            &c[0], &dc[0], &ddc[0], &dddc[0], aux ) )
    return false;
  if ( !_mbs_TabBezCurveDer3f ( spdimen, degc01, c01, nknu, knu, 3*spdimen,
            &c[1*spdimen], &dc[1*spdimen], &ddc[1*spdimen], &dddc[1*spdimen], aux ) )
    return false;
  if ( !_mbs_TabBezCurveDer3f ( spdimen, degc02, c02, nknu, knu, 3*spdimen,
            &c[2*spdimen], &dc[2*spdimen], &ddc[2*spdimen], &dddc[2*spdimen], aux ) )
    return false;
  if ( !_mbs_TabBezCurveDer3f ( spdimen, degd00, d00, nknv, knv, 3*spdimen,
            &d[0], &dd[0], &ddd[0], &dddd[0], aux ) )
    return false;
  if ( !_mbs_TabBezCurveDer3f ( spdimen, degd01, d01, nknv, knv, 3*spdimen,
            &d[1*spdimen], &dd[1*spdimen], &ddd[1*spdimen], &dddd[1*spdimen], aux ) )
    return false;
  if ( !_mbs_TabBezCurveDer3f ( spdimen, degd02, d02, nknv, knv, 3*spdimen,
            &d[2*spdimen], &dd[2*spdimen], &ddd[2*spdimen], &dddd[2*spdimen], aux ) )
    return false;

  if ( p )
    _mbs_TabBezC2Coons0f ( spdimen, nknu, nknv, c, d, pcorners,
                           hfuncu, hfuncv, p, aux );
  if ( pu )
    _mbs_TabBezC2Coons0f ( spdimen, nknu, nknv, dc, d, pcorners,
                           dhfuncu, hfuncv, pu, aux );
  if ( pv )
    _mbs_TabBezC2Coons0f ( spdimen, nknu, nknv, c, dd, pcorners,
                           hfuncu, dhfuncv, pv, aux );
  if ( puu )
    _mbs_TabBezC2Coons0f ( spdimen, nknu, nknv, ddc, d, pcorners,
                           ddhfuncu, hfuncv, puu, aux );
  if ( puv )
    _mbs_TabBezC2Coons0f ( spdimen, nknu, nknv, dc, dd, pcorners,
                           dhfuncu, dhfuncv, puv, aux );
  if ( pvv )
    _mbs_TabBezC2Coons0f ( spdimen, nknu, nknv, c, ddd, pcorners,
                           hfuncu, ddhfuncv, pvv, aux );
  if ( puuu )
    _mbs_TabBezC2Coons0f ( spdimen, nknu, nknv, dddc, d, pcorners,
                           dddhfuncu, hfuncv, puuu, aux );
  if ( puuv )
    _mbs_TabBezC2Coons0f ( spdimen, nknu, nknv, ddc, dd, pcorners,
                           ddhfuncu, dhfuncv, puuv, aux );
  if ( puvv )
    _mbs_TabBezC2Coons0f ( spdimen, nknu, nknv, dc, ddd, pcorners,
                           dhfuncu, ddhfuncv, puvv, aux );
  if ( pvvv )
    _mbs_TabBezC2Coons0f ( spdimen, nknu, nknv, c, dddd, pcorners,
                           hfuncu, dddhfuncv, pvvv, aux );
  return true;
} /*_mbs_TabBezC2Coons0Der3f*/

