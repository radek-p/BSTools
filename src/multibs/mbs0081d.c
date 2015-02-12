
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

boolean _mbs_TabBezC2CoonsDer3d ( int spdimen,
      int nknu, const double *knu, const double *hfuncu,
      const double *dhfuncu, const double *ddhfuncu, const double *dddhfuncu,
      int nknv, const double *knv, const double *hfuncv,
      const double *dhfuncv, const double *ddhfuncv, const double *dddhfuncv,
      int degc00, const double *c00,
      int degc01, const double *c01,
      int degc02, const double *c02,
      int degc10, const double *c10,
      int degc11, const double *c11,
      int degc12, const double *c12,
      int degd00, const double *d00,
      int degd01, const double *d01,
      int degd02, const double *d02,
      int degd10, const double *d10,
      int degd11, const double *d11,
      int degd12, const double *d12,
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
  aux = &pcorners[36*spdimen];

  if ( !_mbs_BezC2CoonsFindCornersd ( spdimen,
              degc00, c00, degc01, c01, degc02, c02,
              degc10, c10, degc11, c11, degc12, c12,
              pcorners, aux ) )
    return false;

  if ( !_mbs_TabBezCurveDer3d ( spdimen, degc00, c00, nknu, knu, 6*spdimen,
            &c[0], &dc[0], &ddc[0], &dddc[0], aux ) )
    return false;
  if ( !_mbs_TabBezCurveDer3d ( spdimen, degc10, c10, nknu, knu, 6*spdimen,
            &c[spdimen], &dc[spdimen], &ddc[spdimen], &dddc[spdimen], aux ) )
    return false;
  if ( !_mbs_TabBezCurveDer3d ( spdimen, degc01, c01, nknu, knu, 6*spdimen,
            &c[2*spdimen], &dc[2*spdimen], &ddc[2*spdimen], &dddc[2*spdimen], aux ) )
    return false;
  if ( !_mbs_TabBezCurveDer3d ( spdimen, degc11, c11, nknu, knu, 6*spdimen,
            &c[3*spdimen], &dc[3*spdimen], &ddc[3*spdimen], &dddc[3*spdimen], aux ) )
    return false;
  if ( !_mbs_TabBezCurveDer3d ( spdimen, degc02, c02, nknu, knu, 6*spdimen,
            &c[4*spdimen], &dc[4*spdimen], &ddc[4*spdimen], &dddc[4*spdimen], aux ) )
    return false;
  if ( !_mbs_TabBezCurveDer3d ( spdimen, degc12, c12, nknu, knu, 6*spdimen,
            &c[5*spdimen], &dc[5*spdimen], &ddc[5*spdimen], &dddc[5*spdimen], aux ) )
    return false;

  if ( !_mbs_TabBezCurveDer3d ( spdimen, degd00, d00, nknv, knv, 6*spdimen,
            &d[0], &dd[0], &ddd[0], &dddd[0], aux ) )
    return false;
  if ( !_mbs_TabBezCurveDer3d ( spdimen, degd10, d10, nknv, knv, 6*spdimen,
            &d[spdimen], &dd[spdimen], &ddd[spdimen], &dddd[spdimen], aux ) )
    return false;
  if ( !_mbs_TabBezCurveDer3d ( spdimen, degd01, d01, nknv, knv, 6*spdimen,
            &d[2*spdimen], &dd[2*spdimen], &ddd[2*spdimen], &dddd[2*spdimen], aux ) )
    return false;
  if ( !_mbs_TabBezCurveDer3d ( spdimen, degd11, d11, nknv, knv, 6*spdimen,
            &d[3*spdimen], &dd[3*spdimen], &ddd[3*spdimen], &dddd[3*spdimen], aux ) )
    return false;
  if ( !_mbs_TabBezCurveDer3d ( spdimen, degd02, d02, nknv, knv, 6*spdimen,
            &d[4*spdimen], &dd[4*spdimen], &ddd[4*spdimen], &dddd[4*spdimen], aux ) )
    return false;
  if ( !_mbs_TabBezCurveDer3d ( spdimen, degd12, d12, nknv, knv, 6*spdimen,
            &d[5*spdimen], &dd[5*spdimen], &ddd[5*spdimen], &dddd[5*spdimen], aux ) )
    return false;

  if ( p )
    _mbs_TabBezC2Coonsd ( spdimen, nknu, nknv, c, d, pcorners,
                          hfuncu, hfuncv, p, aux );
  if ( pu )
    _mbs_TabBezC2Coonsd ( spdimen, nknu, nknv, dc, d, pcorners,
                          dhfuncu, hfuncv, pu, aux );
  if ( pv )
    _mbs_TabBezC2Coonsd ( spdimen, nknu, nknv, c, dd, pcorners,
                          hfuncu, dhfuncv, pv, aux );
  if ( puu )
    _mbs_TabBezC2Coonsd ( spdimen, nknu, nknv, ddc, d, pcorners,
                          ddhfuncu, hfuncv, puu, aux );
  if ( puv )
    _mbs_TabBezC2Coonsd ( spdimen, nknu, nknv, dc, dd, pcorners,
                          dhfuncu, dhfuncv, puv, aux );
  if ( pvv )
    _mbs_TabBezC2Coonsd ( spdimen, nknu, nknv, c, ddd, pcorners,
                          hfuncu, ddhfuncv, pvv, aux );
  if ( puuu )
    _mbs_TabBezC2Coonsd ( spdimen, nknu, nknv, dddc, d, pcorners,
                          dddhfuncu, hfuncv, puuu, aux );
  if ( puuv )
    _mbs_TabBezC2Coonsd ( spdimen, nknu, nknv, ddc, dd, pcorners,
                          ddhfuncu, dhfuncv, puuv, aux );
  if ( puvv )
    _mbs_TabBezC2Coonsd ( spdimen, nknu, nknv, dc, ddd, pcorners,
                          dhfuncu, ddhfuncv, puvv, aux );
  if ( pvvv )
    _mbs_TabBezC2Coonsd ( spdimen, nknu, nknv, c, dddd, pcorners,
                          hfuncu, dddhfuncv, pvvv, aux );
  return true;
} /*_mbs_TabBezC2CoonsDer3d*/

