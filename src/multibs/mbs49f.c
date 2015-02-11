
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

boolean _mbs_BezC1CoonsFindCornersf ( int spdimen,
                                     int degc00, const float *c00,
                                     int degc01, const float *c01,
                                     int degc10, const float *c10,
                                     int degc11, const float *c11,
                                     float *pcorners,
                                     float *workspace )
{
  if ( !_mbs_multiBCHornerDerf ( degc00, 1, spdimen, 0, c00, 0.0,
                  &pcorners[0], &pcorners[spdimen*8], workspace ) )
    return false;
  if ( !_mbs_multiBCHornerDerf ( degc00, 1, spdimen, 0, c00, 1.0,
                  &pcorners[spdimen*4], &pcorners[spdimen*12], workspace ) )
    return false;
  if ( !_mbs_multiBCHornerDerf ( degc10, 1, spdimen, 0, c10, 0.0,
                  &pcorners[spdimen*1], &pcorners[spdimen*9], workspace ) )
    return false;
  if ( !_mbs_multiBCHornerDerf ( degc10, 1, spdimen, 0, c10, 1.0,
                  &pcorners[spdimen*5], &pcorners[spdimen*13], workspace ) )
    return false;
  if ( !_mbs_multiBCHornerDerf ( degc01, 1, spdimen, 0, c01, 0.0,
                  &pcorners[spdimen*2], &pcorners[spdimen*10], workspace ) )
    return false;
  if ( !_mbs_multiBCHornerDerf ( degc01, 1, spdimen, 0, c01, 1.0,
                  &pcorners[spdimen*6], &pcorners[spdimen*14], workspace ) )
    return false;
  if ( !_mbs_multiBCHornerDerf ( degc11, 1, spdimen, 0, c11, 0.0,
                  &pcorners[spdimen*3], &pcorners[spdimen*11], workspace ) )
    return false;
  if ( !_mbs_multiBCHornerDerf ( degc11, 1, spdimen, 0, c11, 1.0,
                  &pcorners[spdimen*7], &pcorners[spdimen*15], workspace ) )
    return false;
  return true;
} /*_mbs_BezC1CoonsFindCornersf*/

boolean mbs_BezC1CoonsFindCornersf ( int spdimen,
                                     int degc00, const float *c00,
                                     int degc01, const float *c01,
                                     int degc10, const float *c10,
                                     int degc11, const float *c11,
                                     float *pcorners )
{
  void    *sp;
  float   *workspace;
  boolean result;
  
  sp = workspace = pkv_GetScratchMemf ( 2*spdimen );
  if ( !workspace )
    return false;
  result = _mbs_BezC1CoonsFindCornersf ( spdimen,
                   degc00, c00, degc01, c01, degc10, c10, degc11, c11,
                   pcorners, workspace );
  pkv_SetScratchMemTop ( sp );
  return result;
} /*mbs_BezC1CoonsFindCornersf*/

static int FindMaxInt ( int a, int b, int c, int d, int e )
{
  if ( b > a ) a = b;  if ( c > a ) a = c;  if ( d > a ) a = d;
  return e > a ? e : a;
} /*FindMaxInt*/

boolean _mbs_BezC1CoonsToBezf ( int spdimen,
                                int degc00, const float *c00,
                                int degc01, const float *c01,
                                int degc10, const float *c10,
                                int degc11, const float *c11,
                                int degd00, const float *d00,
                                int degd01, const float *d01,
                                int degd10, const float *d10,
                                int degd11, const float *d11,
                                int *n, int *m, float *p,
                                float *workspace )
{
  int   degu, degv, d, du, dv;
  int   pitch;
  float *pc, *bc;
  float *p1, *p2, *p3; 
  int   i;

  *n = degu = FindMaxInt ( degc00, degc01, degc10, degc11, 3 );
  *m = degv = FindMaxInt ( degd00, degd01, degd10, degd11, 3 );
  pc = &workspace[2*spdimen];
  p1 = &pc[16*spdimen];
  p2 = &p1[(degu+1)*(degv+1)*spdimen];
  p3 = &p2[(degu+1)*(degv+1)*spdimen];
  bc = &p3[(degu+1)*(degv+1)*spdimen];

        /* construct the patch p1 */
  if ( !mbs_multiBCDegElevf ( 1, spdimen, 0, degc00, c00, degu-degc00,
                              0, &d, &bc[0] ) )
    return false;
  if ( !mbs_multiBCDegElevf ( 1, spdimen, 0, degc01, c01, degu-degc01,
                              0, &d, &bc[(degu+1)*spdimen] ) )
    return false;
  if ( !mbs_multiBCDegElevf ( 1, spdimen, 0, degc10, c10, degu-degc10,
                              0, &d, &bc[2*(degu+1)*spdimen] ) )
    return false;
  if ( !mbs_multiBCDegElevf ( 1, spdimen, 0, degc11, c11, degu-degc11,
                              0, &d, &bc[3*(degu+1)*spdimen] ) )
    return false;
  pitch = (degu+1)*spdimen;
  if ( !mbs_multiInterp2knHermiteBezf ( 1, pitch, 3, 2, pitch, bc,
                       2, pitch, &bc[2*(degu+1)*spdimen], pitch, p1 ) )
    return false;
  pkv_TransposeMatrixc ( 4, degu+1, spdimen*sizeof(float),
                         (degu+1)*spdimen*sizeof(float), (char*)p1,
                         4*spdimen*sizeof(float), (char*)p );
  if ( !mbs_BCDegElevPf ( spdimen, degu, 3, p, 0, degv-3, &du, &dv, p ) )
    return false;

        /* construct the patch p2 */
  if ( !mbs_multiBCDegElevf ( 1, spdimen, 0, degd00, d00, degv-degd00,
                              0, &d, &bc[0] ) )
    return false;
  if ( !mbs_multiBCDegElevf ( 1, spdimen, 0, degd01, d01, degv-degd01,
                              0, &d, &bc[(degv+1)*spdimen] ) )
    return false;
  if ( !mbs_multiBCDegElevf ( 1, spdimen, 0, degd10, d10, degv-degd10,
                              0, &d, &bc[2*(degv+1)*spdimen] ) )
    return false;
  if ( !mbs_multiBCDegElevf ( 1, spdimen, 0, degd11, d11, degv-degd11,
                              0, &d, &bc[3*(degv+1)*spdimen] ) )
    return false;
  pitch = (degv+1)*spdimen;
  if ( !mbs_multiInterp2knHermiteBezf ( 1, pitch, 3, 2, pitch, bc,
                       2, pitch, &bc[2*(degv+1)*spdimen], pitch, p2 ) )
    return false;
  if ( !mbs_BCDegElevPf ( spdimen, 3, degv, p2, degu-3, 0, &du, &dv, p2 ) )
    return false;
  pkn_AddMatrixf ( 1, spdimen*(degu+1)*(degv+1), 0, p, 0, p2, 0, p );

        /* construct the patch p3 */
  if ( !_mbs_BezC1CoonsFindCornersf ( spdimen, degc00, c00, degc01, c01,
                       degc10, c10, degc11, c11, pc, workspace ) )
    return false;

  pkv_Selectf ( 2, 4*spdimen, 4*2*spdimen, 4*spdimen, pc, p3 );
  pkv_Selectf ( 2, 4*spdimen, 4*2*spdimen, 4*spdimen, &pc[4*spdimen], &p3[4*2*spdimen] );
  for ( i = 0; i < 4; i++ ) {
    pkv_Selectf ( 2, spdimen, 2*spdimen, spdimen, &p3[4*i*spdimen], &bc[4*i*spdimen] );
    pkv_Selectf ( 2, spdimen, 2*spdimen, spdimen, &p3[(4*i+1)*spdimen], &bc[(4*i+2)*spdimen] );
  }
  if ( !mbs_multiInterp2knHermiteBezf ( 1, 4*spdimen, 3, 2, 0, bc,
                                        2, 0, &bc[4*2*spdimen], 0, p3 ) )
    return false;
  if ( !mbs_multiInterp2knHermiteBezf ( 4, spdimen, 3, 2, 4*spdimen, p3,
                                  2, 4*spdimen, &p3[2*spdimen], 4*spdimen, p2 ) )
    return false;
  if ( !mbs_BCDegElevPf ( spdimen, 3, 3, p2, degu-3, degv-3, &du, &dv, p3 ) )
    return false;
  pkn_SubtractMatrixf ( 1, spdimen*(degu+1)*(degv+1), 0, p, 0, p3, 0, p );
  return true;
} /*_mbs_BezC1CoonsToBezf*/

boolean mbs_BezC1CoonsToBezf ( int spdimen,
                               int degc00, const float *c00,
                               int degc01, const float *c01,
                               int degc10, const float *c10,
                               int degc11, const float *c11,
                               int degd00, const float *d00,
                               int degd01, const float *d01,
                               int degd10, const float *d10,
                               int degd11, const float *d11,
                               int *n, int *m, float *p )
{
  void    *sp;
  float   *workspace;
  int     degu, degv, size;
  boolean result;

  degu = FindMaxInt ( degc00, degc01, degc10, degc11, 3 );
  degv = FindMaxInt ( degd00, degd01, degd10, degd11, 3 );
  size = (18 + 3*(degu+1)*(degv+1) + 4*(max(degu,degv)+1))*spdimen;
  sp = workspace = pkv_GetScratchMemf ( size );
  if ( !workspace )
    return false;
  result = _mbs_BezC1CoonsToBezf ( spdimen,
                   degc00, c00, degc01, c01, degc10, c10, degc11, c11,
                   degd00, d00, degd01, d01, degd10, d10, degd11, d11,
                   n, m, p, workspace );
  pkv_SetScratchMemTop ( sp );
  return result;
} /*mbs_BezC1CoonsToBezf*/

/* ///////////////////////////////////////////////////////////////////////// */
boolean mbs_TabCubicHFuncDer2f ( float a, float b, int nkn, const float *kn,
                                 float *hfunc, float *dhfunc, float *ddhfunc )
{
  float HFunc[16] = {1.0, 1.0, 0.0, 0.0,                  /* h00 */
                     0.0, 0.0, 1.0, 1.0,                  /* h10 */
                     0.0, (float)(1.0/3.0), 0.0, 0.0,     /* h01 */
                     0.0, 0.0, (float)(-1.0/3.0), 0.0};   /* h11 */
  int   i, j;
  float h, h_1, h_2, wsp[3];

  if ( a == 0.0 && b == 1.0 ) {
     for ( i = j = 0;  i < nkn;  i++, j += 4 )
       if ( !_mbs_multiBCHornerDer2f ( 3, 4, 1, 4, HFunc, kn[i],
                       &hfunc[j], &dhfunc[j], &ddhfunc[j], wsp ) )
         return false;
  }
  else {
    h = b - a;
    h_1 = 1.0/h;
    h_2 = h_1*h_1;
    for ( i = j = 0;  i < nkn;  i++, j += 4 ) {
      if ( !_mbs_multiBCHornerDer2f ( 3, 4, 1, 4, HFunc, (kn[i]-a)/h,
                      &hfunc[j], &dhfunc[j], &ddhfunc[j], wsp ) )
        return false;
      hfunc[j+2]   *= h;      hfunc[j+3]   *= h;
      dhfunc[j]    *= h_1;    dhfunc[j+1]  *= h_1;  
      ddhfunc[j]   *= h_2;    ddhfunc[j+1] *= h_2;  
      ddhfunc[j+2] *= h_1;    ddhfunc[j+3] *= h_1;  
    }
  }
  return true;
} /*mbs_TabCubicHFuncDer2f*/

boolean mbs_TabCubicHFuncDer3f ( float a, float b, int nkn, const float *kn,
                                 float *hfunc, float *dhfunc, float *ddhfunc,
                                 float *dddhfunc )
{
  float HFunc[16] = {1.0, 1.0, 0.0, 0.0,                  /* h00 */
                     0.0, 0.0, 1.0, 1.0,                  /* h10 */
                     0.0, (float)(1.0/3.0), 0.0, 0.0,     /* h01 */
                     0.0, 0.0, (float)(-1.0/3.0), 0.0};   /* h11 */
  int   i, j;
  float h, h_1, h_2, h_3, wsp[4];

  if ( a == 0.0 && b == 1.0 ) {
     for ( i = j = 0;  i < nkn;  i++, j += 4 )
       if ( !_mbs_multiBCHornerDer3f ( 3, 4, 1, 4, HFunc, kn[i],
                            &hfunc[j], &dhfunc[j], &ddhfunc[j], &dddhfunc[j], wsp ) )
         return false;
  }
  else {
    h = b - a;
    h_1 = (float)1.0/h;
    h_2 = h_1*h_1;
    h_3 = h_1*h_2;
    for ( i = j = 0;  i < nkn;  i++, j += 4 ) {
      if ( !_mbs_multiBCHornerDer3f ( 3, 4, 1, 4, HFunc, (kn[i]-a)/h,
                           &hfunc[j], &dhfunc[j], &ddhfunc[j], &dddhfunc[j], wsp ) )
        return false;
      hfunc[j+2]    *= h;      hfunc[j+3]   *= h;
      dhfunc[j]     *= h_1;    dhfunc[j+1]  *= h_1;
      ddhfunc[j]    *= h_2;    ddhfunc[j+1] *= h_2;
      ddhfunc[j+2]  *= h_1;    ddhfunc[j+3] *= h_1;
      dddhfunc[j]   *= h_3;    dddhfunc[j+1] *= h_3;
      dddhfunc[j+2] *= h_2;    dddhfunc[j+3] *= h_2;
    }
  }
  return true;
} /*mbs_TabCubicHFuncDer3f*/

boolean _mbs_TabBezCurveDer2f ( int spdimen, int degree, const float *cp,
                                int nkn, const float *kn,
                                int ppitch,
                                float *p, float *dp, float *ddp,
                                float *workspace )
{
  int i, j;

  for ( i = j = 0;  i < nkn;  i++, j += ppitch )
    if ( !_mbs_multiBCHornerDer2f ( degree, 1, spdimen, 0, cp, kn[i],
                                    &p[j], &dp[j], &ddp[j], workspace ) )
      return false;
  return true;
} /*_mbs_TabBezCurveDer2f*/

boolean mbs_TabBezCurveDer2f ( int spdimen, int degree, const float *cp,
                               int nkn, const float *kn,
                               int ppitch,
                               float *p, float *dp, float *ddp )
{
  void    *sp;
  float   *workspace;
  boolean result;

  sp = workspace = pkv_GetScratchMem ( 3*spdimen );
  if ( !workspace )
    return false;
  result = _mbs_TabBezCurveDer2f ( spdimen, degree, cp, nkn, kn,
                                   ppitch, p, dp, ddp, workspace );
  pkv_SetScratchMemTop ( sp );
  return result;
} /*mbs_TabBezCurveDer2f*/

void _mbs_TabBezC1Coonsf ( int spdimen, int nknu, int nknv,
                   const float *c, const float *d, const float *p,
                   const float *hu, const float *hv, float *pp,
                   float *workspace )
{
  int i, j, k, l;

  memcpy ( workspace, d, 4*nknv*spdimen*sizeof(float) );
  for ( i = 0; i < nknv; i++ )
    for ( j = 0; j < 4; j++ )
      for ( k = 0; k < 4; k++ )
        for ( l = 0; l < spdimen; l++ )
          workspace[(4*i+j)*spdimen+l] -= hv[4*i+k]*p[(4*j+k)*spdimen+l];

  memset ( pp, 0, nknu*nknv*spdimen*sizeof(float) );
  for ( i = 0; i < nknu; i++ )
    for ( j = 0; j < nknv; j++ )
      for ( k = 0; k < 4; k++ )
        for ( l = 0; l < spdimen; l++ )
          pp[(nknv*i+j)*spdimen+l] += c[(4*i+k)*spdimen+l]*hv[4*j+k] +
                                      workspace[(4*j+k)*spdimen+l]*hu[4*i+k];
} /*_mbs_TabBezC1Coonsf*/

boolean _mbs_TabBezC1CoonsDer2f ( int spdimen,
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
      float *p, float *pu, float *pv, float *puu, float *puv, float *pvv,
      float *workspace )
{
  int   ku, kv;
  float *c, *dc, *ddc, *d, *dd, *ddd, *pcorners, *aux;

  ku = 4*spdimen*nknu;
  kv = 4*spdimen*nknv;
  c = workspace;  dc = &c[ku];  ddc = &dc[ku];
  d = &ddc[ku];   dd = &d[kv];  ddd = &dd[kv];
  pcorners = &ddd[kv];
  aux = &pcorners[16*spdimen];

  if ( !_mbs_BezC1CoonsFindCornersf ( spdimen, degc00, c00, degc01, c01,
                               degc10, c10, degc11, c11, pcorners, aux ) )
    return false;

  if ( !_mbs_TabBezCurveDer2f ( spdimen, degc00, c00, nknu, knu, 4*spdimen,
            &c[0], &dc[0], &ddc[0], aux ) )
    return false;
  if ( !_mbs_TabBezCurveDer2f ( spdimen, degc10, c10, nknu, knu, 4*spdimen,
            &c[spdimen], &dc[spdimen], &ddc[spdimen], aux ) )
    return false;
  if ( !_mbs_TabBezCurveDer2f ( spdimen, degc01, c01, nknu, knu, 4*spdimen,
            &c[2*spdimen], &dc[2*spdimen], &ddc[2*spdimen], aux ) )
    return false;
  if ( !_mbs_TabBezCurveDer2f ( spdimen, degc11, c11, nknu, knu, 4*spdimen,
            &c[3*spdimen], &dc[3*spdimen], &ddc[3*spdimen], aux ) )
    return false;

  if ( !_mbs_TabBezCurveDer2f ( spdimen, degd00, d00, nknv, knv, 4*spdimen,
            &d[0], &dd[0], &ddd[0], aux ) )
    return false;
  if ( !_mbs_TabBezCurveDer2f ( spdimen, degd10, d10, nknv, knv, 4*spdimen,
            &d[spdimen], &dd[spdimen], &ddd[spdimen], aux ) )
    return false;
  if ( !_mbs_TabBezCurveDer2f ( spdimen, degd01, d01, nknv, knv, 4*spdimen,
            &d[2*spdimen], &dd[2*spdimen], &ddd[2*spdimen], aux ) )
    return false;
  if ( !_mbs_TabBezCurveDer2f ( spdimen, degd11, d11, nknv, knv, 4*spdimen,
            &d[3*spdimen], &dd[3*spdimen], &ddd[3*spdimen], aux ) )
    return false;

  if ( p )
    _mbs_TabBezC1Coonsf ( spdimen, nknu, nknv, c, d, pcorners,
                                hfuncu, hfuncv, p, aux );
  if ( pu )
    _mbs_TabBezC1Coonsf ( spdimen, nknu, nknv, dc, d, pcorners,
                                dhfuncu, hfuncv, pu, aux );
  if ( pv )
    _mbs_TabBezC1Coonsf ( spdimen, nknu, nknv, c, dd, pcorners,
                                hfuncu, dhfuncv, pv, aux );
  if ( puu )
    _mbs_TabBezC1Coonsf ( spdimen, nknu, nknv, ddc, d, pcorners,
                                ddhfuncu, hfuncv, puu, aux );
  if ( puv )
    _mbs_TabBezC1Coonsf ( spdimen, nknu, nknv, dc, dd, pcorners,
                                dhfuncu, dhfuncv, puv, aux );
  if ( pvv )
    _mbs_TabBezC1Coonsf ( spdimen, nknu, nknv, c, ddd, pcorners,
                                hfuncu, ddhfuncv, pvv, aux );
  return true;
} /*_mbs_TabBezC1CoonsDer2f*/

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

boolean _mbs_TabBezC1CoonsDer3f ( int spdimen,
      int nknu, const float *knu, const float *hfuncu,
      const float *dhfuncu, const float *ddhfuncu, const float *dddhfuncu,
      int nknv, const float *knv, const float *hfuncv,
      const float *dhfuncv, const float *ddhfuncv, const float *dddhfuncv,
      int degc00, const float *c00,
      int degc01, const float *c01,
      int degc10, const float *c10,
      int degc11, const float *c11,
      int degd00, const float *d00,
      int degd01, const float *d01,
      int degd10, const float *d10,
      int degd11, const float *d11,
      float *p, float *pu, float *pv, float *puu, float *puv, float *pvv,
      float *puuu, float *puuv, float *puvv, float *pvvv,
      float *workspace )
{
  int   ku, kv;
  float *c, *dc, *ddc, *dddc, *d, *dd, *ddd, *dddd, *pcorners, *aux;

  ku = 6*spdimen*nknu;
  kv = 6*spdimen*nknv;
  c = workspace;  dc = &c[ku];  ddc = &dc[ku];  dddc = &ddc[ku];
  d = &dddc[ku];  dd = &d[kv];  ddd = &dd[kv];  dddd = &ddd[kv];
  pcorners = &dddd[kv];
  aux = &pcorners[16*spdimen];

  if ( !_mbs_BezC1CoonsFindCornersf ( spdimen, degc00, c00, degc01, c01,
                               degc10, c10, degc11, c11, pcorners, aux ) )
    return false;

  if ( !_mbs_TabBezCurveDer3f ( spdimen, degc00, c00, nknu, knu, 4*spdimen,
            &c[0], &dc[0], &ddc[0], &dddc[0], aux ) )
    return false;
  if ( !_mbs_TabBezCurveDer3f ( spdimen, degc10, c10, nknu, knu, 4*spdimen,
            &c[spdimen], &dc[spdimen], &ddc[spdimen], &dddc[spdimen], aux ) )
    return false;
  if ( !_mbs_TabBezCurveDer3f ( spdimen, degc01, c01, nknu, knu, 4*spdimen,
            &c[2*spdimen], &dc[2*spdimen], &ddc[2*spdimen], &dddc[2*spdimen], aux ) )
    return false;
  if ( !_mbs_TabBezCurveDer3f ( spdimen, degc11, c11, nknu, knu, 4*spdimen,
            &c[3*spdimen], &dc[3*spdimen], &ddc[3*spdimen], &dddc[3*spdimen], aux ) )
    return false;

  if ( !_mbs_TabBezCurveDer3f ( spdimen, degd00, d00, nknv, knv, 4*spdimen,
            &d[0], &dd[0], &ddd[0], &dddd[0], aux ) )
    return false;
  if ( !_mbs_TabBezCurveDer3f ( spdimen, degd10, d10, nknv, knv, 4*spdimen,
            &d[spdimen], &dd[spdimen], &ddd[spdimen], &dddd[spdimen], aux ) )
    return false;
  if ( !_mbs_TabBezCurveDer3f ( spdimen, degd01, d01, nknv, knv, 4*spdimen,
            &d[2*spdimen], &dd[2*spdimen], &ddd[2*spdimen], &dddd[2*spdimen], aux ) )
    return false;
  if ( !_mbs_TabBezCurveDer3f ( spdimen, degd11, d11, nknv, knv, 4*spdimen,
            &d[3*spdimen], &dd[3*spdimen], &ddd[3*spdimen], &dddd[3*spdimen], aux ) )
    return false;

  if ( p )
    _mbs_TabBezC1Coonsf ( spdimen, nknu, nknv, c, d, pcorners,
                                hfuncu, hfuncv, p, aux );
  if ( pu )
    _mbs_TabBezC1Coonsf ( spdimen, nknu, nknv, dc, d, pcorners,
                                dhfuncu, hfuncv, pu, aux );
  if ( pv )
    _mbs_TabBezC1Coonsf ( spdimen, nknu, nknv, c, dd, pcorners,
                                hfuncu, dhfuncv, pv, aux );
  if ( puu )
    _mbs_TabBezC1Coonsf ( spdimen, nknu, nknv, ddc, d, pcorners,
                                ddhfuncu, hfuncv, puu, aux );
  if ( puv )
    _mbs_TabBezC1Coonsf ( spdimen, nknu, nknv, dc, dd, pcorners,
                                dhfuncu, dhfuncv, puv, aux );
  if ( pvv )
    _mbs_TabBezC1Coonsf ( spdimen, nknu, nknv, c, ddd, pcorners,
                                hfuncu, ddhfuncv, pvv, aux );
  if ( puuu )
    _mbs_TabBezC1Coonsf ( spdimen, nknu, nknv, dddc, d, pcorners,
                                dddhfuncu, hfuncv, puuu, aux );
  if ( puuv )
    _mbs_TabBezC1Coonsf ( spdimen, nknu, nknv, ddc, dd, pcorners,
                                ddhfuncu, dhfuncv, puuv, aux );
  if ( puvv )
    _mbs_TabBezC1Coonsf ( spdimen, nknu, nknv, dc, ddd, pcorners,
                                dhfuncu, ddhfuncv, puvv, aux );
  if ( pvvv )
    _mbs_TabBezC1Coonsf ( spdimen, nknu, nknv, c, dddd, pcorners,
                                hfuncu, dddhfuncv, pvvv, aux );
  return true;
} /*_mbs_TabBezC1CoonsDer3f*/

boolean mbs_TabBezC1CoonsDer3f ( int spdimen,
      int nknu, const float *knu, const float *hfuncu,
      const float *dhfuncu, const float *ddhfuncu, const float *dddhfuncu,
      int nknv, const float *knv, const float *hfuncv,
      const float *dhfuncv, const float *ddhfuncv, const float *dddhfuncv,
      int degc00, const float *c00,
      int degc01, const float *c01,
      int degc10, const float *c10,
      int degc11, const float *c11,
      int degd00, const float *d00,
      int degd01, const float *d01,
      int degd10, const float *d10,
      int degd11, const float *d11,
      float *p, float *pu, float *pv, float *puu, float *puv, float *pvv,
      float *puuu, float *puuv, float *puvv, float *pvvv )
{
  void    *sp;
  int     ku, kv;
  float   *workspace;
  boolean result;

  ku = 6*spdimen*nknu;
  kv = 6*spdimen*nknv;
  sp = workspace = pkv_GetScratchMemf ( (24*(nknu+nknv)+20)*spdimen + max(ku,kv) );
  if ( !workspace )
    return false;
  result = _mbs_TabBezC1CoonsDer3f ( spdimen,
                nknu, knu, hfuncu, dhfuncu, ddhfuncu, dddhfuncu,
                nknv, knv, hfuncv, dhfuncv, ddhfuncv, dddhfuncv,
                degc00, c00, degc01, c01, degc10, c10, degc11, c11,
                degd00, d00, degd01, d01, degd10, d10, degd11, d11,
                p, pu, pv, puu, puv, pvv, puuu, puuv, puvv, pvvv, workspace );
  pkv_SetScratchMemTop ( sp );
  return result;
} /*mbs_TabBezC1CoonsDer3f*/

void _mbs_TabBezC1Coons0f ( int spdimen, int nknu, int nknv,
                   const float *c, const float *d, const float *p,
                   const float *hu, const float *hv, float *pp,
                   float *workspace )
{
  int i, j, k, l;

  memcpy ( workspace, d, 2*nknv*spdimen*sizeof(float) );
  for ( i = 0; i < nknv; i++ )
    for ( j = 0; j < 2; j++ )
      for ( k = 0; k < 2; k++ )
        for ( l = 0; l < spdimen; l++ )
          workspace[(2*i+j)*spdimen+l] -= hv[4*i+2*k]*p[(2*j+k)*spdimen+l];

  memset ( pp, 0, nknu*nknv*spdimen*sizeof(float) );
  for ( i = 0; i < nknu; i++ )
    for ( j = 0; j < nknv; j++ )
      for ( k = 0; k < 2; k++ )
        for ( l = 0; l < spdimen; l++ )
          pp[(nknv*i+j)*spdimen+l] += c[(2*i+k)*spdimen+l]*hv[4*j+2*k] +
                                      workspace[(2*j+k)*spdimen+l]*hu[4*i+2*k];
} /*_mbs_TabBezC1Coons0f*/

boolean _mbs_TabBezC1Coons0Der2f ( int spdimen,
      int nknu, const float *knu, const float *hfuncu,
      const float *dhfuncu, const float *ddhfuncu,
      int nknv, const float *knv, const float *hfuncv,
      const float *dhfuncv, const float *ddhfuncv,
      int degc00, const float *c00,
      int degc01, const float *c01,
      int degd00, const float *d00,
      int degd01, const float *d01,
      float *p, float *pu, float *pv, float *puu, float *puv, float *pvv,
      float *workspace )
{
  int   ku, kv;
  float *c, *dc, *ddc, *d, *dd, *ddd, *pcorners, *aux;

  ku = 2*spdimen*nknu;
  kv = 2*spdimen*nknv;
  c = workspace;  dc = &c[ku];  ddc = &dc[ku];
  d = &ddc[ku];   dd = &d[kv];  ddd = &dd[kv];
  pcorners = &ddd[kv];
  aux = &pcorners[16*spdimen];

  if ( !_mbs_multiBCHornerDerf ( degc00, 1, spdimen, 0, c00, 0.0,
      &pcorners[0], &pcorners[spdimen*2], aux ) )
    return false;
  if ( !_mbs_multiBCHornerDerf ( degc01, 1, spdimen, 0, c01, 0.0,
      &pcorners[spdimen*1], &pcorners[spdimen*3], aux) )
    return false;

  if ( !_mbs_TabBezCurveDer2f ( spdimen, degc00, c00, nknu, knu, 2*spdimen,
            &c[0], &dc[0], &ddc[0], aux ) )
    return false;
  if ( !_mbs_TabBezCurveDer2f ( spdimen, degc01, c01, nknu, knu, 2*spdimen,
            &c[spdimen], &dc[spdimen], &ddc[spdimen], aux ) )
    return false;
  if ( !_mbs_TabBezCurveDer2f ( spdimen, degd00, d00, nknv, knv, 2*spdimen,
            &d[0], &dd[0], &ddd[0], aux ) )
    return false;
  if ( !_mbs_TabBezCurveDer2f ( spdimen, degd01, d01, nknv, knv, 2*spdimen,
            &d[spdimen], &dd[spdimen], &ddd[spdimen], aux ) )
    return false;

  if ( p )
    _mbs_TabBezC1Coons0f ( spdimen, nknu, nknv, c, d, pcorners,
                                 hfuncu, hfuncv, p, aux );
  if ( pu )
    _mbs_TabBezC1Coons0f ( spdimen, nknu, nknv, dc, d, pcorners,
                                 dhfuncu, hfuncv, pu, aux );
  if ( pv )
    _mbs_TabBezC1Coons0f ( spdimen, nknu, nknv, c, dd, pcorners,
                                 hfuncu, dhfuncv, pv, aux );
  if ( puu )
    _mbs_TabBezC1Coons0f ( spdimen, nknu, nknv, ddc, d, pcorners,
                                 ddhfuncu, hfuncv, puu, aux );
  if ( puv )
    _mbs_TabBezC1Coons0f ( spdimen, nknu, nknv, dc, dd, pcorners,
                                 dhfuncu, dhfuncv, puv, aux );
  if ( pvv )
    _mbs_TabBezC1Coons0f ( spdimen, nknu, nknv, c, ddd, pcorners,
                                 hfuncu, ddhfuncv, pvv, aux );
  return true;
} /*_mbs_TabBezC1Coons0Der2f*/

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

boolean _mbs_TabBezC1Coons0Der3f ( int spdimen,
      int nknu, const float *knu, const float *hfuncu,
      const float *dhfuncu, const float *ddhfuncu, const float *dddhfuncu,
      int nknv, const float *knv, const float *hfuncv,
      const float *dhfuncv, const float *ddhfuncv, const float *dddhfuncv,
      int degc00, const float *c00,
      int degc01, const float *c01,
      int degd00, const float *d00,
      int degd01, const float *d01,
      float *p, float *pu, float *pv, float *puu, float *puv, float *pvv,
      float *puuu, float *puuv, float *puvv, float *pvvv,
      float *workspace )
{
  int   ku, kv;
  float *c, *dc, *ddc, *dddc, *d, *dd, *ddd, *dddd, *pcorners, *aux;

  ku = 6*spdimen*nknu;
  kv = 6*spdimen*nknv;
  c = workspace;  dc = &c[ku];  ddc = &dc[ku];  dddc = &ddc[ku];
  d = &dddc[ku];  dd = &d[kv];  ddd = &dd[kv];  dddd = &ddd[kv];
  pcorners = &dddd[kv];
  aux = &pcorners[16*spdimen];

  if ( !_mbs_multiBCHornerDerf ( degc00, 1, spdimen, 0, c00, 0.0,
      &pcorners[0], &pcorners[spdimen*2], aux ) )
    return false;
  if ( !_mbs_multiBCHornerDerf ( degc01, 1, spdimen, 0, c01, 0.0,
      &pcorners[spdimen*1], &pcorners[spdimen*3], aux ) )
    return false;

  if ( !_mbs_TabBezCurveDer3f ( spdimen, degc00, c00, nknu, knu, 2*spdimen,
            &c[0], &dc[0], &ddc[0], &dddc[0], aux ) )
    return false;
  if ( !_mbs_TabBezCurveDer3f ( spdimen, degc01, c01, nknu, knu, 2*spdimen,
            &c[spdimen], &dc[spdimen], &ddc[spdimen], &dddc[spdimen], aux ) )
    return false;
  if ( !_mbs_TabBezCurveDer3f ( spdimen, degd00, d00, nknv, knv, 2*spdimen,
            &d[0], &dd[0], &ddd[0], &dddd[0], aux ) )
    return false;
  if ( !_mbs_TabBezCurveDer3f ( spdimen, degd01, d01, nknv, knv, 2*spdimen,
            &d[spdimen], &dd[spdimen], &ddd[spdimen], &dddd[spdimen], aux ) )
    return false;

  if ( p )
    _mbs_TabBezC1Coons0f ( spdimen, nknu, nknv, c, d, pcorners,
                                 hfuncu, hfuncv, p, aux );
  if ( pu )
    _mbs_TabBezC1Coons0f ( spdimen, nknu, nknv, dc, d, pcorners,
                                 dhfuncu, hfuncv, pu, aux );
  if ( pv )
    _mbs_TabBezC1Coons0f ( spdimen, nknu, nknv, c, dd, pcorners,
                                 hfuncu, dhfuncv, pv, aux );
  if ( puu )
    _mbs_TabBezC1Coons0f ( spdimen, nknu, nknv, ddc, d, pcorners,
                                 ddhfuncu, hfuncv, puu, aux );
  if ( puv )
    _mbs_TabBezC1Coons0f ( spdimen, nknu, nknv, dc, dd, pcorners,
                                 dhfuncu, dhfuncv, puv, aux );
  if ( pvv )
    _mbs_TabBezC1Coons0f ( spdimen, nknu, nknv, c, ddd, pcorners,
                                 hfuncu, ddhfuncv, pvv, aux );
  if ( puuu )
    _mbs_TabBezC1Coons0f ( spdimen, nknu, nknv, dddc, d, pcorners,
                                 dddhfuncu, hfuncv, puuu, aux );
  if ( puuv )
    _mbs_TabBezC1Coons0f ( spdimen, nknu, nknv, ddc, dd, pcorners,
                                 ddhfuncu, dhfuncv, puuv, aux );
  if ( puvv )
    _mbs_TabBezC1Coons0f ( spdimen, nknu, nknv, dc, ddd, pcorners,
                                 dhfuncu, ddhfuncv, puvv, aux );
  if ( pvvv )
    _mbs_TabBezC1Coons0f ( spdimen, nknu, nknv, c, dddd, pcorners,
                                 hfuncu, dddhfuncv, pvvv, aux );
  return true;
} /*_mbs_TabBezC1Coons0Der3f*/

boolean mbs_TabBezC1Coons0Der3f ( int spdimen,
      int nknu, const float *knu, const float *hfuncu,
      const float *dhfuncu, const float *ddhfuncu, const float *dddhfuncu,
      int nknv, const float *knv, const float *hfuncv,
      const float *dhfuncv, const float *ddhfuncv, const float *dddhfuncv,
      int degc00, const float *c00,
      int degc01, const float *c01,
      int degd00, const float *d00,
      int degd01, const float *d01,
      float *p, float *pu, float *pv, float *puu, float *puv, float *pvv,
      float *puuu, float *puuv, float *puvv, float *pvvv )
{
  void   *sp;
  float  *workspace;
  boolean result;

  sp = workspace = pkv_GetScratchMemf ( (24*(nknu+nknv)+20)*spdimen );
  if ( !workspace )
    return false;
  result = _mbs_TabBezC1Coons0Der3f ( spdimen,
                nknu, knu, hfuncu, dhfuncu, ddhfuncu, dddhfuncu,
                nknv, knv, hfuncv, dhfuncv, ddhfuncv, dddhfuncv,
                degc00, c00, degc01, c01, degd00, d00, degd01, d01,
                p, pu, pv, puu, puv, pvv, puuu, puuv, puvv, pvvv, workspace );
  pkv_SetScratchMemTop ( sp );
  return result;
} /*mbs_TabBezC1Coons0Der3f*/

