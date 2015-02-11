
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

boolean _mbs_BezC2CoonsFindCornersf ( int spdimen,
                                      int degc00, const float *c00,
                                      int degc01, const float *c01,
                                      int degc02, const float *c02,
                                      int degc10, const float *c10,
                                      int degc11, const float *c11,
                                      int degc12, const float *c12,
                                      float *pcorners,
                                      float *workspace )
{
  if ( !_mbs_multiBCHornerDer2f ( degc00, 1, spdimen, 0, c00, 0.0,
        &pcorners[0], &pcorners[spdimen*12], &pcorners[spdimen*24], workspace ) )
    return false;
  if ( !_mbs_multiBCHornerDer2f ( degc00, 1, spdimen, 0, c00, 1.0,
        &pcorners[spdimen*6], &pcorners[spdimen*18], &pcorners[spdimen*30], workspace ) )
    return false;
  if ( !_mbs_multiBCHornerDer2f ( degc10, 1, spdimen, 0, c10, 0.0,
        &pcorners[spdimen*1], &pcorners[spdimen*13], &pcorners[spdimen*25], workspace ) )
    return false;
  if ( !_mbs_multiBCHornerDer2f ( degc10, 1, spdimen, 0, c10, 1.0,
        &pcorners[spdimen*7], &pcorners[spdimen*19], &pcorners[spdimen*31], workspace ) )
    return false;
  if ( !_mbs_multiBCHornerDer2f ( degc01, 1, spdimen, 0, c01, 0.0,
        &pcorners[spdimen*2], &pcorners[spdimen*14], &pcorners[spdimen*26], workspace ) )
    return false;
  if ( !_mbs_multiBCHornerDer2f ( degc01, 1, spdimen, 0, c01, 1.0,
        &pcorners[spdimen*8], &pcorners[spdimen*20], &pcorners[spdimen*32], workspace ) )
    return false;
  if ( !_mbs_multiBCHornerDer2f ( degc11, 1, spdimen, 0, c11, 0.0,
        &pcorners[spdimen*3], &pcorners[spdimen*15], &pcorners[spdimen*27], workspace ) )
    return false;
  if ( !_mbs_multiBCHornerDer2f ( degc11, 1, spdimen, 0, c11, 1.0,
        &pcorners[spdimen*9], &pcorners[spdimen*21], &pcorners[spdimen*33], workspace ) )
    return false;
  if ( !_mbs_multiBCHornerDer2f ( degc02, 1, spdimen, 0, c02, 0.0,
        &pcorners[spdimen*4], &pcorners[spdimen*16], &pcorners[spdimen*28], workspace ) )
    return false;
  if ( !_mbs_multiBCHornerDer2f ( degc02, 1, spdimen, 0, c02, 1.0,
        &pcorners[spdimen*10], &pcorners[spdimen*22], &pcorners[spdimen*34], workspace ) )
    return false;
  if ( !_mbs_multiBCHornerDer2f ( degc12, 1, spdimen, 0, c12, 0.0,
        &pcorners[spdimen*5], &pcorners[spdimen*17], &pcorners[spdimen*29], workspace ) )
    return false;
  if ( !_mbs_multiBCHornerDer2f ( degc12, 1, spdimen, 0, c12, 1.0,
        &pcorners[spdimen*11], &pcorners[spdimen*23], &pcorners[spdimen*35], workspace ) )
    return false;
  return true;
} /*_mbs_BezC2CoonsFindCornersf*/

boolean mbs_BezC2CoonsFindCornersf ( int spdimen,
                                     int degc00, const float *c00,
                                     int degc01, const float *c01,
                                     int degc02, const float *c02,
                                     int degc10, const float *c10,
                                     int degc11, const float *c11,
                                     int degc12, const float *c12,
                                     float *pcorners )
{
  void    *sp;
  float   *workspace;
  boolean result;

  sp = workspace = pkv_GetScratchMemf ( 3*spdimen );
  if ( !workspace )
    return false;
  result = _mbs_BezC2CoonsFindCornersf ( spdimen,
                   degc00, c00, degc01, c01, degc02, c02,
                   degc10, c10, degc11, c11, degc12, c12,
                   pcorners, workspace );
  pkv_SetScratchMemTop ( sp );
  return result;
} /*mbs_BezC2CoonsFindCornersf*/

static int FindMaxInt ( int a, int b, int c, int d, int e, int f, int g )
{
  a = max ( a, b );  a = max ( a, c );  a = max ( a, d );
  a = max ( a, e );  a = max ( a, f );  a = max ( a, g );
  return a;
} /*FindMaxInt*/

boolean _mbs_BezC2CoonsToBezf ( int spdimen,
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
                                int *n, int *m, float *p,
                                float *workspace )
{
  int   degu, degv, d, du, dv;
  int   pitch;
  float *pc, *bc, *aux;
  float *p1, *p2, *p3; 
  int   i;

  *n = degu = FindMaxInt ( degc00, degc01, degc02, degc10, degc11, degc12, 5 );
  *m = degv = FindMaxInt ( degd00, degd01, degd02, degd10, degd11, degd12, 5 );
  pc = workspace;
  p1 = &pc[36*spdimen];
  p2 = &p1[(degu+1)*(degv+1)*spdimen];
  p3 = &p2[(degu+1)*(degv+1)*spdimen];
  bc = &p3[(degu+1)*(degv+1)*spdimen];
  aux = &bc[6*(max(degu,degv)+1)*spdimen];
        /* construct the patch p1 */
  if ( !mbs_multiBCDegElevf ( 1, spdimen, 0, degc00, c00, degu-degc00,
                              0, &d, &bc[0] ) )
    return false;
  if ( !mbs_multiBCDegElevf ( 1, spdimen, 0, degc01, c01, degu-degc01,
                              0, &d, &bc[(degu+1)*spdimen] ) )
    return false;
  if ( !mbs_multiBCDegElevf ( 1, spdimen, 0, degc02, c02, degu-degc02,
                              0, &d, &bc[2*(degu+1)*spdimen] ) )
    return false;
  if ( !mbs_multiBCDegElevf ( 1, spdimen, 0, degc10, c10, degu-degc10,
                              0, &d, &bc[3*(degu+1)*spdimen] ) )
    return false;
  if ( !mbs_multiBCDegElevf ( 1, spdimen, 0, degc11, c11, degu-degc11,
                              0, &d, &bc[4*(degu+1)*spdimen] ) )
    return false;
  if ( !mbs_multiBCDegElevf ( 1, spdimen, 0, degc12, c12, degu-degc12,
                              0, &d, &bc[5*(degu+1)*spdimen] ) )
    return false;
  pitch = (degu+1)*spdimen;
  if ( !mbs_multiInterp2knHermiteBezf ( 1, pitch, 5, 3, pitch, bc,
           3, pitch, &bc[3*(degu+1)*spdimen], pitch, p1 ) )
    return false;
  pkv_TransposeMatrixc ( 6, degu+1, spdimen*sizeof(float),
                         (degu+1)*spdimen*sizeof(float), (char*)p1,
                         6*spdimen*sizeof(float), (char*)p );
  if ( !mbs_BCDegElevPf ( spdimen, degu, 5, p, 0, degv-5, &du, &dv, p ) )
    return false;

        /* construct the patch p2 */
  if ( !mbs_multiBCDegElevf ( 1, spdimen, 0, degd00, d00, degv-degd00,
                              0, &d, &bc[0] ) )
    return false;
  if ( !mbs_multiBCDegElevf ( 1, spdimen, 0, degd01, d01, degv-degd01,
                              0, &d, &bc[(degv+1)*spdimen] ) )
    return false;
  if ( !mbs_multiBCDegElevf ( 1, spdimen, 0, degd02, d02, degv-degd02,
                              0, &d, &bc[2*(degv+1)*spdimen] ) )
    return false;
  if ( !mbs_multiBCDegElevf ( 1, spdimen, 0, degd10, d10, degv-degd10,
                              0, &d, &bc[3*(degv+1)*spdimen] ) )
    return false;
  if ( !mbs_multiBCDegElevf ( 1, spdimen, 0, degd11, d11, degv-degd11,
                              0, &d, &bc[4*(degv+1)*spdimen] ) )
    return false;
  if ( !mbs_multiBCDegElevf ( 1, spdimen, 0, degd12, d12, degv-degd12,
                              0, &d, &bc[5*(degv+1)*spdimen] ) )
    return false;
  pitch = (degv+1)*spdimen;
  if ( !mbs_multiInterp2knHermiteBezf ( 1, pitch, 5, 3, pitch, bc,
           3, pitch, &bc[3*(degv+1)*spdimen], pitch, p2 ) )
    return false;
  if ( !mbs_BCDegElevPf ( spdimen, 5, degv, p2, degu-5, 0, &du, &dv, p2 ) )
    return false;
  pkn_AddMatrixf ( 1, spdimen*(degu+1)*(degv+1), 0, p, 0, p2, 0, p );

        /* construct the patch p3 */
  if ( !_mbs_BezC2CoonsFindCornersf ( spdimen,
                       degc00, c00, degc01, c01, degc02, c02,
                       degc10, c10, degc11, c11, degc12, c12, pc, aux ) )
    return false;

  pkv_Selectf ( 3, 6*spdimen, 6*2*spdimen, 6*spdimen, pc, p3 );
  pkv_Selectf ( 3, 6*spdimen, 6*2*spdimen, 6*spdimen, &pc[6*spdimen], &p3[6*3*spdimen] );
  for ( i = 0; i < 6; i++ ) {
    pkv_Selectf ( 3, spdimen, 2*spdimen, spdimen, &p3[6*i*spdimen], &bc[6*i*spdimen] );
    pkv_Selectf ( 3, spdimen, 2*spdimen, spdimen, &p3[(6*i+1)*spdimen], &bc[(6*i+3)*spdimen] );
  }
  if ( !mbs_multiInterp2knHermiteBezf ( 1, 6*spdimen, 5, 3, 0, bc,
                                        3, 0, &bc[6*3*spdimen], 0, p3 ) )
    return false;
  if ( !mbs_multiInterp2knHermiteBezf ( 6, spdimen, 5, 3, 6*spdimen, p3,
                                  3, 6*spdimen, &p3[3*spdimen], 6*spdimen, p2 ) )
    return false;
  if ( !mbs_BCDegElevPf ( spdimen, 5, 5, p2, degu-5, degv-5, &du, &dv, p3 ) )
    return false;
  pkn_SubtractMatrixf ( 1, spdimen*(degu+1)*(degv+1), 0, p, 0, p3, 0, p );
  return true;
} /*_mbs_BezC2CoonsToBezf*/

boolean mbs_BezC2CoonsToBezf ( int spdimen,
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
                               int *n, int *m, float *p )
{
  void    *sp;
  int     degu, degv;
  float   *workspace;
  boolean result;


  degu = FindMaxInt ( degc00, degc01, degc02, degc10, degc11, degc12, 5 );
  degv = FindMaxInt ( degd00, degd01, degd02, degd10, degd11, degd12, 5 );
  sp = workspace = pkv_GetScratchMemf ( (39 + 3*(degu+1)*(degv+1) +
                                         6*(max(degu,degv)+1))*spdimen );
  if ( !workspace )
    return false;
  result = _mbs_BezC2CoonsToBezf ( spdimen,
                   degc00, c00, degc01, c01, degc02, c02,
                   degc10, c10, degc11, c11, degc12, c12,
                   degd00, d00, degd01, d01, degd02, d02,
                   degd10, d10, degd11, d11, degd12, d12,
                   n, m, p, workspace );
  pkv_SetScratchMemTop ( sp );
  return result;
} /*mbs_BezC2CoonsToBezf*/

/* ///////////////////////////////////////////////////////////////////////// */
boolean mbs_TabQuinticHFuncDer3f ( float a, float b, int nkn, const float *kn,
                                   float *hfunc, float *dhfunc,
                                   float *ddhfunc, float *dddhfunc )
{
  float HFunc[36] = {1.0, 1.0, 1.0, 0.0, 0.0, 0.0,   /* h00 */
                     0.0, 0.0, 0.0, 1.0, 1.0, 1.0,   /* h10 */
                     0.0, 0.2, 0.4, 0.0, 0.0, 0.0,   /* h01 */
                     0.0, 0.0, 0.0,-0.4,-0.2, 0.0,   /* h11 */
                     0.0, 0.0,0.05, 0.0, 0.0, 0.0,   /* h02 */
                     0.0, 0.0, 0.0,0.05, 0.0, 0.0};  /* h12 */
  int   i, j;
  float h, h2, h_1, h_2, h_3;

  if ( a == 0.0 && b == 1.0 ) {
    for ( i = j = 0;  i < nkn;  i++, j += 6 )
      if ( !mbs_multiBCHornerDer3f ( 5, 6, 1, 6, HFunc, kn[i],
              &hfunc[j], &dhfunc[j], &ddhfunc[j], &dddhfunc[j] ) )
        return false;
  }
  else {
    h = b - a;    
    h2 = h*h;     
    h_1 = 1.0/h;  
    h_2 = 1.0/h2;
    h_3 = h_1*h_2;
    for ( i = j = 0;  i < nkn;  i++, j += 6 ) {
      if ( !mbs_multiBCHornerDer3f ( 5, 6, 1, 6, HFunc, kn[i],
               &hfunc[j], &dhfunc[j], &ddhfunc[j], &dddhfunc[j] ) )
        return false;
      hfunc[j+2] *= h;       hfunc[j+3] *= h;
      hfunc[j+4] *= h2;      hfunc[j+5] *= h2;
      dhfunc[j] *= h_1;      dhfunc[j+1] *= h_1;
      dhfunc[j+4] *= h;      dhfunc[j+5] *= h;  
      ddhfunc[j] *= h_2;     ddhfunc[j+1] *= h_2;
      ddhfunc[j+2] *= h_1;   ddhfunc[j+3] *= h_1;
      dddhfunc[j] *= h_3;    dddhfunc[j+1] *= h_3;
      dddhfunc[j+2] *= h_2;  dddhfunc[j+3] *= h_2;
      dddhfunc[j+4] *= h_1;  dddhfunc[j+5] *= h_1;
    }
  }
  return true;
} /*mbs_TabQuinticHFuncDer3f*/

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

boolean mbs_TabBezCurveDer3f ( int spdimen, int degree, const float *cp,
                               int nkn, const float *kn,
                               int ppitch,
                               float *p, float *dp, float *ddp, float *dddp )
{
  void    *sp;
  float   *workspace;
  boolean result;

  sp = workspace = pkv_GetScratchMemf ( 4*spdimen );
  if ( !workspace )
    return false;
  result = _mbs_TabBezCurveDer3f ( spdimen, degree, cp, nkn, kn,
                                   ppitch, p, dp, ddp, dddp, workspace );
  pkv_SetScratchMemTop ( sp );
  return result;
} /*mbs_TabBezCurveDer3f*/

void _mbs_TabBezC2Coonsf ( int spdimen, int nknu, int nknv,
                           const float *c, const float *d, const float *p,
                           const float *hu, const float *hv, float *pp,
                           float *workspace )
{
  int i, j, k, l;

  memcpy ( workspace, d, 6*nknv*spdimen*sizeof(float) );
  for ( i = 0; i < nknv; i++ )
    for ( j = 0; j < 6; j++ )
      for ( k = 0; k < 6; k++ )
        for ( l = 0; l < spdimen; l++ )
          workspace[(6*i+j)*spdimen+l] -= hv[6*i+k]*p[(6*j+k)*spdimen+l];

  memset ( pp, 0, nknu*nknv*spdimen*sizeof(float) );
  for ( i = 0; i < nknu; i++ )
    for ( j = 0; j < nknv; j++ )
      for ( k = 0; k < 6; k++ )
        for ( l = 0; l < spdimen; l++ )
          pp[(nknv*i+j)*spdimen+l] += c[(6*i+k)*spdimen+l]*hv[6*j+k] +
                                      workspace[(6*j+k)*spdimen+l]*hu[6*i+k];
} /*_mbs_TabBezC2Coonsf*/

boolean _mbs_TabBezC2CoonsDer3f ( int spdimen,
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
  aux = &pcorners[36*spdimen];

  if ( !_mbs_BezC2CoonsFindCornersf ( spdimen,
              degc00, c00, degc01, c01, degc02, c02,
              degc10, c10, degc11, c11, degc12, c12,
              pcorners, aux ) )
    return false;

  if ( !_mbs_TabBezCurveDer3f ( spdimen, degc00, c00, nknu, knu, 6*spdimen,
            &c[0], &dc[0], &ddc[0], &dddc[0], aux ) )
    return false;
  if ( !_mbs_TabBezCurveDer3f ( spdimen, degc10, c10, nknu, knu, 6*spdimen,
            &c[spdimen], &dc[spdimen], &ddc[spdimen], &dddc[spdimen], aux ) )
    return false;
  if ( !_mbs_TabBezCurveDer3f ( spdimen, degc01, c01, nknu, knu, 6*spdimen,
            &c[2*spdimen], &dc[2*spdimen], &ddc[2*spdimen], &dddc[2*spdimen], aux ) )
    return false;
  if ( !_mbs_TabBezCurveDer3f ( spdimen, degc11, c11, nknu, knu, 6*spdimen,
            &c[3*spdimen], &dc[3*spdimen], &ddc[3*spdimen], &dddc[3*spdimen], aux ) )
    return false;
  if ( !_mbs_TabBezCurveDer3f ( spdimen, degc02, c02, nknu, knu, 6*spdimen,
            &c[4*spdimen], &dc[4*spdimen], &ddc[4*spdimen], &dddc[4*spdimen], aux ) )
    return false;
  if ( !_mbs_TabBezCurveDer3f ( spdimen, degc12, c12, nknu, knu, 6*spdimen,
            &c[5*spdimen], &dc[5*spdimen], &ddc[5*spdimen], &dddc[5*spdimen], aux ) )
    return false;

  if ( !_mbs_TabBezCurveDer3f ( spdimen, degd00, d00, nknv, knv, 6*spdimen,
            &d[0], &dd[0], &ddd[0], &dddd[0], aux ) )
    return false;
  if ( !_mbs_TabBezCurveDer3f ( spdimen, degd10, d10, nknv, knv, 6*spdimen,
            &d[spdimen], &dd[spdimen], &ddd[spdimen], &dddd[spdimen], aux ) )
    return false;
  if ( !_mbs_TabBezCurveDer3f ( spdimen, degd01, d01, nknv, knv, 6*spdimen,
            &d[2*spdimen], &dd[2*spdimen], &ddd[2*spdimen], &dddd[2*spdimen], aux ) )
    return false;
  if ( !_mbs_TabBezCurveDer3f ( spdimen, degd11, d11, nknv, knv, 6*spdimen,
            &d[3*spdimen], &dd[3*spdimen], &ddd[3*spdimen], &dddd[3*spdimen], aux ) )
    return false;
  if ( !_mbs_TabBezCurveDer3f ( spdimen, degd02, d02, nknv, knv, 6*spdimen,
            &d[4*spdimen], &dd[4*spdimen], &ddd[4*spdimen], &dddd[4*spdimen], aux ) )
    return false;
  if ( !_mbs_TabBezCurveDer3f ( spdimen, degd12, d12, nknv, knv, 6*spdimen,
            &d[5*spdimen], &dd[5*spdimen], &ddd[5*spdimen], &dddd[5*spdimen], aux ) )
    return false;

  if ( p )
    _mbs_TabBezC2Coonsf ( spdimen, nknu, nknv, c, d, pcorners,
                          hfuncu, hfuncv, p, aux );
  if ( pu )
    _mbs_TabBezC2Coonsf ( spdimen, nknu, nknv, dc, d, pcorners,
                          dhfuncu, hfuncv, pu, aux );
  if ( pv )
    _mbs_TabBezC2Coonsf ( spdimen, nknu, nknv, c, dd, pcorners,
                          hfuncu, dhfuncv, pv, aux );
  if ( puu )
    _mbs_TabBezC2Coonsf ( spdimen, nknu, nknv, ddc, d, pcorners,
                          ddhfuncu, hfuncv, puu, aux );
  if ( puv )
    _mbs_TabBezC2Coonsf ( spdimen, nknu, nknv, dc, dd, pcorners,
                          dhfuncu, dhfuncv, puv, aux );
  if ( pvv )
    _mbs_TabBezC2Coonsf ( spdimen, nknu, nknv, c, ddd, pcorners,
                          hfuncu, ddhfuncv, pvv, aux );
  if ( puuu )
    _mbs_TabBezC2Coonsf ( spdimen, nknu, nknv, dddc, d, pcorners,
                          dddhfuncu, hfuncv, puuu, aux );
  if ( puuv )
    _mbs_TabBezC2Coonsf ( spdimen, nknu, nknv, ddc, dd, pcorners,
                          ddhfuncu, dhfuncv, puuv, aux );
  if ( puvv )
    _mbs_TabBezC2Coonsf ( spdimen, nknu, nknv, dc, ddd, pcorners,
                          dhfuncu, ddhfuncv, puvv, aux );
  if ( pvvv )
    _mbs_TabBezC2Coonsf ( spdimen, nknu, nknv, c, dddd, pcorners,
                          hfuncu, dddhfuncv, pvvv, aux );
  return true;
} /*_mbs_TabBezC2CoonsDer3f*/

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

void _mbs_TabBezC2Coons0f ( int spdimen, int nknu, int nknv,
                            const float *c, const float *d, const float *p,
                            const float *hu, const float *hv, float *pp,
                            float *workspace )
{
  int    i, j, k, l;

  memcpy ( workspace, d, 3*nknv*spdimen*sizeof(float) );
  for ( i = 0; i < nknv; i++ )
    for ( j = 0; j < 3; j++ )
      for ( k = 0; k < 3; k++ )
        for ( l = 0; l < spdimen; l++ )
          workspace[(3*i+j)*spdimen+l] -= hv[6*i+2*k]*p[(3*j+k)*spdimen+l];

  memset ( pp, 0, nknu*nknv*spdimen*sizeof(float) );
  for ( i = 0; i < nknu; i++ )
    for ( j = 0; j < nknv; j++ )
      for ( k = 0; k < 3; k++ )
        for ( l = 0; l < spdimen; l++ )
          pp[(nknv*i+j)*spdimen+l] += c[(3*i+k)*spdimen+l]*hv[6*j+2*k] +
                                      workspace[(3*j+k)*spdimen+l]*hu[6*i+2*k];
} /*_mbs_TabBezC2Coons0f*/

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

boolean mbs_TabBezC2Coons0Der3f ( int spdimen,
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
      float *puuu, float *puuv, float *puvv, float *pvvv )
{
  void   *sp;
  float *workspace;
  boolean result;

  sp = workspace = pkv_GetScratchMemf ( (12*(nknu+nknv) +
                                         3*max(nknu,nknv) + 9)*spdimen );
  if ( !workspace )
    return false;
  result = _mbs_TabBezC2Coons0Der3f ( spdimen,
                nknu, knu, hfuncu, dhfuncu, ddhfuncu, dddhfuncu,
                nknv, knv, hfuncv, dhfuncv, ddhfuncv, dddhfuncv,
                degc00, c00, degc01, c01, degc02, c02,
                degd00, d00, degd01, d01, degd02, d02,
                p, pu, pv, puu, puv, pvv, puuu, puuv, puvv, pvvv, workspace );
  pkv_SetScratchMemTop ( sp );
  return result;
} /*mbs_TabBezC2Coons0Der3f*/

