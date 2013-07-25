
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2005, 2009                            */
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

void mbs_BezC1CoonsFindCornersf ( int spdimen,
                                  int degc00, const float *c00,
                                  int degc01, const float *c01,
                                  int degc10, const float *c10,
                                  int degc11, const float *c11,
                                  float *pcorners )
{
  mbs_multiBCHornerDerf ( degc00, 1, spdimen, 0, c00, 0.0,
                          &pcorners[0], &pcorners[spdimen*8] );
  mbs_multiBCHornerDerf ( degc00, 1, spdimen, 0, c00, 1.0,
                          &pcorners[spdimen*4], &pcorners[spdimen*12] );
  mbs_multiBCHornerDerf ( degc10, 1, spdimen, 0, c10, 0.0,
                          &pcorners[spdimen*1], &pcorners[spdimen*9] );
  mbs_multiBCHornerDerf ( degc10, 1, spdimen, 0, c10, 1.0,
                          &pcorners[spdimen*5], &pcorners[spdimen*13] );
  mbs_multiBCHornerDerf ( degc01, 1, spdimen, 0, c01, 0.0,
                          &pcorners[spdimen*2], &pcorners[spdimen*10] );
  mbs_multiBCHornerDerf ( degc01, 1, spdimen, 0, c01, 1.0,
                          &pcorners[spdimen*6], &pcorners[spdimen*14] );
  mbs_multiBCHornerDerf ( degc11, 1, spdimen, 0, c11, 0.0,
                          &pcorners[spdimen*3], &pcorners[spdimen*11] );
  mbs_multiBCHornerDerf ( degc11, 1, spdimen, 0, c11, 1.0,
                          &pcorners[spdimen*7], &pcorners[spdimen*15] );
} /*mbs_BezC1CoonsFindCornersf*/

static int FindMaxInt ( int a, int b, int c, int d, int e )
{
  a = max ( a, b );  a = max ( a, c );  a = max ( a, d );  a = max ( a, e );
  return a;
} /*FindMaxInt*/

/*
static void Verify ( int spdimen,
                     int degc00, const float *c00,
                     int degc01, const float *c01,
                     int degc10, const float *c10,
                     int degc11, const float *c11,
                     int degd00, const float *d00,
                     int degd01, const float *d01,
                     int degd10, const float *d10,
                     int degd11, const float *d11 )
{
  void  *sp;
  float *cc, *cd, *dd;
  int   i, j, k;

  sp = pkv_GetScratchMemTop ();
  cc = pkv_GetScratchMemf ( 16*spdimen );
  cd = pkv_GetScratchMemf ( 16*spdimen );
  dd = pkv_GetScratchMemf ( 16*spdimen );

  mbs_BezC1CoonsFindCornersf ( spdimen, degc00, c00, degc01, c01,
                                 degc10, c10, degc11, c11, cc );
  mbs_BezC1CoonsFindCornersf ( spdimen, degd00, d00, degd01, d01,
                                 degd10, d10, degd11, d11, cd );
  pkv_TransposeMatrixc ( 4, 4, spdimen*sizeof(float),
                         4*spdimen*sizeof(float), (char*)cd,
                         4*spdimen*sizeof(float), (char*)dd );
  for ( i = 0; i < 4; i++ ) {
    for ( j = 0; j < 4; j++ )
      for ( k = 0; k < spdimen; k++ )
        printf ( "%7.4f ", cc[(4*i+j)*spdimen+k] );
    printf ( "\n" );
  }
  printf ( "\n" );
  for ( i = 0; i < 4; i++ ) {
    for ( j = 0; j < 4; j++ )
      for ( k = 0; k < spdimen; k++ )
        printf ( "%7.4f ", dd[(4*i+j)*spdimen+k] );
    printf ( "\n" );
  }
  printf ( "\n" );
  for ( i = 0; i < 4; i++ ) {
    for ( j = 0; j < 4; j++ )
      for ( k = 0; k < spdimen; k++ )
        printf ( "%7.4f ", dd[(4*i+j)*spdimen+k]-cc[(4*i+j)*spdimen+k] );
    printf ( "\n" );
  }

  pkv_SetScratchMemTop ( sp );
} / *Verify*/

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
  void  *sp;
  int   degu, degv, d, du, dv;
  int   pitch;
  float *pc, *bc;
  float *p1, *p2, *p3; 
  int   i;

  sp = pkv_GetScratchMemTop ();
  *n = degu = FindMaxInt ( degc00, degc01, degc10, degc11, 3 );
  *m = degv = FindMaxInt ( degd00, degd01, degd10, degd11, 3 );
  pc = pkv_GetScratchMemf ( 16*spdimen );
  p1 = pkv_GetScratchMemf ( 3*(degu+1)*(degv+1)*spdimen );
  bc = pkv_GetScratchMemf ( 4*(max(degu,degv)+1)*spdimen );
  if ( !pc || !p1 || !bc )
    goto failure;
  p2 = &p1[(degu+1)*(degv+1)*spdimen];
  p3 = &p2[(degu+1)*(degv+1)*spdimen];

        /* construct the patch p1 */
  mbs_multiBCDegElevf ( 1, spdimen, 0, degc00, c00, degu-degc00,
                        0, &d, &bc[0] );
  mbs_multiBCDegElevf ( 1, spdimen, 0, degc01, c01, degu-degc01,
                        0, &d, &bc[(degu+1)*spdimen] );
  mbs_multiBCDegElevf ( 1, spdimen, 0, degc10, c10, degu-degc10,
                        0, &d, &bc[2*(degu+1)*spdimen] );
  mbs_multiBCDegElevf ( 1, spdimen, 0, degc11, c11, degu-degc11,
                        0, &d, &bc[3*(degu+1)*spdimen] );
  pitch = (degu+1)*spdimen;
  mbs_multiInterp2knHermiteBezf ( 1, pitch, 3, 2, pitch, bc,
           2, pitch, &bc[2*(degu+1)*spdimen], pitch, p1 );
  pkv_TransposeMatrixc ( 4, degu+1, spdimen*sizeof(float),
                         (degu+1)*spdimen*sizeof(float), (char*)p1,
                         4*spdimen*sizeof(float), (char*)p );
  mbs_BCDegElevPf ( spdimen, degu, 3, p, 0, degv-3, &du, &dv, p );

        /* construct the patch p2 */
  mbs_multiBCDegElevf ( 1, spdimen, 0, degd00, d00, degv-degd00,
                        0, &d, &bc[0] );
  mbs_multiBCDegElevf ( 1, spdimen, 0, degd01, d01, degv-degd01,
                        0, &d, &bc[(degv+1)*spdimen] );
  mbs_multiBCDegElevf ( 1, spdimen, 0, degd10, d10, degv-degd10,
                        0, &d, &bc[2*(degv+1)*spdimen] );
  mbs_multiBCDegElevf ( 1, spdimen, 0, degd11, d11, degv-degd11,
                        0, &d, &bc[3*(degv+1)*spdimen] );
  pitch = (degv+1)*spdimen;
  mbs_multiInterp2knHermiteBezf ( 1, pitch, 3, 2, pitch, bc,
           2, pitch, &bc[2*(degv+1)*spdimen], pitch, p2 );
  mbs_BCDegElevPf ( spdimen, 3, degv, p2, degu-3, 0, &du, &dv, p2 );
  pkn_AddMatrixf ( 1, spdimen*(degu+1)*(degv+1), 0, p, 0, p2, 0, p );

        /* construct the patch p3 */
  mbs_BezC1CoonsFindCornersf ( spdimen, degc00, c00, degc01, c01,
                                 degc10, c10, degc11, c11, pc );

  pkv_Selectf ( 2, 4*spdimen, 4*2*spdimen, 4*spdimen, pc, p3 );
  pkv_Selectf ( 2, 4*spdimen, 4*2*spdimen, 4*spdimen, &pc[4*spdimen], &p3[4*2*spdimen] );
  for ( i = 0; i < 4; i++ ) {
    pkv_Selectf ( 2, spdimen, 2*spdimen, spdimen, &p3[4*i*spdimen], &bc[4*i*spdimen] );
    pkv_Selectf ( 2, spdimen, 2*spdimen, spdimen, &p3[(4*i+1)*spdimen], &bc[(4*i+2)*spdimen] );
  }
  mbs_multiInterp2knHermiteBezf ( 1, 4*spdimen, 3, 2, 0, bc,
                                  2, 0, &bc[4*2*spdimen], 0, p3 );
  mbs_multiInterp2knHermiteBezf ( 4, spdimen, 3, 2, 4*spdimen, p3,
                                  2, 4*spdimen, &p3[2*spdimen], 4*spdimen, p2 );
  mbs_BCDegElevPf ( spdimen, 3, 3, p2, degu-3, degv-3, &du, &dv, p3 );
  pkn_SubtractMatrixf ( 1, spdimen*(degu+1)*(degv+1), 0, p, 0, p3, 0, p );

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*mbs_BezC1CoonsToBezf*/

/* ///////////////////////////////////////////////////////////////////////// */
void mbs_TabCubicHFuncDer2f ( float a, float b, int nkn, const float *kn,
                              float *hfunc, float *dhfunc, float *ddhfunc )
{
  float HFunc[16] = {1.0, 1.0, 0.0, 0.0,                  /* h00 */
                     0.0, 0.0, 1.0, 1.0,                  /* h10 */
                     0.0, (float)(1.0/3.0), 0.0, 0.0,     /* h01 */
                     0.0, 0.0, (float)(-1.0/3.0), 0.0};   /* h11 */
  int   i, j;
  float h, h_1, h_2;

  if ( a == 0.0 && b == 1.0 ) {
     for ( i = j = 0;  i < nkn;  i++, j += 4 )
       mbs_multiBCHornerDer2f ( 3, 4, 1, 4, HFunc, kn[i],
                                &hfunc[j], &dhfunc[j], &ddhfunc[j] );
  }
  else {
    h = b - a;
    h_1 = (float)1.0/h;
    h_2 = h_1*h_1;
    for ( i = j = 0;  i < nkn;  i++, j += 4 ) {
      mbs_multiBCHornerDer2f ( 3, 4, 1, 4, HFunc, (kn[i]-a)/h,
                               &hfunc[j], &dhfunc[j], &ddhfunc[j] );
      hfunc[j+2]   *= h;      hfunc[j+3]   *= h;
      dhfunc[j]    *= h_1;    dhfunc[j+1]  *= h_1;
      ddhfunc[j]   *= h_2;    ddhfunc[j+1] *= h_2;
      ddhfunc[j+2] *= h_1;    ddhfunc[j+3] *= h_1;
    }
  }
} /*mbs_TabCubicHFuncDer2f*/

void mbs_TabCubicHFuncDer3f ( float a, float b, int nkn, const float *kn,
                              float *hfunc, float *dhfunc, float *ddhfunc,
                              float *dddhfunc )
{
  float HFunc[16] = {1.0, 1.0, 0.0, 0.0,                  /* h00 */
                     0.0, 0.0, 1.0, 1.0,                  /* h10 */
                     0.0, (float)(1.0/3.0), 0.0, 0.0,     /* h01 */
                     0.0, 0.0, (float)(-1.0/3.0), 0.0};   /* h11 */
  int   i, j;
  float h, h_1, h_2, h_3;

  if ( a == 0.0 && b == 1.0 ) {
     for ( i = j = 0;  i < nkn;  i++, j += 4 )
       mbs_multiBCHornerDer3f ( 3, 4, 1, 4, HFunc, kn[i],
                            &hfunc[j], &dhfunc[j], &ddhfunc[j], &dddhfunc[j] );
  }
  else {
    h = b - a;
    h_1 = (float)1.0/h;
    h_2 = h_1*h_1;
    h_3 = h_1*h_2;
    for ( i = j = 0;  i < nkn;  i++, j += 4 ) {
      mbs_multiBCHornerDer3f ( 3, 4, 1, 4, HFunc, (kn[i]-a)/h,
                           &hfunc[j], &dhfunc[j], &ddhfunc[j], &dddhfunc[j] );
      hfunc[j+2]    *= h;      hfunc[j+3]   *= h;
      dhfunc[j]     *= h_1;    dhfunc[j+1]  *= h_1;
      ddhfunc[j]    *= h_2;    ddhfunc[j+1] *= h_2;
      ddhfunc[j+2]  *= h_1;    ddhfunc[j+3] *= h_1;
      dddhfunc[j]   *= h_3;    dddhfunc[j+1] *= h_3;
      dddhfunc[j+2] *= h_2;    dddhfunc[j+3] *= h_2;
    }
  }
} /*mbs_TabCubicHFuncDer3f*/

void mbs_TabBezCurveDer2f ( int spdimen, int degree, const float *cp,
                            int nkn, const float *kn,
                            int ppitch,
                            float *p, float *dp, float *ddp )
{
  int i, j;

  for ( i = j = 0;  i < nkn;  i++, j += ppitch )
    mbs_multiBCHornerDer2f ( degree, 1, spdimen, 0, cp, kn[i],
                             &p[j], &dp[j], &ddp[j] );
} /*mbs_TabBezCurveDer2f*/

boolean _mbs_TabBezC1Coonsf ( int spdimen, int nknu, int nknv,
                   const float *c, const float *d, const float *p,
                   const float *hu, const float *hv, float *pp )
{
  void  *sp;
  int   i, j, k, l;
  float *a;

  sp = pkv_GetScratchMemTop ();
  a = pkv_GetScratchMemf ( 4*nknv*spdimen );
  if ( !a )
    goto failure;

  memcpy ( a, d, 4*nknv*spdimen*sizeof(float) );
  for ( i = 0; i < nknv; i++ )
    for ( j = 0; j < 4; j++ )
      for ( k = 0; k < 4; k++ )
        for ( l = 0; l < spdimen; l++ )
          a[(4*i+j)*spdimen+l] -= hv[4*i+k]*p[(4*j+k)*spdimen+l];

  memset ( pp, 0, nknu*nknv*spdimen*sizeof(float) );
  for ( i = 0; i < nknu; i++ )
    for ( j = 0; j < nknv; j++ )
      for ( k = 0; k < 4; k++ )
        for ( l = 0; l < spdimen; l++ )
          pp[(nknv*i+j)*spdimen+l] += c[(4*i+k)*spdimen+l]*hv[4*j+k] +
                                      a[(4*j+k)*spdimen+l]*hu[4*i+k];

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*_mbs_TabBezC1Coonsf*/

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
  void  *sp;
  int   ku, kv;
  float *c, *dc, *ddc, *d, *dd, *ddd, *pcorners;

  sp = pkv_GetScratchMemTop ();
  ku = 4*spdimen*nknu;
  kv = 4*spdimen*nknv;
  c = pkv_GetScratchMemf ( (18*(nknu+nknv)+16)*spdimen );
  if ( !c )
    goto failure;
                 dc = &c[ku];  ddc = &dc[ku];
  d = &ddc[ku];  dd = &d[kv];  ddd = &dd[kv];
  pcorners = &ddd[kv];

  mbs_BezC1CoonsFindCornersf ( spdimen, degc00, c00, degc01, c01,
                               degc10, c10, degc11, c11, pcorners );

  mbs_TabBezCurveDer2f ( spdimen, degc00, c00, nknu, knu, 4*spdimen,
            &c[0], &dc[0], &ddc[0] );
  mbs_TabBezCurveDer2f ( spdimen, degc10, c10, nknu, knu, 4*spdimen,
            &c[spdimen], &dc[spdimen], &ddc[spdimen] );
  mbs_TabBezCurveDer2f ( spdimen, degc01, c01, nknu, knu, 4*spdimen,
            &c[2*spdimen], &dc[2*spdimen], &ddc[2*spdimen] );
  mbs_TabBezCurveDer2f ( spdimen, degc11, c11, nknu, knu, 4*spdimen,
            &c[3*spdimen], &dc[3*spdimen], &ddc[3*spdimen] );

  mbs_TabBezCurveDer2f ( spdimen, degd00, d00, nknv, knv, 4*spdimen,
            &d[0], &dd[0], &ddd[0] );
  mbs_TabBezCurveDer2f ( spdimen, degd10, d10, nknv, knv, 4*spdimen,
            &d[spdimen], &dd[spdimen], &ddd[spdimen] );
  mbs_TabBezCurveDer2f ( spdimen, degd01, d01, nknv, knv, 4*spdimen,
            &d[2*spdimen], &dd[2*spdimen], &ddd[2*spdimen] );
  mbs_TabBezCurveDer2f ( spdimen, degd11, d11, nknv, knv, 4*spdimen,
            &d[3*spdimen], &dd[3*spdimen], &ddd[3*spdimen] );

  if ( p )
    if ( !_mbs_TabBezC1Coonsf ( spdimen, nknu, nknv, c, d, pcorners,
                                hfuncu, hfuncv, p ) )
      goto failure;
  if ( pu )
    if ( !_mbs_TabBezC1Coonsf ( spdimen, nknu, nknv, dc, d, pcorners,
                                dhfuncu, hfuncv, pu ) )
      goto failure;
  if ( pv )
    if ( !_mbs_TabBezC1Coonsf ( spdimen, nknu, nknv, c, dd, pcorners,
                                hfuncu, dhfuncv, pv ) )
      goto failure;
  if ( puu )
    if ( !_mbs_TabBezC1Coonsf ( spdimen, nknu, nknv, ddc, d, pcorners,
                                ddhfuncu, hfuncv, puu ) )
      goto failure;
  if ( puv )
    if ( !_mbs_TabBezC1Coonsf ( spdimen, nknu, nknv, dc, dd, pcorners,
                                dhfuncu, dhfuncv, puv ) )
      goto failure;
  if ( pvv )
    if ( !_mbs_TabBezC1Coonsf ( spdimen, nknu, nknv, c, ddd, pcorners,
                                hfuncu, ddhfuncv, pvv ) )
      goto failure;

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*mbs_TabBezC1CoonsDer2f*/

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
  void  *sp;
  int   ku, kv;
  float *c, *dc, *ddc, *dddc, *d, *dd, *ddd, *dddd, *pcorners;

  sp = pkv_GetScratchMemTop ();
  ku = 6*spdimen*nknu;
  kv = 6*spdimen*nknv;
  c = pkv_GetScratchMemf ( (24*(nknu+nknv)+16)*spdimen );
  if ( !c )
    goto failure;
                  dc = &c[ku];  ddc = &dc[ku];  dddc = &ddc[ku];
  d = &dddc[ku];  dd = &d[kv];  ddd = &dd[kv];  dddd = &ddd[kv];
  pcorners = &dddd[kv];

  mbs_BezC1CoonsFindCornersf ( spdimen, degc00, c00, degc01, c01,
                               degc10, c10, degc11, c11, pcorners );

  mbs_TabBezCurveDer3f ( spdimen, degc00, c00, nknu, knu, 4*spdimen,
            &c[0], &dc[0], &ddc[0], &dddc[0] );
  mbs_TabBezCurveDer3f ( spdimen, degc10, c10, nknu, knu, 4*spdimen,
            &c[spdimen], &dc[spdimen], &ddc[spdimen], &dddc[spdimen] );
  mbs_TabBezCurveDer3f ( spdimen, degc01, c01, nknu, knu, 4*spdimen,
            &c[2*spdimen], &dc[2*spdimen], &ddc[2*spdimen], &dddc[2*spdimen] );
  mbs_TabBezCurveDer3f ( spdimen, degc11, c11, nknu, knu, 4*spdimen,
            &c[3*spdimen], &dc[3*spdimen], &ddc[3*spdimen], &dddc[3*spdimen] );

  mbs_TabBezCurveDer3f ( spdimen, degd00, d00, nknv, knv, 4*spdimen,
            &d[0], &dd[0], &ddd[0], &dddd[0] );
  mbs_TabBezCurveDer3f ( spdimen, degd10, d10, nknv, knv, 4*spdimen,
            &d[spdimen], &dd[spdimen], &ddd[spdimen], &dddd[spdimen] );
  mbs_TabBezCurveDer3f ( spdimen, degd01, d01, nknv, knv, 4*spdimen,
            &d[2*spdimen], &dd[2*spdimen], &ddd[2*spdimen], &dddd[2*spdimen] );
  mbs_TabBezCurveDer3f ( spdimen, degd11, d11, nknv, knv, 4*spdimen,
            &d[3*spdimen], &dd[3*spdimen], &ddd[3*spdimen], &dddd[3*spdimen] );

  if ( p )
    if ( !_mbs_TabBezC1Coonsf ( spdimen, nknu, nknv, c, d, pcorners,
                                hfuncu, hfuncv, p ) )
      goto failure;
  if ( pu )
    if ( !_mbs_TabBezC1Coonsf ( spdimen, nknu, nknv, dc, d, pcorners,
                                dhfuncu, hfuncv, pu ) )
      goto failure;
  if ( pv )
    if ( !_mbs_TabBezC1Coonsf ( spdimen, nknu, nknv, c, dd, pcorners,
                                hfuncu, dhfuncv, pv ) )
      goto failure;
  if ( puu )
    if ( !_mbs_TabBezC1Coonsf ( spdimen, nknu, nknv, ddc, d, pcorners,
                                ddhfuncu, hfuncv, puu ) )
      goto failure;
  if ( puv )
    if ( !_mbs_TabBezC1Coonsf ( spdimen, nknu, nknv, dc, dd, pcorners,
                                dhfuncu, dhfuncv, puv ) )
      goto failure;
  if ( pvv )
    if ( !_mbs_TabBezC1Coonsf ( spdimen, nknu, nknv, c, ddd, pcorners,
                                hfuncu, ddhfuncv, pvv ) )
      goto failure;
  if ( puuu )
    if ( !_mbs_TabBezC1Coonsf ( spdimen, nknu, nknv, dddc, d, pcorners,
                                dddhfuncu, hfuncv, puuu ) )
      goto failure;
  if ( puuv )
    if ( !_mbs_TabBezC1Coonsf ( spdimen, nknu, nknv, ddc, dd, pcorners,
                                ddhfuncu, dhfuncv, puuv ) )
      goto failure;
  if ( puvv )
    if ( !_mbs_TabBezC1Coonsf ( spdimen, nknu, nknv, dc, ddd, pcorners,
                                dhfuncu, ddhfuncv, puvv ) )
      goto failure;
  if ( pvvv )
    if ( !_mbs_TabBezC1Coonsf ( spdimen, nknu, nknv, c, dddd, pcorners,
                                hfuncu, dddhfuncv, pvvv ) )
      goto failure;

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*mbs_TabBezC1CoonsDer3f*/

boolean _mbs_TabBezC1Coons0f (
                   int spdimen, int nknu, int nknv,
                   const float *c, const float *d, const float *p,
                   const float *hu, const float *hv, float *pp )
{
  void  *sp;
  int   i, j, k, l;
  float *a;

  sp = pkv_GetScratchMemTop ();
  a = pkv_GetScratchMemf ( 2*nknv*spdimen );
  if ( !a )
    goto failure;

  memcpy ( a, d, 2*nknv*spdimen*sizeof(float) );
  for ( i = 0; i < nknv; i++ )
    for ( j = 0; j < 2; j++ )
      for ( k = 0; k < 2; k++ )
        for ( l = 0; l < spdimen; l++ )
          a[(2*i+j)*spdimen+l] -= hv[4*i+2*k]*p[(2*j+k)*spdimen+l];

  memset ( pp, 0, nknu*nknv*spdimen*sizeof(float) );
  for ( i = 0; i < nknu; i++ )
    for ( j = 0; j < nknv; j++ )
      for ( k = 0; k < 2; k++ )
        for ( l = 0; l < spdimen; l++ )
          pp[(nknv*i+j)*spdimen+l] += c[(2*i+k)*spdimen+l]*hv[4*j+2*k] +
                                     a[(2*j+k)*spdimen+l]*hu[4*i+2*k];

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*_mbs_TabBezC1Coons0f*/

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
  void *sp;
  int  ku, kv;
  float *c, *dc, *ddc, *d, *dd, *ddd, *pcorners;

  sp = pkv_GetScratchMemTop ();
  ku = 2*spdimen*nknu;
  kv = 2*spdimen*nknv;
  c = pkv_GetScratchMemf ( (9*(nknu+nknv)+4)*spdimen );
  if ( !c )
    goto failure;
                 dc = &c[ku];  ddc = &dc[ku];
  d = &ddc[ku];  dd = &d[kv];  ddd = &dd[kv];
  pcorners = &ddd[kv];

  mbs_multiBCHornerDerf ( degc00, 1, spdimen, 0, c00, 0.0,
      &pcorners[0], &pcorners[spdimen*2] );
  mbs_multiBCHornerDerf ( degc01, 1, spdimen, 0, c01, 0.0,
      &pcorners[spdimen*1], &pcorners[spdimen*3] );

  mbs_TabBezCurveDer2f ( spdimen, degc00, c00, nknu, knu, 2*spdimen,
            &c[0], &dc[0], &ddc[0] );
  mbs_TabBezCurveDer2f ( spdimen, degc01, c01, nknu, knu, 2*spdimen,
            &c[spdimen], &dc[spdimen], &ddc[spdimen] );
  mbs_TabBezCurveDer2f ( spdimen, degd00, d00, nknv, knv, 2*spdimen,
            &d[0], &dd[0], &ddd[0] );
  mbs_TabBezCurveDer2f ( spdimen, degd01, d01, nknv, knv, 2*spdimen,
            &d[spdimen], &dd[spdimen], &ddd[spdimen] );

  if ( p )
    if ( !_mbs_TabBezC1Coons0f ( spdimen, nknu, nknv, c, d, pcorners,
                                 hfuncu, hfuncv, p ) )
      goto failure;
  if ( pu )
    if ( !_mbs_TabBezC1Coons0f ( spdimen, nknu, nknv, dc, d, pcorners,
                                 dhfuncu, hfuncv, pu ) )
      goto failure;
  if ( pv )
    if ( !_mbs_TabBezC1Coons0f ( spdimen, nknu, nknv, c, dd, pcorners,
                                 hfuncu, dhfuncv, pv ) )
      goto failure;
  if ( puu )
    if ( !_mbs_TabBezC1Coons0f ( spdimen, nknu, nknv, ddc, d, pcorners,
                                 ddhfuncu, hfuncv, puu ) )
      goto failure;
  if ( puv )
    if ( !_mbs_TabBezC1Coons0f ( spdimen, nknu, nknv, dc, dd, pcorners,
                                 dhfuncu, dhfuncv, puv ) )
      goto failure;
  if ( pvv )
    if ( !_mbs_TabBezC1Coons0f ( spdimen, nknu, nknv, c, ddd, pcorners,
                                 hfuncu, ddhfuncv, pvv ) )
      goto failure;

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*mbs_TabBezC1Coons0Der2f*/

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
  void  *sp;
  int   ku, kv;
  float *c, *dc, *ddc, *dddc, *d, *dd, *ddd, *dddd, *pcorners;

  sp = pkv_GetScratchMemTop ();
  ku = 6*spdimen*nknu;
  kv = 6*spdimen*nknv;
  c = pkv_GetScratchMemf ( (24*(nknu+nknv)+16)*spdimen );
  if ( !c )
    goto failure;
                  dc = &c[ku];  ddc = &dc[ku];  dddc = &ddc[ku];
  d = &dddc[ku];  dd = &d[kv];  ddd = &dd[kv];  dddd = &ddd[kv];
  pcorners = &dddd[kv];

  mbs_multiBCHornerDerf ( degc00, 1, spdimen, 0, c00, 0.0,
      &pcorners[0], &pcorners[spdimen*2] );
  mbs_multiBCHornerDerf ( degc01, 1, spdimen, 0, c01, 0.0,
      &pcorners[spdimen*1], &pcorners[spdimen*3] );

  mbs_TabBezCurveDer3f ( spdimen, degc00, c00, nknu, knu, 2*spdimen,
            &c[0], &dc[0], &ddc[0], &dddc[0] );
  mbs_TabBezCurveDer3f ( spdimen, degc01, c01, nknu, knu, 2*spdimen,
            &c[spdimen], &dc[spdimen], &ddc[spdimen], &dddc[spdimen] );
  mbs_TabBezCurveDer3f ( spdimen, degd00, d00, nknv, knv, 2*spdimen,
            &d[0], &dd[0], &ddd[0], &dddd[0] );
  mbs_TabBezCurveDer3f ( spdimen, degd01, d01, nknv, knv, 2*spdimen,
            &d[spdimen], &dd[spdimen], &ddd[spdimen], &dddd[spdimen] );

  if ( p )
    if ( !_mbs_TabBezC1Coons0f ( spdimen, nknu, nknv, c, d, pcorners,
                                 hfuncu, hfuncv, p ) )
      goto failure;
  if ( pu )
    if ( !_mbs_TabBezC1Coons0f ( spdimen, nknu, nknv, dc, d, pcorners,
                                 dhfuncu, hfuncv, pu ) )
      goto failure;
  if ( pv )
    if ( !_mbs_TabBezC1Coons0f ( spdimen, nknu, nknv, c, dd, pcorners,
                                 hfuncu, dhfuncv, pv ) )
      goto failure;
  if ( puu )
    if ( !_mbs_TabBezC1Coons0f ( spdimen, nknu, nknv, ddc, d, pcorners,
                                 ddhfuncu, hfuncv, puu ) )
      goto failure;
  if ( puv )
    if ( !_mbs_TabBezC1Coons0f ( spdimen, nknu, nknv, dc, dd, pcorners,
                                 dhfuncu, dhfuncv, puv ) )
      goto failure;
  if ( pvv )
    if ( !_mbs_TabBezC1Coons0f ( spdimen, nknu, nknv, c, ddd, pcorners,
                                 hfuncu, ddhfuncv, pvv ) )
      goto failure;
  if ( puuu )
    if ( !_mbs_TabBezC1Coons0f ( spdimen, nknu, nknv, dddc, d, pcorners,
                                 dddhfuncu, hfuncv, puuu ) )
      goto failure;
  if ( puuv )
    if ( !_mbs_TabBezC1Coons0f ( spdimen, nknu, nknv, ddc, dd, pcorners,
                                 ddhfuncu, dhfuncv, puuv ) )
      goto failure;
  if ( puvv )
    if ( !_mbs_TabBezC1Coons0f ( spdimen, nknu, nknv, dc, ddd, pcorners,
                                 dhfuncu, ddhfuncv, puvv ) )
      goto failure;
  if ( pvvv )
    if ( !_mbs_TabBezC1Coons0f ( spdimen, nknu, nknv, c, dddd, pcorners,
                                 hfuncu, dddhfuncv, pvvv ) )
      goto failure;

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*mbs_TabBezC1Coons0Der3f*/

