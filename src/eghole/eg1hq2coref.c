
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2008, 2013                            */
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

#undef CONST_
#define CONST_

#include "eg1holef.h"
#include "eg1hprivatef.h"
#include "eg1herror.h"

/*#define DEBUG*/

/* ///////////////////////////////////////////////////////////////////////// */
void g1h_DestroyQ2PrivateDataf ( GHoleDomainf *domain )
{
  G1HolePrivateRecf *privateG1;

  privateG1 = domain->privateG1;
  if ( privateG1->Q2AMat )
    { free ( privateG1->Q2AMat );  privateG1->Q2AMat = NULL; }
  if ( privateG1->Q2BMat )
    { free ( privateG1->Q2BMat );  privateG1->Q2BMat = NULL; }
  if ( privateG1->Q2LMat )
    { free ( privateG1->Q2LMat );  privateG1->Q2LMat = NULL; }
  if ( privateG1->Q2EAMat )
    { free ( privateG1->Q2EAMat );  privateG1->Q2EAMat = NULL; }
  if ( privateG1->Q2EBMat )
    { free ( privateG1->Q2EBMat );  privateG1->Q2EBMat = NULL; }
  if ( privateG1->Q2ELMat )
    { free ( privateG1->Q2ELMat );  privateG1->Q2ELMat = NULL; }
} /*g1h_DestroyQ2PrivateDataf*/

/* ///////////////////////////////////////////////////////////////////////// */
boolean _g1h_Q2TabDiPatchJac3f ( int nkn, const float *kn, const float *hfunc,
             const float *dhfunc, const float *ddhfunc, const float *dddhfunc,
             const vector2f *c00, const vector2f *c01,
             const vector2f *c10, const vector2f *c11,
             const vector2f *d00, const vector2f *d01,
             const vector2f *d10, const vector2f *d11,
             float *jac, float *trd )
{
  void     *sp;
  int      i, k;
  vector2f *tabpu, *tabpv, *tabpuu, *tabpuv, *tabpvv,
           *tabpuuu, *tabpuuv, *tabpuvv, *tabpvvv;

  sp = pkv_GetScratchMemTop ();
  k = nkn*nkn;
  tabpu = (vector2f*)pkv_GetScratchMem ( 9*k*sizeof(vector2f) );
  if ( !tabpu )
    goto failure;
  tabpv = &tabpu[k];      tabpuu = &tabpv[k];    tabpuv = &tabpuu[k];
  tabpvv = &tabpuv[k];    tabpuuu = &tabpvv[k];  tabpuuv = &tabpuuu[k];
  tabpuvv = &tabpuuv[k];  tabpvvv = &tabpuvv[k];

  if ( !mbs_TabBezC1CoonsDer3f ( 2, nkn, kn, hfunc, dhfunc, ddhfunc, dddhfunc,
           nkn, kn, hfunc, dhfunc, ddhfunc, dddhfunc,
           G1_CROSS00DEG, (float*)c00, G1_CROSS01DEG, (float*)c01,
           G1_CROSS10DEG, (float*)c10, G1_CROSS11DEG, (float*)c11,
           G1_CROSS00DEG, (float*)d00, G1_CROSS01DEG, (float*)d01,
           G1_CROSS10DEG, (float*)d10, G1_CROSS11DEG, (float*)d11,
           NULL, (float*)tabpu, (float*)tabpv, (float*)tabpuu,
           (float*)tabpuv, (float*)tabpvv, (float*)tabpuuu, (float*)tabpuuv,
           (float*)tabpuvv, (float*)tabpvvv ) )
    goto failure;

  for ( i = 0; i < k; i++ )
    if ( !_g2h_DiJacobian3f ( &tabpu[i], &tabpv[i], &tabpuu[i], &tabpuv[i], &tabpvv[i],
                              &tabpuuu[i], &tabpuuv[i], &tabpuvv[i], &tabpvvv[i],
                              &jac[i], &trd[18*i] ) )
      goto failure;

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*_g1h_Q2TabDiPatchJac3f*/

boolean _g1h_Q2TabLaplacianGradf ( int nkn, const float *tkn,
        const float *hfunc, const float *dhfunc, const float *ddhfunc,
        const float *dddhfunc,
        const float *fc00, const float *fc01,
        const float *fc10, const float *fc11,
        const float *fd00, const float *fd01,
        const float *fd10, const float *fd11,
        const float *trd,
        vector2f *lapgrad )
{
  void  *sp;
  int   i, k;
  float *tabfu, *tabfv, *tabfuu, *tabfuv, *tabfvv,
        *tabfuuu, *tabfuuv, *tabfuvv, *tabfvvv;

  sp = pkv_GetScratchMemTop ();
  k = nkn*nkn;

  tabfu = pkv_GetScratchMemf ( 9*k );
  if ( !tabfu )
    goto failure;
  tabfv = &tabfu[k];      tabfuu = &tabfv[k];    tabfuv = &tabfuu[k];
  tabfvv = &tabfuv[k];    tabfuuu = &tabfvv[k];  tabfuuv = &tabfuuu[k];
  tabfuvv = &tabfuuv[k];  tabfvvv = &tabfuvv[k];

  if ( !mbs_TabBezC1CoonsDer3f ( 1, nkn, tkn, hfunc, dhfunc, ddhfunc, dddhfunc,
                     nkn, tkn, hfunc, dhfunc, ddhfunc, dddhfunc,
                     G1_CROSS00DEG, fc00, G1_CROSS01DEG, fc01,
                     G1_CROSS10DEG, fc10, G1_CROSS11DEG, fc11,
                     G1_CROSS00DEG, fd00, G1_CROSS01DEG, fd01,
                     G1_CROSS10DEG, fd10, G1_CROSS11DEG, fd11,
                     NULL, tabfu, tabfv, tabfuu, tabfuv, tabfvv,
                     tabfuuu, tabfuuv, tabfuvv, tabfvvv ) )
    goto failure;

  for ( i = k = 0;  i < nkn*nkn;  i++, k += 18 ) {
    lapgrad[i].x = trd[k]*tabfu[i] + trd[k+1]*tabfv[i] +
                   trd[k+2]*tabfuu[i] + trd[k+3]*tabfuv[i] +
                   trd[k+4]*tabfvv[i] +
                   trd[k+5]*tabfuuu[i] + trd[k+6]*tabfuuv[i] +
                   trd[k+7]*tabfuvv[i] + trd[k+8]*tabfvvv[i];
    lapgrad[i].y = trd[k+9]*tabfu[i] + trd[k+10]*tabfv[i] +
                   trd[k+11]*tabfuu[i] + trd[k+12]*tabfuv[i] +
                   trd[k+13]*tabfvv[i] +
                   trd[k+14]*tabfuuu[i] + trd[k+15]*tabfuuv[i] +
                   trd[k+16]*tabfuvv[i] + trd[k+17]*tabfvvv[i];
  }

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*_g1h_Q2TabLaplacianGradf*/

boolean _g1h_Q2TabLaplacianGrad0f ( int nkn, const float *tkn,
        const float *hfunc, const float *dhfunc, const float *ddhfunc,
        const float *dddhfunc,
        const float *fc00, const float *fc01,
        const float *fd00, const float *fd01,
        const float *trd,
        vector2f *lapgrad )
{
  void  *sp;
  int   i, k;
  float *tabfu, *tabfv, *tabfuu, *tabfuv, *tabfvv,
        *tabfuuu, *tabfuuv, *tabfuvv, *tabfvvv;

  sp = pkv_GetScratchMemTop ();
  k = nkn*nkn;

  tabfu = pkv_GetScratchMemf ( 9*k );
  if ( !tabfu )
    goto failure;
  tabfv = &tabfu[k];      tabfuu = &tabfv[k];    tabfuv = &tabfuu[k];
  tabfvv = &tabfuv[k];    tabfuuu = &tabfvv[k];  tabfuuv = &tabfuuu[k];
  tabfuvv = &tabfuuv[k];  tabfvvv = &tabfuvv[k];

  if ( !mbs_TabBezC1Coons0Der3f ( 1, nkn, tkn, hfunc, dhfunc, ddhfunc, dddhfunc,
                     nkn, tkn, hfunc, dhfunc, ddhfunc, dddhfunc,
                     G1_CROSS00DEG, fc00, G1_CROSS01DEG, fc01,
                     G1_CROSS00DEG, fd00, G1_CROSS01DEG, fd01,
                     NULL, tabfu, tabfv, tabfuu, tabfuv, tabfvv,
                     tabfuuu, tabfuuv, tabfuvv, tabfvvv ) )
    goto failure;

  for ( i = k = 0;  i < nkn*nkn;  i++, k += 18 ) {
    lapgrad[i].x = trd[k]*tabfu[i] + trd[k+1]*tabfv[i] +
                   trd[k+2]*tabfuu[i] + trd[k+3]*tabfuv[i] +
                   trd[k+4]*tabfvv[i] +
                   trd[k+5]*tabfuuu[i] + trd[k+6]*tabfuuv[i] +
                   trd[k+7]*tabfuvv[i] + trd[k+8]*tabfvvv[i];
    lapgrad[i].y = trd[k+9]*tabfu[i] + trd[k+10]*tabfv[i] +
                   trd[k+11]*tabfuu[i] + trd[k+12]*tabfuv[i] +
                   trd[k+13]*tabfvv[i] +
                   trd[k+14]*tabfuuu[i] + trd[k+15]*tabfuuv[i] +
                   trd[k+16]*tabfuvv[i] + trd[k+17]*tabfvvv[i];
  }

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*_g1h_Q2TabLaplacianGrad0f*/

boolean _g1h_TabCurveJacobianf ( int deg, const point2f *cp,
                                 int nkn, const float *kn, float *jac )
{
  int      i;
  point2f  p;
  vector2f dp;

  for ( i = 0; i < nkn; i++ ) {
    if ( !mbs_BCHornerDerC2f ( deg, cp, kn[i], &p, &dp ) )
      return false;
    jac[i] = (float)sqrt ( dp.x*dp.x+dp.y*dp.y );
  }
  return true;
} /*_g1h_TabCurveJacobianf*/

boolean _g1h_LapCoefff ( const vector2f *du, const vector2f *dv,
                         const vector2f *duu, const vector2f *duv,
                         const vector2f *dvv, float *trd )
{
  vector2f gx, gy, gxx, gxy, gyy;
  float    A21[6], A22[9];

  if ( !pkn_f2iDerivatives2f ( du->x, du->y, dv->x, dv->y,
          duu->x, duu->y, duv->x, duv->y, dvv->x, dvv->y,
          (float*)&gx, (float*)&gy, (float*)&gxx, (float*)&gxy, (float*)&gyy ) )
    return false;

  pkn_Setup2DerA21Matrixf ( gxx.x, gxx.y, gxy.x, gxy.y, gyy.x, gyy.y, A21 );
  pkn_Setup2DerA22Matrixf ( gx.x, gx.y, gy.x, gy.y, A22 );

  trd[0]  = A21[0]+A21[4];   trd[1]  = A21[1]+A21[5];
  trd[2]  = A22[0]+A22[6];   trd[3]  = A22[1]+A22[7];   trd[4] = A22[2]+A22[8];
  return true;
} /*_g1h_LapCoefff*/

boolean _g1h_TabCurveLapCoeff0f ( const point2f *c00, const vector2f *c01,
                                  const point2f *c10, const vector2f *c11,
                                  const point2f *d00, const vector2f *d01,
                                  const point2f *d10, const vector2f *d11,
                                  int nkn, const float *tkn, const float *hfunc,
                                  const float *dhfunc, const float *ddhfunc,
                                  const float *atkn, const float *ahfunc,
                                  const float *adhfunc, const float *addhfunc,
                                  float *trdc00, float *trdc10,
                                  float *trdd00, float *trdd10 )
{
  void     *sp;
  vector2f *tabpu, *tabpv, *tabpuu, *tabpuv, *tabpvv;
  int      i, ii;

  sp = pkv_GetScratchMemTop ();
  tabpu = (vector2f*)pkv_GetScratchMem ( 10*nkn*sizeof(vector2f) );
  if ( !tabpu )
    goto failure;
  tabpv = &tabpu[2*nkn];    tabpuu = &tabpv[2*nkn];
  tabpuv = &tabpuu[2*nkn];  tabpvv = &tabpuv[2*nkn];

  if ( !mbs_TabBezC1CoonsDer2f ( 2, nkn, tkn, hfunc, dhfunc, ddhfunc,
                         2, atkn, ahfunc, adhfunc, addhfunc,
                         G1_CROSS00DEG, (float*)c00, G1_CROSS01DEG, (float*)c01,
                         G1_CROSS10DEG, (float*)c10, G1_CROSS11DEG, (float*)c11,
                         G1_CROSS00DEG, (float*)d00, G1_CROSS01DEG, (float*)d01,
                         G1_CROSS10DEG, (float*)d10, G1_CROSS11DEG, (float*)d11,
                         NULL, (float*)tabpu, (float*)tabpv,
                         (float*)tabpuu, (float*)tabpuv, (float*)tabpvv ) )
    goto failure;
  for ( i = ii = 0;  i < nkn;  i++ ) {
    if ( !_g1h_LapCoefff ( &tabpu[ii], &tabpv[ii],
                           &tabpuu[ii], &tabpuv[ii], &tabpvv[ii], &trdc00[5*i] ) )
      goto failure;
    ii++;
    if ( !_g1h_LapCoefff ( &tabpu[ii], &tabpv[ii],
                           &tabpuu[ii], &tabpuv[ii], &tabpvv[ii], &trdc10[5*i] ) )
      goto failure;
    ii++;
  }

  if ( !mbs_TabBezC1CoonsDer2f ( 2, 2, atkn, ahfunc, adhfunc, addhfunc,
                         nkn, tkn, hfunc, dhfunc, ddhfunc,
                         G1_CROSS00DEG, (float*)c00, G1_CROSS01DEG, (float*)c01,
                         G1_CROSS10DEG, (float*)c10, G1_CROSS11DEG, (float*)c11,
                         G1_CROSS00DEG, (float*)d00, G1_CROSS01DEG, (float*)d01,
                         G1_CROSS10DEG, (float*)d10, G1_CROSS11DEG, (float*)d11,
                         NULL, (float*)tabpu, (float*)tabpv,
                         (float*)tabpuu, (float*)tabpuv, (float*)tabpvv ) )
    goto failure;
  for ( i = 0, ii = nkn;  i < nkn;  i++, ii++ ) {
    if ( !_g1h_LapCoefff ( &tabpu[i], &tabpv[i],
                           &tabpuu[i], &tabpuv[i], &tabpvv[i], &trdd00[5*i] ) )
      goto failure;
    if ( !_g1h_LapCoefff ( &tabpu[ii], &tabpv[ii],
                           &tabpuu[ii], &tabpuv[ii], &tabpvv[ii], &trdd10[5*i] ) )
      goto failure;
  }

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*_g1h_TabCurveLapCoeff0f*/

boolean _g1h_TabCurveLapCoeff1f ( const point2f *sicp, int nkn,
                                  const float *tkn, float *trd )
{
  int i;
  point2f d;
  vector2f du, dv, duu, duv, dvv;

  for ( i = 0; i < nkn; i++ ) {
    if ( !mbs_BCHornerDer2Pf ( 3, 3, 2, (float*)sicp, 0.0, tkn[i],
                         &d.x, &du.x, &dv.x, &duu.x, &duv.x, &dvv.x ) )
      return false;
    if ( !_g1h_LapCoefff ( &du, &dv, &duu, &duv, &dvv, &trd[5*i] ) )
      return false;
  }
  return true;
} /*_g1h_TabCurveLapCoeff1f*/

boolean _g1h_Q2TabLaplacianJump0f ( int nkn, const float *tkn,
              const float *hfunc, const float *dhfunc, const float *ddhfunc,
              const float *atkn, const float *ahfunc, const float *adhfunc,
              const float *addhfunc,
              const float *ec00, const float *ec01,
              const float *ed00, const float *ed01, const float *etrdd00,
              const float *fc00, const float *fc01,
              const float *fd00, const float *fd01,
              const float *ftrdc00, const float *ftrdc10, const float *ftrdd10,
              float *lapjc00, float *lapjc10, float *lapjd10 )
{
  void  *sp;
  int   i, j, k;
  float *tabeu, *tabev, *tabeuu, *tabeuv, *tabevv,
        *tabfu, *tabfv, *tabfuu, *tabfuv, *tabfvv;
  float lape, lapf;

  sp = pkv_GetScratchMemTop ();
  tabeu = pkv_GetScratchMemf ( 20*nkn );
  if ( !tabeu )
    goto failure;
  tabev = &tabeu[2*nkn];    tabeuu = &tabev[2*nkn];   tabeuv = &tabeuu[2*nkn];
  tabevv = &tabeuv[2*nkn];  tabfu = &tabevv[2*nkn];   tabfv = &tabfu[2*nkn];
  tabfuu = &tabfv[2*nkn];   tabfuv = &tabfuu[2*nkn];  tabfvv = &tabfuv[2*nkn];
  if ( !mbs_TabBezC1Coons0Der2f ( 1, nkn, tkn, hfunc, dhfunc, ddhfunc,
              2, atkn, ahfunc, adhfunc, addhfunc,
              G1_CROSS00DEG, fc00, G1_CROSS01DEG, fc01,
              G1_CROSS00DEG, fd00, G1_CROSS01DEG, fd01,
              NULL, tabfu, tabfv, tabfuu, tabfuv, tabfvv ) )
    goto failure;
  if ( !mbs_TabBezC1Coons0Der2f ( 1, 1, atkn, ahfunc, adhfunc, addhfunc,
              nkn, tkn, hfunc, dhfunc, ddhfunc,
              G1_CROSS00DEG, ec00, G1_CROSS01DEG, ec01,
              G1_CROSS00DEG, ed00, G1_CROSS01DEG, ed01,
              NULL, tabeu, tabev, tabeuu, tabeuv, tabevv ) )
    goto failure;
  for ( i = j = k = 0;  i < nkn;  i++, j += 5, k += 2 ) {
    lape = etrdd00[j]*tabeu[i] + etrdd00[j+1]*tabev[i] + etrdd00[j+2]*tabeuu[i] +
           etrdd00[j+3]*tabeuv[i] + etrdd00[j+4]*tabevv[i];
    lapf = ftrdc00[j]*tabfu[k] + ftrdc00[j+1]*tabfv[k] + ftrdc00[j+2]*tabfuu[k] +
           ftrdc00[j+3]*tabfuv[k] + ftrdc00[j+4]*tabfvv[k];
    lapjc00[i] = lapf-lape;
  }
  for ( i = j = 0, k = 1;  i < nkn;  i++, j += 5, k += 2 ) {
    lapf = ftrdc10[j]*tabfu[k] + ftrdc10[j+1]*tabfv[k] + ftrdc10[j+2]*tabfuu[k] +
           ftrdc10[j+3]*tabfuv[k] + ftrdc10[j+4]*tabfvv[k];
    lapjc10[i] = -lapf;
  }

  if ( !mbs_TabBezC1Coons0Der2f ( 1,
              1, &atkn[1], &ahfunc[4], &adhfunc[4], &addhfunc[4],
              nkn, tkn, hfunc, dhfunc, ddhfunc,
              G1_CROSS00DEG, fc00, G1_CROSS01DEG, fc01,
              G1_CROSS00DEG, fd00, G1_CROSS01DEG, fd01,
              NULL, tabfu, tabfv, tabfuu, tabfuv, tabfvv ) )
    goto failure;
  for ( i = j = 0;  i < nkn;  i++, j += 5 ) {
    lapf = ftrdd10[j]*tabfu[i] + ftrdd10[j+1]*tabfv[i] + ftrdd10[j+2]*tabfuu[i] +
           ftrdd10[j+3]*tabfuv[i] + ftrdd10[j+4]*tabfvv[i];
    lapjd10[i] = -lapf;
  }

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*_g1h_Q2TabLaplacianJump0f*/

boolean _g1h_Q2TabLaplacianJumpf ( int nkn, const float *tkn,
              const float *hfunc, const float *dhfunc, const float *ddhfunc,
              const float *atkn, const float *ahfunc, const float *adhfunc,
              const float *addhfunc,
              const float *ec00, const float *ec01,
              const float *ec10, const float *ec11,
              const float *ed00, const float *ed01,
              const float *ed10, const float *ed11, const float *etrdd00,
              const float *fc00, const float *fc01,
              const float *fc10, const float *fc11,
              const float *fd00, const float *fd01,
              const float *fd10, const float *fd11,
              const float *ftrdc00, const float *ftrdc10, const float *ftrdd10,
              const float *eicp1, const float *etrdc10,
              const float *eicp2, const float *etrdd10,
              float *lapjc00, float *lapjc10, float *lapjd10 )
{
  void  *sp;
  int   i, j, k;
  float *tabeu, *tabev, *tabeuu, *tabeuv, *tabevv,
        *tabfu, *tabfv, *tabfuu, *tabfuv, *tabfvv;
  float lape, lapf;

  sp = pkv_GetScratchMemTop ();
  tabeu = pkv_GetScratchMemf ( 15*nkn );
  if ( !tabeu )
    goto failure;
  tabev = &tabeu[nkn];      tabeuu = &tabev[nkn];     tabeuv = &tabeuu[nkn];
  tabevv = &tabeuv[nkn];    tabfu = &tabevv[nkn];     tabfv = &tabfu[2*nkn];
  tabfuu = &tabfv[2*nkn];   tabfuv = &tabfuu[2*nkn];  tabfvv = &tabfuv[2*nkn];
  if ( !mbs_TabBezC1CoonsDer2f ( 1, nkn, tkn, hfunc, dhfunc, ddhfunc,
              2, atkn, ahfunc, adhfunc, addhfunc,
              G1_CROSS00DEG, fc00, G1_CROSS01DEG, fc01,
              G1_CROSS10DEG, fc10, G1_CROSS11DEG, fc11,
              G1_CROSS00DEG, fd00, G1_CROSS01DEG, fd01,
              G1_CROSS10DEG, fd10, G1_CROSS11DEG, fd11,
              NULL, tabfu, tabfv, tabfuu, tabfuv, tabfvv ) )
    goto failure;
  if ( !mbs_TabBezC1CoonsDer2f ( 1, 1, atkn, ahfunc, adhfunc, addhfunc,
              nkn, tkn, hfunc, dhfunc, ddhfunc,
              G1_CROSS00DEG, ec00, G1_CROSS01DEG, ec01,
              G1_CROSS10DEG, ec10, G1_CROSS11DEG, ec11,
              G1_CROSS00DEG, ed00, G1_CROSS01DEG, ed01,
              G1_CROSS10DEG, ed10, G1_CROSS11DEG, ed11,
              NULL, tabeu, tabev, tabeuu, tabeuv, tabevv ) )
    goto failure;
  for ( i = j = k = 0;  i < nkn;  i++, j += 5, k += 2 ) {
    lape = etrdd00[j]*tabeu[i] + etrdd00[j+1]*tabev[i] + etrdd00[j+2]*tabeuu[i] +
           etrdd00[j+3]*tabeuv[i] + etrdd00[j+4]*tabevv[i];
    lapf = ftrdc00[j]*tabfu[k] + ftrdc00[j+1]*tabfv[k] + ftrdc00[j+2]*tabfuu[k] +
           ftrdc00[j+3]*tabfuv[k] + ftrdc00[j+4]*tabfvv[k];
    lapjc00[i] = lapf-lape;
  }
  for ( i = j = 0, k = 1;  i < nkn;  i++, j += 5, k += 2 ) {
    if ( !mbs_BCHornerDer2Pf ( 3, 3, 1, eicp1, 0.0, tkn[i],
                         &tabeu[1], tabeu, tabev, tabeuu, tabeuv, tabevv ) )
      goto failure;
    lape = etrdc10[j]*tabeu[0] + etrdc10[j+1]*tabev[0] + etrdc10[j+2]*tabeuu[0] +
           etrdc10[j+3]*tabeuv[0] + etrdc10[j+4]*tabevv[0];
    lapf = ftrdc10[j]*tabfu[k] + ftrdc10[j+1]*tabfv[k] + ftrdc10[j+2]*tabfuu[k] +
           ftrdc10[j+3]*tabfuv[k] + ftrdc10[j+4]*tabfvv[k];
    lapjc10[i] = lape-lapf;
  }

  if ( !mbs_TabBezC1CoonsDer2f ( 1,
              1, &atkn[1], &ahfunc[4], &adhfunc[4], &addhfunc[4],
              nkn, tkn, hfunc, dhfunc, ddhfunc,
              G1_CROSS00DEG, fc00, G1_CROSS01DEG, fc01,
              G1_CROSS10DEG, fc10, G1_CROSS11DEG, fc11,
              G1_CROSS00DEG, fd00, G1_CROSS01DEG, fd01,
              G1_CROSS10DEG, fd10, G1_CROSS11DEG, fd11,
              NULL, tabfu, tabfv, tabfuu, tabfuv, tabfvv ) )
    goto failure;
  for ( i = j = 0;  i < nkn;  i++, j += 5 ) {
    if ( !mbs_BCHornerDer2Pf ( 3, 3, 1, eicp2, 0.0, tkn[i],
                         &tabeu[1], tabeu, tabev, tabeuu, tabeuv, tabevv ) )
      goto failure;
    lape = etrdd10[j]*tabeu[0] + etrdd10[j+1]*tabev[0] + etrdd10[j+2]*tabeuu[0] +
           etrdd10[j+3]*tabeuv[0] + etrdd10[j+4]*tabevv[0];
    lapf = ftrdd10[j]*tabfu[i] + ftrdd10[j+1]*tabfv[i] + ftrdd10[j+2]*tabfuu[i] +
           ftrdd10[j+3]*tabfuv[i] + ftrdd10[j+4]*tabfvv[i];
    lapjd10[i] = lape-lapf;
  }

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*_g1h_Q2TabLaplacianJumpf*/

float _g1h_Q2Integralf ( int hole_k, int nquad, float *jac,
                         unsigned short supp1, float *lapj1,
                         unsigned short supp2, float *lapj2 )
{
  int            i, k;
  double         s;
  unsigned short supp;
  float          *jc, *f1, *f2;

  supp = (unsigned short)(supp1 & supp2);
  s = 0.0;
  for ( k = 0; k < hole_k; k++ )
    if ( supp & (0x0001 << k) ) {
      jc = &jac[k*3*nquad];
      f1 = &lapj1[k*3*nquad];
      f2 = &lapj2[k*3*nquad];
      for ( i = 0; i < 3*nquad; i++ )
        s += f1[i]*f2[i]*jc[i];
    }
  s /= (double)nquad;
  return (float)s;
} /*_g1h_Q2Integralf*/

boolean g1h_Q2ComputeFormMatrixf ( GHoleDomainf *domain )
{
  void     *sp;
  GHolePrivateRecf    *privateG;
  G1HolePrivateRecf   *privateG1;
  int      hole_k, nfunc_a, nfunc_b;
  int      i, j, l, n, fn;
  float    *tkn, *hfunc, *dhfunc, *ddhfunc, *dddhfunc;
  float    *atkn, *ahfunc, *adhfunc, *addhfunc;
  vector2f *c00, *c01, *c10, *c11, *d00, *d01, *d10, *d11;
  float    *ec00, *ec01, *ec10, *ec11, *ed00, *ed01, *ed10, *ed11;
  float    *fc00, *fc01, *fc10, *fc11, *fd00, *fd01, *fd10, *fd11;
  float    *trd, *lapj, *jac, *amat, *bmat;
  vector2f *lgr;
  unsigned short *support_b, supp;
  int      option, ndata, *idata;
  float    *fdata, *eicp1, *eicp2, C;
  point2f  *sicp;
#ifdef DEBUG
FILE *f;
#endif

  sp = pkv_GetScratchMemTop ();
  if ( !domain->basisG1 ) {
    if ( !g1h_ComputeBasisf ( domain ) )
      goto failure;
  }
  hole_k    = domain->hole_k;
  privateG  = domain->privateG;
  privateG1 = domain->privateG1;
  nfunc_a   = privateG1->nfunc_a;
  nfunc_b   = privateG1->nfunc_b;
  support_b = privateG->support_b;
  if ( !(amat = privateG1->Q2AMat) )
    amat = privateG1->Q2AMat = malloc ( (nfunc_a*(nfunc_a+1)/2)*sizeof(float) );
  if ( !(bmat = privateG1->Q2BMat) )
    bmat = privateG1->Q2BMat = malloc ( nfunc_a*nfunc_b*sizeof(float) );
  if ( !amat || !bmat )
    goto failure;
  if ( !privateG->diam )
    privateG->diam = gh_DomainDiamf ( domain );

  n = hole_k*G1_NQUADSQ;

  tkn = pkv_GetScratchMemf ( 17*G1_NQUAD );
  jac = pkv_GetScratchMemf ( n );
  trd = pkv_GetScratchMemf ( 18*n );
  lgr = (vector2f*)pkv_GetScratchMem ( (nfunc_a+1)*n*sizeof(vector2f) );
  if ( !tkn || !jac || !trd || !lgr )
    goto failure;
  hfunc = &tkn[G1_NQUAD];         dhfunc = &hfunc[4*G1_NQUAD];
  ddhfunc = &dhfunc[4*G1_NQUAD];  dddhfunc = &ddhfunc[4*G1_NQUAD];
  _gh_PrepareTabKnotsf ( G1_NQUAD, privateG1->opt_quad, tkn );
  if ( !mbs_TabCubicHFuncDer3f ( 0.0, 1.0, G1_NQUAD, tkn,
                                 hfunc, dhfunc, ddhfunc, dddhfunc ) )
    goto failure;

        /* integrate the Laplacian gradients */
  for ( i = 0; i < hole_k; i++ ) {
    _g1h_GetDiPatchCurvesf ( domain, i,
                             &c00, &c01, &c10, &c11, &d00, &d01, &d10, &d11 );
    _g1h_Q2TabDiPatchJac3f ( G1_NQUAD, tkn, hfunc, dhfunc, ddhfunc, dddhfunc,
                             c00, c01, c10, c11, d00, d01, d10, d11,
                             &jac[i*G1_NQUADSQ], &trd[i*G1_NQUADSQ*18] );
  }

  for ( fn = 0; fn < nfunc_a; fn++ )
    for ( i = 0; i < hole_k; i++ ) {
      _g1h_GetBFAPatchCurvesf ( domain, fn, i, &fc00, &fc01, &fd00, &fd01 );
      _g1h_Q2TabLaplacianGrad0f ( G1_NQUAD, tkn, hfunc, dhfunc, ddhfunc, dddhfunc,
                                  fc00, fc01, fd00, fd01,
                                  &trd[i*G1_NQUADSQ*18],
                                  &lgr[(fn*hole_k+i)*G1_NQUADSQ] );
    }
  for ( i = l = 0;  i < nfunc_a;  i++ ) {
    for ( j = 0;  j <= i;  j++, l++ )
      amat[l] = _g2h_Integralf ( hole_k, G1_NQUADSQ, jac,
                                 0xFFFF, &lgr[i*n], 0xFFFF, &lgr[j*n] );
  }

  for ( fn = 0; fn < nfunc_b; fn++ ) {
    for ( i = 0; i < hole_k; i++ )
      if ( support_b[fn] & (0x0001 << i) ) {
        _g1h_GetBFBPatchCurvesf ( domain, fn, i, &fc00, &fc01, &fc10, &fc11,
                                                 &fd00, &fd01, &fd10, &fd11 );
        _g1h_Q2TabLaplacianGradf ( G1_NQUAD, tkn, hfunc, dhfunc, ddhfunc, dddhfunc,
                                   fc00, fc01, fc10, fc11, fd00, fd01, fd10, fd11,
                                   &trd[i*G1_NQUADSQ*18],
                                   &lgr[(nfunc_a*hole_k+i)*G1_NQUADSQ] );
      }
    for ( i = 0; i < nfunc_a; i++ )
      bmat[i*nfunc_b+fn] = _g2h_Integralf ( hole_k, G1_NQUADSQ, jac,
                             0xFFFF, &lgr[i*n], support_b[fn], &lgr[nfunc_a*n] );
  }

        /* integrate the Laplacian jumps */
  option = privateG1->GetOption ( domain, G1HQUERY_Q2_FORM_CONSTANT, 0,
                                  &ndata, &idata, &fdata );
  switch ( option ) {
case G1H_DEFAULT:
    privateG1->C1b = C = 1.0;
                       /* no idea whether this default value makes any sense */
    break;
case G1H_Q2_USE_SUPPLIED_CONSTANT:
    privateG1->C1b = C = *fdata;
    break;
default:
    goto failure;
  }
  C *= (float)(5.0/privateG->diam);

        /* the arrays allocated for the previous stage are reused */
        /* their length is more than sufficient */
  lapj = &lgr[0].x;
        /* two additional arrays are however needed */
  atkn = pkv_GetScratchMemf ( 26 );
  sicp = (point2f*)pkv_GetScratchMem ( 16*sizeof(point2f) );
  if ( !atkn || !sicp )
    goto failure;
  ahfunc = &atkn[2];  adhfunc = &ahfunc[8];  addhfunc = &adhfunc[8];
  eicp1 = (float*)sicp;  eicp2 = &eicp1[16];
  atkn[0] = 0.0;  atkn[1] = 1.0;
  if ( !mbs_TabCubicHFuncDer2f ( 0.0, 1.0, 2, atkn, ahfunc, adhfunc, addhfunc ) )
    goto failure;

  for ( i = 0; i < hole_k; i++ ) {
    j = (i+1) % hole_k;
          /* compute the curve Jacobians */
    _g1h_GetDiPatchCurvesf ( domain, i,
                             &c00, &c01, &c10, &c11, &d00, &d01, &d10, &d11 );
    if ( !_g1h_TabCurveJacobianf ( G1_CROSS00DEG, c00, G1_NQUAD, tkn,
                                   &jac[3*i*G1_NQUAD] ) )
      goto failure;
    if ( !_g1h_TabCurveJacobianf ( G1_CROSS10DEG, c10, G1_NQUAD, tkn,
                                   &jac[(3*i+1)*G1_NQUAD] ) )
      goto failure;
    if ( !_g1h_TabCurveJacobianf ( G1_CROSS10DEG, d10, G1_NQUAD, tkn,
                                   &jac[(3*i+2)*G1_NQUAD] ) )
      goto failure;
          /* compute the coefficients for computing the Laplacians */
    if ( !_g1h_TabCurveLapCoeff0f ( c00, c01, c10, c11, d00, d01, d10, d11,
                            G1_NQUAD, tkn, hfunc, dhfunc, ddhfunc,
                            atkn, ahfunc, adhfunc, addhfunc,
                            &trd[30*i*G1_NQUAD], &trd[(30*i+10)*G1_NQUAD],
                            &trd[(30*j+5)*G1_NQUAD], &trd[(30*i+20)*G1_NQUAD] ) )
      goto failure;
    if ( !_gh_FindDomSurrndPatchf ( domain, j, 1, sicp ) )
      goto failure;
    if ( !_g1h_TabCurveLapCoeff1f ( sicp, G1_NQUAD, tkn, &trd[(30*i+15)*G1_NQUAD] ) )
      goto failure;
    if ( !_gh_FindDomSurrndPatchf ( domain, i, 2, sicp ) )
      goto failure;
    if ( !_g1h_TabCurveLapCoeff1f ( sicp, G1_NQUAD, tkn, &trd[(30*i+25)*G1_NQUAD] ) )
      goto failure;
  }

  for ( fn = 0; fn < nfunc_a; fn++ ) {
    _g1h_GetBFAPatchCurvesf ( domain, fn, hole_k-1,
                              &ec00, &ec01, &ed00, &ed01 );
    for ( i = 0; i < hole_k; i++ ) {
      _g1h_GetBFAPatchCurvesf ( domain, fn, i,
                                &fc00, &fc01, &fd00, &fd01 );
      if ( !_g1h_Q2TabLaplacianJump0f ( G1_NQUAD, tkn, hfunc, dhfunc, ddhfunc,
                          atkn, ahfunc, adhfunc, addhfunc,
                          ec00, ec01, ed00, ed01, &trd[(30*i+5)*G1_NQUAD],
                          fc00, fc01, fd00, fd01, &trd[30*i*G1_NQUAD],
                          &trd[(30*i+10)*G1_NQUAD], &trd[(30*i+20)*G1_NQUAD],
                          &lapj[3*(fn*hole_k+i)*G1_NQUAD],
                          &lapj[(3*(fn*hole_k+i)+1)*G1_NQUAD],
                          &lapj[(3*(fn*hole_k+i)+2)*G1_NQUAD] ) )
        goto failure;
      ec00 = fc00;  ec01 = fc01;  ed00 = fd00;  ed01 = fd01;
    }
  }
#ifdef DEBUG
f = fopen ( "lapj.txt", "w+" );
for ( fn = 0; fn < nfunc_a; fn++ ) {
  fprintf ( f, "fn:%d\n", fn );
  for ( i = 0; i < hole_k; i++ ) {
    for ( j = 0; j < 3*G1_NQUAD; j++ )
      fprintf ( f, "%8.5e,", lapj[(fn*hole_k+i)*3*G1_NQUAD+j] );
    fprintf ( f, "\n" );
  }
  fprintf ( f, "\n" );
}
#endif
  for ( i = l = 0;  i < nfunc_a;  i++ ) {
    for ( j = 0;  j <= i;  j++, l++ )
      amat[l] += C*_g1h_Q2Integralf ( hole_k, G1_NQUAD, jac,
                                      0xFFFF, &lapj[i*hole_k*3*G1_NQUAD],
                                      0xFFFF, &lapj[j*hole_k*3*G1_NQUAD] );
  }
  for ( fn = 0; fn < nfunc_b; fn++ ) {
    supp = _g1h_ExtendSupport ( hole_k, support_b[fn] );
    for ( i = 0; i < hole_k; i++ ) {
      j = (i+1) % hole_k;
      if ( supp & (0x0001 << i) ) {
        _g1h_GetBFBPatchCurvesf ( domain, fn, (i+hole_k-1) % hole_k,
                  &ec00, &ec01, &ec10, &ec11, &ed00, &ed01, &ed10, &ed11 );
        _g1h_GetBFBPatchCurvesf ( domain, fn, i,
                  &fc00, &fc01, &fc10, &fc11, &fd00, &fd01, &fd10, &fd11 );
        gh_GetDomSurrndBFuncf ( domain, fn, j, 1, eicp1 );
        gh_GetDomSurrndBFuncf ( domain, fn, i, 2, eicp2 );
        if ( !_g1h_Q2TabLaplacianJumpf ( G1_NQUAD, tkn, hfunc, dhfunc, ddhfunc,
                           atkn, ahfunc, adhfunc, addhfunc,
                           ec00, ec01, ec10, ec11, ed00, ed01, ed10, ed11,
                           &trd[(30*i+5)*G1_NQUAD],
                           fc00, fc01, fc10, fc11, fd00, fd01, fd10, fd11,
                           &trd[30*i*G1_NQUAD], &trd[(30*i+10)*G1_NQUAD],
                           &trd[(30*i+20)*G1_NQUAD],
                           eicp1, &trd[(30*i+15)*G1_NQUAD],
                           eicp2, &trd[(30*i+25)*G1_NQUAD],
                           &lapj[3*(nfunc_a*hole_k+i)*G1_NQUAD],
                           &lapj[(3*(nfunc_a*hole_k+i)+1)*G1_NQUAD],
                           &lapj[(3*(nfunc_a*hole_k+i)+2)*G1_NQUAD] ) )
          goto failure;
      }
    }
#ifdef DEBUG
fprintf ( f, "fn:%d\n", fn );
for ( i = 0; i < hole_k; i++ ) {
  for ( j = 0; j < 3*G1_NQUAD; j++ )
    fprintf ( f, "%8.5e,", lapj[(nfunc_a*hole_k+i)*3*G1_NQUAD+j] );
  fprintf ( f, "\n" );
}
fprintf ( f, "\n" );
#endif
    for ( i = 0; i < nfunc_a; i++ )
      bmat[i*nfunc_b+fn] += C*_g1h_Q2Integralf ( hole_k, G1_NQUAD, jac,
                                   0xFFFF, &lapj[i*hole_k*3*G1_NQUAD],
                                   supp, &lapj[nfunc_a*hole_k*3*G1_NQUAD] );
  }
#ifdef DEBUG
fclose ( f );
#endif

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*g1h_Q2ComputeFormMatrixf*/

boolean g1h_Q2DecomposeMatrixf ( GHoleDomainf *domain )
{
  void  *sp;
  G1HolePrivateRecf *privateG1;
  float *lmat;
  int   nfunc_a;

  sp = pkv_GetScratchMemTop ();
  privateG1 = domain->privateG1;
  if ( !privateG1->Q2AMat )
    if ( !g1h_Q2ComputeFormMatrixf ( domain ) )
      goto failure;
  if ( !privateG1->Q2LMat ) {
    nfunc_a = privateG1->nfunc_a;
    if ( !(lmat = privateG1->Q2LMat) )
      lmat = privateG1->Q2LMat =
                 malloc ( ((nfunc_a*(nfunc_a+1))/2)*sizeof(float) );
    if ( !lmat )
      goto failure;
    memcpy ( lmat, privateG1->Q2AMat, ((nfunc_a*(nfunc_a+1))/2)*sizeof(float) );
    if ( !pkn_CholeskyDecompf ( nfunc_a, lmat ) ) {
      domain->error_code = G1H_ERROR_NONPOSITIVE_MATRIX;
      goto failure;
    }
  }

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*g1h_Q2DecomposeMatrixf*/

boolean g1h_Q2FillHolef ( GHoleDomainf *domain,
                          int spdimen, CONST_ float *hole_cp,
                          float *acoeff, void *usrptr,
                          void (*outpatch) ( int n, int m, const float *cp,
                                             void *usrptr ) )
{
  void  *sp;
  G1HolePrivateRecf *privateG1;
  int   hole_k, nfunc_a, nfunc_b;
  float *lmat, *b, *x, *fc00;

  sp = pkv_GetScratchMemTop ();
  hole_k = domain->hole_k;
  privateG1 = domain->privateG1;
  if ( !privateG1->Q2AMat )
    if ( !g1h_Q2ComputeFormMatrixf ( domain ) )
      goto failure;
  if ( !privateG1->Q2LMat )
    if ( !g1h_Q2DecomposeMatrixf ( domain ) )
      goto failure;

  nfunc_a = privateG1->nfunc_a;
  nfunc_b = privateG1->nfunc_b;
  lmat = privateG1->Q2LMat;
  x = pkv_GetScratchMemf ( spdimen*max(nfunc_a,nfunc_b) );
  b = pkv_GetScratchMemf ( spdimen*nfunc_a );
  fc00 = pkv_GetScratchMemf ( (G1_CROSSDEGSUM+4)*2*hole_k*spdimen );
  if ( !x || !b || !fc00 )
    goto failure;

  if ( !_g1h_SetRightSidef ( domain, privateG1->Q2BMat,
                             spdimen, hole_cp, fc00, b ) )
    goto failure;

  pkn_LowerTrMatrixSolvef ( nfunc_a, lmat, spdimen, spdimen, b, spdimen, x );
  pkn_UpperTrMatrixSolvef ( nfunc_a, lmat, spdimen, spdimen, x, spdimen, x );
  if ( acoeff )
    memcpy ( acoeff, x, spdimen*nfunc_a*sizeof(float) );

  if ( !_g1h_OutputPatchesf ( domain, spdimen, x, fc00, usrptr, outpatch ) )
    goto failure;

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*g1h_Q2FillHolef*/

