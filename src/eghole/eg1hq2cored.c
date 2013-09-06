
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

#include "eg1holed.h"
#include "eg1hprivated.h"
#include "eg1herror.h"

/*#define DEBUG*/

/* ///////////////////////////////////////////////////////////////////////// */
void g1h_DestroyQ2PrivateDatad ( GHoleDomaind *domain )
{
  G1HolePrivateRecd *privateG1;

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
} /*g1h_DestroyQ2PrivateDatad*/

/* ///////////////////////////////////////////////////////////////////////// */
boolean _g1h_Q2TabDiPatchJac3d ( int nkn, const double *kn, const double *hfunc,
             const double *dhfunc, const double *ddhfunc, const double *dddhfunc,
             const vector2d *c00, const vector2d *c01,
             const vector2d *c10, const vector2d *c11,
             const vector2d *d00, const vector2d *d01,
             const vector2d *d10, const vector2d *d11,
             double *jac, double *trd )
{
  void     *sp;
  int      i, k;
  vector2d *tabpu, *tabpv, *tabpuu, *tabpuv, *tabpvv,
           *tabpuuu, *tabpuuv, *tabpuvv, *tabpvvv;

  sp = pkv_GetScratchMemTop ();
  k = nkn*nkn;
  tabpu = (vector2d*)pkv_GetScratchMem ( 9*k*sizeof(vector2d) );
  if ( !tabpu )
    goto failure;
  tabpv = &tabpu[k];      tabpuu = &tabpv[k];    tabpuv = &tabpuu[k];
  tabpvv = &tabpuv[k];    tabpuuu = &tabpvv[k];  tabpuuv = &tabpuuu[k];
  tabpuvv = &tabpuuv[k];  tabpvvv = &tabpuvv[k];

  if ( !mbs_TabBezC1CoonsDer3d ( 2, nkn, kn, hfunc, dhfunc, ddhfunc, dddhfunc,
           nkn, kn, hfunc, dhfunc, ddhfunc, dddhfunc,
           G1_CROSS00DEG, (double*)c00, G1_CROSS01DEG, (double*)c01,
           G1_CROSS10DEG, (double*)c10, G1_CROSS11DEG, (double*)c11,
           G1_CROSS00DEG, (double*)d00, G1_CROSS01DEG, (double*)d01,
           G1_CROSS10DEG, (double*)d10, G1_CROSS11DEG, (double*)d11,
           NULL, (double*)tabpu, (double*)tabpv, (double*)tabpuu,
           (double*)tabpuv, (double*)tabpvv, (double*)tabpuuu, (double*)tabpuuv,
           (double*)tabpuvv, (double*)tabpvvv ) )
    goto failure;

  for ( i = 0; i < k; i++ )
    _g2h_DiJacobian3d ( &tabpu[i], &tabpv[i], &tabpuu[i], &tabpuv[i], &tabpvv[i],
                        &tabpuuu[i], &tabpuuv[i], &tabpuvv[i], &tabpvvv[i],
                        &jac[i], &trd[18*i] );

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*_g1h_Q2TabDiPatchJac3d*/

boolean _g1h_Q2TabLaplacianGradd ( int nkn, const double *tkn,
        const double *hfunc, const double *dhfunc, const double *ddhfunc,
        const double *dddhfunc,
        const double *fc00, const double *fc01,
        const double *fc10, const double *fc11,
        const double *fd00, const double *fd01,
        const double *fd10, const double *fd11,
        const double *trd,
        vector2d *lapgrad )
{
  void   *sp;
  int    i, k;
  double *tabfu, *tabfv, *tabfuu, *tabfuv, *tabfvv,
         *tabfuuu, *tabfuuv, *tabfuvv, *tabfvvv;

  sp = pkv_GetScratchMemTop ();
  k = nkn*nkn;

  tabfu = pkv_GetScratchMemd ( 9*k );
  if ( !tabfu )
    goto failure;
  tabfv = &tabfu[k];      tabfuu = &tabfv[k];    tabfuv = &tabfuu[k];
  tabfvv = &tabfuv[k];    tabfuuu = &tabfvv[k];  tabfuuv = &tabfuuu[k];
  tabfuvv = &tabfuuv[k];  tabfvvv = &tabfuvv[k];

  if ( !mbs_TabBezC1CoonsDer3d ( 1, nkn, tkn, hfunc, dhfunc, ddhfunc, dddhfunc,
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
} /*_g1h_Q2TabLaplacianGradd*/

boolean _g1h_Q2TabLaplacianGrad0d ( int nkn, const double *tkn,
        const double *hfunc, const double *dhfunc, const double *ddhfunc,
        const double *dddhfunc,
        const double *fc00, const double *fc01,
        const double *fd00, const double *fd01,
        const double *trd,
        vector2d *lapgrad )
{
  void   *sp;
  int    i, k;
  double *tabfu, *tabfv, *tabfuu, *tabfuv, *tabfvv,
         *tabfuuu, *tabfuuv, *tabfuvv, *tabfvvv;

  sp = pkv_GetScratchMemTop ();
  k = nkn*nkn;

  tabfu = pkv_GetScratchMemd ( 9*k );
  if ( !tabfu )
    goto failure;
  tabfv = &tabfu[k];      tabfuu = &tabfv[k];    tabfuv = &tabfuu[k];
  tabfvv = &tabfuv[k];    tabfuuu = &tabfvv[k];  tabfuuv = &tabfuuu[k];
  tabfuvv = &tabfuuv[k];  tabfvvv = &tabfuvv[k];

  if ( !mbs_TabBezC1Coons0Der3d ( 1, nkn, tkn, hfunc, dhfunc, ddhfunc, dddhfunc,
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
} /*_g1h_Q2TabLaplacianGrad0d*/

void _g1h_TabCurveJacobiand ( int deg, const point2d *cp,
                              int nkn, const double *kn, double *jac )
{
  int      i;
  point2d  p;
  vector2d dp;

  for ( i = 0; i < nkn; i++ ) {
    mbs_BCHornerDerC2d ( deg, cp, kn[i], &p, &dp );
    jac[i] = sqrt ( dp.x*dp.x+dp.y*dp.y );
  }
} /*_g1h_TabCurveJacobiand*/

void _g1h_LapCoeffd ( const vector2d *du, const vector2d *dv,
                     const vector2d *duu, const vector2d *duv,
                     const vector2d *dvv, double *trd )
{
  vector2d gx, gy, gxx, gxy, gyy;
  double   A21[6], A22[9];

  pkn_f2iDerivatives2d ( du->x, du->y, dv->x, dv->y,
      duu->x, duu->y, duv->x, duv->y, dvv->x, dvv->y,
      (double*)&gx, (double*)&gy, (double*)&gxx, (double*)&gxy, (double*)&gyy );

  pkn_Setup2DerA21Matrixd ( gxx.x, gxx.y, gxy.x, gxy.y, gyy.x, gyy.y, A21 );
  pkn_Setup2DerA22Matrixd ( gx.x, gx.y, gy.x, gy.y, A22 );

  trd[0]  = A21[0]+A21[4];   trd[1]  = A21[1]+A21[5];
  trd[2]  = A22[0]+A22[6];   trd[3]  = A22[1]+A22[7];   trd[4] = A22[2]+A22[8];
} /*_g1h_LapCoeffd*/

boolean _g1h_TabCurveLapCoeff0d ( const point2d *c00, const vector2d *c01,
                                  const point2d *c10, const vector2d *c11,
                                  const point2d *d00, const vector2d *d01,
                                  const point2d *d10, const vector2d *d11,
                                  int nkn, const double *tkn, const double *hfunc,
                                  const double *dhfunc, const double *ddhfunc,
                                  const double *atkn, const double *ahfunc,
                                  const double *adhfunc, const double *addhfunc,
                                  double *trdc00, double *trdc10,
                                  double *trdd00, double *trdd10 )
{
  void     *sp;
  vector2d *tabpu, *tabpv, *tabpuu, *tabpuv, *tabpvv;
  int      i, ii;

  sp = pkv_GetScratchMemTop ();
  tabpu = (vector2d*)pkv_GetScratchMem ( 10*nkn*sizeof(vector2d) );
  if ( !tabpu )
    goto failure;
  tabpv = &tabpu[2*nkn];    tabpuu = &tabpv[2*nkn];
  tabpuv = &tabpuu[2*nkn];  tabpvv = &tabpuv[2*nkn];

  if ( !mbs_TabBezC1CoonsDer2d ( 2, nkn, tkn, hfunc, dhfunc, ddhfunc,
                         2, atkn, ahfunc, adhfunc, addhfunc,
                         G1_CROSS00DEG, (double*)c00, G1_CROSS01DEG, (double*)c01,
                         G1_CROSS10DEG, (double*)c10, G1_CROSS11DEG, (double*)c11,
                         G1_CROSS00DEG, (double*)d00, G1_CROSS01DEG, (double*)d01,
                         G1_CROSS10DEG, (double*)d10, G1_CROSS11DEG, (double*)d11,
                         NULL, (double*)tabpu, (double*)tabpv,
                         (double*)tabpuu, (double*)tabpuv, (double*)tabpvv ) )
    goto failure;
  for ( i = ii = 0;  i < nkn;  i++ ) {
    _g1h_LapCoeffd ( &tabpu[ii], &tabpv[ii],
                     &tabpuu[ii], &tabpuv[ii], &tabpvv[ii], &trdc00[5*i] );
    ii++;
    _g1h_LapCoeffd ( &tabpu[ii], &tabpv[ii],
                     &tabpuu[ii], &tabpuv[ii], &tabpvv[ii], &trdc10[5*i] );
    ii++;
  }

  if ( !mbs_TabBezC1CoonsDer2d ( 2, 2, atkn, ahfunc, adhfunc, addhfunc,
                         nkn, tkn, hfunc, dhfunc, ddhfunc,
                         G1_CROSS00DEG, (double*)c00, G1_CROSS01DEG, (double*)c01,
                         G1_CROSS10DEG, (double*)c10, G1_CROSS11DEG, (double*)c11,
                         G1_CROSS00DEG, (double*)d00, G1_CROSS01DEG, (double*)d01,
                         G1_CROSS10DEG, (double*)d10, G1_CROSS11DEG, (double*)d11,
                         NULL, (double*)tabpu, (double*)tabpv,
                         (double*)tabpuu, (double*)tabpuv, (double*)tabpvv ) )
    goto failure;
  for ( i = 0, ii = nkn;  i < nkn;  i++, ii++ ) {
    _g1h_LapCoeffd ( &tabpu[i], &tabpv[i],
                     &tabpuu[i], &tabpuv[i], &tabpvv[i], &trdd00[5*i] );
    _g1h_LapCoeffd ( &tabpu[ii], &tabpv[ii],
                     &tabpuu[ii], &tabpuv[ii], &tabpvv[ii], &trdd10[5*i] );
  }

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*_g1h_TabCurveLapCoeff0d*/

void _g1h_TabCurveLapCoeff1d ( const point2d *sicp, int nkn,
                               const double *tkn, double *trd )
{
  int i;
  point2d d;
  vector2d du, dv, duu, duv, dvv;

  for ( i = 0; i < nkn; i++ ) {
    mbs_BCHornerDer2Pd ( 3, 3, 2, (double*)sicp, 0.0, tkn[i],
                         &d.x, &du.x, &dv.x, &duu.x, &duv.x, &dvv.x );
    _g1h_LapCoeffd ( &du, &dv, &duu, &duv, &dvv, &trd[5*i] );
  }
} /*_g1h_TabCurveLapCoeff1d*/

boolean _g1h_Q2TabLaplacianJump0d ( int nkn, const double *tkn,
              const double *hfunc, const double *dhfunc, const double *ddhfunc,
              const double *atkn, const double *ahfunc, const double *adhfunc,
              const double *addhfunc,
              const double *ec00, const double *ec01,
              const double *ed00, const double *ed01, const double *etrdd00,
              const double *fc00, const double *fc01,
              const double *fd00, const double *fd01,
              const double *ftrdc00, const double *ftrdc10, const double *ftrdd10,
              double *lapjc00, double *lapjc10, double *lapjd10 )
{
  void   *sp;
  int    i, j, k;
  double *tabeu, *tabev, *tabeuu, *tabeuv, *tabevv,
         *tabfu, *tabfv, *tabfuu, *tabfuv, *tabfvv;
  double lape, lapf;

  sp = pkv_GetScratchMemTop ();
  tabeu = pkv_GetScratchMemd ( 20*nkn );
  if ( !tabeu )
    goto failure;
  tabev = &tabeu[2*nkn];    tabeuu = &tabev[2*nkn];   tabeuv = &tabeuu[2*nkn];
  tabevv = &tabeuv[2*nkn];  tabfu = &tabevv[2*nkn];   tabfv = &tabfu[2*nkn];
  tabfuu = &tabfv[2*nkn];   tabfuv = &tabfuu[2*nkn];  tabfvv = &tabfuv[2*nkn];
  if ( !mbs_TabBezC1Coons0Der2d ( 1, nkn, tkn, hfunc, dhfunc, ddhfunc,
              2, atkn, ahfunc, adhfunc, addhfunc,
              G1_CROSS00DEG, fc00, G1_CROSS01DEG, fc01,
              G1_CROSS00DEG, fd00, G1_CROSS01DEG, fd01,
              NULL, tabfu, tabfv, tabfuu, tabfuv, tabfvv ) )
    goto failure;
  if ( !mbs_TabBezC1Coons0Der2d ( 1, 1, atkn, ahfunc, adhfunc, addhfunc,
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

  if ( !mbs_TabBezC1Coons0Der2d ( 1,
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
} /*_g1h_Q2TabLaplacianJump0d*/

boolean _g1h_Q2TabLaplacianJumpd ( int nkn, const double *tkn,
              const double *hfunc, const double *dhfunc, const double *ddhfunc,
              const double *atkn, const double *ahfunc, const double *adhfunc,
              const double *addhfunc,
              const double *ec00, const double *ec01,
              const double *ec10, const double *ec11,
              const double *ed00, const double *ed01,
              const double *ed10, const double *ed11, const double *etrdd00,
              const double *fc00, const double *fc01,
              const double *fc10, const double *fc11,
              const double *fd00, const double *fd01,
              const double *fd10, const double *fd11,
              const double *ftrdc00, const double *ftrdc10, const double *ftrdd10,
              const double *eicp1, const double *etrdc10,
              const double *eicp2, const double *etrdd10,
              double *lapjc00, double *lapjc10, double *lapjd10 )
{
  void   *sp;
  int    i, j, k;
  double *tabeu, *tabev, *tabeuu, *tabeuv, *tabevv,
         *tabfu, *tabfv, *tabfuu, *tabfuv, *tabfvv;
  double lape, lapf;

  sp = pkv_GetScratchMemTop ();
  tabeu = pkv_GetScratchMemd ( 15*nkn );
  if ( !tabeu )
    goto failure;
  tabev = &tabeu[nkn];      tabeuu = &tabev[nkn];     tabeuv = &tabeuu[nkn];
  tabevv = &tabeuv[nkn];    tabfu = &tabevv[nkn];     tabfv = &tabfu[2*nkn];
  tabfuu = &tabfv[2*nkn];   tabfuv = &tabfuu[2*nkn];  tabfvv = &tabfuv[2*nkn];
  if ( !mbs_TabBezC1CoonsDer2d ( 1, nkn, tkn, hfunc, dhfunc, ddhfunc,
              2, atkn, ahfunc, adhfunc, addhfunc,
              G1_CROSS00DEG, fc00, G1_CROSS01DEG, fc01,
              G1_CROSS10DEG, fc10, G1_CROSS11DEG, fc11,
              G1_CROSS00DEG, fd00, G1_CROSS01DEG, fd01,
              G1_CROSS10DEG, fd10, G1_CROSS11DEG, fd11,
              NULL, tabfu, tabfv, tabfuu, tabfuv, tabfvv ) )
    goto failure;
  if ( !mbs_TabBezC1CoonsDer2d ( 1, 1, atkn, ahfunc, adhfunc, addhfunc,
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
    mbs_BCHornerDer2Pd ( 3, 3, 1, eicp1, 0.0, tkn[i],
                         &tabeu[1], tabeu, tabev, tabeuu, tabeuv, tabevv );
    lape = etrdc10[j]*tabeu[0] + etrdc10[j+1]*tabev[0] + etrdc10[j+2]*tabeuu[0] +
           etrdc10[j+3]*tabeuv[0] + etrdc10[j+4]*tabevv[0];
    lapf = ftrdc10[j]*tabfu[k] + ftrdc10[j+1]*tabfv[k] + ftrdc10[j+2]*tabfuu[k] +
           ftrdc10[j+3]*tabfuv[k] + ftrdc10[j+4]*tabfvv[k];
    lapjc10[i] = lape-lapf;
  }

  if ( !mbs_TabBezC1CoonsDer2d ( 1,
              1, &atkn[1], &ahfunc[4], &adhfunc[4], &addhfunc[4],
              nkn, tkn, hfunc, dhfunc, ddhfunc,
              G1_CROSS00DEG, fc00, G1_CROSS01DEG, fc01,
              G1_CROSS10DEG, fc10, G1_CROSS11DEG, fc11,
              G1_CROSS00DEG, fd00, G1_CROSS01DEG, fd01,
              G1_CROSS10DEG, fd10, G1_CROSS11DEG, fd11,
              NULL, tabfu, tabfv, tabfuu, tabfuv, tabfvv ) )
    goto failure;
  for ( i = j = 0;  i < nkn;  i++, j += 5 ) {
    mbs_BCHornerDer2Pd ( 3, 3, 1, eicp2, 0.0, tkn[i],
                         &tabeu[1], tabeu, tabev, tabeuu, tabeuv, tabevv );
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
} /*_g1h_Q2TabLaplacianJumpd*/

double _g1h_Q2Integrald ( int hole_k, int nquad, double *jac,
                          unsigned short supp1, double *lapj1,
                          unsigned short supp2, double *lapj2 )
{
  int            i, k;
  long double    s;
  unsigned short supp;
  double         *jc, *f1, *f2;

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
  return (double)s;
} /*_g1h_Q2Integrald*/

boolean g1h_Q2ComputeFormMatrixd ( GHoleDomaind *domain )
{
  void     *sp;
  GHolePrivateRecd    *privateG;
  G1HolePrivateRecd   *privateG1;
  int      hole_k, nfunc_a, nfunc_b;
  int      i, j, l, n, fn;
  double   *tkn, *hfunc, *dhfunc, *ddhfunc, *dddhfunc;
  double   *atkn, *ahfunc, *adhfunc, *addhfunc;
  vector2d *c00, *c01, *c10, *c11, *d00, *d01, *d10, *d11;
  double   *ec00, *ec01, *ec10, *ec11, *ed00, *ed01, *ed10, *ed11;
  double   *fc00, *fc01, *fc10, *fc11, *fd00, *fd01, *fd10, *fd11;
  double   *trd, *lapj, *jac, *amat, *bmat;
  vector2d *lgr;
  unsigned short *support_b, supp;
  int      option, ndata, *idata;
  double   *fdata, *eicp1, *eicp2, C;
  point2d  *sicp;
#ifdef DEBUG
FILE *f;
#endif

  sp = pkv_GetScratchMemTop ();
  if ( !domain->basisG1 ) {
    if ( !g1h_ComputeBasisd ( domain ) )
      goto failure;
  }
  hole_k    = domain->hole_k;
  privateG  = domain->privateG;
  privateG1 = domain->privateG1;
  nfunc_a   = privateG1->nfunc_a;
  nfunc_b   = privateG1->nfunc_b;
  support_b = privateG->support_b;
  if ( !(amat = privateG1->Q2AMat) )
    amat = privateG1->Q2AMat = malloc ( (nfunc_a*(nfunc_a+1)/2)*sizeof(double) );
  if ( !(bmat = privateG1->Q2BMat) )
    bmat = privateG1->Q2BMat = malloc ( nfunc_a*nfunc_b*sizeof(double) );
  if ( !amat || !bmat )
    goto failure;
  if ( !privateG->diam )
    privateG->diam = gh_DomainDiamd ( domain );

  n = hole_k*G1_NQUADSQ;

  tkn = pkv_GetScratchMemd ( 17*G1_NQUAD );
  jac = pkv_GetScratchMemd ( n );
  trd = pkv_GetScratchMemd ( 18*n );
  lgr = (vector2d*)pkv_GetScratchMem ( (nfunc_a+1)*n*sizeof(vector2d) );
  if ( !tkn || !jac || !trd || !lgr )
    goto failure;
  hfunc = &tkn[G1_NQUAD];         dhfunc = &hfunc[4*G1_NQUAD];
  ddhfunc = &dhfunc[4*G1_NQUAD];  dddhfunc = &ddhfunc[4*G1_NQUAD];
  _gh_PrepareTabKnotsd ( G1_NQUAD, privateG1->opt_quad, tkn );
  mbs_TabCubicHFuncDer3d ( 0.0, 1.0, G1_NQUAD, tkn,
                           hfunc, dhfunc, ddhfunc, dddhfunc );

        /* integrate the Laplacian gradients */
  for ( i = 0; i < hole_k; i++ ) {
    _g1h_GetDiPatchCurvesd ( domain, i,
                             &c00, &c01, &c10, &c11, &d00, &d01, &d10, &d11 );
    _g1h_Q2TabDiPatchJac3d ( G1_NQUAD, tkn, hfunc, dhfunc, ddhfunc, dddhfunc,
                             c00, c01, c10, c11, d00, d01, d10, d11,
                             &jac[i*G1_NQUADSQ], &trd[i*G1_NQUADSQ*18] );
  }

  for ( fn = 0; fn < nfunc_a; fn++ )
    for ( i = 0; i < hole_k; i++ ) {
      _g1h_GetBFAPatchCurvesd ( domain, fn, i, &fc00, &fc01, &fd00, &fd01 );
      _g1h_Q2TabLaplacianGrad0d ( G1_NQUAD, tkn, hfunc, dhfunc, ddhfunc, dddhfunc,
                                  fc00, fc01, fd00, fd01,
                                  &trd[i*G1_NQUADSQ*18],
                                  &lgr[(fn*hole_k+i)*G1_NQUADSQ] );
    }
  for ( i = l = 0;  i < nfunc_a;  i++ ) {
    for ( j = 0;  j <= i;  j++, l++ )
      amat[l] = _g2h_Integrald ( hole_k, G1_NQUADSQ, jac,
                                 0xFFFF, &lgr[i*n], 0xFFFF, &lgr[j*n] );
  }

  for ( fn = 0; fn < nfunc_b; fn++ ) {
    for ( i = 0; i < hole_k; i++ )
      if ( support_b[fn] & (0x0001 << i) ) {
        _g1h_GetBFBPatchCurvesd ( domain, fn, i, &fc00, &fc01, &fc10, &fc11,
                                                 &fd00, &fd01, &fd10, &fd11 );
        _g1h_Q2TabLaplacianGradd ( G1_NQUAD, tkn, hfunc, dhfunc, ddhfunc, dddhfunc,
                                   fc00, fc01, fc10, fc11, fd00, fd01, fd10, fd11,
                                   &trd[i*G1_NQUADSQ*18],
                                   &lgr[(nfunc_a*hole_k+i)*G1_NQUADSQ] );
      }
    for ( i = 0; i < nfunc_a; i++ )
      bmat[i*nfunc_b+fn] = _g2h_Integrald ( hole_k, G1_NQUADSQ, jac,
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
  C *= 5.0/privateG->diam;

        /* the arrays allocated for the previous stage are reused */
        /* their length is more than sufficient */
  lapj = &lgr[0].x;
        /* two additional arrays are however needed */
  atkn = pkv_GetScratchMemd ( 26 );
  sicp = (point2d*)pkv_GetScratchMem ( 16*sizeof(point2d) );
  if ( !atkn || !sicp )
    goto failure;
  ahfunc = &atkn[2];  adhfunc = &ahfunc[8];  addhfunc = &adhfunc[8];
  eicp1 = (double*)sicp;  eicp2 = &eicp1[16];
  atkn[0] = 0.0;  atkn[1] = 1.0;
  mbs_TabCubicHFuncDer2d ( 0.0, 1.0, 2, atkn, ahfunc, adhfunc, addhfunc );

  for ( i = 0; i < hole_k; i++ ) {
    j = (i+1) % hole_k;
          /* compute the curve Jacobians */
    _g1h_GetDiPatchCurvesd ( domain, i,
                             &c00, &c01, &c10, &c11, &d00, &d01, &d10, &d11 );
    _g1h_TabCurveJacobiand ( G1_CROSS00DEG, c00, G1_NQUAD, tkn, &jac[3*i*G1_NQUAD] );
    _g1h_TabCurveJacobiand ( G1_CROSS10DEG, c10, G1_NQUAD, tkn, &jac[(3*i+1)*G1_NQUAD] );
    _g1h_TabCurveJacobiand ( G1_CROSS10DEG, d10, G1_NQUAD, tkn, &jac[(3*i+2)*G1_NQUAD] );
          /* compute the coefficients for computing the Laplacians */
    if ( !_g1h_TabCurveLapCoeff0d ( c00, c01, c10, c11, d00, d01, d10, d11,
                            G1_NQUAD, tkn, hfunc, dhfunc, ddhfunc,
                            atkn, ahfunc, adhfunc, addhfunc,
                            &trd[30*i*G1_NQUAD], &trd[(30*i+10)*G1_NQUAD],
                            &trd[(30*j+5)*G1_NQUAD], &trd[(30*i+20)*G1_NQUAD] ) )
      goto failure;
    if ( !_gh_FindDomSurrndPatchd ( domain, j, 1, sicp ) )
      goto failure;
    _g1h_TabCurveLapCoeff1d ( sicp, G1_NQUAD, tkn, &trd[(30*i+15)*G1_NQUAD] );
    if ( !_gh_FindDomSurrndPatchd ( domain, i, 2, sicp ) )
      goto failure;
    _g1h_TabCurveLapCoeff1d ( sicp, G1_NQUAD, tkn, &trd[(30*i+25)*G1_NQUAD] );
  }

  for ( fn = 0; fn < nfunc_a; fn++ ) {
    _g1h_GetBFAPatchCurvesd ( domain, fn, hole_k-1,
                              &ec00, &ec01, &ed00, &ed01 );
    for ( i = 0; i < hole_k; i++ ) {
      _g1h_GetBFAPatchCurvesd ( domain, fn, i,
                                &fc00, &fc01, &fd00, &fd01 );
      if ( !_g1h_Q2TabLaplacianJump0d ( G1_NQUAD, tkn, hfunc, dhfunc, ddhfunc,
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
      amat[l] += C*_g1h_Q2Integrald ( hole_k, G1_NQUAD, jac,
                                      0xFFFF, &lapj[i*hole_k*3*G1_NQUAD],
                                      0xFFFF, &lapj[j*hole_k*3*G1_NQUAD] );
  }
  for ( fn = 0; fn < nfunc_b; fn++ ) {
    supp = _g1h_ExtendSupport ( hole_k, support_b[fn] );
    for ( i = 0; i < hole_k; i++ ) {
      j = (i+1) % hole_k;
      if ( supp & (0x0001 << i) ) {
        _g1h_GetBFBPatchCurvesd ( domain, fn, (i+hole_k-1) % hole_k,
                  &ec00, &ec01, &ec10, &ec11, &ed00, &ed01, &ed10, &ed11 );
        _g1h_GetBFBPatchCurvesd ( domain, fn, i,
                  &fc00, &fc01, &fc10, &fc11, &fd00, &fd01, &fd10, &fd11 );
        gh_GetDomSurrndBFuncd ( domain, fn, j, 1, eicp1 );
        gh_GetDomSurrndBFuncd ( domain, fn, i, 2, eicp2 );
        if ( !_g1h_Q2TabLaplacianJumpd ( G1_NQUAD, tkn, hfunc, dhfunc, ddhfunc,
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
      bmat[i*nfunc_b+fn] += C*_g1h_Q2Integrald ( hole_k, G1_NQUAD, jac,
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
} /*g1h_Q2ComputeFormMatrixd*/

boolean g1h_Q2DecomposeMatrixd ( GHoleDomaind *domain )
{
  void   *sp;
  G1HolePrivateRecd *privateG1;
  double *lmat;
  int    nfunc_a;

  sp = pkv_GetScratchMemTop ();
  privateG1 = domain->privateG1;
  if ( !privateG1->Q2AMat )
    if ( !g1h_Q2ComputeFormMatrixd ( domain ) )
      goto failure;
  if ( !privateG1->Q2LMat ) {
    nfunc_a = privateG1->nfunc_a;
    if ( !(lmat = privateG1->Q2LMat) )
      lmat = privateG1->Q2LMat =
                 malloc ( ((nfunc_a*(nfunc_a+1))/2)*sizeof(double) );
    if ( !lmat )
      goto failure;
    memcpy ( lmat, privateG1->Q2AMat, ((nfunc_a*(nfunc_a+1))/2)*sizeof(double) );
    if ( !pkn_CholeskyDecompd ( nfunc_a, lmat ) ) {
      domain->error_code = G1H_ERROR_NONPOSITIVE_MATRIX;
      goto failure;
    }
  }

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*g1h_Q2DecomposeMatrixd*/

boolean g1h_Q2FillHoled ( GHoleDomaind *domain,
                          int spdimen, CONST_ double *hole_cp,
                          double *acoeff, void *usrptr,
                          void (*outpatch) ( int n, int m, const double *cp,
                                             void *usrptr ) )
{
  void   *sp;
  G1HolePrivateRecd *privateG1;
  int    hole_k, nfunc_a, nfunc_b;
  double *lmat, *b, *x, *fc00;

  sp = pkv_GetScratchMemTop ();
  hole_k = domain->hole_k;
  privateG1 = domain->privateG1;
  if ( !privateG1->Q2AMat )
    if ( !g1h_Q2ComputeFormMatrixd ( domain ) )
      goto failure;
  if ( !privateG1->Q2LMat )
    if ( !g1h_Q2DecomposeMatrixd ( domain ) )
      goto failure;

  nfunc_a = privateG1->nfunc_a;
  nfunc_b = privateG1->nfunc_b;
  lmat = privateG1->Q2LMat;
  x = pkv_GetScratchMemd ( spdimen*max(nfunc_a,nfunc_b) );
  b = pkv_GetScratchMemd ( spdimen*nfunc_a );
  fc00 = pkv_GetScratchMemd ( (G1_CROSSDEGSUM+4)*2*hole_k*spdimen );
  if ( !x || !b || !fc00 )
    goto failure;

  if ( !_g1h_SetRightSided ( domain, privateG1->Q2BMat,
                             spdimen, hole_cp, fc00, b ) )
    goto failure;

  pkn_LowerTrMatrixSolved ( nfunc_a, lmat, spdimen, spdimen, b, spdimen, x );
  pkn_UpperTrMatrixSolved ( nfunc_a, lmat, spdimen, spdimen, x, spdimen, x );
  if ( acoeff )
    memcpy ( acoeff, x, spdimen*nfunc_a*sizeof(double) );

  if ( !_g1h_OutputPatchesd ( domain, spdimen, x, fc00, usrptr, outpatch ) )
    goto failure;

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*g1h_Q2FillHoled*/

