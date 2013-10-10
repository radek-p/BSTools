
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

#include "eg1holef.h"
#include "eg2holef.h"
#include "eg1hprivatef.h"
#include "eg2hprivatef.h"
#include "eg1herror.h"


#define _DEBUG_HESSIAN
#define _DEBUG_FVAL

#ifdef DEBUG_FVAL
/* ////////////////////////////////////////////////////////////////////////// */
static boolean ComputePDer ( GHoleDomainf *domain, G1HNLPrivatef *nlpr,
                             const float *acoeff,
                             int k, int cn, float knot,
                             float *p, float *pu, float *pv,
                             float *jpuu, float *jpuv, float *jpvv, vector2f *tang )
{
#define N ((G1H_FINALDEG+1)*(G1H_FINALDEG+1))
  void     *sp;
  G1HolePrivateRecf *privateG1;
  GHolePrivateRecf  *privateG;
  point2f  *dicp;
  float    *fc00, *fc01, *fc10, *fc11, *fd00, *fd01, *fd10, *fd11;
  float    *pc00, *pc01, *pc10, *pc11, *pd00, *pd01, *pd10, *pd11, *pcp;
  vector2f di, diu, div, diuu, diuv, divv;
  vector2f ei, eiu, eiv, eiuu, eiuv, eivv;
  float    e, eu, ev, euu, euv, evv, f, fu, fv, fuu, fuv, fvv,
           qu, qv, puu, puv, pvv, bvz;
  int      hole_k, i, j, n, m, fn, nfunc_a, nfunc_b, nfunc_c, nfunc;
  unsigned char *bfcpn;

  sp = pkv_GetScratchMemTop ();
  privateG = domain->privateG;
  privateG1 = domain->privateG1;
  bfcpn = privateG->bfcpn;
  nfunc_a = privateG1->nfunc_a;
  nfunc_b = privateG1->nfunc_b;
  hole_k = domain->hole_k;
  nfunc_c = hole_k*G1_DBDIM;
  nfunc = nfunc_a+nfunc_c;
  dicp = pkv_GetScratchMem ( N*sizeof(point2f) );
  pc00 = pkv_GetScratchMemf ( 2*(G1_CROSSDEGSUM+4) );
  pcp = pkv_GetScratchMemf ( N );
  if ( !dicp || !pc00 || !pcp ) {
    PKV_SIGNALERROR ( LIB_EGHOLE, ERRCODE_2, ERRMSG_2 );
    goto failure;
  }
  pc01 = &pc00[G1_CROSS00DEG+1];  pc10 = &pc01[G1_CROSS01DEG+1];
  pc11 = &pc10[G1_CROSS10DEG+1];  pd00 = &pc11[G1_CROSS11DEG+1];
  pd01 = &pd00[G1_CROSS00DEG+1];  pd10 = &pd01[G1_CROSS01DEG+1];
  pd11 = &pd10[G1_CROSS10DEG+1];

  pkv_Selectf ( N, 2, 3, 2, &nlpr->nldi[k*N], dicp );
  memset ( pc00, 0, 2*(G1_CROSSDEGSUM+4)*sizeof(float) );
  for ( fn = 0; fn < nfunc_b; fn++ ) {
    bvz = nlpr->rhole_cp[bfcpn[fn]].z;
    _g1h_GetBFBPatchCurvesf ( domain, fn, k, &fc00, &fc01, &fc10, &fc11,
                                             &fd00, &fd01, &fd10, &fd11 );
    pkn_AddMatrixMf ( 1, G1_CROSS00DEG+1, 0, pc00, 0, fc00, bvz, 0, pc00 );
    pkn_AddMatrixMf ( 1, G1_CROSS01DEG+1, 0, pc01, 0, fc01, bvz, 0, pc01 );
    pkn_AddMatrixMf ( 1, G1_CROSS10DEG+1, 0, pc10, 0, fc10, bvz, 0, pc10 );
    pkn_AddMatrixMf ( 1, G1_CROSS11DEG+1, 0, pc11, 0, fc11, bvz, 0, pc11 );
    pkn_AddMatrixMf ( 1, G1_CROSS00DEG+1, 0, pd00, 0, fd00, bvz, 0, pd00 );
    pkn_AddMatrixMf ( 1, G1_CROSS01DEG+1, 0, pd01, 0, fd01, bvz, 0, pd01 );
    pkn_AddMatrixMf ( 1, G1_CROSS10DEG+1, 0, pd10, 0, fd10, bvz, 0, pd10 );
    pkn_AddMatrixMf ( 1, G1_CROSS11DEG+1, 0, pd11, 0, fd11, bvz, 0, pd11 );
  }
  for ( fn = 0; fn < nfunc_a; fn++ ) {
    bvz = -acoeff[nfunc_c+fn];
    _g1h_GetBFAPatchCurvesf ( domain, fn, k, &fc00, &fc01, &fd00, &fd01 );
    pkn_AddMatrixMf ( 1, G1_CROSS00DEG+1, 0, pc00, 0, fc00, bvz, 0, pc00 );
    pkn_AddMatrixMf ( 1, G1_CROSS01DEG+1, 0, pc01, 0, fc01, bvz, 0, pc01 );
    pkn_AddMatrixMf ( 1, G1_CROSS00DEG+1, 0, pd00, 0, fd00, bvz, 0, pd00 );
    pkn_AddMatrixMf ( 1, G1_CROSS01DEG+1, 0, pd01, 0, fd01, bvz, 0, pd01 );
  }
  mbs_BezC1CoonsToBezf ( 1, G1_CROSS00DEG, pc00, G1_CROSS01DEG, pc01,
                         G1_CROSS10DEG, pc10, G1_CROSS11DEG, pc11,
                         G1_CROSS00DEG, pd00, G1_CROSS01DEG, pd01,
                         G1_CROSS10DEG, pd10, G1_CROSS11DEG, pd11,
                         &n, &m, pcp );
  for ( i = 0; i < G1H_FINALDEG-3; i++ )
    for ( j = 0; j < G1H_FINALDEG-3; j++ ) {
      fn = k*G1_DBDIM+i*(G1H_FINALDEG-3)+j;
      pcp[(i+2)*(G1H_FINALDEG+1)+j+2] -= acoeff[fn];
    }
  switch ( cn ) {
case 0:
    if ( !mbs_BCHornerDer2Pf ( G1H_FINALDEG, G1H_FINALDEG, 2, &dicp[0].x, knot, 0.0,
                               &di.x, &diu.x, &div.x, &diuu.x, &diuv.x, &divv.x ) )
      goto failure;
    *tang = diu;
    if ( !mbs_BCHornerDer2Pf ( n, m, 1, pcp, knot, 0.0,
                               &f, &fu, &fv, &fuu, &fuv, &fvv ) )
      goto failure;
    if ( !pkn_Comp2iDerivatives2f ( diu.x, diu.y, div.x, div.y, diuu.x, diuu.y,
                                    diuv.x, diuv.y, divv.x, divv.y,
                                    1, &fu, &fv, &fuu, &fuv, &fvv,
                                    pu, pv, jpuu, jpuv, jpvv ) )
      goto failure;
    k = (k+hole_k-1) % hole_k;
    pkv_Selectf ( N, 2, 3, 2, &nlpr->nldi[k*N], dicp );
    memset ( pc00, 0, 2*(G1_CROSSDEGSUM+4)*sizeof(float) );
    for ( fn = 0; fn < nfunc_b; fn++ ) {
      bvz = nlpr->rhole_cp[bfcpn[fn]].z;
      _g1h_GetBFBPatchCurvesf ( domain, fn, k, &fc00, &fc01, &fc10, &fc11,
                                               &fd00, &fd01, &fd10, &fd11 );
      pkn_AddMatrixMf ( 1, G1_CROSS00DEG+1, 0, pc00, 0, fc00, bvz, 0, pc00 );
      pkn_AddMatrixMf ( 1, G1_CROSS01DEG+1, 0, pc01, 0, fc01, bvz, 0, pc01 );
      pkn_AddMatrixMf ( 1, G1_CROSS10DEG+1, 0, pc10, 0, fc10, bvz, 0, pc10 );
      pkn_AddMatrixMf ( 1, G1_CROSS11DEG+1, 0, pc11, 0, fc11, bvz, 0, pc11 );
      pkn_AddMatrixMf ( 1, G1_CROSS00DEG+1, 0, pd00, 0, fd00, bvz, 0, pd00 );
      pkn_AddMatrixMf ( 1, G1_CROSS01DEG+1, 0, pd01, 0, fd01, bvz, 0, pd01 );
      pkn_AddMatrixMf ( 1, G1_CROSS10DEG+1, 0, pd10, 0, fd10, bvz, 0, pd10 );
      pkn_AddMatrixMf ( 1, G1_CROSS11DEG+1, 0, pd11, 0, fd11, bvz, 0, pd11 );
    }
    for ( fn = 0; fn < nfunc_a; fn++ ) {
      bvz = -acoeff[nfunc_c+fn];
      _g1h_GetBFAPatchCurvesf ( domain, fn, k, &fc00, &fc01, &fd00, &fd01 );
      pkn_AddMatrixMf ( 1, G1_CROSS00DEG+1, 0, pc00, 0, fc00, bvz, 0, pc00 );
      pkn_AddMatrixMf ( 1, G1_CROSS01DEG+1, 0, pc01, 0, fc01, bvz, 0, pc01 );
      pkn_AddMatrixMf ( 1, G1_CROSS00DEG+1, 0, pd00, 0, fd00, bvz, 0, pd00 );
      pkn_AddMatrixMf ( 1, G1_CROSS01DEG+1, 0, pd01, 0, fd01, bvz, 0, pd01 );
    }
    mbs_BezC1CoonsToBezf ( 1, G1_CROSS00DEG, pc00, G1_CROSS01DEG, pc01,
                           G1_CROSS10DEG, pc10, G1_CROSS11DEG, pc11,
                           G1_CROSS00DEG, pd00, G1_CROSS01DEG, pd01,
                           G1_CROSS10DEG, pd10, G1_CROSS11DEG, pd11,
                           &n, &m, pcp );
    for ( i = 0; i < G1H_FINALDEG-3; i++ )
      for ( j = 0; j < G1H_FINALDEG-3; j++ ) {
        fn = k*G1_DBDIM+i*(G1H_FINALDEG-3)+j;
        pcp[(i+2)*(G1H_FINALDEG+1)+j+2] -= acoeff[fn];
      }
    if ( !mbs_BCHornerDer2Pf ( G1H_FINALDEG, G1H_FINALDEG, 2, &dicp[0].x, 0.0, knot,
                               &ei.x, &eiu.x, &eiv.x, &eiuu.x, &eiuv.x, &eivv.x ) )
      goto failure;
    if ( !mbs_BCHornerDer2Pf ( n, m, 1, pcp, 0.0, knot,
                               &e, &eu, &ev, &euu, &euv, &evv ) )
      goto failure;
    if ( !pkn_Comp2iDerivatives2f ( eiu.x, eiu.y, eiv.x, eiv.y, eiuu.x, eiuu.y,
                                    eiuv.x, eiuv.y, eivv.x, eivv.y,
                                    1, &eu, &ev, &euu, &euv, &evv,
                                    &qu, &qv, &puu, &puv, &pvv ) )
      goto failure;
    *jpuu -= puu;
    *jpuv -= puv;
    *jpvv -= pvv;
    break;

case 1:
    if ( !mbs_BCHornerDer2Pf ( G1H_FINALDEG, G1H_FINALDEG, 2, &dicp[0].x, knot, 1.0,
                               &di.x, &diu.x, &div.x, &diuu.x, &diuv.x, &divv.x ) )
      goto failure;
    *tang = diu;
    if ( !mbs_BCHornerDer2Pf ( n, m, 1, pcp, knot, 1.0,
                               &f, &fu, &fv, &fuu, &fuv, &fvv ) )
      goto failure;
    if ( !pkn_Comp2iDerivatives2f ( diu.x, diu.y, div.x, div.y, diuu.x, diuu.y,
                                    diuv.x, diuv.y, divv.x, divv.y,
                                    1, &fu, &fv, &fuu, &fuv, &fvv,
                                    pu, pv, &puu, &puv, &pvv ) )
      goto failure;
    *jpuu = -puu;
    *jpuv = -puv;
    *jpvv = -pvv;
    if ( !_gh_FindDomSurrndPatchf ( domain, (k+1) % hole_k, 1, dicp ) )
      goto failure;
    if ( !mbs_BCHornerDer2Pf ( 3, 3, 2, &dicp[0].x, 0.0, knot, &ei.x, &eiu.x, &eiv.x,
                               &eiuu.x, &eiuv.x, &eivv.x ) )
      goto failure;
    memset ( pcp, 0, 16*sizeof(float) );
    for ( fn = 0; fn < nfunc_b; fn++ ) {
      bvz = nlpr->rhole_cp[bfcpn[fn]].z;
      gh_GetDomSurrndBFuncf ( domain, fn, (k+1) % hole_k, 1, pc00 );
      pkn_AddMatrixMf ( 1, 16, 0, pcp, 0, pc00, bvz, 0, pcp );
    }
    if ( !mbs_BCHornerDer2Pf ( 3, 3, 1, pcp, 0.0, knot,
                               &e, &eu, &ev, &euu, &euv, &evv ) )
      goto failure;
    if ( !pkn_Comp2iDerivatives2f ( eiu.x, eiu.y, eiv.x, eiv.y, eiuu.x, eiuu.y,
                                    eiuv.x, eiuv.y, eivv.x, eivv.y,
                                    1, &eu, &ev, &euu, &euv, &evv,
                                    &qu, &qv, &puu, &puv, &pvv ) )
      goto failure;
    *jpuu += puu;
    *jpuv += puv;
    *jpvv += pvv;
    break;

case 2:
    if ( !mbs_BCHornerDer2Pf ( G1H_FINALDEG, G1H_FINALDEG, 2, &dicp[0].x, 1.0, knot,
                               &di.x, &diu.x, &div.x, &diuu.x, &diuv.x, &divv.x ) )
      goto failure;
    *tang = div;
    if ( !mbs_BCHornerDer2Pf ( n, m, 1, pcp, 1.0, knot,
                               &f, &fu, &fv, &fuu, &fuv, &fvv ) )
      goto failure;
    if ( !pkn_Comp2iDerivatives2f ( diu.x, diu.y, div.x, div.y, diuu.x, diuu.y,
                                    diuv.x, diuv.y, divv.x, divv.y,
                                    1, &fu, &fv, &fuu, &fuv, &fvv,
                                    pu, pv, &puu, &puv, &pvv ) )
      goto failure;
    *jpuu = -puu;
    *jpuv = -puv;
    *jpvv = -pvv;
    if ( !_gh_FindDomSurrndPatchf ( domain, k, 2, dicp ) )
      goto failure;
    if ( !mbs_BCHornerDer2Pf ( 3, 3, 2, &dicp[0].x, 0.0, knot, &ei.x, &eiu.x, &eiv.x,
                               &eiuu.x, &eiuv.x, &eivv.x ) )
      goto failure;
    memset ( pcp, 0, 16*sizeof(float) );
    for ( fn = 0; fn < nfunc_b; fn++ ) {
      bvz = nlpr->rhole_cp[bfcpn[fn]].z;
      gh_GetDomSurrndBFuncf ( domain, fn, k, 2, pc00 );
      pkn_AddMatrixMf ( 1, 16, 0, pcp, 0, pc00, bvz, 0, pcp );
    }
    if ( !mbs_BCHornerDer2Pf ( 3, 3, 1, pcp, 0.0, knot,
                               &e, &eu, &ev, &euu, &euv, &evv ) )
      goto failure;
    if ( !pkn_Comp2iDerivatives2f ( eiu.x, eiu.y, eiv.x, eiv.y, eiuu.x, eiuu.y,
                                    eiuv.x, eiuv.y, eivv.x, eivv.y,
                                    1, &eu, &ev, &euu, &euv, &evv,
                                    &qu, &qv, &puu, &puv, &pvv ) )
      goto failure;
    *jpuu += puu;
    *jpuv += puv;
    *jpvv += pvv;
    break;

default:
    goto failure;
  }

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
#undef N
} /*ComputePDer*/

static void DComputePDer ( GHoleDomainf *domain, G1HNLPrivatef *nlpr,
                           const float *acoeff, int knotno, G1Q2HNLFuncf *cf )
{
  G1HolePrivateRecf *privateG1;
  G1Q2HNLFuncf      acf;
  int   hole_k, k, cn;
  float p, tkn[G1_NQUAD];
  
  privateG1 = domain->privateG1;
  _gh_PrepareTabKnotsf ( G1_NQUAD, privateG1->opt_quad, tkn );
  hole_k = domain->hole_k;
  k = knotno / (3*G1_NQUAD);
  knotno %= 3*G1_NQUAD;
  cn = knotno / G1_NQUAD;
  knotno %= G1_NQUAD;
  ComputePDer ( domain, nlpr, acoeff, k, cn, tkn[knotno],
                &p, &acf.pu, &acf.pv, &acf.jpuu, &acf.jpuv, &acf.jpvv, &acf.tang );
} /*DComputePDer*/

/* ////////////////////////////////////////////////////////////////////////// */
#endif

static G1HNLPrivatef *_g1hq2_InitExtNLprf ( int hole_k, int nfunc_a,
                                            int nconstr )
{
#define N (G1H_FINALDEG+1)*(G1H_FINALDEG+1)
  G1HNLPrivatef *nlpr;
  int  nfunc_c;

  if ( (nlpr = pkv_GetScratchMem ( sizeof(G1HNLPrivatef) )) ) {
    nfunc_c = hole_k*G1_DBDIM;
    nlpr->auxc = 0;
    nlpr->nldi = pkv_GetScratchMem ( N*hole_k*sizeof(point3f) );
    nlpr->acoeff = pkv_GetScratchMem ( (nfunc_a+nfunc_c)*sizeof(vector3f) );
    nlpr->rhole_cp = pkv_GetScratchMem ( (12*hole_k+1+nconstr)*sizeof(vector3f) );
    nlpr->diu = pkv_GetScratchMem ( 9*G1_NQUADSQ*hole_k*sizeof(vector2f) );
    nlpr->jac = pkv_GetScratchMemf ( G1_NQUADSQ*hole_k );
    nlpr->psiu = pkv_GetScratchMemf ( 9*((nfunc_a+1)*hole_k+nfunc_c)*G1_NQUADSQ );
    nlpr->ctang = pkv_GetScratchMem ( 3*G1_NQUAD*hole_k*sizeof(vector2f) );
    nlpr->cpsiu = pkv_GetScratchMemf ( 5*(3*(nfunc_a+1)*hole_k+4*nfunc_c)*G1_NQUAD );

    if ( !nlpr->nldi || !nlpr->acoeff || !nlpr->rhole_cp ||
         !nlpr->diu || !nlpr->jac || !nlpr->psiu ||
         !nlpr->ctang || !nlpr->cpsiu )
      return NULL;

    nlpr->div = &nlpr->diu[G1_NQUADSQ*hole_k];
    nlpr->diuu = &nlpr->div[G1_NQUADSQ*hole_k];
    nlpr->diuv = &nlpr->diuu[G1_NQUADSQ*hole_k];
    nlpr->divv = &nlpr->diuv[G1_NQUADSQ*hole_k];
    nlpr->diuuu = &nlpr->divv[G1_NQUADSQ*hole_k];
    nlpr->diuuv = &nlpr->diuuu[G1_NQUADSQ*hole_k];
    nlpr->diuvv = &nlpr->diuuv[G1_NQUADSQ*hole_k];
    nlpr->divvv = &nlpr->diuvv[G1_NQUADSQ*hole_k];
    nlpr->psiv = &nlpr->psiu[((nfunc_a+1)*hole_k+nfunc_c)*G1_NQUADSQ];
    nlpr->psiuu = &nlpr->psiv[((nfunc_a+1)*hole_k+nfunc_c)*G1_NQUADSQ];
    nlpr->psiuv = &nlpr->psiuu[((nfunc_a+1)*hole_k+nfunc_c)*G1_NQUADSQ];
    nlpr->psivv = &nlpr->psiuv[((nfunc_a+1)*hole_k+nfunc_c)*G1_NQUADSQ];
    nlpr->psiuuu = &nlpr->psivv[((nfunc_a+1)*hole_k+nfunc_c)*G1_NQUADSQ];
    nlpr->psiuuv = &nlpr->psiuuu[((nfunc_a+1)*hole_k+nfunc_c)*G1_NQUADSQ];
    nlpr->psiuvv = &nlpr->psiuuv[((nfunc_a+1)*hole_k+nfunc_c)*G1_NQUADSQ];
    nlpr->psivvv = &nlpr->psiuvv[((nfunc_a+1)*hole_k+nfunc_c)*G1_NQUADSQ];
    nlpr->cpsiv = &nlpr->cpsiu[(3*(nfunc_a+1)*hole_k+4*nfunc_c)*G1_NQUAD];
    nlpr->cpsiuu = &nlpr->cpsiv[(3*(nfunc_a+1)*hole_k+4*nfunc_c)*G1_NQUAD];
    nlpr->cpsiuv = &nlpr->cpsiuu[(3*(nfunc_a+1)*hole_k+4*nfunc_c)*G1_NQUAD];
    nlpr->cpsivv = &nlpr->cpsiuv[(3*(nfunc_a+1)*hole_k+4*nfunc_c)*G1_NQUAD];
  }
  return nlpr;
#undef N
} /*_g1hq2_InitExtNLprf*/

static boolean _g1hq2_TabExtNLBasisFunctionsf ( GHoleDomainf *domain,
                                                G1HNLPrivatef *nlpr )
{
#define N ((G1H_FINALDEG+1)*(G1H_FINALDEG+1))
  void     *sp;
  G1HolePrivateRecf *privateG1;
  int      hole_k, nfunc_a, nfunc_c, i, j, k, kk, l, f, fN, bN, s;
  float    *tkn, *tbez, *tbezu, *tbezv, *tbezuu, *tbezuv, *tbezvv,
           *tbezuuu, *tbezuuv, *tbezuvv, *tbezvvv;
  float    *psiu, *psiv, *psiuu, *psiuv, *psivv,
           *psiuuu, *psiuuv, *psiuvv, *psivvv;
  float    *b, *ctr, *ctrd, *ctrdd, *A11, *A21, *A22, *bc00;
  vector2f *diu, *div, *diuu, *diuv, *divv, *diuuu, *diuuv, *diuvv, *divvv;

  sp = pkv_GetScratchMemTop ();
  hole_k = domain->hole_k;
  bc00 = pkv_GetScratchMemf ( 2*(G1_CROSSDEGSUM+4)*hole_k );
  if ( !bc00 ) {
    domain->error_code = G1H_ERROR_NO_SCRATCH_MEMORY;
    goto failure;
  }
  if ( !_g1hq2_TabNLBasisFunctionsOmegaf ( domain, G1_NQUAD, nlpr, bc00 ) )
    goto failure;

  privateG1 = domain->privateG1;
  nfunc_a = privateG1->nfunc_a;
  nfunc_c = hole_k*G1_DBDIM;
  tkn = pkv_GetScratchMemf ( G1_NQUAD );
  tbez = pkv_GetScratchMemf ( 10*G1_NQUADSQ*G1_DBDIM );
  if ( !tkn || !tbez ) {
    domain->error_code = G1H_ERROR_NO_SCRATCH_MEMORY;
    goto failure;
  }
  tbezu = &tbez[G1_NQUADSQ*G1_DBDIM];       tbezv = &tbezu[G1_NQUADSQ*G1_DBDIM];
  tbezuu = &tbezv[G1_NQUADSQ*G1_DBDIM];     tbezuv = &tbezuu[G1_NQUADSQ*G1_DBDIM];
  tbezvv = &tbezuv[G1_NQUADSQ*G1_DBDIM];    tbezuuu = &tbezvv[G1_NQUADSQ*G1_DBDIM];
  tbezuuv = &tbezuuu[G1_NQUADSQ*G1_DBDIM];  tbezuvv = &tbezuuv[G1_NQUADSQ*G1_DBDIM];
  tbezvvv = &tbezuvv[G1_NQUADSQ*G1_DBDIM];
  _gh_PrepareTabKnotsf ( G1_NQUAD, privateG1->opt_quad, tkn );

        /* evaluate the derivatives of the C block functions in Omega_i */
  if ( !_g1h_TabTensBezPolyDer3f ( G1_NQUAD, tkn, tbez, tbezu, tbezv, tbezuu,
                                   tbezuv, tbezvv, tbezuuu, tbezuuv,
                                   tbezuvv, tbezvvv ) )
    goto failure;
  for ( k = f = 0;  k < hole_k;  k++ ) {
    diu = &nlpr->diu[k*G1_NQUADSQ];      div = &nlpr->div[k*G1_NQUADSQ];
    diuu = &nlpr->diuu[k*G1_NQUADSQ];    diuv = &nlpr->diuv[k*G1_NQUADSQ];
    divv = &nlpr->divv[k*G1_NQUADSQ];    diuuu = &nlpr->diuuu[k*G1_NQUADSQ];
    diuuv = &nlpr->diuuv[k*G1_NQUADSQ];  diuvv = &nlpr->diuvv[k*G1_NQUADSQ];
    divvv = &nlpr->divvv[k*G1_NQUADSQ];
    for ( i = bN = 0;  i < G1H_FINALDEG-3; i++ )
      for ( j = 0;  j < G1H_FINALDEG-3;  j++, f++, bN += G1_NQUADSQ ) {
        fN = ((nfunc_a+1)*hole_k+f)*G1_NQUADSQ;
        psiu = &nlpr->psiu[fN];      psiv = &nlpr->psiv[fN];
        psiuu = &nlpr->psiuu[fN];    psiuv = &nlpr->psiuv[fN];
        psivv = &nlpr->psivv[fN];
        psiuuu = &nlpr->psiuuu[fN];  psiuuv = &nlpr->psiuuv[fN];
        psiuvv = &nlpr->psiuvv[fN];  psivvv = &nlpr->psivvv[fN];
        for ( l = 0; l < G1_NQUADSQ; l++ )
          if ( !pkn_Comp2iDerivatives3f ( diu[l].x, diu[l].y, div[l].x, div[l].y,
                  diuu[l].x, diuu[l].y, diuv[l].x, diuv[l].y,
                  divv[l].x, divv[l].y, diuuu[l].x, diuuu[l].y,
                  diuuv[l].x, diuuv[l].y, diuvv[l].x, diuvv[l].y,
                  divvv[l].x, divvv[l].y, 1, &tbezu[bN+l], &tbezv[bN+l],
                  &tbezuu[bN+l], &tbezuv[bN+l], &tbezvv[bN+l],
                  &tbezuuu[bN+l], &tbezuuv[bN+l], &tbezuvv[bN+l], &tbezvvv[bN+l],
                  &psiu[l], &psiv[l], &psiuu[l], &psiuv[l], &psivv[l],
                  &psiuuu[l], &psiuuv[l], &psiuvv[l], &psivvv[l] ) )
          goto failure;
      }
  }
  if ( !(ctr = pkv_GetScratchMemf ( 3*G1_NQUAD*38*hole_k)) ) {
    domain->error_code = G1H_ERROR_NO_SCRATCH_MEMORY;
    goto failure;
  }
  if ( !_g1hq2_TabNLBasisFunctionsGammaf ( domain, G1_NQUAD, nlpr, ctr, bc00 ) )
    goto failure;
        /* evaluate the derivatives and jumps of the C block functions */
        /* at the boundary curves */
  psiu = &nlpr->cpsiu[(nfunc_a+1)*hole_k*3*G1_NQUAD];
  psiv = &nlpr->cpsiv[(nfunc_a+1)*hole_k*3*G1_NQUAD];
  psiuu = &nlpr->cpsiuu[(nfunc_a+1)*hole_k*3*G1_NQUAD];
  psiuv = &nlpr->cpsiuv[(nfunc_a+1)*hole_k*3*G1_NQUAD];
  psivv = &nlpr->cpsivv[(nfunc_a+1)*hole_k*3*G1_NQUAD];
  s = nfunc_c*4*G1_NQUAD;
  memset ( psiu, 0, s*sizeof(float) );
  memset ( psiv, 0, s*sizeof(float) );
  memset ( psiuu, 0, s*sizeof(float) );
  memset ( psiuv, 0, s*sizeof(float) );
  memset ( psivv, 0, s*sizeof(float) );
  b = pkv_GetScratchMemf ( N );
  if ( !b ) {
    domain->error_code = G1H_ERROR_NO_SCRATCH_MEMORY;
    goto failure;
  }
  memset ( b, 0, N*sizeof(float) );
  for ( k = f = bN = 0, kk = 1;  k < hole_k;  k++, kk = (k+1) % hole_k ) {
    ctrd = &ctr[38*k*3*G1_NQUAD];
    ctrdd = &ctr[38*kk*3*G1_NQUAD+19];
    for ( i = 0; i < G1H_FINALDEG-3; i++ )
      for ( j = 0;  j < G1H_FINALDEG-3;  j++, f++, bN += 4*G1_NQUAD ) {
        b[(i+2)*(G1H_FINALDEG+1)+j+2] = 1.0;
        if ( i == 0 ) {
          for ( l = 0; l < G1_NQUAD; l++ ) {
            A11 = &ctrdd[38*l];  A21 = &A11[4];  A22 = &A21[6];
            if ( !mbs_BCHornerDer2Pf ( G1H_FINALDEG, G1H_FINALDEG, 1, b, 0.0, tkn[l],
                             &tbezu[1], tbezu, tbezv, tbezuu, tbezuv, tbezvv ) )
              goto failure;
            psiuu[bN+3*G1_NQUAD+l] = -(A21[0]*tbezu[0] + A21[1]*tbezv[0] +
                          A22[0]*tbezuu[0] + A22[1]*tbezuv[0] + A22[2]*tbezvv[0]);
            psiuv[bN+3*G1_NQUAD+l] = -(A21[2]*tbezu[0] + A21[3]*tbezv[0] +
                          A22[3]*tbezuu[0] + A22[4]*tbezuv[0] + A22[5]*tbezvv[0]);
            psivv[bN+3*G1_NQUAD+l] = -(A21[4]*tbezu[0] + A21[5]*tbezv[0] +
                          A22[6]*tbezuu[0] + A22[7]*tbezuv[0] + A22[8]*tbezvv[0]);
          }
        }
        else if ( i == G1H_FINALDEG-4 ) {
          for ( l = 0; l < G1_NQUAD; l++ ) {
            A11 = &ctrd[38*(2*G1_NQUAD+l)];  A21 = &A11[4];  A22 = &A21[6];
            if ( !mbs_BCHornerDer2Pf ( G1H_FINALDEG, G1H_FINALDEG, 1, b, 1.0, tkn[l],
                             &tbezu[1], tbezu, tbezv, tbezuu, tbezuv, tbezvv ) )
              goto failure;
            psiuu[bN+2*G1_NQUAD+l] = -(A21[0]*tbezu[0] + A21[1]*tbezv[0] +
                          A22[0]*tbezuu[0] + A22[1]*tbezuv[0] + A22[2]*tbezvv[0]);
            psiuv[bN+2*G1_NQUAD+l] = -(A21[2]*tbezu[0] + A21[3]*tbezv[0] +
                          A22[3]*tbezuu[0] + A22[4]*tbezuv[0] + A22[5]*tbezvv[0]);
            psivv[bN+2*G1_NQUAD+l] = -(A21[4]*tbezu[0] + A21[5]*tbezv[0] +
                          A22[6]*tbezuu[0] + A22[7]*tbezuv[0] + A22[8]*tbezvv[0]);
          }
        }
        if ( j == 0 ) {
          for ( l = 0; l < G1_NQUAD; l++  ) {
            A11 = &ctrd[38*l];  A21 = &A11[4];  A22 = &A21[6];
            if ( !mbs_BCHornerDer2Pf ( G1H_FINALDEG, G1H_FINALDEG, 1, b, tkn[l], 0.0,
                             &tbezu[1], tbezu, tbezv, tbezuu, tbezuv, tbezvv ) )
              goto failure;
            psiuu[bN+l] = A21[0]*tbezu[0] + A21[1]*tbezv[0] +
                          A22[0]*tbezuu[0] + A22[1]*tbezuv[0] + A22[2]*tbezvv[0];
            psiuv[bN+l] = A21[2]*tbezu[0] + A21[3]*tbezv[0] +
                          A22[3]*tbezuu[0] + A22[4]*tbezuv[0] + A22[5]*tbezvv[0];
            psivv[bN+l] = A21[4]*tbezu[0] + A21[5]*tbezv[0] +
                          A22[6]*tbezuu[0] + A22[7]*tbezuv[0] + A22[8]*tbezvv[0];
          }
        }
        else if ( j == G1H_FINALDEG-4 ) {
          for ( l = 0; l < G1_NQUAD; l++ ) {
            A11 = &ctrd[38*(G1_NQUAD+l)];  A21 = &A11[4];  A22 = &A21[6];
            if ( !mbs_BCHornerDer2Pf ( G1H_FINALDEG, G1H_FINALDEG, 1, b, tkn[l], 1.0,
                             &tbezu[1], tbezu, tbezv, tbezuu, tbezuv, tbezvv ) )
              goto failure;
            psiuu[bN+G1_NQUAD+l] = -(A21[0]*tbezu[0] + A21[1]*tbezv[0] +
                          A22[0]*tbezuu[0] + A22[1]*tbezuv[0] + A22[2]*tbezvv[0]);
            psiuv[bN+G1_NQUAD+l] = -(A21[2]*tbezu[0] + A21[3]*tbezv[0] +
                          A22[3]*tbezuu[0] + A22[4]*tbezuv[0] + A22[5]*tbezvv[0]);
            psivv[bN+G1_NQUAD+l] = -(A21[4]*tbezu[0] + A21[5]*tbezv[0] +
                          A22[6]*tbezuu[0] + A22[7]*tbezuv[0] + A22[8]*tbezvv[0]);
          }
        }
        b[(i+2)*(G1H_FINALDEG+1)+j+2] = 0.0;
      }
  }

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
#undef N
} /*_g1hq2_TabExtNLBasisFunctionsf*/

/* ///////////////////////////////////////////////////////////////////////// */
static float _g1hq2_ComputeExtNLFuncf ( GHoleDomainf *domain,
                                        G1HNLPrivatef *nlprivate,
                                        const float *coeff,
                                        float C )
{
  G1HolePrivateRecf *privateG1;
  G2HNLFuncf   f;
  G1Q2HNLFuncf cf;
  int   hole_k, nfunc_a, nfunc_c;
  int   i, ii, k, ki, kn, knot, fi;
  float c, funct, cfunct;

  privateG1 = domain->privateG1;
  hole_k = domain->hole_k;
  nfunc_a = privateG1->nfunc_a;
  nfunc_c = hole_k*G1_DBDIM;

  funct = cfunct = 0.0;
        /* integration in the area Omega */
  for ( k = knot = 0;  k < hole_k;  k++ ) {
    for ( kn = 0;  kn < G1_NQUADSQ;  kn++, knot++ ) {
          /* evaluate the function and its derivatives */
      fi = nfunc_a*G1_NQUADSQ*hole_k + knot;
      f.pu = nlprivate->psiu[fi];      f.pv = nlprivate->psiv[fi];
      f.puu = nlprivate->psiuu[fi];    f.puv = nlprivate->psiuv[fi];
      f.pvv = nlprivate->psivv[fi];    f.puuu = nlprivate->psiuuu[fi];
      f.puuv = nlprivate->psiuuv[fi];  f.puvv = nlprivate->psiuvv[fi];
      f.pvvv = nlprivate->psivvv[fi];
      for ( i = 0; i < nfunc_a; i++ ) {
        fi = i*G1_NQUADSQ*hole_k + knot;
        c = coeff[nfunc_c+i];
        f.pu -= c*nlprivate->psiu[fi];      f.pv -= c*nlprivate->psiv[fi];
        f.puu -= c*nlprivate->psiuu[fi];    f.puv -= c*nlprivate->psiuv[fi];
        f.pvv -= c*nlprivate->psivv[fi];    f.puuu -= c*nlprivate->psiuuu[fi];
        f.puuv -= c*nlprivate->psiuuv[fi];  f.puvv -= c*nlprivate->psiuvv[fi];
        f.pvvv -= c*nlprivate->psivvv[fi];
      }
      for ( i = 0; i < G1_DBDIM; i++ ) {
        fi = ((nfunc_a+1)*hole_k+k*G1_DBDIM+i)*G1_NQUADSQ+kn;
        c = coeff[k*G1_DBDIM+i];
        f.pu -= c*nlprivate->psiu[fi];      f.pv -= c*nlprivate->psiv[fi];
        f.puu -= c*nlprivate->psiuu[fi];    f.puv -= c*nlprivate->psiuv[fi];
        f.pvv -= c*nlprivate->psivv[fi];    f.puuu -= c*nlprivate->psiuuu[fi];
        f.puuv -= c*nlprivate->psiuuv[fi];  f.puvv -= c*nlprivate->psiuvv[fi];
        f.pvvv -= c*nlprivate->psivvv[fi];
      }
      f.jac = nlprivate->jac[knot];

        /* integrate the functional */
      _g2h_IntFunc1af ( &f, &funct );
    }
  }
  funct /= (float)G1_NQUADSQ;

        /* integration along the curves */
  for ( k = knot = 0, ki = hole_k-1;
        k < hole_k;
        ki = k++ ) {
    for ( kn = 0;  kn < G1_NQUAD*3;  kn++, knot++ ) {
          /* evaluate the function and its derivatives */
      fi = nfunc_a*3*G1_NQUAD*hole_k + knot;
      cf.pu = nlprivate->cpsiu[fi];     cf.pv = nlprivate->cpsiv[fi];
      cf.jpuu = nlprivate->cpsiuu[fi];  cf.jpuv = nlprivate->cpsiuv[fi];
      cf.jpvv = nlprivate->cpsivv[fi];
      for ( i = 0; i < nfunc_a; i++ ) {
        fi = i*3*G1_NQUAD*hole_k + knot;
        c = coeff[nfunc_c+i];
        cf.pu -= c*nlprivate->cpsiu[fi];     cf.pv -= c*nlprivate->cpsiv[fi];
        cf.jpuu -= c*nlprivate->cpsiuu[fi];  cf.jpuv -= c*nlprivate->cpsiuv[fi];
        cf.jpvv -= c*nlprivate->cpsivv[fi];
      }
      for ( i = 0; i < G1_DBDIM; i++ ) {
        ii = k*G1_DBDIM+i;
        fi = ((nfunc_a+1)*3*hole_k + ii*4)*G1_NQUAD + kn;
        c = coeff[ii];
        cf.pu -= c*nlprivate->cpsiu[fi];     cf.pv -= c*nlprivate->cpsiv[fi];
        cf.jpuu -= c*nlprivate->cpsiuu[fi];  cf.jpuv -= c*nlprivate->cpsiuv[fi];
        cf.jpvv -= c*nlprivate->cpsivv[fi];
      }
      if ( kn < G1_NQUAD )
        for ( i = 0; i < G1_DBDIM; i++ ) {
          ii = ki*G1_DBDIM+i;
          fi = ((nfunc_a+1)*3*hole_k + ii*4 + 3)*G1_NQUAD + kn;
          c = coeff[ii];
          cf.pu -= c*nlprivate->cpsiu[fi];     cf.pv -= c*nlprivate->cpsiv[fi];
          cf.jpuu -= c*nlprivate->cpsiuu[fi];  cf.jpuv -= c*nlprivate->cpsiuv[fi];
          cf.jpvv -= c*nlprivate->cpsivv[fi];
        }
      cf.tang = nlprivate->ctang[knot];

        /* integrate the functional */
      _g1hq2_IntFunc1af ( &cf, &cfunct );
    }
  }
  cfunct /= (float)G1_NQUAD;
  return funct + C*cfunct;
} /*_g1hq2_ComputeExtNLFuncf*/

static boolean _g1hq2_ComputeExtNLFuncGradf ( GHoleDomainf *domain,
                                              G1HNLPrivatef *nlprivate,
                                              const float *coeff,
                                              float C,
                                              float *func, float *grad )
{
  void  *sp;
  G1HolePrivateRecf *privateG1;
  G2HNLFuncf   f;
  G1Q2HNLFuncf cf;
  int   hole_k, nfunc_a, nfunc_c, nfunc;
  int   i, ii, k, ki, kn, knot, fi;
  float c, funct, cfunct, *cgrad;

  sp = pkv_GetScratchMemTop ();
  privateG1 = domain->privateG1;
  hole_k = domain->hole_k;
  nfunc_a = privateG1->nfunc_a;
  nfunc_c = hole_k*G1_DBDIM;
  nfunc = nfunc_a+nfunc_c;
  cgrad = pkv_GetScratchMemf ( nfunc );
  if ( !cgrad ) {
    domain->error_code = G1H_ERROR_NO_SCRATCH_MEMORY;
    goto failure;
  }

  funct = cfunct = 0.0;
  memset ( grad, 0, nfunc*sizeof(float) );
  memset ( cgrad, 0, nfunc*sizeof(float) );
        /* integration in the area Omega */
  for ( k = knot = 0;  k < hole_k;  k++ ) {
    for ( kn = 0;  kn < G1_NQUADSQ;  kn++, knot++ ) {
          /* evaluate the function and its derivatives */
      fi = nfunc_a*G1_NQUADSQ*hole_k + knot;
      f.pu = nlprivate->psiu[fi];      f.pv = nlprivate->psiv[fi];
      f.puu = nlprivate->psiuu[fi];    f.puv = nlprivate->psiuv[fi];
      f.pvv = nlprivate->psivv[fi];    f.puuu = nlprivate->psiuuu[fi];
      f.puuv = nlprivate->psiuuv[fi];  f.puvv = nlprivate->psiuvv[fi];
      f.pvvv = nlprivate->psivvv[fi];
      for ( i = 0; i < nfunc_a; i++ ) {
        fi = i*G1_NQUADSQ*hole_k + knot;
        c = coeff[nfunc_c+i];
        f.pu -= c*nlprivate->psiu[fi];      f.pv -= c*nlprivate->psiv[fi];
        f.puu -= c*nlprivate->psiuu[fi];    f.puv -= c*nlprivate->psiuv[fi];
        f.pvv -= c*nlprivate->psivv[fi];    f.puuu -= c*nlprivate->psiuuu[fi];
        f.puuv -= c*nlprivate->psiuuv[fi];  f.puvv -= c*nlprivate->psiuvv[fi];
        f.pvvv -= c*nlprivate->psivvv[fi];
      }
      for ( i = 0; i < G1_DBDIM; i++ ) {
        fi = ((nfunc_a+1)*hole_k+k*G1_DBDIM+i)*G1_NQUADSQ+kn;
        c = coeff[k*G1_DBDIM+i];
        f.pu -= c*nlprivate->psiu[fi];      f.pv -= c*nlprivate->psiv[fi];
        f.puu -= c*nlprivate->psiuu[fi];    f.puv -= c*nlprivate->psiuv[fi];
        f.pvv -= c*nlprivate->psivv[fi];    f.puuu -= c*nlprivate->psiuuu[fi];
        f.puuv -= c*nlprivate->psiuuv[fi];  f.puvv -= c*nlprivate->psiuvv[fi];
        f.pvvv -= c*nlprivate->psivvv[fi];
      }
      f.jac = nlprivate->jac[knot];

        /* integrate the functional */
      _g2h_IntFunc1bf ( &f, &funct );

        /* integrate the gradient */
      for ( i = 0; i < nfunc_a; i++ ) {
        fi = i*G1_NQUADSQ*hole_k+knot;
        f.psiu = nlprivate->psiu[fi];      f.psiv = nlprivate->psiv[fi];
        f.psiuu = nlprivate->psiuu[fi];    f.psiuv = nlprivate->psiuv[fi];
        f.psivv = nlprivate->psivv[fi];    f.psiuuu = nlprivate->psiuuu[fi];
        f.psiuuv = nlprivate->psiuuv[fi];  f.psiuvv = nlprivate->psiuvv[fi];
        f.psivvv = nlprivate->psivvv[fi];
        _g2h_IntFunc2bf ( &f, &grad[nfunc_c+i] );
      }
      for ( i = 0; i < G1_DBDIM; i++ ) {
        fi = ((nfunc_a+1)*hole_k+k*G1_DBDIM+i)*G1_NQUADSQ+kn;
        f.psiu = nlprivate->psiu[fi];      f.psiv = nlprivate->psiv[fi];
        f.psiuu = nlprivate->psiuu[fi];    f.psiuv = nlprivate->psiuv[fi];
        f.psivv = nlprivate->psivv[fi];    f.psiuuu = nlprivate->psiuuu[fi];
        f.psiuuv = nlprivate->psiuuv[fi];  f.psiuvv = nlprivate->psiuvv[fi];
        f.psivvv = nlprivate->psivvv[fi];
        _g2h_IntFunc2bf ( &f, &grad[k*G1_DBDIM+i] );
      }
    }
  }
  funct /= (float)G1_NQUADSQ;
  pkn_MultMatrixNumf ( 1, nfunc, 0, grad, 1.0/(float)G1_NQUADSQ, 0, grad );

        /* integration along the curves */
  for ( k = knot = 0, ki = hole_k-1;
        k < hole_k;
        ki = k++ ) {
    for ( kn = 0;  kn < G1_NQUAD*3;  kn++, knot++ ) {
          /* evaluate the function and its derivatives */
      fi = nfunc_a*3*G1_NQUAD*hole_k + knot;
      cf.pu = nlprivate->cpsiu[fi];     cf.pv = nlprivate->cpsiv[fi];
      cf.jpuu = nlprivate->cpsiuu[fi];  cf.jpuv = nlprivate->cpsiuv[fi];
      cf.jpvv = nlprivate->cpsivv[fi];
      for ( i = 0; i < nfunc_a; i++ ) {
        fi = i*3*G1_NQUAD*hole_k + knot;
        c = coeff[nfunc_c+i];
        cf.pu -= c*nlprivate->cpsiu[fi];     cf.pv -= c*nlprivate->cpsiv[fi];
        cf.jpuu -= c*nlprivate->cpsiuu[fi];  cf.jpuv -= c*nlprivate->cpsiuv[fi];
        cf.jpvv -= c*nlprivate->cpsivv[fi];
      }
      for ( i = 0; i < G1_DBDIM; i++ ) {
        ii = k*G1_DBDIM+i;
        fi = ((nfunc_a+1)*3*hole_k + ii*4)*G1_NQUAD + kn;
        c = coeff[ii];
        cf.pu -= c*nlprivate->cpsiu[fi];     cf.pv -= c*nlprivate->cpsiv[fi];
        cf.jpuu -= c*nlprivate->cpsiuu[fi];  cf.jpuv -= c*nlprivate->cpsiuv[fi];
        cf.jpvv -= c*nlprivate->cpsivv[fi];
      }
      if ( kn < G1_NQUAD )
        for ( i = 0; i < G1_DBDIM; i++ ) {
          ii = ki*G1_DBDIM+i;
          fi = ((nfunc_a+1)*3*hole_k + ii*4 + 3)*G1_NQUAD + kn;
          c = coeff[ii];
          cf.pu -= c*nlprivate->cpsiu[fi];     cf.pv -= c*nlprivate->cpsiv[fi];
          cf.jpuu -= c*nlprivate->cpsiuu[fi];  cf.jpuv -= c*nlprivate->cpsiuv[fi];
          cf.jpvv -= c*nlprivate->cpsivv[fi];
        }
      cf.tang = nlprivate->ctang[knot];

        /* integrate the functional */
      _g1hq2_IntFunc1bf ( &cf, &cfunct );

        /* integrate the gradient */
      for ( i = 0; i < nfunc_a; i++ ) {
        fi = i*3*G1_NQUAD*hole_k + knot;
        cf.psiu = nlprivate->cpsiu[fi];     cf.psiv = nlprivate->cpsiv[fi];
        cf.jpsiuu = nlprivate->cpsiuu[fi];  cf.jpsiuv = nlprivate->cpsiuv[fi];
        cf.jpsivv = nlprivate->cpsivv[fi];
        _g1hq2_IntFunc2bf ( &cf, &cgrad[nfunc_c+i] );
      }
      for ( i = 0; i < G1_DBDIM; i++ ) {
        ii = k*G1_DBDIM+i;
        fi = ((nfunc_a+1)*3*hole_k + ii*4)*G1_NQUAD + kn;
        cf.psiu = nlprivate->cpsiu[fi];     cf.psiv = nlprivate->cpsiv[fi];
        cf.jpsiuu = nlprivate->cpsiuu[fi];  cf.jpsiuv = nlprivate->cpsiuv[fi];
        cf.jpsivv = nlprivate->cpsivv[fi];
        _g1hq2_IntFunc2bf ( &cf, &cgrad[ii] );
      }
      if ( kn < G1_NQUAD )
        for ( i = 0; i < G1_DBDIM; i++ ) {
          ii = ki*G1_DBDIM+i;
          fi = ((nfunc_a+1)*3*hole_k + ii*4 + 3)*G1_NQUAD + kn;
          cf.psiu = nlprivate->cpsiu[fi];     cf.psiv = nlprivate->cpsiv[fi];
          cf.jpsiuu = nlprivate->cpsiuu[fi];  cf.jpsiuv = nlprivate->cpsiuv[fi];
          cf.jpsivv = nlprivate->cpsivv[fi];
          _g1hq2_IntFunc2bf ( &cf, &cgrad[ii] );
        }
    }
  }
  cfunct /= (float)G1_NQUAD;
  *func = funct+C*cfunct;
  pkn_AddMatrixMf ( 1, nfunc, 0, grad, 0, cgrad, C/(float)G1_NQUAD, 0, grad );

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*_g1hq2_ComputeExtNLFuncGradf*/

static boolean _g1hq2_ComputeExtNLFuncGradHessianf ( GHoleDomainf *domain,
                          G1HNLPrivatef *nlprivate,
                          const float *coeff,
                          float C,
                          float *func, float *grad, float *hessian,
                          float *chessian )
{
  void  *sp;
  G1HolePrivateRecf *privateG1;
  G2HNLFuncf   f;
  G1Q2HNLFuncf cf;
  int   hole_k, nfunc_a, nfunc_c, nfunc, asize;
  int   i, ii, j, jj, k, ki, kn, knot, fi, fj;
  float *Li, *Bi, *BiLT, *Di, *hii, *hij, *hjj, *hki, *hkj, *hkk, *hp;
  float c, funct, cfunct, *cgrad;

  sp = pkv_GetScratchMemTop ();
  privateG1 = domain->privateG1;
  hole_k = domain->hole_k;
  nfunc_a = privateG1->nfunc_a;
  nfunc_c = hole_k*G1_DBDIM;
  nfunc = nfunc_a+nfunc_c;
  asize = pkn_Block3ArraySize ( hole_k-1, G1_DBDIM, G1_DBDIM+nfunc_a );
  cgrad = pkv_GetScratchMemf ( nfunc );
  Li = pkv_GetScratchMemf ( 8*nfunc );
  if ( !cgrad || !Li ) {
    domain->error_code = G1H_ERROR_NO_SCRATCH_MEMORY;
    goto failure;
  }
  Bi = &Li[2*nfunc];  BiLT = &Bi[3*nfunc];  Di = &BiLT[2*nfunc];

  funct = cfunct = 0.0;
  memset ( grad, 0, nfunc*sizeof(float) );
  memset ( cgrad, 0, nfunc*sizeof(float) );
  memset ( hessian, 0, asize*sizeof(float) );
  memset ( chessian, 0, asize*sizeof(float) );
        /* integration in the area Omega */
  hkk = &hessian[pkn_Block3FindBlockPos (hole_k-1, G1_DBDIM, G1_DBDIM+nfunc_a,
                                         hole_k-1, hole_k-1)];
  fj = 0;
  for ( k = knot = 0;  k < hole_k;  k++ ) {
    hii = &hessian[pkn_Block3FindBlockPos(hole_k-1, G1_DBDIM, G1_DBDIM+nfunc_a,
                                          k, k )];
    hki = &hessian[pkn_Block3FindBlockPos(hole_k-1, G1_DBDIM, G1_DBDIM+nfunc_a,
                                          hole_k-1, k)];
    for ( kn = 0;  kn < G1_NQUADSQ;  kn++, knot++ ) {
          /* evaluate the function and its derivatives */
      fi = nfunc_a*G1_NQUADSQ*hole_k + knot;
      f.pu = nlprivate->psiu[fi];      f.pv = nlprivate->psiv[fi];
      f.puu = nlprivate->psiuu[fi];    f.puv = nlprivate->psiuv[fi];
      f.pvv = nlprivate->psivv[fi];    f.puuu = nlprivate->psiuuu[fi];
      f.puuv = nlprivate->psiuuv[fi];  f.puvv = nlprivate->psiuvv[fi];
      f.pvvv = nlprivate->psivvv[fi];
      for ( i = 0; i < nfunc_a; i++ ) {
        fi = i*G1_NQUADSQ*hole_k + knot;
        c = coeff[nfunc_c+i];
        f.pu -= c*nlprivate->psiu[fi];      f.pv -= c*nlprivate->psiv[fi];
        f.puu -= c*nlprivate->psiuu[fi];    f.puv -= c*nlprivate->psiuv[fi];
        f.pvv -= c*nlprivate->psivv[fi];    f.puuu -= c*nlprivate->psiuuu[fi];
        f.puuv -= c*nlprivate->psiuuv[fi];  f.puvv -= c*nlprivate->psiuvv[fi];
        f.pvvv -= c*nlprivate->psivvv[fi];
      }
      for ( i = 0; i < G1_DBDIM; i++ ) {
        fi = ((nfunc_a+1)*hole_k+k*G1_DBDIM+i)*G1_NQUADSQ+kn;
        c = coeff[k*G1_DBDIM+i];
        f.pu -= c*nlprivate->psiu[fi];      f.pv -= c*nlprivate->psiv[fi];
        f.puu -= c*nlprivate->psiuu[fi];    f.puv -= c*nlprivate->psiuv[fi];
        f.pvv -= c*nlprivate->psivv[fi];    f.puuu -= c*nlprivate->psiuuu[fi];
        f.puuv -= c*nlprivate->psiuuv[fi];  f.puvv -= c*nlprivate->psiuvv[fi];
        f.pvvv -= c*nlprivate->psivvv[fi];
      }
      f.jac = nlprivate->jac[knot];

        /* integrate the functional */
      _g2h_IntFunc1cf ( &f, &funct );

        /* integrate the gradient */
      for ( i = 0; i < nfunc_a; i++ ) {
        fi = i*G1_NQUADSQ*hole_k+knot;
        f.psiu = nlprivate->psiu[fi];      f.psiv = nlprivate->psiv[fi];
        f.psiuu = nlprivate->psiuu[fi];    f.psiuv = nlprivate->psiuv[fi];
        f.psivv = nlprivate->psivv[fi];    f.psiuuu = nlprivate->psiuuu[fi];
        f.psiuuv = nlprivate->psiuuv[fi];  f.psiuvv = nlprivate->psiuvv[fi];
        f.psivvv = nlprivate->psivvv[fi];
        _g2h_IntFunc2cf ( &f, (vector2f*)(&Li[2*(nfunc_c+i)]),
                          &Bi[3*(nfunc_c+i)], (vector2f*)(&BiLT[2*(nfunc_c+i)]),
                          &Di[nfunc_c+i], &grad[nfunc_c+i] );
      }
      for ( i = 0; i < G1_DBDIM; i++ ) {
        fi = ((nfunc_a+1)*hole_k+k*G1_DBDIM+i)*G1_NQUADSQ+kn;
        f.psiu = nlprivate->psiu[fi];      f.psiv = nlprivate->psiv[fi];
        f.psiuu = nlprivate->psiuu[fi];    f.psiuv = nlprivate->psiuv[fi];
        f.psivv = nlprivate->psivv[fi];    f.psiuuu = nlprivate->psiuuu[fi];
        f.psiuuv = nlprivate->psiuuv[fi];  f.psiuvv = nlprivate->psiuvv[fi];
        f.psivvv = nlprivate->psivvv[fi];
        _g2h_IntFunc2cf ( &f, (vector2f*)(&Li[2*(k*G1_DBDIM+i)]),
                          &Bi[3*(k*G1_DBDIM+i)], (vector2f*)(&BiLT[2*(k*G1_DBDIM+i)]),
                          &Di[k*G1_DBDIM+i], &grad[k*G1_DBDIM+i] );
      }

        /* integrate the Hessian */
      for ( i = 0; i < nfunc_a; i++ ) {
        fi = i*G1_NQUADSQ*hole_k+knot;
        f.psiu = nlprivate->psiu[fi];      f.psiv = nlprivate->psiv[fi];
        f.psiuu = nlprivate->psiuu[fi];    f.psiuv = nlprivate->psiuv[fi];
        f.psivv = nlprivate->psivv[fi];    f.psiuuu = nlprivate->psiuuu[fi];
        f.psiuuv = nlprivate->psiuuv[fi];  f.psiuvv = nlprivate->psiuvv[fi];
        f.psivvv = nlprivate->psivvv[fi];
        for ( j = 0; j <= i; j++ ) {
          fj = j*G1_NQUADSQ*hole_k+knot;
          f.psju = nlprivate->psiu[fj];      f.psjv = nlprivate->psiv[fj];
          f.psjuu = nlprivate->psiuu[fj];    f.psjuv = nlprivate->psiuv[fj];
          f.psjvv = nlprivate->psivv[fj];    f.psjuuu = nlprivate->psiuuu[fj];
          f.psjuuv = nlprivate->psiuuv[fj];  f.psjuvv = nlprivate->psiuvv[fj];
          f.psjvvv = nlprivate->psivvv[fj];
          _g2h_IntFunc3cf ( &f, (vector2f*)(&Li[2*(nfunc_c+i)]),
                                (vector2f*)(&Li[2*(nfunc_c+j)]),
                                (vector2f*)(&BiLT[2*(nfunc_c+i)]),
                                (vector2f*)(&BiLT[2*(nfunc_c+j)]),
                                Di[nfunc_c+i], Di[nfunc_c+j],
                                &hkk[pkn_SymMatIndex(G1_DBDIM+i,G1_DBDIM+j)] );
        }
        for ( j = 0; j < G1_DBDIM; j++ ) {
          fj = ((nfunc_a+1)*hole_k+k*G1_DBDIM+j)*G1_NQUADSQ+kn;
          f.psju = nlprivate->psiu[fj];      f.psjv = nlprivate->psiv[fj];
          f.psjuu = nlprivate->psiuu[fj];    f.psjuv = nlprivate->psiuv[fj];
          f.psjvv = nlprivate->psivv[fj];    f.psjuuu = nlprivate->psiuuu[fj];
          f.psjuuv = nlprivate->psiuuv[fj];  f.psjuvv = nlprivate->psiuvv[fj];
          f.psjvvv = nlprivate->psivvv[fj];
          if ( k < hole_k-1 )
            hp = &hki[(G1_DBDIM+i)*G1_DBDIM+j];
          else
            hp = &hkk[pkn_SymMatIndex(G1_DBDIM+i,j)];
          _g2h_IntFunc3cf ( &f, (vector2f*)(&Li[2*(nfunc_c+i)]),
                                (vector2f*)(&Li[2*(k*G1_DBDIM+j)]),
                                (vector2f*)(&BiLT[2*(nfunc_c+i)]),
                                (vector2f*)(&BiLT[2*(k*G1_DBDIM+j)]),
                                Di[nfunc_c+i], Di[k*G1_DBDIM+j],
                                hp );
        }
      }
      for ( i = 0; i < G1_DBDIM; i++ ) {
        fi = ((nfunc_a+1)*hole_k+k*G1_DBDIM+i)*G1_NQUADSQ+kn;
        f.psiu = nlprivate->psiu[fi];      f.psiv = nlprivate->psiv[fi];
        f.psiuu = nlprivate->psiuu[fi];    f.psiuv = nlprivate->psiuv[fi];
        f.psivv = nlprivate->psivv[fi];    f.psiuuu = nlprivate->psiuuu[fi];
        f.psiuuv = nlprivate->psiuuv[fi];  f.psiuvv = nlprivate->psiuvv[fi];
        f.psivvv = nlprivate->psivvv[fi];
        for ( j = 0; j <= i; j++ ) {
          fj = ((nfunc_a+1)*hole_k+k*G1_DBDIM+j)*G1_NQUADSQ+kn;
          f.psju = nlprivate->psiu[fj];      f.psjv = nlprivate->psiv[fj];
          f.psjuu = nlprivate->psiuu[fj];    f.psjuv = nlprivate->psiuv[fj];
          f.psjvv = nlprivate->psivv[fj];    f.psjuuu = nlprivate->psiuuu[fj];
          f.psjuuv = nlprivate->psiuuv[fj];  f.psjuvv = nlprivate->psiuvv[fj];
          f.psjvvv = nlprivate->psivvv[fj];
          _g2h_IntFunc3cf ( &f, (vector2f*)(&Li[2*(k*G1_DBDIM+i)]),
                                (vector2f*)(&Li[2*(k*G1_DBDIM+j)]),
                                (vector2f*)(&BiLT[2*(k*G1_DBDIM+i)]),
                                (vector2f*)(&BiLT[2*(k*G1_DBDIM+j)]),
                                Di[k*G1_DBDIM+i], Di[k*G1_DBDIM+j],
                                &hii[pkn_SymMatIndex(i,j)] );
        }
      }
    }
  }
  funct /= (float)G1_NQUADSQ;
  pkn_MultMatrixNumf ( 1, nfunc, 0, grad, 1.0/(float)G1_NQUADSQ, 0, grad );
  pkn_MultMatrixNumf ( 1, asize, 0, hessian, 1.0/(float)G1_NQUADSQ, 0, hessian );

        /* integration along the curves */
  hkk = &chessian[pkn_Block3FindBlockPos (hole_k-1, G1_DBDIM, G1_DBDIM+nfunc_a,
                                          hole_k-1, hole_k-1)];
  for ( k = knot = 0, ki = hole_k-1;
        k < hole_k;
        ki = k++ ) {
    hii = &chessian[pkn_Block3FindBlockPos(hole_k-1, G1_DBDIM, G1_DBDIM+nfunc_a,
                                           k, k )];
    hjj = &chessian[pkn_Block3FindBlockPos(hole_k-1, G1_DBDIM, G1_DBDIM+nfunc_a,
                                           ki, ki )];
    if ( k > 0 ) {
      hij = &chessian[pkn_Block3FindBlockPos(hole_k-1, G1_DBDIM, G1_DBDIM+nfunc_a,
                                           k, ki )];
      hkj = &chessian[pkn_Block3FindBlockPos(hole_k-1, G1_DBDIM, G1_DBDIM+nfunc_a,
                                           hole_k-1, ki )];
    }
    else {
      hij = hkj = &chessian[pkn_Block3FindBlockPos(hole_k-1,
                                  G1_DBDIM, G1_DBDIM+nfunc_a, hole_k-1, 0 )];
    }
    hki = &chessian[pkn_Block3FindBlockPos(hole_k-1, G1_DBDIM, G1_DBDIM+nfunc_a,
                                           hole_k-1, k)];
    for ( kn = 0;  kn < G1_NQUAD*3;  kn++, knot++ ) {
          /* evaluate the function and its derivatives */
      fi = nfunc_a*3*G1_NQUAD*hole_k + knot;
      cf.pu = nlprivate->cpsiu[fi];     cf.pv = nlprivate->cpsiv[fi];
      cf.jpuu = nlprivate->cpsiuu[fi];  cf.jpuv = nlprivate->cpsiuv[fi];
      cf.jpvv = nlprivate->cpsivv[fi];
      for ( i = 0; i < nfunc_a; i++ ) {
        fi = i*3*G1_NQUAD*hole_k + knot;
        c = coeff[nfunc_c+i];
        cf.pu -= c*nlprivate->cpsiu[fi];     cf.pv -= c*nlprivate->cpsiv[fi];
        cf.jpuu -= c*nlprivate->cpsiuu[fi];  cf.jpuv -= c*nlprivate->cpsiuv[fi];
        cf.jpvv -= c*nlprivate->cpsivv[fi];
      }
      for ( i = 0; i < G1_DBDIM; i++ ) {
        ii = k*G1_DBDIM+i;
        fi = ((nfunc_a+1)*3*hole_k + ii*4)*G1_NQUAD + kn;
        c = coeff[ii];
        cf.pu -= c*nlprivate->cpsiu[fi];     cf.pv -= c*nlprivate->cpsiv[fi];
        cf.jpuu -= c*nlprivate->cpsiuu[fi];  cf.jpuv -= c*nlprivate->cpsiuv[fi];
        cf.jpvv -= c*nlprivate->cpsivv[fi];
      }
      if ( kn < G1_NQUAD )
        for ( i = 0; i < G1_DBDIM; i++ ) {
          ii = ki*G1_DBDIM+i;
          fi = ((nfunc_a+1)*3*hole_k + ii*4 + 3)*G1_NQUAD + kn;
          c = coeff[ii];
          cf.pu -= c*nlprivate->cpsiu[fi];     cf.pv -= c*nlprivate->cpsiv[fi];
          cf.jpuu -= c*nlprivate->cpsiuu[fi];  cf.jpuv -= c*nlprivate->cpsiuv[fi];
          cf.jpvv -= c*nlprivate->cpsivv[fi];
        }
      cf.tang = nlprivate->ctang[knot];

#ifdef DEBUG_FVAL
DComputePDer ( domain, nlprivate, coeff, knot, &cf );
#endif

        /* integrate the functional */
      _g1hq2_IntFunc1cf ( &cf, &cfunct );

        /* integrate the gradient */
      for ( i = 0; i < nfunc_a; i++ ) {
        fi = i*3*G1_NQUAD*hole_k + knot;
        cf.psiu = nlprivate->cpsiu[fi];     cf.psiv = nlprivate->cpsiv[fi];
        cf.jpsiuu = nlprivate->cpsiuu[fi];  cf.jpsiuv = nlprivate->cpsiuv[fi];
        cf.jpsivv = nlprivate->cpsivv[fi];
        _g1hq2_IntFunc2cf ( &cf, &Di[nfunc_c+i], &Bi[nfunc_c+i],
                            &cgrad[nfunc_c+i] );
      }
      for ( i = 0; i < G1_DBDIM; i++ ) {
        ii = k*G1_DBDIM+i;
        fi = ((nfunc_a+1)*3*hole_k + ii*4)*G1_NQUAD + kn;
        cf.psiu = nlprivate->cpsiu[fi];     cf.psiv = nlprivate->cpsiv[fi];
        cf.jpsiuu = nlprivate->cpsiuu[fi];  cf.jpsiuv = nlprivate->cpsiuv[fi];
        cf.jpsivv = nlprivate->cpsivv[fi];
        _g1hq2_IntFunc2cf ( &cf, &Di[ii], &Bi[ii], &cgrad[ii] );
      }
      if ( kn < G1_NQUAD )
        for ( i = 0; i < G1_DBDIM; i++ ) {
          ii = ki*G1_DBDIM+i;
          fi = ((nfunc_a+1)*3*hole_k + ii*4 + 3)*G1_NQUAD + kn;
          cf.psiu = nlprivate->cpsiu[fi];     cf.psiv = nlprivate->cpsiv[fi];
          cf.jpsiuu = nlprivate->cpsiuu[fi];  cf.jpsiuv = nlprivate->cpsiuv[fi];
          cf.jpsivv = nlprivate->cpsivv[fi];
          _g1hq2_IntFunc2cf ( &cf, &Di[ii], &Bi[ii], &cgrad[ii] );
        }

        /* integrate the Hessian */
          /* A x A */
      for ( i = 0; i < nfunc_a; i++ ) {
        fi = i*3*G1_NQUAD*hole_k+k;
        ii = nfunc_c+i;
        cf.psiu = nlprivate->cpsiu[fi];     cf.psiv = nlprivate->cpsiv[fi];
        cf.jpsiuu = nlprivate->cpsiuu[fi];  cf.jpsiuv = nlprivate->cpsiuv[fi];
        cf.jpsivv = nlprivate->cpsivv[fi];
        for ( j = 0; j <= i; j++ ) {
          fj = j*3*G1_NQUAD*hole_k+k;
          jj = nfunc_c+j;
          cf.psju = nlprivate->cpsiu[fj];     cf.psjv = nlprivate->cpsiv[fj];
          cf.jpsjuu = nlprivate->cpsiuu[fj];  cf.jpsjuv = nlprivate->cpsiuv[fj];
          cf.jpsjvv = nlprivate->cpsivv[fj];
          _g1hq2_IntFunc3cf ( &cf, Di[ii], Bi[ii], Di[jj], Bi[jj],
                              &hkk[pkn_SymMatIndex(G1_DBDIM+i,G1_DBDIM+j)] );
        }
          /* A x C */
        for ( j = 0; j < G1_DBDIM; j++ ) {
          jj = k*G1_DBDIM+j;
          fj = ((nfunc_a+1)*3*hole_k + jj*4)*G1_NQUAD + kn;
          cf.psju = nlprivate->cpsiu[fj];     cf.psjv = nlprivate->cpsiv[fj];
          cf.jpsjuu = nlprivate->cpsiuu[fj];  cf.jpsjuv = nlprivate->cpsiuv[fj];
          cf.jpsjvv = nlprivate->cpsivv[fj];
          if ( k < hole_k-1 )
            hp = &hki[(G1_DBDIM+i)*G1_DBDIM+j];
          else
            hp = &hki[pkn_SymMatIndex(G1_DBDIM+i,j)];
          _g1hq2_IntFunc3cf ( &cf, Di[ii], Bi[ii], Di[jj], Bi[jj], hp );
        }
        if ( kn < G1_NQUAD ) {
          for ( j = 0; j < G1_DBDIM; j++ ) {
            jj = ki*G1_DBDIM+j;
            fi = ((nfunc_a+1)*3*hole_k + jj*4 + 3)*G1_NQUAD + kn;
            cf.psju = nlprivate->cpsiu[fj];     cf.psjv = nlprivate->cpsiv[fj];
            cf.jpsjuu = nlprivate->cpsiuu[fj];  cf.jpsjuv = nlprivate->cpsiuv[fj];
            cf.jpsjvv = nlprivate->cpsivv[fj];
            if ( k > 0 )
              hp = &hkj[(i+G1_DBDIM)*G1_DBDIM+j];
            else
              hp = &hkk[pkn_SymMatIndex(i+G1_DBDIM,j)];
            _g1hq2_IntFunc3cf ( &cf, Di[ii], Bi[ii], Di[jj], Bi[jj], hp );
          }
        }
      }
          /* C x C */
      for ( i = 0; i < G1_DBDIM; i++ ) {
        ii = k*G1_DBDIM+i;
        fi = ((nfunc_a+1)*3*hole_k + ii*4)*G1_NQUAD + kn;
        cf.psiu = nlprivate->cpsiu[fi];     cf.psiv = nlprivate->cpsiv[fi];
        cf.jpsiuu = nlprivate->cpsiuu[fi];  cf.jpsiuv = nlprivate->cpsiuv[fi];
        cf.jpsivv = nlprivate->cpsivv[fi];
        for ( j = 0; j <= i; j++ ) {
          jj = k*G1_DBDIM+j;
          fj = ((nfunc_a+1)*3*hole_k + jj*4)*G1_NQUAD + kn;
          cf.psju = nlprivate->cpsiu[fj];     cf.psjv = nlprivate->cpsiv[fj];
          cf.jpsjuu = nlprivate->cpsiuu[fj];  cf.jpsjuv = nlprivate->cpsiuv[fj];
          cf.jpsjvv = nlprivate->cpsivv[fj];
          _g1hq2_IntFunc3cf ( &cf, Di[ii], Bi[ii], Di[jj], Bi[jj],
                              &hii[pkn_SymMatIndex(i,j)] );
        }
      }
      if ( kn < G1_NQUAD ) {
        for ( i = 0; i < G1_DBDIM; i++ ) {
          ii = k*G1_DBDIM+i;
          fi = ((nfunc_a+1)*3*hole_k + ii*4)*G1_NQUAD + kn;
          cf.psiu = nlprivate->cpsiu[fi];     cf.psiv = nlprivate->cpsiv[fi];
          cf.jpsiuu = nlprivate->cpsiuu[fi];  cf.jpsiuv = nlprivate->cpsiuv[fi];
          cf.jpsivv = nlprivate->cpsivv[fi];
          for ( j = 0; j < G1_DBDIM; j++ ) {
            jj = ki*G1_DBDIM+j;
            fi = ((nfunc_a+1)*3*hole_k + jj*4 + 3)*G1_NQUAD + kn;
            cf.psju = nlprivate->cpsiu[fj];     cf.psjv = nlprivate->cpsiv[fj];
            cf.jpsjuu = nlprivate->cpsiuu[fj];  cf.jpsjuv = nlprivate->cpsiuv[fj];
            cf.jpsjvv = nlprivate->cpsivv[fj];
            if ( k > 0 )
              hp = &hij[i*G1_DBDIM+j];
            else
              hp = &hij[j*G1_DBDIM+i];
            _g1hq2_IntFunc3cf ( &cf, Di[ii], Bi[ii], Di[jj], Bi[jj], hp );
          }
        }
        for ( i = 0; i < G1_DBDIM; i++ ) {
          ii = ki*G1_DBDIM+i;
          fi = ((nfunc_a+1)*3*hole_k + ii*4 + 3)*G1_NQUAD + kn;
          cf.psiu = nlprivate->cpsiu[fi];     cf.psiv = nlprivate->cpsiv[fi];
          cf.jpsiuu = nlprivate->cpsiuu[fi];  cf.jpsiuv = nlprivate->cpsiuv[fi];
          cf.jpsivv = nlprivate->cpsivv[fi];
          for ( j = 0; j <= i; j++ ) {
            jj = ki*G1_DBDIM+j;
            fj = ((nfunc_a+1)*3*hole_k + jj*4 + 3)*G1_NQUAD + kn;
            cf.psju = nlprivate->cpsiu[fj];     cf.psjv = nlprivate->cpsiv[fj];
            cf.jpsjuu = nlprivate->cpsiuu[fj];  cf.jpsjuv = nlprivate->cpsiuv[fj];
            cf.jpsjvv = nlprivate->cpsivv[fj];
            _g1hq2_IntFunc3cf ( &cf, Di[ii], Bi[ii], Di[jj], Bi[jj],
                                &hjj[pkn_SymMatIndex(i,j)] );
          }
        }
      }
    }
  }
  cfunct /= (float)G1_NQUAD;

printf ( "f1 = %10.6f, f2 = %10.6f, C = %10.6f, f = %10.6f\n",
         funct, cfunct, C, funct+C*cfunct ); 

  *func = funct+C*cfunct;
  pkn_AddMatrixMf ( 1, nfunc, 0, grad, 0, cgrad, C/(float)G1_NQUAD, 0, grad );
  pkn_AddMatrixMf ( 1, asize, 0, hessian, 0, chessian, C/(float)G1_NQUAD,
                    0, hessian );

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*_g1hq2_ComputeExtNLFuncGradHessianf*/

/* ///////////////////////////////////////////////////////////////////////// */
#ifdef DEBUG_HESSIAN
static boolean _TestHessian ( GHoleDomainf *domain, G1HNLPrivatef *nlpr, float C )
{
  void   *sp;
  G1HolePrivateRecf *privateG1;
  GHolePrivateRecf  *privateG;
  int    hole_k, nfunc_a, nfunc_b, nfunc_c, nfunc, asize;
  float  func, *grad, *hessian, *coeff;
  float  *amat;
  int    i, j, pos;
  FILE   *f;

  sp = pkv_GetScratchMemTop ();
  privateG1 = domain->privateG1;
  privateG = domain->privateG;
  hole_k = domain->hole_k;
  nfunc_a = privateG1->nfunc_a;
  nfunc_b = privateG1->nfunc_b;
  nfunc_c = hole_k*G1_DBDIM;
  nfunc = nfunc_a+nfunc_c;
  coeff = pkv_GetScratchMemf ( nfunc );
  grad = pkv_GetScratchMemf ( nfunc );
  asize = pkn_Block3ArraySize ( hole_k-1, G1_DBDIM, G1_DBDIM+nfunc_a );
  hessian = pkv_GetScratchMemf ( 2*asize );
  if ( !coeff || !grad || !hessian ) {
    PKV_SIGNALERROR ( LIB_EGHOLE, ERRCODE_2, ERRMSG_2 );
    goto failure;
  }
  memset ( coeff, 0, nfunc*sizeof(float) );
  if ( !_g1hq2_ComputeExtNLFuncGradHessianf ( domain, nlpr, coeff,
                                5.0*C/nlpr->ddiam,
                                &func, grad, hessian, &hessian[asize] ) )
    goto failure;
  f = fopen ( "g1q2ehessianf.txt", "w+" );
  amat = privateG1->Q2EAMat;
  for ( i = 0; i < nfunc; i++ )
    for ( j = 0; j <= i; j++ ) {
      pos = pkn_Block3FindElemPos ( hole_k-1, G1_DBDIM, G1_DBDIM+nfunc_a, i, j );
      if ( pos >= 0 )
        if ( amat[pos] || hessian[pos] )
          fprintf ( f, "%4d,%4d: %12.5f %12.5f %14.7f\n", i, j,
                    0.5*amat[pos], hessian[pos], 0.5*amat[pos]-hessian[pos] );
    }
  fclose ( f );
  printf ( "%s\n", "g1q2ehessianf.txt" );
  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*_TestHessian*/
#endif
/* ///////////////////////////////////////////////////////////////////////// */
static boolean g1hq2_ExtNLNewtonf ( GHoleDomainf *domain,
                                    G1HNLPrivatef *nlprivate )
{
#define EPSF 2.0e-4
  void *sp;
  G1HolePrivateRecf *privateG1;
  int     itn, jtn, ktn, hole_k, nfunc_a, nfunc_c, nfunc, asize;
  float   *coeff, *dcoeff, *grad, *hii, *chii;
  float   func, func0, func1, aux, gn, gn0 = 0.0, gn1, dco, dyn;
  float   C;
  boolean positive;

  sp = pkv_GetScratchMemTop ();
  hole_k  = domain->hole_k;
  privateG1 = domain->privateG1;
        /* the penalty constant C1e is set by the procedure */
        /* computing the initial approximation of the solution */
  C = (float)(5.0*privateG1->C1e/nlprivate->ddiam);

printf ( "C = %f\n", privateG1->C1e );

  nfunc_a = privateG1->nfunc_a;
  nfunc_c = hole_k*G1_DBDIM;
  nfunc   = nfunc_a+nfunc_c;

  coeff = pkv_GetScratchMemf ( 3*nfunc );
  asize = pkn_Block3ArraySize ( hole_k-1, G1_DBDIM, G1_DBDIM+nfunc_a );
  hii    = pkv_GetScratchMemf ( 2*asize );
  if ( !coeff || !hii ) {
    domain->error_code = G1H_ERROR_NO_SCRATCH_MEMORY;
    goto failure;
  }
  dcoeff = &coeff[nfunc];
  grad   = &dcoeff[nfunc];
  chii   = &hii[asize];   

        /* setup the initial point */
  pkv_Selectf ( nfunc, 1, 3, 1, &nlprivate->acoeff[0].z, coeff );

        /* Newton iterations */
  for ( itn = ktn = 0; ; itn++ ) {
    if ( !_g1hq2_ComputeExtNLFuncGradHessianf ( domain, nlprivate, coeff, C,
                             &func, grad, hii, chii ) )
      goto failure;
    gn = (float)sqrt ( pkn_ScalarProductf ( nfunc, grad, grad ) );
    if ( itn == 0 ) {

printf ( "func = %f, gn0 = %f\n", func, gn );

      gn0 = gn;
    }
    memcpy ( chii, hii, asize*sizeof(float) );

    if ( (positive = pkn_Block3CholeskyDecompMf ( hole_k-1, G1_DBDIM,
                           G1_DBDIM+nfunc_a, hii )) ) {
      pkn_Block3LowerTrMSolvef ( hole_k-1, G1_DBDIM, G1_DBDIM+nfunc_a,
                           hii, 1, 1, grad );
      pkn_Block3UpperTrMSolvef ( hole_k-1, G1_DBDIM, G1_DBDIM+nfunc_a,
                           hii, 1, 1, grad );
    }
    else {

printf ( "! " );

      pkn_Block3SymMatrixMultf ( hole_k-1, G1_DBDIM, G1_DBDIM+nfunc_a, chii,
                                 1, 1, grad, 1, dcoeff );
      aux = (float)pkn_ScalarProductf ( nfunc, grad, dcoeff );
      if ( gn < 0.0 || aux < EPSF*gn ) {
        domain->error_code = G1H_ERROR_NL_MINIMIZATION;
        goto failure;
      }
      pkn_MultMatrixNumf ( 1, nfunc, 0, grad, gn/aux, 0, grad );
    }
    dco = (float)sqrt ( pkn_ScalarProductf(nfunc, coeff, coeff) );
    dyn = (float)sqrt ( pkn_ScalarProductf(nfunc, grad, grad) );

    for ( aux = 1.0; aux > EPSF; aux *= 0.5 ) {
      pkn_AddMatrixMf ( 1, nfunc, 0, coeff, 0, grad, aux, 0, dcoeff );
      func1 = _g1hq2_ComputeExtNLFuncf ( domain, nlprivate, dcoeff, C );
      if ( func1 < func )
        break;
    }

    memcpy ( coeff, dcoeff, nfunc*sizeof(float) );
    func = func1;
    ktn ++;
    if ( positive && aux > 0.1 ) {
        /* Now the Hessian is positive-definite; */
        /* as it is expensive to compute, we try to make some */
        /* extra iterations with the same Hessian. */

      for ( jtn = 0; jtn < 10; jtn++ ) {

printf ( "+" );

        _g1hq2_ComputeExtNLFuncGradf ( domain, nlprivate, coeff, C,
                                       &func0, grad );
        gn1 = (float)sqrt ( pkn_ScalarProductf ( nfunc, grad, grad ) );
        pkn_Block3LowerTrMSolvef ( hole_k-1, G1_DBDIM, G1_DBDIM+nfunc_a,
                                   hii, 1, 1, grad );
        pkn_Block3UpperTrMSolvef ( hole_k-1, G1_DBDIM, G1_DBDIM+nfunc_a,
                                   hii, 1, 1, grad );
        pkn_AddMatrixf ( 1, nfunc, 0, coeff, 0, grad, 0, dcoeff );
        func1 = _g1hq2_ComputeExtNLFuncf ( domain, nlprivate, dcoeff, C );
        if ( func1 >= func0 || gn1 >= 0.5*gn )
          break;

printf ( "    func = %f, gn = %f\n", func1, gn );

        memcpy ( coeff, dcoeff, nfunc*sizeof(float) );
        func = func1;
        ktn ++;
        gn = gn1;
        aux = 1.0;
        dco = (float)sqrt ( pkn_ScalarProductf(nfunc, coeff, coeff) );
        dyn = (float)sqrt ( pkn_ScalarProductf(nfunc, grad, grad) );
      }
    }

    if ( _g1h_StopItf ( itn, gn0, gn, dco, dyn, aux ) )
      break;
  }

printf ( "itn = %d, ktn = %d\n", itn+1, ktn );
printf ( "func = %f\n", func1 );

  pkv_Selectf ( nfunc, 1, 1, 3, coeff, &nlprivate->acoeff[0].z );

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
#undef EPSF
} /*g1hq2_ExtNLNewtonf*/

boolean g1h_Q2NLExtFillHolef ( GHoleDomainf *domain,
                    const point3f *hole_cp,
                    float *acoeff, void *usrptr,
                    void (*outpatch) ( int n, int m, const point3f *cp,
                                       void *usrptr ) )
{
  typedef void outscf ( int n, int m, const float *cp, void *usrptr );

  void   *sp;
  G1HNLPrivatef     *nlprivate;
  G1HolePrivateRecf *privateG1;
  float  *fc00, *Bi, *Bk;
  int    hole_k, nfunc_a, nfunc_b, nfunc_c;

  sp = pkv_GetScratchMemTop ();
  privateG1 = domain->privateG1;
  hole_k = domain->hole_k;
  nfunc_a = privateG1->nfunc_a;
  nfunc_b = privateG1->nfunc_b;
  nfunc_c = hole_k*G1_DBDIM;
  _g1h_nlprivf = nlprivate = _g1hq2_InitExtNLprf ( hole_k, nfunc_a, 0 );
  if ( !nlprivate )
    goto failure;

  if ( !_g1h_ComputeNLNormalf ( domain, nlprivate, hole_cp ) )
    goto failure;
  g1h_ReflectVectorsf ( 12*hole_k+1, hole_cp, nlprivate->rhole_cp );

  nlprivate->auxc = 0;
  if ( !g1h_Q2ExtFillHolef ( domain, 3, (float*)nlprivate->rhole_cp,
                             (float*)nlprivate->acoeff, NULL,
                             g1h_nonlinoutpatchf ) )
    goto failure;

  if ( !_g1hq2_TabExtNLBasisFunctionsf ( domain, nlprivate ) )
    goto failure;

#ifdef DEBUG_HESSIAN
_TestHessian ( domain, nlprivate, privateG1->C1e );
#endif

  if ( !g1hq2_ExtNLNewtonf ( domain, nlprivate ) )
    goto failure;

  g1h_ReflectVectorsf ( nfunc_a+nfunc_c, nlprivate->acoeff,
                        nlprivate->acoeff );
  if ( acoeff )
    memcpy ( acoeff, nlprivate->acoeff,
            (nfunc_a+nfunc_c)*sizeof(vector3f) );

  fc00 = pkv_GetScratchMem ( (G1_CROSSDEGSUM+4)*2*hole_k*sizeof(vector3f) );
  if ( !fc00 ) {
    domain->error_code = G1H_ERROR_NO_SCRATCH_MEMORY;
    goto failure;
  }
  Bi = privateG1->Q2EBMat;
  Bk = &Bi[hole_k*G1_DBDIM*nfunc_b];
  if ( !_g1h_SetExtRightSidef ( domain, Bi, Bk, 3, (float*)hole_cp, fc00, NULL ) )
    goto failure;

  if ( !_g1h_OutputExtPatchesf ( domain, 3, (float*)nlprivate->acoeff, fc00,
                                 usrptr, (outscf*)outpatch ) )
    goto failure;

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*g1h_Q2NLExtFillHolef*/

/* ///////////////////////////////////////////////////////////////////////// */
static boolean g1hq2_NLExtConstrNewtonf ( GHoleDomainf *domain,
                                          G1HNLPrivatef *nlprivate, int nconstr,
                                          float *ECmat )
{
#define EPSF 2.0e-4
  void    *sp;
  G1HolePrivateRecf *privateG1;
  int     itn, jtn, ktn, i, j, hole_k, nfunc_a, nfunc_c, nfunc_ac, nfunc;
  int     diagblsize, subdiagblsize, sideblsize, esideblsize, asize, esize;
  int     mS, ms;
  float   *coeff, *grad, *hessian, *hij, *hki, *hkk;
  float   *cT, *E22ii, *E22ij, *E22kk, *E22ki, *cE22ii;
  float   *aa, *D1, *y, *y1, *M, *f;
  float   func, func0, func1, aux, gn, gn0 = 0.0, gn1, dco, dyn, dyn1;
  float   C;
  boolean positive;

  sp = pkv_GetScratchMemTop ();
  hole_k    = domain->hole_k;
  privateG1 = domain->privateG1;
        /* the penalty constant C1e is set by the procedure */
        /* computing the initial approximation of the solution */
  C = (float)(5.0*privateG1->C1e/nlprivate->ddiam);
  nfunc_a   = privateG1->nfunc_a;
  nfunc_c   = hole_k*G1_DBDIM;
  nfunc_ac  = nfunc_a+nfunc_c;
  nfunc     = nfunc_ac-nconstr;
  mS        = nfunc_a+G1_DBDIM;
  ms        = mS-nconstr;
  diagblsize    = G1_DBDIM*(G1_DBDIM+1)/2;
  subdiagblsize = G1_DBDIM*G1_DBDIM;
  sideblsize    = G1_DBDIM*mS;
  esideblsize   = G1_DBDIM*ms;

  coeff    = pkv_GetScratchMemf ( nfunc_ac );
  grad     = pkv_GetScratchMemf ( nfunc_ac );
  asize    = pkn_Block3ArraySize ( hole_k-1, G1_DBDIM, mS );
  hessian  = pkv_GetScratchMemf ( asize );
  cT       = pkv_GetScratchMemf ( mS*nconstr );
  aa       = pkv_GetScratchMemf ( 2*nconstr );
  D1       = pkv_GetScratchMemf ( (nconstr*(nconstr+1))/2);
  esize    = pkn_Block3ArraySize ( hole_k-1, G1_DBDIM, ms );
  E22ii    = pkv_GetScratchMemf ( asize );
  cE22ii   = pkv_GetScratchMemf ( esize );
  y        = pkv_GetScratchMemf ( nfunc_ac );
  y1       = pkv_GetScratchMemf ( nfunc_ac );
  M        = pkv_GetScratchMemf ( (mS*(mS+1))/2 );
  f        = pkv_GetScratchMemf ( nfunc_ac );
  if ( !coeff || !grad || !hessian || !cT || !aa ||
       !D1 || !E22ii || !cE22ii || !y || !y1 || !M || !f )
    goto failure;
  hij = &hessian[pkn_Block3FindBlockPos ( hole_k-1, G1_DBDIM, mS, 1, 0 )];
  hkk = &hessian[pkn_Block3FindBlockPos ( hole_k-1, G1_DBDIM, mS,
                                          hole_k-1, hole_k-1 )];
  hki = &hessian[pkn_Block3FindBlockPos ( hole_k-1, G1_DBDIM, mS, hole_k-1, 0 )];
  E22ij = &E22ii[pkn_Block3FindBlockPos ( hole_k-1, G1_DBDIM, ms, 1, 0 )];
  E22kk = &E22ii[pkn_Block3FindBlockPos ( hole_k-1, G1_DBDIM, ms,
                                          hole_k-1, hole_k-1 )];
  E22ki = &E22ii[pkn_Block3FindBlockPos ( hole_k-1, G1_DBDIM, ms,
                                          hole_k-1, 0 )];

        /* setup the initial point */
  pkv_Selectf ( nfunc_ac, 1, 3, 1, &nlprivate->acoeff[0].z, coeff );

        /* step 1: decompose the constraint equations matrix */
  pkv_TransposeMatrixf ( nconstr, mS, nfunc_ac, &ECmat[nfunc_ac-mS],
                         nconstr, cT );
  pkn_QRDecomposeMatrixf ( mS, nconstr, cT, aa );
  for ( i = 0; i < nconstr; i++ )
    for ( j = i; j < nconstr; j++ )
      D1[pkn_SymMatIndex(i,j)] = cT[i*nconstr+j];

        /* Newton iterations */
  for ( itn = ktn = 0; ; itn++ ) {
          /* step 2 */
    memcpy ( y, coeff, nfunc_ac*sizeof(float) );
    pkn_multiReflectVectorf ( mS, nconstr, cT, aa, 1, 1, &y[nfunc_ac-mS] );
    pkn_LowerTrMatrixSolvef ( nconstr, D1, 1,
                              3, &nlprivate->rhole_cp[12*hole_k+1].z,
                              1, &y[nfunc_ac-mS] );
    for ( i = nfunc_ac-mS; i < nfunc_ac-mS+nconstr; i++ )
      y[i] = -y[i];

          /* step 3 */
    if ( !_g1hq2_ComputeExtNLFuncGradHessianf ( domain, nlprivate, coeff, C,
                                                &func, grad, hessian, E22ii ) )
      goto failure;

    if ( !pkn_ComputeQTSQf ( mS, hkk, nconstr, cT, aa, M ) )
      goto failure;
    for ( i = 0; i < ms; i++ )
      for ( j = i; j < ms; j++ )
        E22kk[pkn_SymMatIndex(i,j)] = M[pkn_SymMatIndex(nconstr+i,nconstr+j)];
    for ( i = 0; i < hole_k-1; i++ )
      pkn_multiReflectVectorf ( mS, nconstr, cT, aa,
                                G1_DBDIM, G1_DBDIM, &hki[i*sideblsize] );
    pkv_Selectf ( hole_k-1, esideblsize, sideblsize, esideblsize,
                  &hki[nconstr*G1_DBDIM], E22ki );
    memcpy ( E22ii, hessian, (hole_k-1)*diagblsize*sizeof(float) );
    memcpy ( E22ij, hij, (hole_k-2)*subdiagblsize*sizeof(float) );
    memcpy ( cE22ii, E22ii, esize*sizeof(float) );
          /* step 4 */
    memcpy ( f, grad, nfunc_ac*sizeof(float) );
    pkn_multiReflectVectorf ( mS, nconstr, cT, aa, 1, 1, &f[nfunc_ac-mS] );
    memmove ( &f[nfunc_ac-mS], &f[nfunc_ac-mS+nconstr], ms*sizeof(float) );
    gn = (float)sqrt ( pkn_ScalarProductf ( nfunc, f, f ) );
    if ( itn == 0 ) {
/*
printf ( "func = %f, gn0 = %f\n", func, gn );
*/
      gn0 = gn;
    }
          /* step 5 */
    if ( (positive = pkn_Block3CholeskyDecompMf ( hole_k-1, G1_DBDIM,
                                      ms, E22ii )) ) {
      pkn_Block3LowerTrMSolvef ( hole_k-1, G1_DBDIM, ms, E22ii, 1, 1, f );
      pkn_Block3UpperTrMSolvef ( hole_k-1, G1_DBDIM, ms, E22ii, 1, 1, f );
    }
    else {

printf ( "! " );

      pkn_Block3SymMatrixMultf ( hole_k-1, G1_DBDIM, ms,
                                 cE22ii, 1, 1, f, 1, y1 );
      aux = (float)pkn_ScalarProductf ( nfunc, f, y1 );
      if ( aux <= 0.0 || aux < EPSF*gn )
        goto failure;
      pkn_MultMatrixNumf ( 1, nfunc, 1, f, gn/aux, 1, f );
    }

          /* step 6 */
    dyn = (float)sqrt ( pkn_ScalarProductf ( nfunc, f, f ) );
    memmove ( &f[nfunc_ac-mS+nconstr], &f[nfunc_ac-mS],
              ms*sizeof(float) );
    dco = (float)sqrt ( pkn_ScalarProductf ( nfunc_ac, coeff, coeff ) );
    for ( aux = 1.0; aux >= EPSF; aux *= 0.5 ) {
      for ( i = 0; i < nfunc_c; i++ )
        y1[i] = y[i]+aux*f[i];
      memcpy ( &y1[nfunc_ac-mS], &y[nfunc_ac-mS], nconstr*sizeof(float) );
      for ( i = nfunc_ac-mS+nconstr; i < nfunc_ac; i++ )
        y1[i] = y[i]+aux*f[i];
      pkn_multiInvReflectVectorf ( mS, nconstr, cT, aa, 1, 1, &y1[nfunc_ac-mS] );
      func1 = _g1hq2_ComputeExtNLFuncf ( domain, nlprivate, y1, C );
      if ( func1 < func )
        break;
    }
    memcpy ( coeff, y1, nfunc_ac*sizeof(float) );
    func = func1;
    ktn ++;

    if ( positive && aux > 0.1 ) {
        /* Now the Hessian is positive-definite; */
        /* as it is expensive to compute, we try to make some */
        /* extra iterations with the same Hessian. */

      for ( jtn = 0; jtn < 10; jtn++ ) {
/*
printf ( "+" );
*/
              /* step 2' */
        memcpy ( y, coeff, nfunc_ac*sizeof(float) );
        pkn_multiReflectVectorf ( mS, nconstr, cT, aa, 1, 1, &y[nfunc_ac-mS] );
        pkn_LowerTrMatrixSolvef ( nconstr, D1, 1,
                                  3, &nlprivate->rhole_cp[12*hole_k+1].z,
                                  1, &y[nfunc_ac-mS] );
        for ( i = nfunc_ac-mS; i < nfunc_ac-mS+nconstr; i++ )
          y[i] = -y[i];
              /* step 3' */
        _g1hq2_ComputeExtNLFuncGradf ( domain, nlprivate, coeff, C, &func0, grad );
              /* step 4' */
        memcpy ( f, grad, nfunc_ac*sizeof(float) );
        pkn_multiReflectVectorf ( mS, nconstr, cT, aa, 1, 1, &f[nfunc_ac-mS] );
        memmove ( &f[nfunc_ac-mS], &f[nfunc_ac-mS+nconstr],
                  ms*sizeof(float) );
        gn1 = (float)sqrt ( pkn_ScalarProductf ( nfunc, f, f ) );
              /* step 5' */
        pkn_Block3LowerTrMSolvef ( hole_k-1, G1_DBDIM, ms, E22ii, 1, 1, f );
        pkn_Block3UpperTrMSolvef ( hole_k-1, G1_DBDIM, ms, E22ii, 1, 1, f );
              /* step 6' */
        dyn1 = (float)sqrt ( pkn_ScalarProductf ( nfunc, f, f ) );
        memmove ( &f[nfunc_ac-mS+nconstr], &f[nfunc_ac-mS],
                  ms*sizeof(float) );
        for ( i = 0; i < nfunc_ac-mS; i++ )
          y1[i] = y[i]+f[i];
        memcpy ( &y1[nfunc_ac-mS], &y[nfunc_ac-mS], nconstr*sizeof(float) );
        for ( i = nfunc_ac-mS+nconstr; i < nfunc_ac; i++ )
          y1[i] = y[i]+f[i];
        pkn_multiInvReflectVectorf ( mS, nconstr, cT, aa, 1, 1, &y1[nfunc_ac-mS] );

        func1 = _g1hq2_ComputeExtNLFuncf ( domain, nlprivate, y1, C );
        if ( func1 >= func0 || gn1 >= 0.5*gn )
          break;

printf ( "    func = %f, gn = %f\n", func1, gn );

        memcpy ( coeff, y1, nfunc_ac*sizeof(float) );
        func = func1;
        ktn ++;
        gn = gn1;
        aux = 1.0;
        dco = (float)sqrt ( pkn_ScalarProductf(nfunc_ac, coeff, coeff) );
        dyn = dyn1;
      }
    }
    if ( _g1h_StopItf ( itn, gn0, gn, dco, dyn, aux ) )
      break;
  }

printf ( "itn = %d, ktn = %d\n", itn+1, ktn );
printf ( "func = %f\n", func1 );

  pkv_Selectf ( nfunc_ac, 1, 1, 3, coeff, &nlprivate->acoeff[0].z );

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
#undef EPSF
} /*g1hq2_NLExtConstrNewtonf*/

boolean g1h_Q2NLExtFillHoleConstrf ( GHoleDomainf *domain, CONST_ point3f *hole_cp,
                                  int nconstr, CONST_ vector3f *constr,
                                  float *acoeff, void *usrptr,
                                  void (*outpatch) ( int n, int m,
                                                     const point3f *cp,
                                                     void *usrptr ) )
{
  typedef void outscf ( int n, int m, const float *cp, void *usrptr );

  void   *sp;
  G1HNLPrivatef     *nlprivate;
  G1HolePrivateRecf *privateG1;
  float  *fc00, *Bi, *Bk;
  int    hole_k, nfunc_a, nfunc_b, nfunc_c, nfunc_ac;

  sp = pkv_GetScratchMemTop ();
  privateG1 = domain->privateG1;
  hole_k  = domain->hole_k;
  nfunc_a = privateG1->nfunc_a;
  nfunc_b = privateG1->nfunc_b;
  nfunc_c = hole_k*G1_DBDIM;
  nfunc_ac = nfunc_a+nfunc_c;
  _g1h_nlprivf = nlprivate = _g1hq2_InitExtNLprf ( hole_k, nfunc_a, nconstr );
  if ( !nlprivate )
    goto failure;

  if ( !_g1h_ComputeNLNormalf ( domain, nlprivate, hole_cp ) )
    goto failure;
  g1h_ReflectVectorsf ( 12*hole_k+1, hole_cp, nlprivate->rhole_cp );
  g1h_ReflectVectorsf ( nconstr, constr, &nlprivate->rhole_cp[12*hole_k+1] );

  nlprivate->auxc = 0;
  if ( !g1h_Q2ExtFillHoleConstrf ( domain, 3, (float*)nlprivate->rhole_cp,
                    nconstr, (float*)&nlprivate->rhole_cp[12*hole_k+1],
                    (float*)nlprivate->acoeff, NULL, g1h_nonlinoutpatchf ) )
    goto failure;

  if ( !_g1hq2_TabExtNLBasisFunctionsf ( domain, nlprivate ) )
    goto failure;

  if ( !g1hq2_NLExtConstrNewtonf ( domain, nlprivate, nconstr, privateG1->ECmat ) )
    goto failure;
  g1h_ReflectVectorsf ( nfunc_ac, nlprivate->acoeff, nlprivate->acoeff );
  if ( acoeff )
    memcpy ( acoeff, nlprivate->acoeff, nfunc_ac*sizeof(vector3f) );

  fc00 = pkv_GetScratchMem ( (G1_CROSSDEGSUM+4)*2*hole_k*sizeof(vector3f) );
  if ( !fc00 )
    goto failure;
  Bi = privateG1->Q2EBMat;
  Bk = &Bi[hole_k*G1_DBDIM*nfunc_b];
  if ( !_g1h_SetExtRightSidef ( domain, Bi, Bk, 3, (float*)hole_cp, fc00, NULL ) )
    goto failure;

  if ( !_g1h_OutputExtPatchesf ( domain, 3, (float*)nlprivate->acoeff, fc00,
                                 usrptr, (outscf*)outpatch ) )
    goto failure;

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*g1h_Q2NLExtFillHoleConstrf*/

/* ///////////////////////////////////////////////////////////////////////// */
boolean g1h_Q2NLExtFillHoleAltConstrf ( GHoleDomainf *domain, CONST_ point3f *hole_cp,
                                     int naconstr, CONST_ float *constr,
                                     float *acoeff, void *usrptr,
                                     void (*outpatch) ( int n, int m,
                                                        const point3f *cp,
                                                        void *usrptr ) )
{
  typedef void outscf ( int n, int m, const float *cp, void *usrptr );

  void     *sp;
  G1HNLPrivatef     *nlprivate;
  G1HolePrivateRecf *privateG1;
  float   *fc00, *Bi, *Bk;
  float   *saveconstr = NULL, *newconstr, *ECmat, *rsconstr;
  int      hole_k, nfunc_a, nfunc_b, nfunc_c, nfunc_ac, i, j;
  boolean  restore;

  sp = pkv_GetScratchMemTop ();
  restore = false;
  privateG1 = domain->privateG1;
  if ( naconstr <= 0 || privateG1->extacdim != 3 ||
       privateG1->extnaconstr != naconstr ) {
    domain->error_code = G1H_ERROR_INCONSISTENT_CONSTR;
    goto failure;
  }
  hole_k  = domain->hole_k;
  nfunc_a = privateG1->nfunc_a;
  nfunc_b = privateG1->nfunc_b;
  nfunc_c = hole_k*G1_DBDIM;
  nfunc_ac = nfunc_a+nfunc_c;
  saveconstr = pkv_GetScratchMemf ( 7*nfunc_ac*naconstr );
  _g1h_nlprivf = nlprivate = _g1hq2_InitExtNLprf ( hole_k, nfunc_a, naconstr );
  if ( !saveconstr || !nlprivate )
    goto failure;

  newconstr = &saveconstr[3*nfunc_ac*naconstr];
  ECmat = &newconstr[3*nfunc_ac*naconstr];
  memcpy ( saveconstr, privateG1->AECmat, nfunc_ac*naconstr*3*sizeof(float) );
  memcpy ( newconstr, privateG1->AECmat, nfunc_ac*naconstr*3*sizeof(float) );
  if ( !_g1h_ComputeNLNormalf ( domain, nlprivate, hole_cp ) )
    goto failure;
  g1h_ReflectVectorsf ( 12*hole_k+1, hole_cp, nlprivate->rhole_cp );
  if ( !_g1h_ReflectExtAltConstrMatrixf ( domain, &nlprivate->reflv, newconstr ) )
    goto failure;

  restore = true;
  if ( !g1h_SetExtAltConstraintMatrixf ( domain, 3, naconstr, newconstr ) )
    goto failure;
  nlprivate->auxc = 0;
  if ( !g1h_Q2ExtFillHoleAltConstrf ( domain, 3, (float*)nlprivate->rhole_cp,
                    naconstr, constr,
                    (float*)nlprivate->acoeff, NULL, g1h_nonlinoutpatchf ) )
    goto failure;

  pkv_Selectf ( naconstr, nfunc_ac, 3*nfunc_ac, nfunc_ac,
                &newconstr[2*nfunc_ac], ECmat );
  rsconstr = &nlprivate->rhole_cp[12*hole_k+1].z;
  pkv_Selectf ( naconstr, 1, 1, 3, constr, rsconstr );
  for ( j = 0; j < naconstr; j++ )
    for ( i = 0; i < nfunc_ac; i++ )
      rsconstr[3*j] += newconstr[3*nfunc_ac*j+i]*nlprivate->acoeff[i].x +
                       newconstr[(3*j+1)*nfunc_ac+i]*nlprivate->acoeff[i].y;

  if ( !_g1hq2_TabExtNLBasisFunctionsf ( domain, nlprivate ) )
    goto failure;

  if ( !g1hq2_NLExtConstrNewtonf ( domain, nlprivate, naconstr, ECmat ) )
    goto failure;

  g1h_ReflectVectorsf ( nfunc_ac, nlprivate->acoeff, nlprivate->acoeff );
  if ( acoeff )
    memcpy ( acoeff, nlprivate->acoeff, nfunc_ac*sizeof(vector3f) );

  fc00 = pkv_GetScratchMem ( (G1_CROSSDEGSUM+4)*2*hole_k*sizeof(vector3f) );
  if ( !fc00 )
    goto failure;
  Bi = privateG1->Q2EBMat;
  Bk = &Bi[nfunc_c*nfunc_b];
  if ( !_g1h_SetExtRightSidef ( domain, Bi, Bk, 3, (float*)hole_cp, fc00, NULL ) )
    goto failure;

  if ( !_g1h_OutputExtPatchesf ( domain, 3, (float*)nlprivate->acoeff, fc00,
                                 usrptr, (outscf*)outpatch ) )
    goto failure;

  g1h_SetExtAltConstraintMatrixf ( domain, 3, naconstr, saveconstr );
  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  if ( restore )
    g1h_SetExtAltConstraintMatrixf ( domain, 3, naconstr, saveconstr );
  pkv_SetScratchMemTop ( sp );
  return false;
} /*g1h_Q2NLExtFillHoleAltConstrf*/


