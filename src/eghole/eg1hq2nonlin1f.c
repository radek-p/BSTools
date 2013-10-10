
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


/* ///////////////////////////////////////////////////////////////////////// */
boolean _g1hq2_FindDomSurrndPatchf ( GHoleDomainf *domain,
                                     G1HNLPrivatef *nlpr,
                                     int i, int j, point2f *bezcp )
{
  void    *sp;
  int     hole_k, k;
  int     *ind;
  point2f *q;
  float   *ukn, *vkn;

  sp = pkv_GetScratchMemTop ();
  ind = (int*)pkv_GetScratchMem ( 16*sizeof(int) );
  q = (point2f*)pkv_GetScratchMem ( 16*sizeof(point2f) );
  if ( !ind || !q ) {
    domain->error_code = G1H_ERROR_NO_SCRATCH_MEMORY;
    goto failure;
  }

  hole_k = domain->hole_k;
  ukn = _gh_GetKnotSequencef ( domain, i-1 );  ukn += 3;
  vkn = _gh_GetKnotSequencef ( domain, i );    vkn += j;
  gh_GetBspInd ( hole_k, i, j, ind );
  for ( k = 0; k < 16; k++ )
    memcpy ( &q[k], &nlpr->rhole_cp[ind[k]], sizeof(point2f) );
  mbs_BSPatchToBezf ( 2, 3, 7, ukn, 3, 7, vkn, 8, (float*)q,
                      NULL, NULL, NULL, NULL, NULL, NULL, 8, (float*)bezcp );
  if ( j == 1 )
    mbs_multiReverseBSCurvef ( 3, 0, NULL, 4, 2, 8, (float*)bezcp );

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*_g1hq2_FindDomSurrndPatchf*/

boolean _g1hq2_FindNLDomainDiameterf ( GHoleDomainf *domain,
                                       G1HNLPrivatef *nlpr )
{
  void     *sp;
  int      hole_k, i, j, k;
  point2f  *p, *q;
  vector2f v;
  float    ddiam, d;

  sp = pkv_GetScratchMemTop ();
  hole_k = domain->hole_k;
  p = pkv_GetScratchMem ( 6*hole_k*sizeof(point2f) );
  q = pkv_GetScratchMem ( 16*sizeof(point2f) );
  if ( !p || !q )
    goto failure;
  for ( i = k = 0;  i < hole_k;  i++ ) {
    _g1hq2_FindDomSurrndPatchf ( domain, nlpr, i, 1, q );
    p[k++] = q[0];  p[k++] = q[1];  p[k++] = q[2];
    _g1hq2_FindDomSurrndPatchf ( domain, nlpr, i, 2, q );
    p[k++] = q[1];  p[k++] = q[2];  p[k++] = q[3];
  }
  ddiam = 0.0;
  for ( i = 1; i < k; i++ )
    for ( j = 0; j < i; j++ ) {
      SubtractPoints2f ( &p[i], &p[j], &v );
      d = v.x*v.x+v.y*v.y;
      ddiam = max ( ddiam, d );
    }
  nlpr->ddiam = (float)sqrt(ddiam);

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*_g1hq2_FindNLDomainDiameterf*/

/* ///////////////////////////////////////////////////////////////////////// */
boolean _g1hq2_TabNLDer0f ( int nkn, const float *tkn,
             const float *hfunc, const float *dhfunc, const float *ddhfunc,
             const float *dddhfunc,
             vector2f *diu, vector2f *div,
             vector2f *diuu, vector2f *diuv, vector2f *divv,
             vector2f *diuuu, vector2f *diuuv, vector2f *diuvv, vector2f *divvv,
             float *fc00, float *fc01, float *fd00, float *fd01,
             float *psiu, float *psiv,
             float *psiuu, float *psiuv, float *psivv,
             float *psiuuu, float *psiuuv, float *psiuvv, float *psivvv )
{
  void  *sp;
  int   i;
  float *hu, *hv, *huu, *huv, *hvv, *huuu, *huuv, *huvv, *hvvv;

  sp = pkv_GetScratchMemTop ();
  hu = pkv_GetScratchMemf ( 9*nkn*nkn );
  if ( !hu )
    goto failure;
  hv = &hu[nkn*nkn];      huu = &hv[nkn*nkn];
  huv = &huu[nkn*nkn];    hvv = &huv[nkn*nkn];
  huuu = &hvv[nkn*nkn];   huuv = &huuu[nkn*nkn];
  huvv = &huuv[nkn*nkn];  hvvv = &huvv[nkn*nkn];

  if ( !mbs_TabBezC1Coons0Der3f ( 1, nkn, tkn, hfunc, dhfunc, ddhfunc, dddhfunc,
                     nkn, tkn, hfunc, dhfunc, ddhfunc, dddhfunc,
                     G1_CROSS00DEG, fc00, G1_CROSS01DEG, fc01,
                     G1_CROSS00DEG, fd00, G1_CROSS01DEG, fd01,
                     NULL, hu, hv, huu, huv, hvv, huuu, huuv, huvv, hvvv ) )
    goto failure;

  for ( i = 0; i < nkn*nkn; i++ )
    if ( !pkn_Comp2iDerivatives3f ( diu[i].x, diu[i].y, div[i].x, div[i].y,
              diuu[i].x, diuu[i].y, diuv[i].x, diuv[i].y, divv[i].x, divv[i].y,
              diuuu[i].x, diuuu[i].y, diuuv[i].x, diuuv[i].y,
              diuvv[i].x, diuvv[i].y, divvv[i].x, divvv[i].y,
              1, &hu[i], &hv[i], &huu[i], &huv[i], &hvv[i],
              &huuu[i], &huuv[i], &huvv[i], &hvvv[i],
              &psiu[i], &psiv[i], &psiuu[i], &psiuv[i], &psivv[i],
              &psiuuu[i], &psiuuv[i], &psiuvv[i], &psivvv[i] ) )
      goto failure;

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*_g1hq2_TabNLDer0f*/

boolean _g1hq2_TabNLDerf ( int nkn, float *tkn,
             const float *hfunc, const float *dhfunc, const float *ddhfunc,
             const float *dddhfunc,
             vector2f *diu, vector2f *div,
             vector2f *diuu, vector2f *diuv, vector2f *divv,
             vector2f *diuuu, vector2f *diuuv, vector2f *diuvv, vector2f *divvv,
             float *fc00, float *fc01, float *fc10, float *fc11,
             float *fd00, float *fd01, float *fd10, float *fd11,
             float *psiu, float *psiv,
             float *psiuu, float *psiuv, float *psivv,
             float *psiuuu, float *psiuuv, float *psiuvv, float *psivvv )
{
  void  *sp;
  int   i;
  float *hu, *hv, *huu, *huv, *hvv, *huuu, *huuv, *huvv, *hvvv;

  sp = pkv_GetScratchMemTop ();
  hu = pkv_GetScratchMemf ( 9*nkn*nkn );
  if ( !hu )
    goto failure;
  hv = &hu[nkn*nkn];      huu = &hv[nkn*nkn];
  huv = &huu[nkn*nkn];    hvv = &huv[nkn*nkn];
  huuu = &hvv[nkn*nkn];   huuv = &huuu[nkn*nkn];
  huvv = &huuv[nkn*nkn];  hvvv = &huvv[nkn*nkn];

  if ( !mbs_TabBezC1CoonsDer3f ( 1, nkn, tkn, hfunc, dhfunc, ddhfunc, dddhfunc,
                     nkn, tkn, hfunc, dhfunc, ddhfunc, dddhfunc,
                     G1_CROSS00DEG, fc00, G1_CROSS01DEG, fc01,
                     G1_CROSS10DEG, fc10, G1_CROSS11DEG, fc11,
                     G1_CROSS00DEG, fd00, G1_CROSS01DEG, fd01,
                     G1_CROSS10DEG, fd10, G1_CROSS11DEG, fd11,
                     NULL, hu, hv, huu, huv, hvv, huuu, huuv, huvv, hvvv ) )
    goto failure;

  for ( i = 0; i < nkn*nkn; i++ )
    if ( !pkn_Comp2iDerivatives3f ( diu[i].x, diu[i].y, div[i].x, div[i].y,
              diuu[i].x, diuu[i].y, diuv[i].x, diuv[i].y, divv[i].x, divv[i].y,
              diuuu[i].x, diuuu[i].y, diuuv[i].x, diuuv[i].y,
              diuvv[i].x, diuvv[i].y, divvv[i].x, divvv[i].y,
              1, &hu[i], &hv[i], &huu[i], &huv[i], &hvv[i],
              &huuu[i], &huuv[i], &huvv[i], &hvvv[i],
              &psiu[i], &psiv[i], &psiuu[i], &psiuv[i], &psivv[i],
              &psiuuu[i], &psiuuv[i], &psiuvv[i], &psivvv[i] ) )
      goto failure;

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*_g1hq2_TabNLDerf*/

boolean _g1hq2_SetupCTrdf ( const vector2f *cdiu, const vector2f *cdiv,
           const vector2f *cdiuu, const vector2f *cdiuv, const vector2f *cdivv,
           float *ctrd )
{
  vector2f gx, gy, gxx, gxy, gyy;

  if ( !pkn_f2iDerivatives2f ( cdiu->x, cdiu->y, cdiv->x, cdiv->y,
            cdiuu->x, cdiuu->y, cdiuv->x, cdiuv->y, cdivv->x, cdivv->y,
            &gx.x, &gy.x, &gxx.x, &gxy.x, &gyy.x ) )
    return false;
  pkn_Setup2DerA11Matrixf ( gx.x, gx.y, gy.x, gy.y, ctrd );
  pkn_Setup2DerA21Matrixf ( gxx.x, gxx.y, gxy.x, gxy.y, gyy.x, gyy.y, &ctrd[4] );
  pkn_Setup2DerA22Matrixf ( gx.x, gx.y, gy.x, gy.y, &ctrd[10] );
  return true;
} /*_g1hq2_SetupCTrdf*/

boolean _g1hq2_TabNLBasisFunctionsOmegaf ( GHoleDomainf *domain, int nkn,
                                           G1HNLPrivatef *nlpr, float *bc00 )
{
#define N ((G1H_FINALDEG+1)*(G1H_FINALDEG+1))
  void     *sp;
  GHolePrivateRecf  *privateG;
  G1HolePrivateRecf *privateG1;
  float    *tkn;
  float    *fc00, *fc01, *fc10, *fc11, *fd00, *fd01, *fd10, *fd11,
           *bc01, *bc10, *bc11, *bd00, *bd01, *bd10, *bd11;
  float    *hfunc, *dhfunc, *ddhfunc, *dddhfunc;
  int      i, j, k, l, kN, kNQ2, hole_k, f, fN, nfunc_a, nfunc_b, nkn2;
  point2f  *di, ddi;
  float    bvz;
  unsigned char *bfcpn;

  sp      = pkv_GetScratchMemTop ();
  if ( !_g1hq2_FindNLDomainDiameterf ( domain, nlpr ) )
    goto failure;
  hole_k  = domain->hole_k;
  privateG  = domain->privateG;
  privateG1 = domain->privateG1;
  nfunc_a = privateG1->nfunc_a;
  nfunc_b = privateG1->nfunc_b;
  di      = pkv_GetScratchMem ( N*sizeof(point2f) );
  tkn     = pkv_GetScratchMemf ( nkn );
  hfunc = pkv_GetScratchMemf ( 16*nkn );
  if ( !tkn || !di || !hfunc ) {
    domain->error_code = G1H_ERROR_NO_SCRATCH_MEMORY;
    goto failure;
  }

  dhfunc = &hfunc[4*nkn];
  ddhfunc = &dhfunc[4*nkn];
  dddhfunc = &ddhfunc[4*nkn];

  nkn2 = nkn*nkn;
  _gh_PrepareTabKnotsf ( nkn, privateG1->opt_quad, tkn );
  mbs_TabCubicHFuncDer3f ( 0.0, 1.0, nkn, tkn,
                           hfunc, dhfunc, ddhfunc, dddhfunc );

  for ( k = kN = kNQ2 = 0;
        k < hole_k;
        k++, kN += N, kNQ2 += nkn2 ) {
    pkv_Selectf ( N, 2, 3, 2, &nlpr->nldi[kN], di );
    for ( i = l = 0;  i < nkn;  i++ )
         /* ***** this is not an optimal algorithm, to be improved */
      for ( j = 0;  j < nkn;  j++, l++ ) {
        if ( !mbs_BCHornerDer3Pf ( G1H_FINALDEG, G1H_FINALDEG, 2, (float*)di,
                    tkn[i], tkn[j], (float*)&ddi,
                    (float*)&nlpr->diu[kNQ2+l], (float*)&nlpr->div[kNQ2+l],
                    (float*)&nlpr->diuu[kNQ2+l], (float*)&nlpr->diuv[kNQ2+l],
                    (float*)&nlpr->divv[kNQ2+l],
                    (float*)&nlpr->diuuu[kNQ2+l], (float*)&nlpr->diuuv[kNQ2+l],
                    (float*)&nlpr->diuvv[kNQ2+l], (float*)&nlpr->divvv[kNQ2+l] ) )
          goto failure;
        nlpr->jac[kNQ2+l] = (float)det2f ( &nlpr->diu[kNQ2+l], &nlpr->div[kNQ2+l] );
      }
  }

  /* verify the Jacobian */
  if ( nlpr->jac[0] > 0.0 ) {
    for ( i = 1; i < hole_k*nkn2; i++ )
      if ( nlpr->jac[i] <= 0.0 ) {
        domain->error_code = G1H_ERROR_NL_JACOBIAN;
        goto failure;
      }
  }
  else {
    for ( i = 0; i < hole_k*nkn2; i++ ) {
      if ( nlpr->jac[i] >= 0.0 ) {
        domain->error_code = G1H_ERROR_NL_JACOBIAN;
        goto failure;
      }
      nlpr->jac[i] = -nlpr->jac[i];
    }
  }

  for ( f = fN = 0;  f < nfunc_a;  f++ ) {
    for ( k = kNQ2 = 0;  k < hole_k;  k++, fN += nkn2, kNQ2 += nkn2 ) {
      _g1h_GetBFAPatchCurvesf ( domain, f, k, &fc00, &fc01, &fd00, &fd01 );
      if ( !_g1hq2_TabNLDer0f ( nkn, tkn,
                        hfunc, dhfunc, ddhfunc, dddhfunc,
                        &nlpr->diu[kNQ2], &nlpr->div[kNQ2],
                        &nlpr->diuu[kNQ2], &nlpr->diuv[kNQ2], &nlpr->divv[kNQ2],
                        &nlpr->diuuu[kNQ2], &nlpr->diuuv[kNQ2],
                        &nlpr->diuvv[kNQ2], &nlpr->divvv[kNQ2],
                        fc00, fc01, fd00, fd01,
                        &nlpr->psiu[fN], &nlpr->psiv[fN],
                        &nlpr->psiuu[fN], &nlpr->psiuv[fN], &nlpr->psivv[fN],
                        &nlpr->psiuuu[fN], &nlpr->psiuuv[fN],
                        &nlpr->psiuvv[fN], &nlpr->psivvv[fN] ) ) {
        domain->error_code = G1H_ERROR_NO_SCRATCH_MEMORY;
        goto failure;
      }
    }
  }

  bc01 = &bc00[(G1_CROSS00DEG+1)*hole_k];  bc10 = &bc01[(G1_CROSS01DEG+1)*hole_k];
  bc11 = &bc10[(G1_CROSS10DEG+1)*hole_k];  bd00 = &bc11[(G1_CROSS11DEG+1)*hole_k];
  bd01 = &bd00[(G1_CROSS00DEG+1)*hole_k];  bd10 = &bd01[(G1_CROSS01DEG+1)*hole_k];
  bd11 = &bd10[(G1_CROSS10DEG+1)*hole_k];

  memset ( bc00, 0, 2*(G1_CROSSDEGSUM+4)*hole_k*sizeof(float) );
  bfcpn = privateG->bfcpn;
  for ( f = 0; f < nfunc_b; f++ ) {
        /* find Coons representation of the constant part of the solution */
    bvz = nlpr->rhole_cp[bfcpn[f]].z;
    for ( k = 0;  k < hole_k;  k++ ) {
      _g1h_GetBFBPatchCurvesf ( domain, f, k, &fc00, &fc01, &fc10, &fc11,
                                              &fd00, &fd01, &fd10, &fd11 );
      pkn_AddMatrixMf ( 1, G1_CROSS00DEG+1, 0, &bc00[(G1_CROSS00DEG+1)*k], 0, fc00,
                        bvz, 0, &bc00[(G1_CROSS00DEG+1)*k] );
      pkn_AddMatrixMf ( 1, G1_CROSS01DEG+1, 0, &bc01[(G1_CROSS01DEG+1)*k], 0, fc01,
                        bvz, 0, &bc01[(G1_CROSS01DEG+1)*k] );
      pkn_AddMatrixMf ( 1, G1_CROSS10DEG+1, 0, &bc10[(G1_CROSS10DEG+1)*k], 0, fc10,
                        bvz, 0, &bc10[(G1_CROSS10DEG+1)*k] );
      pkn_AddMatrixMf ( 1, G1_CROSS11DEG+1, 0, &bc11[(G1_CROSS11DEG+1)*k], 0, fc11,
                        bvz, 0, &bc11[(G1_CROSS11DEG+1)*k] );
      pkn_AddMatrixMf ( 1, G1_CROSS00DEG+1, 0, &bd00[(G1_CROSS00DEG+1)*k], 0, fd00,
                        bvz, 0, &bd00[(G1_CROSS00DEG+1)*k] );
      pkn_AddMatrixMf ( 1, G1_CROSS01DEG+1, 0, &bd01[(G1_CROSS01DEG+1)*k], 0, fd01,
                        bvz, 0, &bd01[(G1_CROSS01DEG+1)*k] );
      pkn_AddMatrixMf ( 1, G1_CROSS10DEG+1, 0, &bd10[(G1_CROSS10DEG+1)*k], 0, fd10,
                        bvz, 0, &bd10[(G1_CROSS10DEG+1)*k] );
      pkn_AddMatrixMf ( 1, G1_CROSS11DEG+1, 0, &bd11[(G1_CROSS11DEG+1)*k], 0, fd11,
                        bvz, 0, &bd11[(G1_CROSS11DEG+1)*k] );
    }
  }
  for ( k = kNQ2 = 0, fN = nfunc_a*nkn2*hole_k;
        k < hole_k;
        k++, kNQ2 += nkn2, fN += nkn2 ) {
    if ( !_g1hq2_TabNLDerf ( nkn, tkn, hfunc, dhfunc, ddhfunc, dddhfunc,
             &nlpr->diu[kNQ2], &nlpr->div[kNQ2],
             &nlpr->diuu[kNQ2], &nlpr->diuv[kNQ2], &nlpr->divv[kNQ2],
             &nlpr->diuuu[kNQ2], &nlpr->diuuv[kNQ2],
             &nlpr->diuvv[kNQ2], &nlpr->divvv[kNQ2],
             &bc00[(G1_CROSS00DEG+1)*k], &bc01[(G1_CROSS01DEG+1)*k],
             &bc10[(G1_CROSS10DEG+1)*k], &bc11[(G1_CROSS11DEG+1)*k],
             &bd00[(G1_CROSS00DEG+1)*k], &bd01[(G1_CROSS01DEG+1)*k],
             &bd10[(G1_CROSS10DEG+1)*k], &bd11[(G1_CROSS11DEG+1)*k],
             &nlpr->psiu[fN], &nlpr->psiv[fN],
             &nlpr->psiuu[fN], &nlpr->psiuv[fN], &nlpr->psivv[fN],
             &nlpr->psiuuu[fN], &nlpr->psiuuv[fN],
             &nlpr->psiuvv[fN], &nlpr->psivvv[fN] ) ) {
      domain->error_code = G1H_ERROR_NO_SCRATCH_MEMORY;
      goto failure;
    }
  }

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
#undef N
} /*_g1hq2_TabNLBasisFunctionsOmegaf*/

boolean _g1hq2_TabNLBasisFunctionsGammaf ( GHoleDomainf *domain, int nkn,
                                           G1HNLPrivatef *nlpr, float *ctrd,
                                           float *bc00 )
{
#define N ((G1H_FINALDEG+1)*(G1H_FINALDEG+1))
  void     *sp;
  GHolePrivateRecf  *privateG;
  G1HolePrivateRecf *privateG1;
  float    *tkn, *hfunc, *dhfunc, *ddhfunc, *dddhfunc,
           *atkn, *ahfunc, *adhfunc, *addhfunc;
  float    *ec00, *ec01, *ed00, *ed01, *fc00, *fc01, *fd00, *fd01,
           *bc01, *bc10, *bc11, *bd00, *bd01, *bd10, *bd11;
  int      i, j, k, kk, kN, kNQ2, hole_k, f, fN, nfunc_a, nfunc_b;
  point2f  *di;
  float    bvz;
  unsigned char *bfcpn;
  vector2f cdi, cdiu, cdiv, cdiuu, cdiuv, cdivv, *sicp;
  float    *tabeu, *tabev, *tabeuu, *tabeuv, *tabevv,
           *tabfu, *tabfv, *tabfuu, *tabfuv, *tabfvv;
  float    *A11, *A21, *A22;

  sp      = pkv_GetScratchMemTop ();
  hole_k  = domain->hole_k;
  privateG  = domain->privateG;
  privateG1 = domain->privateG1;
  nfunc_a = privateG1->nfunc_a;
  nfunc_b = privateG1->nfunc_b;
  di      = pkv_GetScratchMem ( N*sizeof(point2f) );
  tkn     = pkv_GetScratchMemf ( nkn );
  hfunc = pkv_GetScratchMemf ( 16*nkn );
  if ( !tkn || !di || !hfunc ) {
    domain->error_code = G1H_ERROR_NO_SCRATCH_MEMORY;
    goto failure;
  }

  dhfunc = &hfunc[4*nkn];
  ddhfunc = &dhfunc[4*nkn];
  dddhfunc = &ddhfunc[4*nkn];

  _gh_PrepareTabKnotsf ( nkn, privateG1->opt_quad, tkn );
  mbs_TabCubicHFuncDer3f ( 0.0, 1.0, nkn, tkn,
                           hfunc, dhfunc, ddhfunc, dddhfunc );

        /* now evaluate the 1st order derivatives and jumps of 2nd order */
        /* derivatives at the boundary curves of the areas Omega_i */
  sicp = (point2f*)pkv_GetScratchMem ( 16*sizeof(point2f) );
  tabeu = pkv_GetScratchMemf ( 20*nkn );
  atkn = pkv_GetScratchMemf ( 26 );
  if ( !sicp || !tabeu || !atkn )
    goto failure;

  tabev = &tabeu[2*nkn];    tabeuu = &tabev[2*nkn];   tabeuv = &tabeuu[2*nkn];
  tabevv = &tabeuv[2*nkn];  tabfu = &tabevv[2*nkn];   tabfv = &tabfu[2*nkn];
  tabfuu = &tabfv[2*nkn];   tabfuv = &tabfuu[2*nkn];  tabfvv = &tabfuv[2*nkn];

  ahfunc = &atkn[2];  adhfunc = &ahfunc[8];  addhfunc = &adhfunc[8];
  atkn[0] = 0.0;  atkn[1] = 1.0;
  mbs_TabCubicHFuncDer2f ( 0.0, 1.0, 2, atkn, ahfunc, adhfunc, addhfunc );

  memset ( ctrd, 0, 3*nkn*38*hole_k*sizeof(float) );
  for ( k = kN = kNQ2 = 0, kk = 1;
        k < hole_k;
        k++, kk = (k+1) % hole_k, kN += N, kNQ2 += 3*nkn ) {
    pkv_Selectf ( N, 2, 3, 2, &nlpr->nldi[kN], di );
    for ( i = 0; i < nkn; i++ ) {
      if ( !mbs_BCHornerDer2Pf ( G1H_FINALDEG, G1H_FINALDEG, 2, (float*)di,
                                 tkn[i], 0.0, &cdi.x, &cdiu.x, &cdiv.x,
                                 &cdiuu.x, &cdiuv.x, &cdivv.x ) )
        goto failure;
      nlpr->ctang[kNQ2+i] = cdiu;
      _g1hq2_SetupCTrdf ( &cdiu, &cdiv, &cdiuu, &cdiuv, &cdivv,
                          &ctrd[38*(kNQ2+i)] );
      if ( !mbs_BCHornerDer2Pf ( G1H_FINALDEG, G1H_FINALDEG, 2, (float*)di,
                                 tkn[i], 1.0, &cdi.x, &cdiu.x, &cdiv.x,
                                 &cdiuu.x, &cdiuv.x, &cdivv.x ) )
        goto failure;
      nlpr->ctang[kNQ2+nkn+i] = cdiu;
      _g1hq2_SetupCTrdf ( &cdiu, &cdiv, &cdiuu, &cdiuv, &cdivv,
                          &ctrd[38*(kNQ2+nkn+i)] );
      if ( !mbs_BCHornerDer2Pf ( G1H_FINALDEG, G1H_FINALDEG, 2, (float*)di,
                                 1.0, tkn[i], &cdi.x, &cdiu.x, &cdiv.x,
                                 &cdiuu.x, &cdiuv.x, &cdivv.x ) )
        goto failure;
      nlpr->ctang[kNQ2+2*nkn+i] = cdiv;
      _g1hq2_SetupCTrdf ( &cdiu, &cdiv, &cdiuu, &cdiuv, &cdivv,
                          &ctrd[38*(kNQ2+2*nkn+i)] );
      if ( !mbs_BCHornerDer2Pf ( G1H_FINALDEG, G1H_FINALDEG, 2, (float*)di,
                                 0.0, tkn[i], &cdi.x, &cdiu.x, &cdiv.x,
                                 &cdiuu.x, &cdiuv.x, &cdivv.x ) )
        goto failure;
      _g1hq2_SetupCTrdf ( &cdiu, &cdiv, &cdiuu, &cdiuv, &cdivv,
                          &ctrd[38*(kk*3*nkn+i)+19] );
    }

    if ( !_g1hq2_FindDomSurrndPatchf ( domain, nlpr, kk, 1, sicp ) )
      goto failure;
    for ( i = 0; i < nkn; i++ ) {
      if ( !mbs_BCHornerDer2Pf ( 3, 3, 2, (float*)sicp,
                                 0.0, tkn[i], &cdi.x, &cdiu.x, &cdiv.x,
                                 &cdiuu.x, &cdiuv.x, &cdivv.x ) )
        goto failure;
      _g1hq2_SetupCTrdf ( &cdiu, &cdiv, &cdiuu, &cdiuv, &cdivv,
                          &ctrd[38*(kNQ2+nkn+i)+19] );
    }
    if ( !_g1hq2_FindDomSurrndPatchf ( domain, nlpr, k, 2, sicp ) )
      goto failure;
    for ( i = 0; i < nkn; i++ ) {
      if ( !mbs_BCHornerDer2Pf ( 3, 3, 2, (float*)sicp,
                                 0.0, tkn[i], &cdi.x, &cdiu.x, &cdiv.x,
                                 &cdiuu.x, &cdiuv.x, &cdivv.x ) )
        goto failure;
      _g1hq2_SetupCTrdf ( &cdiu, &cdiv, &cdiuu, &cdiuv, &cdivv,
                          &ctrd[38*(kNQ2+2*nkn+i)+19] );
    }
  }

  memset ( nlpr->cpsiu, 0, 3*nkn*(nfunc_a+1)*hole_k*sizeof(float) );
  memset ( nlpr->cpsiv, 0, 3*nkn*(nfunc_a+1)*hole_k*sizeof(float) );
  memset ( nlpr->cpsiuu, 0, 3*nkn*(nfunc_a+1)*hole_k*sizeof(float) );
  memset ( nlpr->cpsiuv, 0, 3*nkn*(nfunc_a+1)*hole_k*sizeof(float) );
  memset ( nlpr->cpsivv, 0, 3*nkn*(nfunc_a+1)*hole_k*sizeof(float) );
        /* block A functions */
  for ( f = 0; f < nfunc_a; f++ ) {
    _g1h_GetBFAPatchCurvesf ( domain, f, hole_k-1,
                              &ec00, &ec01, &ed00, &ed01 );
    for ( k = 0; k < hole_k; k++ ) {
      _g1h_GetBFAPatchCurvesf ( domain, f, k,
                                &fc00, &fc01, &fd00, &fd01 );
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
      for ( i = j = 0, fN = 3*(f*hole_k+k)*nkn;
            i < nkn;
            i++, j += 2, fN++ ) {
        A11 = &ctrd[38*(3*k*nkn+i)];  A21 = &A11[4];  A22 = &A21[6];
        nlpr->cpsiu[fN] = A11[0]*tabfu[j] + A11[1]*tabfv[j];
        nlpr->cpsiv[fN] = A11[2]*tabfu[j] + A11[3]*tabfv[j];
        nlpr->cpsiuu[fN] = A21[0]*tabfu[j] + A21[1]*tabfv[j] +
                           A22[0]*tabfuu[j] + A22[1]*tabfuv[j] + A22[2]*tabfvv[j];
        nlpr->cpsiuv[fN] = A21[2]*tabfu[j] + A21[3]*tabfv[j] +
                           A22[3]*tabfuu[j] + A22[4]*tabfuv[j] + A22[5]*tabfvv[j];
        nlpr->cpsivv[fN] = A21[4]*tabfu[j] + A21[5]*tabfv[j] +
                           A22[6]*tabfuu[j] + A22[7]*tabfuv[j] + A22[8]*tabfvv[j];
        A11 = &ctrd[38*(3*nkn*k+i)+19];  A21 = &A11[4];  A22 = &A21[6];
        nlpr->cpsiuu[fN] -= A21[0]*tabeu[i] + A21[1]*tabev[i] +
                            A22[0]*tabeuu[i] + A22[1]*tabeuv[i] + A22[2]*tabevv[i];
        nlpr->cpsiuv[fN] -= A21[2]*tabeu[i] + A21[3]*tabev[i] +
                            A22[3]*tabeuu[i] + A22[4]*tabeuv[i] + A22[5]*tabevv[i];
        nlpr->cpsivv[fN] -= A21[4]*tabeu[i] + A21[5]*tabev[i] +
                            A22[6]*tabeuu[i] + A22[7]*tabeuv[i] + A22[8]*tabevv[i];
      }
      for ( i = 0, j = 1, fN = (3*(f*hole_k+k)+1)*nkn;
            i < nkn;
            i++, j += 2, fN++ ) {
        A11 = &ctrd[38*((3*k+1)*nkn+i)];  A21 = &A11[4];  A22 = &A21[6];
        nlpr->cpsiu[fN] = A11[0]*tabfu[j] + A11[1]*tabfv[j];
        nlpr->cpsiv[fN] = A11[2]*tabfu[j] + A11[3]*tabfv[j];
        nlpr->cpsiuu[fN] = -(A21[0]*tabfu[j] + A21[1]*tabfv[j] +
                           A22[0]*tabfuu[j] + A22[1]*tabfuv[j] + A22[2]*tabfvv[j]);
        nlpr->cpsiuv[fN] = -(A21[2]*tabfu[j] + A21[3]*tabfv[j] +
                           A22[3]*tabfuu[j] + A22[4]*tabfuv[j] + A22[5]*tabfvv[j]);
        nlpr->cpsivv[fN] = -(A21[4]*tabfu[j] + A21[5]*tabfv[j] +
                           A22[6]*tabfuu[j] + A22[7]*tabfuv[j] + A22[8]*tabfvv[j]);
      }
      if ( !mbs_TabBezC1Coons0Der2f ( 1,
                  1, &atkn[1], &ahfunc[4], &adhfunc[4], &addhfunc[4],
                  nkn, tkn, hfunc, dhfunc, ddhfunc,
                  G1_CROSS00DEG, fc00, G1_CROSS01DEG, fc01,
                  G1_CROSS00DEG, fd00, G1_CROSS01DEG, fd01,
                  NULL, tabfu, tabfv, tabfuu, tabfuv, tabfvv ) )
        goto failure;
      for ( i = 0, fN = (3*(f*hole_k+k)+2)*nkn;
            i < nkn;
            i++, fN++ ) {
        A11 = &ctrd[38*((3*k+2)*nkn+i)];  A21 = &A11[4];  A22 = &A21[6];
        nlpr->cpsiu[fN] = A11[0]*tabfu[i] + A11[1]*tabfv[i];
        nlpr->cpsiv[fN] = A11[2]*tabfu[i] + A11[3]*tabfv[i];
        nlpr->cpsiuu[fN] = -(A21[0]*tabfu[i] + A21[1]*tabfv[i] +
                           A22[0]*tabfuu[i] + A22[1]*tabfuv[i] + A22[2]*tabfvv[i]);
        nlpr->cpsiuv[fN] = -(A21[2]*tabfu[i] + A21[3]*tabfv[i] +
                           A22[3]*tabfuu[i] + A22[4]*tabfuv[i] + A22[5]*tabfvv[i]);
        nlpr->cpsivv[fN] = -(A21[4]*tabfu[i] + A21[5]*tabfv[i] +
                           A22[6]*tabfuu[i] + A22[7]*tabfuv[i] + A22[8]*tabfvv[i]);
      }
      ec00 = fc00;  ec01 = fc01;  ed00 = fd00;  ed01 = fd01;
    }
  }
        /* fixed linear combination of the block B functions */
          /* the Coons patches representing this function in Omega */
          /* are already found, so the algorithm starts from */
          /* the derivatives inside Omega */
  bc01 = &bc00[(G1_CROSS00DEG+1)*hole_k];  bc10 = &bc01[(G1_CROSS01DEG+1)*hole_k];
  bc11 = &bc10[(G1_CROSS10DEG+1)*hole_k];  bd00 = &bc11[(G1_CROSS11DEG+1)*hole_k];
  bd01 = &bd00[(G1_CROSS00DEG+1)*hole_k];  bd10 = &bd01[(G1_CROSS01DEG+1)*hole_k];
  bd11 = &bd10[(G1_CROSS10DEG+1)*hole_k];
  for ( k = 0, kk = 1;  k < hole_k;  k++, kk = (k+1) % hole_k ) {
    if ( !mbs_TabBezC1CoonsDer2f ( 1, nkn, tkn, hfunc, dhfunc, ddhfunc,
                2, atkn, ahfunc, adhfunc, addhfunc,
                G1_CROSS00DEG, &bc00[(G1_CROSS00DEG+1)*k],
                G1_CROSS01DEG, &bc01[(G1_CROSS01DEG+1)*k],
                G1_CROSS10DEG, &bc10[(G1_CROSS10DEG+1)*k],
                G1_CROSS11DEG, &bc11[(G1_CROSS11DEG+1)*k],
                G1_CROSS00DEG, &bd00[(G1_CROSS00DEG+1)*k],
                G1_CROSS01DEG, &bd01[(G1_CROSS01DEG+1)*k],
                G1_CROSS10DEG, &bd10[(G1_CROSS10DEG+1)*k],
                G1_CROSS11DEG, &bd11[(G1_CROSS11DEG+1)*k],
                NULL, tabfu, tabfv, tabfuu, tabfuv, tabfvv ) )
      goto failure;
    for ( i = j = 0, fN = 3*(nfunc_a*hole_k+k)*nkn;
          i < nkn;
          i++, j += 2, fN++ ) {
      A11 = &ctrd[38*(3*k*nkn+i)];  A21 = &A11[4];  A22 = &A21[6];
      nlpr->cpsiu[fN] = A11[0]*tabfu[j] + A11[1]*tabfv[j];
      nlpr->cpsiv[fN] = A11[2]*tabfu[j] + A11[3]*tabfv[j];
      nlpr->cpsiuu[fN] += A21[0]*tabfu[j] + A21[1]*tabfv[j] +
                          A22[0]*tabfuu[j] + A22[1]*tabfuv[j] + A22[2]*tabfvv[j];
      nlpr->cpsiuv[fN] += A21[2]*tabfu[j] + A21[3]*tabfv[j] +
                          A22[3]*tabfuu[j] + A22[4]*tabfuv[j] + A22[5]*tabfvv[j];
      nlpr->cpsivv[fN] += A21[4]*tabfu[j] + A21[5]*tabfv[j] +
                          A22[6]*tabfuu[j] + A22[7]*tabfuv[j] + A22[8]*tabfvv[j];
    }
    for ( i = 0, j = 1, fN = (3*(nfunc_a*hole_k+k)+1)*nkn;
          i < nkn;
          i++, j += 2, fN++ ) {
      A11 = &ctrd[38*((3*k+1)*nkn+i)];  A21 = &A11[4];  A22 = &A21[6];
      nlpr->cpsiu[fN] = A11[0]*tabfu[j] + A11[1]*tabfv[j];
      nlpr->cpsiv[fN] = A11[2]*tabfu[j] + A11[3]*tabfv[j];
      nlpr->cpsiuu[fN] = -(A21[0]*tabfu[j] + A21[1]*tabfv[j] +
                          A22[0]*tabfuu[j] + A22[1]*tabfuv[j] + A22[2]*tabfvv[j]);
      nlpr->cpsiuv[fN] = -(A21[2]*tabfu[j] + A21[3]*tabfv[j] +
                          A22[3]*tabfuu[j] + A22[4]*tabfuv[j] + A22[5]*tabfvv[j]);
      nlpr->cpsivv[fN] = -(A21[4]*tabfu[j] + A21[5]*tabfv[j] +
                          A22[6]*tabfuu[j] + A22[7]*tabfuv[j] + A22[8]*tabfvv[j]);
    }
    if ( !mbs_TabBezC1CoonsDer2f ( 1, 2, atkn, ahfunc, adhfunc, addhfunc,
                nkn, tkn, hfunc, dhfunc, ddhfunc,
                G1_CROSS00DEG, &bc00[(G1_CROSS00DEG+1)*k],
                G1_CROSS01DEG, &bc01[(G1_CROSS01DEG+1)*k],
                G1_CROSS10DEG, &bc10[(G1_CROSS10DEG+1)*k],
                G1_CROSS11DEG, &bc11[(G1_CROSS11DEG+1)*k],
                G1_CROSS00DEG, &bd00[(G1_CROSS00DEG+1)*k],
                G1_CROSS01DEG, &bd01[(G1_CROSS01DEG+1)*k],
                G1_CROSS10DEG, &bd10[(G1_CROSS10DEG+1)*k],
                G1_CROSS11DEG, &bd11[(G1_CROSS11DEG+1)*k],
                NULL, tabfu, tabfv, tabfuu, tabfuv, tabfvv ) )
      goto failure;
    for ( i = 0, fN = 3*(nfunc_a*hole_k+kk)*nkn;
          i < nkn;
          i++, fN++ ) {
      A11 = &ctrd[38*(3*kk*nkn+i)+19];  A21 = &A11[4];  A22 = &A21[6];
      nlpr->cpsiuu[fN] -= A21[0]*tabfu[i] + A21[1]*tabfv[i] +
                          A22[0]*tabfuu[i] + A22[1]*tabfuv[i] + A22[2]*tabfvv[i];
      nlpr->cpsiuv[fN] -= A21[2]*tabfu[i] + A21[3]*tabfv[i] +
                          A22[3]*tabfuu[i] + A22[4]*tabfuv[i] + A22[5]*tabfvv[i];
      nlpr->cpsivv[fN] -= A21[4]*tabfu[i] + A21[5]*tabfv[i] +
                          A22[6]*tabfuu[i] + A22[7]*tabfuv[i] + A22[8]*tabfvv[i];
    }
    for ( i = 0, j = nkn, fN = (3*(nfunc_a*hole_k+k)+2)*nkn;
          i < nkn;
          i++, j++, fN++ ) {
      A11 = &ctrd[38*((3*k+2)*nkn+i)];  A21 = &A11[4];  A22 = &A21[6];
      nlpr->cpsiu[fN] = A11[0]*tabfu[j] + A11[1]*tabfv[j];
      nlpr->cpsiv[fN] = A11[2]*tabfu[j] + A11[3]*tabfv[j];
      nlpr->cpsiuu[fN] = -(A21[0]*tabfu[j] + A21[1]*tabfv[j] +
                          A22[0]*tabfuu[j] + A22[1]*tabfuv[j] + A22[2]*tabfvv[j]);
      nlpr->cpsiuv[fN] = -(A21[2]*tabfu[j] + A21[3]*tabfv[j] +
                          A22[3]*tabfuu[j] + A22[4]*tabfuv[j] + A22[5]*tabfvv[j]);
      nlpr->cpsivv[fN] = -(A21[4]*tabfu[j] + A21[5]*tabfv[j] +
                          A22[6]*tabfuu[j] + A22[7]*tabfuv[j] + A22[8]*tabfvv[j]);
    }
  }
          /* now deal with the derivatives outside Omega, */
          /* the arrays bc00 .. bd11 may be reused */
  bfcpn = privateG->bfcpn;
  memset ( bc00, 0, 2*(G1_CROSSDEGSUM+4)*hole_k*sizeof(float) );
  for ( f = 0; f < nfunc_b; f++ ) {
    bvz = nlpr->rhole_cp[bfcpn[f]].z;
    for ( k = 0, kk = 1;  k < hole_k;  k++, kk = (k+1) % hole_k ) {
      gh_GetDomSurrndBFuncf ( domain, f, kk, 1, (float*)di );
      pkn_AddMatrixMf ( 1, 16, 0, &bc00[k*2*16], 0, (float*)di, bvz,
                        0, &bc00[k*2*16] );
      gh_GetDomSurrndBFuncf ( domain, f, k, 2, (float*)di );
      pkn_AddMatrixMf ( 1, 16, 0, &bc00[(k*2+1)*16], 0, (float*)di, bvz,
                        0, &bc00[(k*2+1)*16] );
    }
  }
  for ( k = 0; k < hole_k; k++ ) {
    for ( i = 0, fN = (3*(nfunc_a*hole_k+k)+1)*nkn;  i < nkn;  i++, fN++ ) {
      if ( !mbs_BCHornerDer2Pf ( 3, 3, 1, &bc00[k*2*16], 0.0, tkn[i],
                                 &tabeu[1], tabeu, tabev, tabeuu, tabeuv, tabevv ) )
        goto failure;
      A11 = &ctrd[38*((3*k+1)*nkn+i)+19];  A21 = &A11[4];  A22 = &A21[6];
      nlpr->cpsiuu[fN] += A21[0]*tabeu[0] + A21[1]*tabev[0] +
                          A22[0]*tabeuu[0] + A22[1]*tabeuv[0] + A22[2]*tabevv[0];
      nlpr->cpsiuv[fN] += A21[2]*tabeu[0] + A21[3]*tabev[0] +
                          A22[3]*tabeuu[0] + A22[4]*tabeuv[0] + A22[5]*tabevv[0];
      nlpr->cpsivv[fN] += A21[4]*tabeu[0] + A21[5]*tabev[0] +
                          A22[6]*tabeuu[0] + A22[7]*tabeuv[0] + A22[8]*tabevv[0];
    }
    for ( i = 0, fN = (3*(nfunc_a*hole_k+k)+2)*nkn;  i < nkn;  i++, fN++ ) {
      if ( !mbs_BCHornerDer2Pf ( 3, 3, 1, &bc00[(k*2+1)*16], 0.0, tkn[i],
                                 &tabeu[1], tabeu, tabev, tabeuu, tabeuv, tabevv ) )
        goto failure;
      A11 = &ctrd[38*((3*k+2)*nkn+i)+19];  A21 = &A11[4];  A22 = &A21[6];
      nlpr->cpsiuu[fN] += A21[0]*tabeu[0] + A21[1]*tabev[0] +
                          A22[0]*tabeuu[0] + A22[1]*tabeuv[0] + A22[2]*tabevv[0];
      nlpr->cpsiuv[fN] += A21[2]*tabeu[0] + A21[3]*tabev[0] +
                          A22[3]*tabeuu[0] + A22[4]*tabeuv[0] + A22[5]*tabevv[0];
      nlpr->cpsivv[fN] += A21[4]*tabeu[0] + A21[5]*tabev[0] +
                          A22[6]*tabeuu[0] + A22[7]*tabeuv[0] + A22[8]*tabevv[0];
    }
  }

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
#undef N
} /*_g1hq2_TabNLBasisFunctionsGammaf*/

/* ///////////////////////////////////////////////////////////////////////// */
void _g1hq2_IntFunc1af ( G1Q2HNLFuncf *f, float *funct )
{
  float e0, e1, e2, e9, e0p2, pup2, pvp2;
  float A, B, x, y, xx, xy, yy, tGt;

  pup2 = f->pu*f->pu;
  pvp2 = f->pv*f->pv;
  f->e1 = e1 = (float)(1.0 + pup2);
  f->e2 = e2 = (float)(1.0 + pvp2);
  f->e0 = e0 = e1 + pvp2;
  e0p2 = e0*e0;
  f->e9 = e9 = f->pu*f->pv;
  x = f->tang.x;
  y = f->tang.y;
  xx = x*x;  xy = x*y;  yy = y*y;
  f->tGt = tGt = (float)(e1*xx + 2.0*e9*xy + e2*yy);
  f->A = A = (float)(0.5*(e2*f->jpuu+e1*f->jpvv) - e9*f->jpuv);
  f->B = B = (float)(sqrt(tGt)/(e0*e0p2));
  *funct += A*A*B;
} /*_g1hq2_IntFunc1af*/

void _g1hq2_IntFunc1bf ( G1Q2HNLFuncf *f, float *funct )
{
  float e0, e1, e2, e9, e0p2, e0p4, pup2, pvp2;
  float A, B, x, y, xx, xy, yy, tGt;
  float b11, b12, b22, dd;

  pup2 = f->pu*f->pu;
  pvp2 = f->pv*f->pv;
  f->e1 = e1 = (float)(1.0 + pup2);
  f->e2 = e2 = (float)(1.0 + pvp2);
  f->e0 = e0 = e1 + pvp2;
  e0p2 = e0*e0;
  f->e9 = e9 = f->pu*f->pv;
  x = f->tang.x;
  y = f->tang.y;
  xx = x*x;  xy = x*y;  yy = y*y;
  f->tGt = tGt = (float)(e1*xx + 2.0*e9*xy + e2*yy);
  f->A = A = (float)(0.5*(e2*f->jpuu+e1*f->jpvv) - e9*f->jpuv);
  f->B = B = (float)(sqrt(tGt)/(e0*e0p2));
  *funct += A*A*B;
        /* data needed for the computation of derivatives of B */
  e0p4 = e0p2*e0p2;
  dd = (float)(1.0/(e0p4*sqrt(tGt)));
  b11 = (float)(f->pu*(pvp2-5.0*e1));
  b12 = (float)(0.5*f->pv*(e2-11.0*pup2));
  b22 = (float)(-6.0*f->pu*e2);
  f->b1 = (float)((b11*xx + 2.0*b12*xy + b22*yy)*dd);
  b11 = (float)(-6.0*f->pv*e1);
  b12 = (float)(0.5*f->pu*(e1-11.0*pvp2));
  b22 = (float)(f->pv*(pup2-5.0*e2));
  f->b2 = (float)((b11*xx + 2.0*b12*xy + b22*yy)*dd);
} /*_g1hq2_IntFunc1bf*/

void _g1hq2_IntFunc1cf ( G1Q2HNLFuncf *f, float *funct )
{
  float e0, e1, e2, e9, e0p2, e0p4, e1p2, e2p2, e9p2, pup2, pvp2;
  float A, B, x, y, xx, xy, yy, xxxx, xxxy, xxyy, xyyy, yyyy, tGt;
  float b11, b12, b22, b23, b33, dd;

  pup2 = f->pu*f->pu;
  pvp2 = f->pv*f->pv;
  f->e1 = e1 = (float)(1.0 + pup2);
  f->e2 = e2 = (float)(1.0 + pvp2);
  f->e0 = e0 = e1 + pvp2;
  e0p2 = e0*e0;
  f->e9 = e9 = f->pu*f->pv;
  x = f->tang.x;
  y = f->tang.y;
  xx = x*x;  xy = x*y;  yy = y*y;
  f->tGt = tGt = (float)(e1*xx + 2.0*e9*xy + e2*yy);
  f->A = A = (float)(0.5*(e2*f->jpuu+e1*f->jpvv) - e9*f->jpuv);
  f->B = B = (float)(sqrt(tGt)/(e0*e0p2));
  *funct += A*A*B;
        /* data needed for the computation of derivatives of B */
  e0p4 = e0p2*e0p2;
  dd = (float)(1.0/(e0p4*sqrt(tGt)));
  b11 = (float)(f->pu*(pvp2-5.0*e1));
  b12 = (float)(0.5*f->pv*(e2-11.0*pup2));
  b22 = (float)(-6.0*f->pu*e2);
  f->b1 = (float)((b11*xx + 2.0*b12*xy + b22*yy)*dd);
  b11 = (float)(-6.0*f->pv*e1);
  b12 = (float)(0.5*f->pu*(e1-11.0*pvp2));
  b22 = (float)(f->pv*(pup2-5.0*e2));
  f->b2 = (float)((b11*xx + 2.0*b12*xy + b22*yy)*dd);
        /* data for the derivatives of B of the second order */  
  e1p2 = e1*e1;
  e2p2 = e2*e2;
  e9p2 = e9*e9;
  dd = (float)(1.0/(e0p4*e0*tGt*sqrt(tGt)));
  xxxx = xx*xx;
  xxxy = xx*xy;
  xxyy = xx*yy;
  xyyy = xy*yy;
  yyyy = yy*yy;
  b11 = (float)(pup2*((30.0*pup2+55.0)*pup2-18.0*e9p2-22.0*pvp2+20.0) + (pvp2-5.0)*e2);
  b12 = (float)(e9*(pup2*(66.0*pup2+48.0-30.0*pvp2)-18.0*e2));
  b22 = (float)(e9p2*(216.0*pup2-72.0*pvp2-10.0)+(73.0*pup2+62.0)*pup2-11.0*e2p2);
  b23 = (float)(e9*e2*(78.0*pup2-18.0*e2));
  b33 = (float)(6.0*e2p2*(7.0*pup2-e2));
  f->b3 = (float)((b11*xxxx + 2.0*b12*xxxy + b22*xxyy + 2.0*b23*xyyy + b33*yyyy)*dd);
  b11 = (float)(6.0*e9*e1*(7.0*e1-pvp2));
  b12 = (float)(0.5*(pup2*(pvp2*(168.0*pup2-18.0*pvp2+164.0)-(6.0*pup2+11.0)*pup2-4.0)+
             (1.0-5.0*pvp2)*e2));
  b22 = (float)(e9*(252.0*e9p2-18.0*(pup2*pup2+pvp2*pvp2)+66.0*(pup2+pvp2)+84.0));
  b23 = (float)(0.5*(pvp2*(pup2*(168.0*pvp2-18.0*pup2+164.0)-(6.0*pvp2+11.0)*pvp2-4.0)+
             (1.0-5.0*pup2)*e1));
  b33 = (float)(6.0*e9*e2*(7.0*e2-pup2));
  f->b4 = (float)((b11*xxxx + 2.0*b12*xxxy + b22*xxyy + 2.0*b23*xyyy + b33*yyyy)*dd);
  b11 = (float)(6.0*e1p2*(7.0*pvp2-e1));
  b12 = (float)(e9*e1*(78.0*pvp2-18.0*e1));
  b22 = (float)(e9p2*(216.0*pvp2-72.0*pup2-10.0)+(73.0*pvp2+62.0)*pvp2-11.0*e1);
  b23 = (float)(e9*(pvp2*(66.0*pvp2+48.0-30.0*pup2)-18.0*e1));
  b33 = (float)(pvp2*((30.0*pvp2+55.0)*pvp2-18.0*e9p2-22.0*pup2+20.0) + (pup2-5.0)*e1);
  f->b5 = (float)((b11*xxxx + 2.0*b12*xxxy + b22*xxyy + 2.0*b23*xyyy + b33*yyyy)*dd);
} /*_g1hq2_IntFunc1cf*/

void _g1hq2_IntFunc2bf ( G1Q2HNLFuncf *f, float *grad )
{
  float A_i, B_i;
  float e1, e2, e9;

  e1 = f->e1;  e2 = f->e2;  e9 = f->e9;
  A_i = (float)((f->pu*f->jpvv-f->pv*f->jpuv)*f->psiu
       -(f->pu*f->jpuv-f->pv*f->jpuu)*f->psiv
       +0.5*(e2*f->jpsiuu+e1*f->jpsivv)-e9*f->jpsiuv);
  B_i = f->b1*f->psiu + f->b2*f->psiv;
  *grad += (float)(f->A*(2.0*A_i*f->B + f->A*B_i));
} /*_g1hq2_IntFunc2bf*/

void _g1hq2_IntFunc2cf ( G1Q2HNLFuncf *f, float *Ai, float *Bi, float *grad )
{
  float A_i, B_i;
  float e1, e2, e9;

  e1 = f->e1;  e2 = f->e2;  e9 = f->e9;
  *Ai = A_i = (float)((f->pu*f->jpvv-f->pv*f->jpuv)*f->psiu
             -(f->pu*f->jpuv-f->pv*f->jpuu)*f->psiv
             +0.5*(e2*f->jpsiuu+e1*f->jpsivv)-e9*f->jpsiuv);
  *Bi = B_i = f->b1*f->psiu + f->b2*f->psiv;
  *grad += (float)(f->A*(2.0*A_i*f->B + f->A*B_i));
} /*_g1hq2_IntFunc2cf*/

void _g1hq2_IntFunc3cf ( G1Q2HNLFuncf *f, float Ai, float Bi, float Aj, float Bj,
                        float *hessian )
{
  float Aij, Bij;

  Aij = f->pu*(f->psiu*f->jpsjvv+f->jpsivv*f->psju-f->psiv*f->jpsjuv-f->jpsiuv*f->psjv)
       +f->pv*(f->psiv*f->jpsjuu+f->jpsiuu*f->psjv-f->psiu*f->jpsjuv-f->jpsiuv*f->psju)
       +f->jpuu*f->psiv*f->psjv
       -f->jpuv*(f->psiu*f->psjv+f->psiv*f->psju)
       +f->jpvv*f->psiu*f->psju;
  Bij = f->b3*f->psiu*f->psju
       +f->b4*(f->psiu*f->psjv+f->psiv*f->psju)
       +f->b5*f->psiv*f->psjv;
  *hessian += (float)(2.0*(f->A*(Aij*f->B+Ai*Bj+Aj*Bi)+ Ai*Aj*f->B) + f->A*f->A*Bij);
} /*_g1hq2_IntFunc3cf*/

