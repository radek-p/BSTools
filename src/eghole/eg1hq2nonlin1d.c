
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2008, 2015                            */
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

#include "eg1holed.h"
#include "eg2holed.h"
#include "eg1hprivated.h"
#include "eg2hprivated.h"
#include "eg1herror.h"


/* ///////////////////////////////////////////////////////////////////////// */
boolean _g1hq2_FindDomSurrndPatchd ( GHoleDomaind *domain,
                                     G1HNLPrivated *nlpr,
                                     int i, int j, point2d *bezcp )
{
  void    *sp;
  int     hole_k, k;
  int     *ind;
  point2d *q;
  double  *ukn, *vkn;

  sp = pkv_GetScratchMemTop ();
  ind = (int*)pkv_GetScratchMem ( 16*sizeof(int) );
  q = (point2d*)pkv_GetScratchMem ( 16*sizeof(point2d) );
  if ( !ind || !q ) {
    domain->error_code = G1H_ERROR_NO_SCRATCH_MEMORY;
    goto failure;
  }

  hole_k = domain->hole_k;
  ukn = _gh_GetKnotSequenced ( domain, i-1 );  ukn += 3;
  vkn = _gh_GetKnotSequenced ( domain, i );    vkn += j;
  gh_GetBspInd ( hole_k, i, j, ind );
  for ( k = 0; k < 16; k++ )
    memcpy ( &q[k], &nlpr->rhole_cp[ind[k]], sizeof(point2d) );
  if ( !mbs_BSPatchToBezd ( 2, 3, 7, ukn, 3, 7, vkn, 8, (double*)q,
                            NULL, NULL, NULL, NULL, NULL, NULL, 8, (double*)bezcp ) )
    goto failure;
  if ( j == 1 )
    mbs_multiReverseBSCurved ( 3, 0, NULL, 4, 2, 8, (double*)bezcp );

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*_g1hq2_FindDomSurrndPatchd*/

boolean _g1hq2_FindNLDomainDiameterd ( GHoleDomaind *domain,
                                       G1HNLPrivated *nlpr )
{
  void     *sp;
  int      hole_k, i, j, k;
  point2d  *p, *q;
  vector2d v;
  double   ddiam, d;

  sp = pkv_GetScratchMemTop ();
  hole_k = domain->hole_k;
  p = pkv_GetScratchMem ( 6*hole_k*sizeof(point2d) );
  q = pkv_GetScratchMem ( 16*sizeof(point2d) );
  if ( !p || !q )
    goto failure;
  for ( i = k = 0;  i < hole_k;  i++ ) {
    _g1hq2_FindDomSurrndPatchd ( domain, nlpr, i, 1, q );
    p[k++] = q[0];  p[k++] = q[1];  p[k++] = q[2];
    _g1hq2_FindDomSurrndPatchd ( domain, nlpr, i, 2, q );
    p[k++] = q[1];  p[k++] = q[2];  p[k++] = q[3];
  }
  ddiam = 0.0;
  for ( i = 1; i < k; i++ )
    for ( j = 0; j < i; j++ ) {
      SubtractPoints2d ( &p[i], &p[j], &v );
      d = v.x*v.x+v.y*v.y;
      ddiam = max ( ddiam, d );
    }
  nlpr->ddiam = sqrt(ddiam);

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*_g1hq2_FindNLDomainDiameterd*/

/* ///////////////////////////////////////////////////////////////////////// */
boolean _g1hq2_TabNLDer0d ( int nkn, const double *tkn,
             const double *hfunc, const double *dhfunc, const double *ddhfunc,
             const double *dddhfunc,
             vector2d *diu, vector2d *div,
             vector2d *diuu, vector2d *diuv, vector2d *divv,
             vector2d *diuuu, vector2d *diuuv, vector2d *diuvv, vector2d *divvv,
             double *fc00, double *fc01, double *fd00, double *fd01,
             double *psiu, double *psiv,
             double *psiuu, double *psiuv, double *psivv,
             double *psiuuu, double *psiuuv, double *psiuvv, double *psivvv )
{
  void   *sp;
  int    i;
  double *hu, *hv, *huu, *huv, *hvv, *huuu, *huuv, *huvv, *hvvv;
  double wsp[4];

  sp = pkv_GetScratchMemTop ();
  hu = pkv_GetScratchMemd ( 9*nkn*nkn );
  if ( !hu )
    goto failure;
  hv = &hu[nkn*nkn];      huu = &hv[nkn*nkn];
  huv = &huu[nkn*nkn];    hvv = &huv[nkn*nkn];
  huuu = &hvv[nkn*nkn];   huuv = &huuu[nkn*nkn];
  huvv = &huuv[nkn*nkn];  hvvv = &huvv[nkn*nkn];

  if ( !mbs_TabBezC1Coons0Der3d ( 1, nkn, tkn, hfunc, dhfunc, ddhfunc, dddhfunc,
                     nkn, tkn, hfunc, dhfunc, ddhfunc, dddhfunc,
                     G1_CROSS00DEG, fc00, G1_CROSS01DEG, fc01,
                     G1_CROSS00DEG, fd00, G1_CROSS01DEG, fd01,
                     NULL, hu, hv, huu, huv, hvv, huuu, huuv, huvv, hvvv ) )
    goto failure;

  for ( i = 0; i < nkn*nkn; i++ )
    if ( !_pkn_Comp2iDerivatives3d ( diu[i].x, diu[i].y, div[i].x, div[i].y,
              diuu[i].x, diuu[i].y, diuv[i].x, diuv[i].y, divv[i].x, divv[i].y,
              diuuu[i].x, diuuu[i].y, diuuv[i].x, diuuv[i].y,
              diuvv[i].x, diuvv[i].y, divvv[i].x, divvv[i].y,
              1, &hu[i], &hv[i], &huu[i], &huv[i], &hvv[i],
              &huuu[i], &huuv[i], &huvv[i], &hvvv[i],
              &psiu[i], &psiv[i], &psiuu[i], &psiuv[i], &psivv[i],
              &psiuuu[i], &psiuuv[i], &psiuvv[i], &psivvv[i], wsp ) )
      goto failure;

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*_g1hq2_TabNLDer0d*/

boolean _g1hq2_TabNLDerd ( int nkn, double *tkn,
             const double *hfunc, const double *dhfunc, const double *ddhfunc,
             const double *dddhfunc,
             vector2d *diu, vector2d *div,
             vector2d *diuu, vector2d *diuv, vector2d *divv,
             vector2d *diuuu, vector2d *diuuv, vector2d *diuvv, vector2d *divvv,
             double *fc00, double *fc01, double *fc10, double *fc11,
             double *fd00, double *fd01, double *fd10, double *fd11,
             double *psiu, double *psiv,
             double *psiuu, double *psiuv, double *psivv,
             double *psiuuu, double *psiuuv, double *psiuvv, double *psivvv )
{
  void   *sp;
  int    i;
  double *hu, *hv, *huu, *huv, *hvv, *huuu, *huuv, *huvv, *hvvv;
  double wsp[4];

  sp = pkv_GetScratchMemTop ();
  hu = pkv_GetScratchMemd ( 9*nkn*nkn );
  if ( !hu )
    goto failure;
  hv = &hu[nkn*nkn];      huu = &hv[nkn*nkn];
  huv = &huu[nkn*nkn];    hvv = &huv[nkn*nkn];
  huuu = &hvv[nkn*nkn];   huuv = &huuu[nkn*nkn];
  huvv = &huuv[nkn*nkn];  hvvv = &huvv[nkn*nkn];

  if ( !mbs_TabBezC1CoonsDer3d ( 1, nkn, tkn, hfunc, dhfunc, ddhfunc, dddhfunc,
                     nkn, tkn, hfunc, dhfunc, ddhfunc, dddhfunc,
                     G1_CROSS00DEG, fc00, G1_CROSS01DEG, fc01,
                     G1_CROSS10DEG, fc10, G1_CROSS11DEG, fc11,
                     G1_CROSS00DEG, fd00, G1_CROSS01DEG, fd01,
                     G1_CROSS10DEG, fd10, G1_CROSS11DEG, fd11,
                     NULL, hu, hv, huu, huv, hvv, huuu, huuv, huvv, hvvv ) )
    goto failure;

  for ( i = 0; i < nkn*nkn; i++ )
    if ( !_pkn_Comp2iDerivatives3d ( diu[i].x, diu[i].y, div[i].x, div[i].y,
              diuu[i].x, diuu[i].y, diuv[i].x, diuv[i].y, divv[i].x, divv[i].y,
              diuuu[i].x, diuuu[i].y, diuuv[i].x, diuuv[i].y,
              diuvv[i].x, diuvv[i].y, divvv[i].x, divvv[i].y,
              1, &hu[i], &hv[i], &huu[i], &huv[i], &hvv[i],
              &huuu[i], &huuv[i], &huvv[i], &hvvv[i],
              &psiu[i], &psiv[i], &psiuu[i], &psiuv[i], &psivv[i],
              &psiuuu[i], &psiuuv[i], &psiuvv[i], &psivvv[i], wsp ) )
      goto failure;

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*_g1hq2_TabNLDerd*/

boolean _g1hq2_SetupCTrdd ( const vector2d *cdiu, const vector2d *cdiv,
           const vector2d *cdiuu, const vector2d *cdiuv, const vector2d *cdivv,
           double *ctrd )
{
  vector2d gx, gy, gxx, gxy, gyy;

  if ( !pkn_f2iDerivatives2d ( cdiu->x, cdiu->y, cdiv->x, cdiv->y,
            cdiuu->x, cdiuu->y, cdiuv->x, cdiuv->y, cdivv->x, cdivv->y,
            &gx.x, &gy.x, &gxx.x, &gxy.x, &gyy.x ) )
    return false;
  pkn_Setup2DerA11Matrixd ( gx.x, gx.y, gy.x, gy.y, ctrd );
  pkn_Setup2DerA21Matrixd ( gxx.x, gxx.y, gxy.x, gxy.y, gyy.x, gyy.y, &ctrd[4] );
  pkn_Setup2DerA22Matrixd ( gx.x, gx.y, gy.x, gy.y, &ctrd[10] );
  return true;
} /*_g1hq2_SetupCTrdd*/

boolean _g1hq2_TabNLBasisFunctionsOmegad ( GHoleDomaind *domain, int nkn,
                                           G1HNLPrivated *nlpr, double *bc00 )
{
#define N ((G1H_FINALDEG+1)*(G1H_FINALDEG+1))
  void     *sp;
  GHolePrivateRecd  *privateG;
  G1HolePrivateRecd *privateG1;
  double   *tkn;
  double   *fc00, *fc01, *fc10, *fc11, *fd00, *fd01, *fd10, *fd11,
           *bc01, *bc10, *bc11, *bd00, *bd01, *bd10, *bd11;
  double   *hfunc, *dhfunc, *ddhfunc, *dddhfunc;
  int      i, j, k, l, kN, kNQ2, hole_k, f, fN, nfunc_a, nfunc_b, nkn2;
  point2d  *di, ddi;
  double   bvz;
  unsigned char *bfcpn;

  sp      = pkv_GetScratchMemTop ();
  if ( !_g1hq2_FindNLDomainDiameterd ( domain, nlpr ) )
    goto failure;
  hole_k  = domain->hole_k;
  privateG  = domain->privateG;
  privateG1 = domain->privateG1;
  nfunc_a = privateG1->nfunc_a;
  nfunc_b = privateG1->nfunc_b;
  di      = pkv_GetScratchMem ( N*sizeof(point2d) );
  tkn     = pkv_GetScratchMemd ( nkn );
  hfunc = pkv_GetScratchMemd ( 16*nkn );
  if ( !tkn || !di || !hfunc ) {
    domain->error_code = G1H_ERROR_NO_SCRATCH_MEMORY;
    goto failure;
  }

  dhfunc = &hfunc[4*nkn];
  ddhfunc = &dhfunc[4*nkn];
  dddhfunc = &ddhfunc[4*nkn];

  nkn2 = nkn*nkn;
  _gh_PrepareTabKnotsd ( nkn, privateG1->opt_quad, tkn );
  if ( !mbs_TabCubicHFuncDer3d ( 0.0, 1.0, nkn, tkn,
                                 hfunc, dhfunc, ddhfunc, dddhfunc ) )
    goto failure;

  for ( k = kN = kNQ2 = 0;
        k < hole_k;
        k++, kN += N, kNQ2 += nkn2 ) {
    pkv_Selectd ( N, 2, 3, 2, &nlpr->nldi[kN], di );
    for ( i = l = 0;  i < nkn;  i++ )
         /* ***** this is not an optimal algorithm, to be improved */
      for ( j = 0;  j < nkn;  j++, l++ ) {
        if ( !mbs_BCHornerDer3Pd ( G1H_FINALDEG, G1H_FINALDEG, 2, (double*)di,
                    tkn[i], tkn[j], (double*)&ddi,
                    (double*)&nlpr->diu[kNQ2+l], (double*)&nlpr->div[kNQ2+l],
                    (double*)&nlpr->diuu[kNQ2+l], (double*)&nlpr->diuv[kNQ2+l],
                    (double*)&nlpr->divv[kNQ2+l],
                    (double*)&nlpr->diuuu[kNQ2+l], (double*)&nlpr->diuuv[kNQ2+l],
                    (double*)&nlpr->diuvv[kNQ2+l], (double*)&nlpr->divvv[kNQ2+l] ) )
          goto failure;
        nlpr->jac[kNQ2+l] = (double)det2d ( &nlpr->diu[kNQ2+l], &nlpr->div[kNQ2+l] );
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
      _g1h_GetBFAPatchCurvesd ( domain, f, k, &fc00, &fc01, &fd00, &fd01 );
      if ( !_g1hq2_TabNLDer0d ( nkn, tkn,
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

  memset ( bc00, 0, 2*(G1_CROSSDEGSUM+4)*hole_k*sizeof(double) );
  bfcpn = privateG->bfcpn;
  for ( f = 0; f < nfunc_b; f++ ) {
        /* find Coons representation of the constant part of the solution */
    bvz = nlpr->rhole_cp[bfcpn[f]].z;
    for ( k = 0;  k < hole_k;  k++ ) {
      _g1h_GetBFBPatchCurvesd ( domain, f, k, &fc00, &fc01, &fc10, &fc11,
                                              &fd00, &fd01, &fd10, &fd11 );
      pkn_AddMatrixMd ( 1, G1_CROSS00DEG+1, 0, &bc00[(G1_CROSS00DEG+1)*k], 0, fc00,
                        bvz, 0, &bc00[(G1_CROSS00DEG+1)*k] );
      pkn_AddMatrixMd ( 1, G1_CROSS01DEG+1, 0, &bc01[(G1_CROSS01DEG+1)*k], 0, fc01,
                        bvz, 0, &bc01[(G1_CROSS01DEG+1)*k] );
      pkn_AddMatrixMd ( 1, G1_CROSS10DEG+1, 0, &bc10[(G1_CROSS10DEG+1)*k], 0, fc10,
                        bvz, 0, &bc10[(G1_CROSS10DEG+1)*k] );
      pkn_AddMatrixMd ( 1, G1_CROSS11DEG+1, 0, &bc11[(G1_CROSS11DEG+1)*k], 0, fc11,
                        bvz, 0, &bc11[(G1_CROSS11DEG+1)*k] );
      pkn_AddMatrixMd ( 1, G1_CROSS00DEG+1, 0, &bd00[(G1_CROSS00DEG+1)*k], 0, fd00,
                        bvz, 0, &bd00[(G1_CROSS00DEG+1)*k] );
      pkn_AddMatrixMd ( 1, G1_CROSS01DEG+1, 0, &bd01[(G1_CROSS01DEG+1)*k], 0, fd01,
                        bvz, 0, &bd01[(G1_CROSS01DEG+1)*k] );
      pkn_AddMatrixMd ( 1, G1_CROSS10DEG+1, 0, &bd10[(G1_CROSS10DEG+1)*k], 0, fd10,
                        bvz, 0, &bd10[(G1_CROSS10DEG+1)*k] );
      pkn_AddMatrixMd ( 1, G1_CROSS11DEG+1, 0, &bd11[(G1_CROSS11DEG+1)*k], 0, fd11,
                        bvz, 0, &bd11[(G1_CROSS11DEG+1)*k] );
    }
  }
  for ( k = kNQ2 = 0, fN = nfunc_a*nkn2*hole_k;
        k < hole_k;
        k++, kNQ2 += nkn2, fN += nkn2 ) {
    if ( !_g1hq2_TabNLDerd ( nkn, tkn, hfunc, dhfunc, ddhfunc, dddhfunc,
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
} /*_g1hq2_TabNLBasisFunctionsOmegad*/

boolean _g1hq2_TabNLBasisFunctionsGammad ( GHoleDomaind *domain, int nkn,
                                           G1HNLPrivated *nlpr, double *ctrd,
                                           double *bc00 )
{
#define N ((G1H_FINALDEG+1)*(G1H_FINALDEG+1))
  void     *sp;
  GHolePrivateRecd  *privateG;
  G1HolePrivateRecd *privateG1;
  double   *tkn, *hfunc, *dhfunc, *ddhfunc, *dddhfunc,
           *atkn, *ahfunc, *adhfunc, *addhfunc;
  double   *ec00, *ec01, *ed00, *ed01, *fc00, *fc01, *fd00, *fd01,
           *bc01, *bc10, *bc11, *bd00, *bd01, *bd10, *bd11;
  int      i, j, k, kk, kN, kNQ2, hole_k, f, fN, nfunc_a, nfunc_b;
  point2d  *di;
  double   bvz;
  unsigned char *bfcpn;
  vector2d cdi, cdiu, cdiv, cdiuu, cdiuv, cdivv, *sicp;
  double   *tabeu, *tabev, *tabeuu, *tabeuv, *tabevv,
           *tabfu, *tabfv, *tabfuu, *tabfuv, *tabfvv;
  double   *A11, *A21, *A22;

  sp      = pkv_GetScratchMemTop ();
  hole_k  = domain->hole_k;
  privateG  = domain->privateG;
  privateG1 = domain->privateG1;
  nfunc_a = privateG1->nfunc_a;
  nfunc_b = privateG1->nfunc_b;
  di      = pkv_GetScratchMem ( N*sizeof(point2d) );
  tkn     = pkv_GetScratchMemd ( nkn );
  hfunc = pkv_GetScratchMemd ( 16*nkn );
  if ( !tkn || !di || !hfunc ) {
    domain->error_code = G1H_ERROR_NO_SCRATCH_MEMORY;
    goto failure;
  }

  dhfunc = &hfunc[4*nkn];
  ddhfunc = &dhfunc[4*nkn];
  dddhfunc = &ddhfunc[4*nkn];

  _gh_PrepareTabKnotsd ( nkn, privateG1->opt_quad, tkn );
  if ( !mbs_TabCubicHFuncDer3d ( 0.0, 1.0, nkn, tkn,
                                 hfunc, dhfunc, ddhfunc, dddhfunc ) )
    goto failure;

        /* now evaluate the 1st order derivatives and jumps of 2nd order */
        /* derivatives at the boundary curves of the areas Omega_i */
  sicp = (point2d*)pkv_GetScratchMem ( 16*sizeof(point2d) );
  tabeu = pkv_GetScratchMemd ( 20*nkn );
  atkn = pkv_GetScratchMemd ( 26 );
  if ( !sicp || !tabeu || !atkn )
    goto failure;

  tabev = &tabeu[2*nkn];    tabeuu = &tabev[2*nkn];   tabeuv = &tabeuu[2*nkn];
  tabevv = &tabeuv[2*nkn];  tabfu = &tabevv[2*nkn];   tabfv = &tabfu[2*nkn];
  tabfuu = &tabfv[2*nkn];   tabfuv = &tabfuu[2*nkn];  tabfvv = &tabfuv[2*nkn];

  ahfunc = &atkn[2];  adhfunc = &ahfunc[8];  addhfunc = &adhfunc[8];
  atkn[0] = 0.0;  atkn[1] = 1.0;
  if ( !mbs_TabCubicHFuncDer2d ( 0.0, 1.0, 2, atkn, ahfunc, adhfunc, addhfunc ) )
    goto failure;

  memset ( ctrd, 0, 3*nkn*38*hole_k*sizeof(double) );
  for ( k = kN = kNQ2 = 0, kk = 1;
        k < hole_k;
        k++, kk = (k+1) % hole_k, kN += N, kNQ2 += 3*nkn ) {
    pkv_Selectd ( N, 2, 3, 2, &nlpr->nldi[kN], di );
    for ( i = 0; i < nkn; i++ ) {
      if ( !mbs_BCHornerDer2Pd ( G1H_FINALDEG, G1H_FINALDEG, 2, (double*)di,
                                 tkn[i], 0.0, &cdi.x, &cdiu.x, &cdiv.x,
                                 &cdiuu.x, &cdiuv.x, &cdivv.x ) )
        goto failure;
      nlpr->ctang[kNQ2+i] = cdiu;
      if ( !_g1hq2_SetupCTrdd ( &cdiu, &cdiv, &cdiuu, &cdiuv, &cdivv,
                                &ctrd[38*(kNQ2+i)] ) )
        goto failure;
      if ( !mbs_BCHornerDer2Pd ( G1H_FINALDEG, G1H_FINALDEG, 2, (double*)di,
                                 tkn[i], 1.0, &cdi.x, &cdiu.x, &cdiv.x,
                                 &cdiuu.x, &cdiuv.x, &cdivv.x ) )
        goto failure;
      nlpr->ctang[kNQ2+nkn+i] = cdiu;
      if ( !_g1hq2_SetupCTrdd ( &cdiu, &cdiv, &cdiuu, &cdiuv, &cdivv,
                                &ctrd[38*(kNQ2+nkn+i)] ) )
        goto failure;
      if ( !mbs_BCHornerDer2Pd ( G1H_FINALDEG, G1H_FINALDEG, 2, (double*)di,
                                 1.0, tkn[i], &cdi.x, &cdiu.x, &cdiv.x,
                                 &cdiuu.x, &cdiuv.x, &cdivv.x ) )
        goto failure;
      nlpr->ctang[kNQ2+2*nkn+i] = cdiv;
      if ( !_g1hq2_SetupCTrdd ( &cdiu, &cdiv, &cdiuu, &cdiuv, &cdivv,
                                &ctrd[38*(kNQ2+2*nkn+i)] ) )
        goto failure;
      if ( !mbs_BCHornerDer2Pd ( G1H_FINALDEG, G1H_FINALDEG, 2, (double*)di,
                                 0.0, tkn[i], &cdi.x, &cdiu.x, &cdiv.x,
                                 &cdiuu.x, &cdiuv.x, &cdivv.x ) )
        goto failure;
      if ( !_g1hq2_SetupCTrdd ( &cdiu, &cdiv, &cdiuu, &cdiuv, &cdivv,
                                &ctrd[38*(kk*3*nkn+i)+19] ) )
        goto failure;
    }

    if ( !_g1hq2_FindDomSurrndPatchd ( domain, nlpr, kk, 1, sicp ) )
      goto failure;
    for ( i = 0; i < nkn; i++ ) {
      if ( !mbs_BCHornerDer2Pd ( 3, 3, 2, (double*)sicp,
                                 0.0, tkn[i], &cdi.x, &cdiu.x, &cdiv.x,
                                 &cdiuu.x, &cdiuv.x, &cdivv.x ) )
        goto failure;
      if ( !_g1hq2_SetupCTrdd ( &cdiu, &cdiv, &cdiuu, &cdiuv, &cdivv,
                                &ctrd[38*(kNQ2+nkn+i)+19] ) )
        goto failure;
    }
    if ( !_g1hq2_FindDomSurrndPatchd ( domain, nlpr, k, 2, sicp ) )
      goto failure;
    for ( i = 0; i < nkn; i++ ) {
      if ( !mbs_BCHornerDer2Pd ( 3, 3, 2, (double*)sicp,
                                 0.0, tkn[i], &cdi.x, &cdiu.x, &cdiv.x,
                                 &cdiuu.x, &cdiuv.x, &cdivv.x ) )
        goto failure;
      if ( !_g1hq2_SetupCTrdd ( &cdiu, &cdiv, &cdiuu, &cdiuv, &cdivv,
                                &ctrd[38*(kNQ2+2*nkn+i)+19] ) )
        goto failure;
    }
  }

  memset ( nlpr->cpsiu, 0, 3*nkn*(nfunc_a+1)*hole_k*sizeof(double) );
  memset ( nlpr->cpsiv, 0, 3*nkn*(nfunc_a+1)*hole_k*sizeof(double) );
  memset ( nlpr->cpsiuu, 0, 3*nkn*(nfunc_a+1)*hole_k*sizeof(double) );
  memset ( nlpr->cpsiuv, 0, 3*nkn*(nfunc_a+1)*hole_k*sizeof(double) );
  memset ( nlpr->cpsivv, 0, 3*nkn*(nfunc_a+1)*hole_k*sizeof(double) );
        /* block A functions */
  for ( f = 0; f < nfunc_a; f++ ) {
    _g1h_GetBFAPatchCurvesd ( domain, f, hole_k-1,
                              &ec00, &ec01, &ed00, &ed01 );
    for ( k = 0; k < hole_k; k++ ) {
      _g1h_GetBFAPatchCurvesd ( domain, f, k,
                                &fc00, &fc01, &fd00, &fd01 );
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
      if ( !mbs_TabBezC1Coons0Der2d ( 1,
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
    if ( !mbs_TabBezC1CoonsDer2d ( 1, nkn, tkn, hfunc, dhfunc, ddhfunc,
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
    if ( !mbs_TabBezC1CoonsDer2d ( 1, 2, atkn, ahfunc, adhfunc, addhfunc,
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
  memset ( bc00, 0, 2*(G1_CROSSDEGSUM+4)*hole_k*sizeof(double) );
  for ( f = 0; f < nfunc_b; f++ ) {
    bvz = nlpr->rhole_cp[bfcpn[f]].z;
    for ( k = 0, kk = 1;  k < hole_k;  k++, kk = (k+1) % hole_k ) {
      gh_GetDomSurrndBFuncd ( domain, f, kk, 1, (double*)di );
      pkn_AddMatrixMd ( 1, 16, 0, &bc00[k*2*16], 0, (double*)di, bvz,
                        0, &bc00[k*2*16] );
      gh_GetDomSurrndBFuncd ( domain, f, k, 2, (double*)di );
      pkn_AddMatrixMd ( 1, 16, 0, &bc00[(k*2+1)*16], 0, (double*)di, bvz,
                        0, &bc00[(k*2+1)*16] );
    }
  }
  for ( k = 0; k < hole_k; k++ ) {
    for ( i = 0, fN = (3*(nfunc_a*hole_k+k)+1)*nkn;  i < nkn;  i++, fN++ ) {
      if ( !mbs_BCHornerDer2Pd ( 3, 3, 1, &bc00[k*2*16], 0.0, tkn[i],
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
      if ( !mbs_BCHornerDer2Pd ( 3, 3, 1, &bc00[(k*2+1)*16], 0.0, tkn[i],
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
} /*_g1hq2_TabNLBasisFunctionsGammad*/

/* ///////////////////////////////////////////////////////////////////////// */
void _g1hq2_IntFunc1ad ( G1Q2HNLFuncd *f, double *funct )
{
  double e0, e1, e2, e9, e0p2, pup2, pvp2;
  double A, B, x, y, xx, xy, yy, tGt;

  pup2 = f->pu*f->pu;
  pvp2 = f->pv*f->pv;
  f->e1 = e1 = 1.0 + pup2;
  f->e2 = e2 = 1.0 + pvp2;
  f->e0 = e0 = e1 + pvp2;
  e0p2 = e0*e0;
  f->e9 = e9 = f->pu*f->pv;
  x = f->tang.x;
  y = f->tang.y;
  xx = x*x;  xy = x*y;  yy = y*y;
  f->tGt = tGt = e1*xx + 2.0*e9*xy + e2*yy;
  f->A = A = 0.5*(e2*f->jpuu+e1*f->jpvv) - e9*f->jpuv;
  f->B = B = sqrt(tGt)/(e0*e0p2);
  *funct += A*A*B;
} /*_g1hq2_IntFunc1ad*/

void _g1hq2_IntFunc1bd ( G1Q2HNLFuncd *f, double *funct )
{
  double e0, e1, e2, e9, e0p2, e0p4, pup2, pvp2;
  double A, B, x, y, xx, xy, yy, tGt;
  double b11, b12, b22, dd;

  pup2 = f->pu*f->pu;
  pvp2 = f->pv*f->pv;
  f->e1 = e1 = 1.0 + pup2;
  f->e2 = e2 = 1.0 + pvp2;
  f->e0 = e0 = e1 + pvp2;
  e0p2 = e0*e0;
  f->e9 = e9 = f->pu*f->pv;
  x = f->tang.x;
  y = f->tang.y;
  xx = x*x;  xy = x*y;  yy = y*y;
  f->tGt = tGt = e1*xx + 2.0*e9*xy + e2*yy;
  f->A = A = 0.5*(e2*f->jpuu+e1*f->jpvv) - e9*f->jpuv;
  f->B = B = sqrt(tGt)/(e0*e0p2);
  *funct += A*A*B;
        /* data needed for the computation of derivatives of B */
  e0p4 = e0p2*e0p2;
  dd = 1.0/(e0p4*sqrt(tGt));
  b11 = f->pu*(pvp2-5.0*e1);
  b12 = 0.5*f->pv*(e2-11.0*pup2);
  b22 = -6.0*f->pu*e2;
  f->b1 = (b11*xx + 2.0*b12*xy + b22*yy)*dd;
  b11 = -6.0*f->pv*e1;
  b12 = 0.5*f->pu*(e1-11.0*pvp2);
  b22 = f->pv*(pup2-5.0*e2);
  f->b2 = (b11*xx + 2.0*b12*xy + b22*yy)*dd;
} /*_g1hq2_IntFunc1bd*/

void _g1hq2_IntFunc1cd ( G1Q2HNLFuncd *f, double *funct )
{
  double e0, e1, e2, e9, e0p2, e0p4, e1p2, e2p2, e9p2, pup2, pvp2;
  double A, B, x, y, xx, xy, yy, xxxx, xxxy, xxyy, xyyy, yyyy, tGt;
  double b11, b12, b22, b23, b33, dd;

  pup2 = f->pu*f->pu;
  pvp2 = f->pv*f->pv;
  f->e1 = e1 = 1.0 + pup2;
  f->e2 = e2 = 1.0 + pvp2;
  f->e0 = e0 = e1 + pvp2;
  e0p2 = e0*e0;
  f->e9 = e9 = f->pu*f->pv;
  x = f->tang.x;
  y = f->tang.y;
  xx = x*x;  xy = x*y;  yy = y*y;
  f->tGt = tGt = e1*xx + 2.0*e9*xy + e2*yy;
  f->A = A = 0.5*(e2*f->jpuu+e1*f->jpvv) - e9*f->jpuv;
  f->B = B = sqrt(tGt)/(e0*e0p2);
  *funct += A*A*B;
        /* data needed for the computation of derivatives of B */
  e0p4 = e0p2*e0p2;
  dd = 1.0/(e0p4*sqrt(tGt));
  b11 = f->pu*(pvp2-5.0*e1);
  b12 = 0.5*f->pv*(e2-11.0*pup2);
  b22 = -6.0*f->pu*e2;
  f->b1 = (b11*xx + 2.0*b12*xy + b22*yy)*dd;
  b11 = -6.0*f->pv*e1;
  b12 = 0.5*f->pu*(e1-11.0*pvp2);
  b22 = f->pv*(pup2-5.0*e2);
  f->b2 = (b11*xx + 2.0*b12*xy + b22*yy)*dd;
        /* data for the derivatives of B of the second order */  
  e1p2 = e1*e1;
  e2p2 = e2*e2;
  e9p2 = e9*e9;
  dd = 1.0/(e0p4*e0*tGt*sqrt(tGt));
  xxxx = xx*xx;
  xxxy = xx*xy;
  xxyy = xx*yy;
  xyyy = xy*yy;
  yyyy = yy*yy;
  b11 = pup2*((30.0*pup2+55.0)*pup2-18.0*e9p2-22.0*pvp2+20.0) + (pvp2-5.0)*e2;
  b12 = e9*(pup2*(66.0*pup2+48.0-30.0*pvp2)-18.0*e2);
  b22 = e9p2*(216.0*pup2-72.0*pvp2-10.0)+(73.0*pup2+62.0)*pup2-11.0*e2p2;
  b23 = e9*e2*(78.0*pup2-18.0*e2);
  b33 = 6.0*e2p2*(7.0*pup2-e2);
  f->b3 = (b11*xxxx + 2.0*b12*xxxy + b22*xxyy + 2.0*b23*xyyy + b33*yyyy)*dd;
  b11 = 6.0*e9*e1*(7.0*e1-pvp2);
  b12 = 0.5*(pup2*(pvp2*(168.0*pup2-18.0*pvp2+164.0)-(6.0*pup2+11.0)*pup2-4.0)+
             (1.0-5.0*pvp2)*e2);
  b22 = e9*(252.0*e9p2-18.0*(pup2*pup2+pvp2*pvp2)+66.0*(pup2+pvp2)+84.0);
  b23 = 0.5*(pvp2*(pup2*(168.0*pvp2-18.0*pup2+164.0)-(6.0*pvp2+11.0)*pvp2-4.0)+
             (1.0-5.0*pup2)*e1);
  b33 = 6.0*e9*e2*(7.0*e2-pup2);
  f->b4 = (b11*xxxx + 2.0*b12*xxxy + b22*xxyy + 2.0*b23*xyyy + b33*yyyy)*dd;
  b11 = 6.0*e1p2*(7.0*pvp2-e1);
  b12 = e9*e1*(78.0*pvp2-18.0*e1);
  b22 = e9p2*(216.0*pvp2-72.0*pup2-10.0)+(73.0*pvp2+62.0)*pvp2-11.0*e1;
  b23 = e9*(pvp2*(66.0*pvp2+48.0-30.0*pup2)-18.0*e1);
  b33 = pvp2*((30.0*pvp2+55.0)*pvp2-18.0*e9p2-22.0*pup2+20.0) + (pup2-5.0)*e1;
  f->b5 = (b11*xxxx + 2.0*b12*xxxy + b22*xxyy + 2.0*b23*xyyy + b33*yyyy)*dd;
} /*_g1hq2_IntFunc1cd*/

void _g1hq2_IntFunc2bd ( G1Q2HNLFuncd *f, double *grad )
{
  double A_i, B_i;
  double e1, e2, e9;

  e1 = f->e1;  e2 = f->e2;  e9 = f->e9;
  A_i = (f->pu*f->jpvv-f->pv*f->jpuv)*f->psiu
       -(f->pu*f->jpuv-f->pv*f->jpuu)*f->psiv
       +0.5*(e2*f->jpsiuu+e1*f->jpsivv)-e9*f->jpsiuv;
  B_i = f->b1*f->psiu + f->b2*f->psiv;
  *grad += f->A*(2.0*A_i*f->B + f->A*B_i);
} /*_g1hq2_IntFunc2bd*/

void _g1hq2_IntFunc2cd ( G1Q2HNLFuncd *f, double *Ai, double *Bi, double *grad )
{
  double A_i, B_i;
  double e1, e2, e9;

  e1 = f->e1;  e2 = f->e2;  e9 = f->e9;
  *Ai = A_i = (f->pu*f->jpvv-f->pv*f->jpuv)*f->psiu
             -(f->pu*f->jpuv-f->pv*f->jpuu)*f->psiv
             +0.5*(e2*f->jpsiuu+e1*f->jpsivv)-e9*f->jpsiuv;
  *Bi = B_i = f->b1*f->psiu + f->b2*f->psiv;
  *grad += f->A*(2.0*A_i*f->B + f->A*B_i);
} /*_g1hq2_IntFunc2cd*/

void _g1hq2_IntFunc3cd ( G1Q2HNLFuncd *f, double Ai, double Bi, double Aj, double Bj,
                         double *hessian )
{
  double Aij, Bij;

  Aij = f->pu*(f->psiu*f->jpsjvv+f->jpsivv*f->psju-f->psiv*f->jpsjuv-f->jpsiuv*f->psjv)
       +f->pv*(f->psiv*f->jpsjuu+f->jpsiuu*f->psjv-f->psiu*f->jpsjuv-f->jpsiuv*f->psju)
       +f->jpuu*f->psiv*f->psjv
       -f->jpuv*(f->psiu*f->psjv+f->psiv*f->psju)
       +f->jpvv*f->psiu*f->psju;
  Bij = f->b3*f->psiu*f->psju
       +f->b4*(f->psiu*f->psjv+f->psiv*f->psju)
       +f->b5*f->psiv*f->psjv;
  *hessian += 2.0*(f->A*(Aij*f->B+Ai*Bj+Aj*Bi)+ Ai*Aj*f->B) + f->A*f->A*Bij;
} /*_g1hq2_IntFunc3cd*/

