
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

#include "eg2holef.h"
#include "eg2hprivatef.h"
#include "eg2herror.h"


G2HNLPrivatef *_g2h_nlprivf;

/* ///////////////////////////////////////////////////////////////////////// */
void g2h_ReflectVectorsf ( int n, const vector3f *v, vector3f *w )
{
  vector3f r;
  double   a;
  int      i;

  r = _g2h_nlprivf->reflv;
  for ( i = 0; i < n; i++ ) {
    a = DotProduct3f ( &v[i], &r );
    AddVector3Mf ( &v[i], &r, -2.0*a, &w[i] );
  }
} /*g2h_ReflectVectorsf*/

void g2h_nonlinoutpatchf ( int n, int m, const float *cp, void *usrptr )
{
#define N (G2H_FINALDEG+1)*(G2H_FINALDEG+1)
  memcpy ( &_g2h_nlprivf->nldi[_g2h_nlprivf->auxc*N], cp, N*sizeof(vector3f) );
  _g2h_nlprivf->auxc ++;
#undef N
} /*g2h_nonlinoutpatchf*/

boolean _g2h_StopItf ( int itn, float gn0, float gn,
                       float cn, float dcn, float scf )
{
#define MAXITER 20
#define EPS0    5.0e-6
#define EPS1    1.0e-4

/*
printf ( "%2d: func = %f, gn = %f, cn = %f, dcn = %f, scf = %f\n",
         itn, func, gn, cn, dcn, scf );
*/
  if ( itn >= MAXITER )
    return true;
  if ( gn <= EPS0*gn0 || (gn <= EPS1*gn0 && scf < 0.75 ) )
    return true;
  if ( scf*dcn <= EPS0*cn )
    return true;
  return false;
#undef EPS1
#undef EPS0
#undef MAXITER
} /*_g2h_StopItf*/

/* ///////////////////////////////////////////////////////////////////////// */
boolean g2h_GetHoleSurrndPatchf ( GHoleDomainf *domain,
                                  const point3f *hole_cp,
                                  int i, int j, point3f *bcp )
{
  void    *sp;
  int     *ind;
  float   *ukn, *vkn;
  point3f *q;
  int     hole_k, k;

  sp  = pkv_GetScratchMemTop ();
  ind = pkv_GetScratchMem ( 16*sizeof(int) );
  q   = pkv_GetScratchMem ( 16*sizeof(point3f) );
  if ( !ind || !q ) {
    domain->error_code = G2H_ERROR_NO_SCRATCH_MEMORY;
    goto failure;
  }

  hole_k = domain->hole_k;
  gh_GetBspInd ( hole_k, i, j, ind );
  for ( k = 0; k < 16; k++ )
    q[k] = hole_cp[ind[k]];
  ukn = &domain->hole_knots[11*((i+hole_k-1) % hole_k)+3 ];
  vkn = &domain->hole_knots[11*i+j];
  if ( !mbs_BSPatchToBezf ( 3, 3, 7, ukn, 3, 7, vkn, 12, (float*)q,
                            NULL, NULL, NULL, NULL, NULL, NULL, 12, (float*)bcp ) )
    goto failure;

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*g2h_GetHoleSurrndPatchf*/

boolean g2h_ComputeNLNormalf ( GHoleDomainf *domain,
                               const point3f *hole_cp,
                               vector3f *anv )
{
#define DENSITY 16
  void     *sp;
  int      hole_k;
  point3f  *bcp, p;
  int      i, j, k;
  vector3f pu, pv, nu, nv;
  float    t;

  sp = pkv_GetScratchMemTop ();
  hole_k = domain->hole_k;
  bcp = pkv_GetScratchMem ( 16*sizeof(point3f) );
  if ( !bcp ) {
    domain->error_code = G2H_ERROR_NO_SCRATCH_MEMORY;
    goto failure;
  }

  SetVector3f ( &nv, 0.0, 0.0, 0.0 );
  for ( i = 0; i < hole_k; i++ )
    for ( j = 1; j < 3; j++ ) {
      g2h_GetHoleSurrndPatchf ( domain, hole_cp, i, j, bcp );
          /* integrate the unit normal vector at the boundary */
      for ( k = 0; k < DENSITY; k++ ) {
        t = (float)(k+k+1)/(float)(2*DENSITY);
        if ( !mbs_BCHornerDerP3f ( 3, 3, bcp, 0.0, t, &p, &pu, &pv ) )
          goto failure;
        NormalizeVector3f ( &pv );
        CrossProduct3f ( &pu, &pv, &nu );
        AddVector3f ( &nv, &nu, &nv );
      }
    }
  NormalizeVector3f ( &nv );

        /* now a verification */
  for ( i = 0; i < hole_k; i++ )
    for ( j = 1; j < 3; j++ ) {
      g2h_GetHoleSurrndPatchf ( domain, hole_cp, i, j, bcp );
          /* integrate the unit normal vector at the boundary */
      for ( k = 0; k < DENSITY; k++ ) {
        t = (float)(k+k+1)/(float)(2*DENSITY);
        if ( !mbs_BCHornerDerP3f ( 3, 3, bcp, 0.0, t, &p, &pu, &pv ) )
          goto failure;
        CrossProduct3f ( &pu, &pv, &nu );
        if ( DotProduct3f ( &nu, &nv ) <= 0.0 ) {
          domain->error_code = G2H_ERROR_NL_CANNOT_PROJECT;
          goto failure;
        }
      }
    }
  *anv = nv;

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
#undef DENSITY
} /*g2h_ComputeNLNormalf*/

boolean _g2h_ComputeNLNormalf ( GHoleDomainf *domain,
                                G2HNLPrivatef *nlprivate,
                                const point3f *hole_cp )
{
  vector3f nv;

  if ( g2h_ComputeNLNormalf ( domain, hole_cp, &nv ) ) {
    nlprivate->nlnv = nv;
    if ( nv.z > 0.0 )
      nv.z += 1.0;
    else
      nv.z -= 1.0;
    NormalizeVector3f ( &nv );
    nlprivate->reflv = nv;
    return true;
  }
  else
    return false;
} /*_g2h_ComputeNLNormalf*/

boolean _g2h_TabNLDer0f ( GHoleDomainf *domain,
             int nkn, const float *tkn,
             const float *hfunc, const float *dhfunc,
             const float *ddhfunc, const float *dddhfunc,
             vector2f *diu, vector2f *div,
             vector2f *diuu, vector2f *diuv, vector2f *divv,
             vector2f *diuuu, vector2f *diuuv, vector2f *diuvv, vector2f *divvv,
             float *fc00, float *fc01, float *fc02,
             float *fd00, float *fd01, float *fd02,
             float *psiu, float *psiv,
             float *psiuu, float *psiuv, float *psivv,
             float *psiuuu, float *psiuuv, float *psiuvv, float *psivvv )
{
  void  *sp;
  int   i;
  float *hu, *hv, *huu, *huv, *hvv, *huuu, *huuv, *huvv, *hvvv;
  float wsp[4];

  sp = pkv_GetScratchMemTop ();
  hu = pkv_GetScratchMemf ( 9*nkn*nkn );
  if ( !hu ) {
    domain->error_code = G2H_ERROR_NO_SCRATCH_MEMORY;
    goto failure;
  }
  hv = &hu[nkn*nkn];      huu = &hv[nkn*nkn];
  huv = &huu[nkn*nkn];    hvv = &huv[nkn*nkn];
  huuu = &hvv[nkn*nkn];   huuv = &huuu[nkn*nkn];
  huvv = &huuv[nkn*nkn];  hvvv = &huvv[nkn*nkn];

  if ( !mbs_TabBezC2Coons0Der3f ( 1, nkn, tkn, hfunc, dhfunc, ddhfunc, dddhfunc,
                     nkn, tkn, hfunc, dhfunc, ddhfunc, dddhfunc,
                     G2_CROSS00DEG, fc00, G2_CROSS01DEG, fc01, G2_CROSS02DEG, fc02,
                     G2_CROSS00DEG, fd00, G2_CROSS01DEG, fd01, G2_CROSS02DEG, fd02,
                     NULL, hu, hv, huu, huv, hvv, huuu, huuv, huvv, hvvv ) )
    goto failure;

  for ( i = 0; i < nkn*nkn; i++ )
    if ( !_pkn_Comp2iDerivatives3f ( diu[i].x, diu[i].y, div[i].x, div[i].y,
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
} /*_g2h_TabNLDer0f*/

boolean _g2h_TabNLDerf ( GHoleDomainf *domain,
             int nkn, float *tkn,
             const float *hfunc, const float *dhfunc,
             const float *ddhfunc, const float *dddhfunc,
             vector2f *diu, vector2f *div,
             vector2f *diuu, vector2f *diuv, vector2f *divv,
             vector2f *diuuu, vector2f *diuuv, vector2f *diuvv, vector2f *divvv,
             float *fc00, float *fc01, float *fc02,
             float *fc10, float *fc11, float *fc12,
             float *fd00, float *fd01, float *fd02,
             float *fd10, float *fd11, float *fd12,
             float *psiu, float *psiv,
             float *psiuu, float *psiuv, float *psivv,
             float *psiuuu, float *psiuuv, float *psiuvv, float *psivvv )
{
  void  *sp;
  int   i;
  float *hu, *hv, *huu, *huv, *hvv, *huuu, *huuv, *huvv, *hvvv;
  float wsp[4];

  sp = pkv_GetScratchMemTop ();
  hu = pkv_GetScratchMemf ( 9*nkn*nkn );
  if ( !hu ) {
    domain->error_code = G2H_ERROR_NO_SCRATCH_MEMORY;
    goto failure;
  }
  hv = &hu[nkn*nkn];      huu = &hv[nkn*nkn];
  huv = &huu[nkn*nkn];    hvv = &huv[nkn*nkn];
  huuu = &hvv[nkn*nkn];   huuv = &huuu[nkn*nkn];
  huvv = &huuv[nkn*nkn];  hvvv = &huvv[nkn*nkn];

  if ( !mbs_TabBezC2CoonsDer3f ( 1, nkn, tkn, hfunc, dhfunc, ddhfunc, dddhfunc,
            nkn, tkn, hfunc, dhfunc, ddhfunc, dddhfunc,
            G2_CROSS00DEG, fc00, G2_CROSS01DEG, fc01, G2_CROSS02DEG, fc02,
            G2_CROSS10DEG, fc10, G2_CROSS11DEG, fc11, G2_CROSS12DEG, fc12,
            G2_CROSS00DEG, fd00, G2_CROSS01DEG, fd01, G2_CROSS02DEG, fd02,
            G2_CROSS10DEG, fd10, G2_CROSS11DEG, fd11, G2_CROSS12DEG, fd12,
            NULL, hu, hv, huu, huv, hvv, huuu, huuv, huvv, hvvv ) )
    goto failure;

  for ( i = 0; i < nkn*nkn; i++ )
    if ( !_pkn_Comp2iDerivatives3f ( diu[i].x, diu[i].y, div[i].x, div[i].y,
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
} /*_g2h_TabNLDerf*/

boolean _g2h_TabNLBasisFunctionsf ( GHoleDomainf *domain,
                                    G2HNLPrivatef *nlpr )
{
#define N ((G2H_FINALDEG+1)*(G2H_FINALDEG+1))
  void     *sp;
  GHolePrivateRecf  *privateG;
  G2HolePrivateRecf *privateG2;
  float    *tkn;
  float    *fc00, *fc01, *fc02, *fc10, *fc11, *fc12,
           *fd00, *fd01, *fd02, *fd10, *fd11, *fd12,
           *bc00, *bc01, *bc02, *bc10, *bc11, *bc12,
           *bd00, *bd01, *bd02, *bd10, *bd11, *bd12;
  float    *hfunc, *dhfunc, *ddhfunc, *dddhfunc;
  int      i, j, k, l, kN, kNQ2, hole_k, f, fN, nfunc_a, nfunc_b;
  point2f  *di, ddi;
  float    bvz;
  unsigned char *bfcpn;

  sp      = pkv_GetScratchMemTop ();
  hole_k  = domain->hole_k;
  privateG  = domain->privateG;
  privateG2 = domain->privateG2;
  nfunc_a = privateG2->nfunc_a;
  nfunc_b = privateG2->nfunc_b;
  di      = pkv_GetScratchMem ( N*sizeof(point2f) );
  tkn     = pkv_GetScratchMemf ( G2_NQUAD );
  hfunc   = pkv_GetScratchMemf ( 24*G2_NQUAD );
  if ( !tkn || !di || !hfunc ) {
    domain->error_code = G2H_ERROR_NO_SCRATCH_MEMORY;
    goto failure;
  }

  dhfunc = &hfunc[6*G2_NQUAD];
  ddhfunc = &dhfunc[6*G2_NQUAD];
  dddhfunc = &ddhfunc[6*G2_NQUAD];

  _gh_PrepareTabKnotsf ( G2_NQUAD, privateG2->opt_quad, tkn );
  if ( !mbs_TabQuinticHFuncDer3f ( 0.0, 1.0, G2_NQUAD, tkn,
                                   hfunc, dhfunc, ddhfunc, dddhfunc ) )
    goto failure;

  for ( k = kN = kNQ2 = 0;
        k < hole_k;
        k++, kN += N, kNQ2 += G2_NQUADSQ ) {
    pkv_Selectf ( N, 2, 3, 2, &nlpr->nldi[kN], di );
    for ( i = l = 0;  i < G2_NQUAD;  i++ )
         /* ***** this is not an optimal algorithm, to be improved */
      for ( j = 0;  j < G2_NQUAD;  j++, l++ ) {
        if ( !mbs_BCHornerDer3Pf ( G2H_FINALDEG, G2H_FINALDEG, 2, (float*)di,
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
    for ( i = 1; i < hole_k*G2_NQUADSQ; i++ )
      if ( nlpr->jac[i] <= 0.0 ) {
        domain->error_code = G2H_ERROR_NL_JACOBIAN;
        goto failure;
      }
  }
  else {
    for ( i = 0; i < hole_k*G2_NQUADSQ; i++ ) {
      if ( nlpr->jac[i] >= 0.0 ) {
        domain->error_code = G2H_ERROR_NL_JACOBIAN;
        goto failure;
      }
      nlpr->jac[i] = -nlpr->jac[i];
    }
  }

  for ( f = fN = 0;  f < nfunc_a;  f++ ) {
    for ( k = kNQ2 = 0;  k < hole_k;  k++, fN += G2_NQUADSQ, kNQ2 += G2_NQUADSQ ) {
      _g2h_GetBFAPatchCurvesf ( domain, f, k,
                                &fc00, &fc01, &fc02, &fd00, &fd01, &fd02 );
      if ( !_g2h_TabNLDer0f ( domain, G2_NQUAD, tkn,
                        hfunc, dhfunc, ddhfunc, dddhfunc,
                        &nlpr->diu[kNQ2], &nlpr->div[kNQ2],
                        &nlpr->diuu[kNQ2], &nlpr->diuv[kNQ2], &nlpr->divv[kNQ2],
                        &nlpr->diuuu[kNQ2], &nlpr->diuuv[kNQ2],
                        &nlpr->diuvv[kNQ2], &nlpr->divvv[kNQ2],
                        fc00, fc01, fc02, fd00, fd01, fd02,
                        &nlpr->psiu[fN], &nlpr->psiv[fN],
                        &nlpr->psiuu[fN], &nlpr->psiuv[fN], &nlpr->psivv[fN], 
                        &nlpr->psiuuu[fN], &nlpr->psiuuv[fN],
                        &nlpr->psiuvv[fN], &nlpr->psivvv[fN] ) ) {
        domain->error_code = G2H_ERROR_NO_SCRATCH_MEMORY;
        goto failure;
      }
    }
  }

  bc00 = pkv_GetScratchMemf ( 2*(G2_CROSSDEGSUM+6)*hole_k );
  if ( !bc00 ) {
    domain->error_code = G2H_ERROR_NO_SCRATCH_MEMORY;
    goto failure;
  }
  bc01 = &bc00[(G2_CROSS00DEG+1)*hole_k];  bc02 = &bc01[(G2_CROSS01DEG+1)*hole_k];
  bc10 = &bc02[(G2_CROSS02DEG+1)*hole_k];  bc11 = &bc10[(G2_CROSS10DEG+1)*hole_k];
  bc12 = &bc11[(G2_CROSS11DEG+1)*hole_k];  bd00 = &bc12[(G2_CROSS12DEG+1)*hole_k];
  bd01 = &bd00[(G2_CROSS00DEG+1)*hole_k];  bd02 = &bd01[(G2_CROSS01DEG+1)*hole_k];
  bd10 = &bd02[(G2_CROSS02DEG+1)*hole_k];  bd11 = &bd10[(G2_CROSS10DEG+1)*hole_k];
  bd12 = &bd11[(G2_CROSS11DEG+1)*hole_k];

  memset ( bc00, 0, 2*(G2_CROSSDEGSUM+6)*hole_k*sizeof(float) );
  bfcpn = privateG->bfcpn;
  for ( f = 0; f < nfunc_b; f++ ) {
        /* find Coons representation of the constant part of the solution */
    bvz = nlpr->rhole_cp[bfcpn[f]].z;
    for ( k = 0;  k < hole_k;  k++ ) {
      _g2h_GetBFBPatchCurvesf ( domain, f, k,
                                &fc00, &fc01, &fc02, &fc10, &fc11, &fc12,
                                &fd00, &fd01, &fd02, &fd10, &fd11, &fd12 );
      pkn_AddMatrixMf ( 1, G2_CROSS00DEG+1, 0, &bc00[(G2_CROSS00DEG+1)*k], 0, fc00,
                        bvz, 0, &bc00[(G2_CROSS00DEG+1)*k] );
      pkn_AddMatrixMf ( 1, G2_CROSS01DEG+1, 0, &bc01[(G2_CROSS01DEG+1)*k], 0, fc01,
                        bvz, 0, &bc01[(G2_CROSS01DEG+1)*k] );
      pkn_AddMatrixMf ( 1, G2_CROSS02DEG+1, 0, &bc02[(G2_CROSS02DEG+1)*k], 0, fc02,
                        bvz, 0, &bc02[(G2_CROSS02DEG+1)*k] );
      pkn_AddMatrixMf ( 1, G2_CROSS10DEG+1, 0, &bc10[(G2_CROSS10DEG+1)*k], 0, fc10,
                        bvz, 0, &bc10[(G2_CROSS10DEG+1)*k] );
      pkn_AddMatrixMf ( 1, G2_CROSS11DEG+1, 0, &bc11[(G2_CROSS11DEG+1)*k], 0, fc11,
                        bvz, 0, &bc11[(G2_CROSS11DEG+1)*k] );
      pkn_AddMatrixMf ( 1, G2_CROSS12DEG+1, 0, &bc12[(G2_CROSS12DEG+1)*k], 0, fc12,
                        bvz, 0, &bc12[(G2_CROSS12DEG+1)*k] );
      pkn_AddMatrixMf ( 1, G2_CROSS00DEG+1, 0, &bd00[(G2_CROSS00DEG+1)*k], 0, fd00,
                        bvz, 0, &bd00[(G2_CROSS00DEG+1)*k] );
      pkn_AddMatrixMf ( 1, G2_CROSS01DEG+1, 0, &bd01[(G2_CROSS01DEG+1)*k], 0, fd01,
                        bvz, 0, &bd01[(G2_CROSS01DEG+1)*k] );
      pkn_AddMatrixMf ( 1, G2_CROSS02DEG+1, 0, &bd02[(G2_CROSS02DEG+1)*k], 0, fd02,
                        bvz, 0, &bd02[(G2_CROSS02DEG+1)*k] );
      pkn_AddMatrixMf ( 1, G2_CROSS10DEG+1, 0, &bd10[(G2_CROSS10DEG+1)*k], 0, fd10,
                        bvz, 0, &bd10[(G2_CROSS10DEG+1)*k] );
      pkn_AddMatrixMf ( 1, G2_CROSS11DEG+1, 0, &bd11[(G2_CROSS11DEG+1)*k], 0, fd11,
                        bvz, 0, &bd11[(G2_CROSS11DEG+1)*k] );
      pkn_AddMatrixMf ( 1, G2_CROSS12DEG+1, 0, &bd12[(G2_CROSS12DEG+1)*k], 0, fd12,
                        bvz, 0, &bd12[(G2_CROSS12DEG+1)*k] );
    }
  }
  for ( k = kNQ2 = 0, fN = nfunc_a*G2_NQUADSQ*hole_k;
        k < hole_k;
        k++, kNQ2 += G2_NQUADSQ, fN += G2_NQUADSQ ) {
    if ( !_g2h_TabNLDerf ( domain, G2_NQUAD, tkn, hfunc, dhfunc, ddhfunc, dddhfunc,
             &nlpr->diu[kNQ2], &nlpr->div[kNQ2],
             &nlpr->diuu[kNQ2], &nlpr->diuv[kNQ2], &nlpr->divv[kNQ2],
             &nlpr->diuuu[kNQ2], &nlpr->diuuv[kNQ2],
             &nlpr->diuvv[kNQ2], &nlpr->divvv[kNQ2],
             &bc00[(G2_CROSS00DEG+1)*k], &bc01[(G2_CROSS01DEG+1)*k],
             &bc02[(G2_CROSS02DEG+1)*k], &bc10[(G2_CROSS10DEG+1)*k],
             &bc11[(G2_CROSS11DEG+1)*k], &bc12[(G2_CROSS12DEG+1)*k],
             &bd00[(G2_CROSS00DEG+1)*k], &bd01[(G2_CROSS01DEG+1)*k],
             &bd02[(G2_CROSS02DEG+1)*k], &bd10[(G2_CROSS10DEG+1)*k],
             &bd11[(G2_CROSS11DEG+1)*k], &bd12[(G2_CROSS12DEG+1)*k],
             &nlpr->psiu[fN], &nlpr->psiv[fN],
             &nlpr->psiuu[fN], &nlpr->psiuv[fN], &nlpr->psivv[fN],
             &nlpr->psiuuu[fN], &nlpr->psiuuv[fN],
             &nlpr->psiuvv[fN], &nlpr->psivvv[fN] ) ) {
      domain->error_code = G2H_ERROR_NO_SCRATCH_MEMORY;
      goto failure;
    }
  }
  pkv_SetScratchMemTop ( sp );

/*
DrawFGraph ( hole_k, G2_NQUAD, nlpr->di, &nlpr->psiuu[3*hole_k*G2_NQUADSQ] );
*/
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
#undef N
} /*_g2h_TabNLBasisFunctionsf*/

/* ////////////////////////////////////////////////////////////////////////// */
void _g2h_IntFunc1af ( G2HNLFuncf *f, float *funct )
{
  float a, d, e0, e1, e2, e7, e8, pu2, pv2, pupv;

  f->B[0] = e2 = (float)(1.0 + (pv2 = f->pv*f->pv));
  f->B[1] = -(pupv = f->pu*f->pv);
  f->B[2] = e1 = (float)(1.0 + (pu2 = f->pu*f->pu));

  e7 = (float)(1.0 - 2.0*pu2 + pv2);
  e8 = (float)(1.0 + pu2 - 2.0*pv2);
  e0 = (float)(1.0 + pu2 + pv2);
  a = e0*e0;
  f->D = d = (float)(1.0/(a*a*e0*sqrt(e0)));  /* e0^(-11/5) */

  f->L.x = (float)((f->puuu*e2 -2.0*f->puuv*pupv +f->puvv*e1)*e0
         +f->pu*f->puu*(6.0*pupv*f->puv-f->pvv*e8/*(1.0+pu2-2.0*pv2)*/-3.0*f->puu*e2)
         -f->pu*f->puv*(2.0*f->puv*e8/*(1.0+pu2-2.0*pv2)*/+3.0*pupv*f->pvv)
         -3.0*f->pv*f->puv*(f->puu*e2+f->pvv));
  f->L.y = (float)((f->puuv*e2 -2.0*f->puvv*pupv +f->pvvv*e1)*e0
         +f->pv*f->pvv*(6.0*pupv*f->puv-f->puu*e7/*(1.0+pv2-2.0*pu2)*/-3.0*f->pvv*e1)
         -f->pv*f->puv*(2.0*f->puv*e7/*(1.0+pv2-2.0*pu2)*/+3.0*pupv*f->puu)
         -3.0*f->pu*f->puv*(f->pvv*e1+f->puu));

  f->BLT.x = f->B[0]*f->L.x + f->B[1]*f->L.y;
  f->BLT.y = f->B[1]*f->L.x + f->B[2]*f->L.y;
  f->LBLT = f->L.x*f->BLT.x + f->L.y*f->BLT.y;
/*
  *funct += 0.25*(f->L.x*(f->L.x*f->B[0] + 2.0*f->L.y*f->B[1]) + f->L.y*f->L.y*f->B[2])*d*f->jac;
*/
  *funct += (float)(0.25*f->LBLT*d*f->jac);
} /*_g2h_IntFunc1af*/

void _g2h_IntFunc1bf ( G2HNLFuncf *f, float *funct )
{
  double a, d, e0, e1, e2, e3, e4, e5, e6, e7, e8, pu2, pv2, pupv;

  f->B[0] = (float)(e2 = 1.0 + (pv2 = f->pv*f->pv));
  f->B[1] = (float)(-(pupv = f->pu*f->pv));
  f->B[2] = (float)(e1 = 1.0 + (pu2 = f->pu*f->pu));

  e3 = 1.0 + 3.0*pu2 + pv2;
  e4 = 1.0 + pu2 + 3.0*pv2;
  e5 = 1.0 + 3.0*pu2 - 2.0*pv2;
  e6 = 1.0 - 2.0*pu2 + 3.0*pv2;
  e7 = 1.0 - 2.0*pu2 + pv2;
  e8 = 1.0 + pu2 - 2.0*pv2;
  e0 = 1.0 + pu2 + pv2;
  a = e0*e0;
  f->D = (float)(d = 1.0/(a*a*e0*sqrt(e0)));  /* e0^(-11/5) */

  f->L.x = (float)((f->puuu*e2 -2.0*f->puuv*pupv +f->puvv*e1)*e0
         +f->pu*f->puu*(6.0*pupv*f->puv-f->pvv*e8-3.0*f->puu*e2)
         -f->pu*f->puv*(2.0*f->puv*e8+3.0*pupv*f->pvv)
         -3.0*f->pv*f->puv*(f->puu*e2+f->pvv));
  f->L.y = (float)((f->puuv*e2 -2.0*f->puvv*pupv +f->pvvv*e1)*e0
         +f->pv*f->pvv*(6.0*pupv*f->puv-f->puu*e7-3.0*f->pvv*e1)
         -f->pv*f->puv*(2.0*f->puv*e7+3.0*pupv*f->puu)
         -3.0*f->pu*f->puv*(f->pvv*e1+f->puu));

  f->Bpu[0] = 0.0;     f->Bpu[1] = -f->pv;  f->Bpu[2] = (float)(2.0*f->pu);
  f->Bpv[0] = (float)(2.0*f->pv);  f->Bpv[1] = -f->pu;  f->Bpv[2] = 0.0;
  f->Dpu = (float)(-11.0*f->pu*d/e0);
  f->Dpv = (float)(-11.0*f->pv*d/e0);
  f->Lpu.x = (float)(2.0*(f->puuu*f->pu*e2-f->puuv*f->pv*e3
                +f->puvv*f->pu*(e0+e1))
           -3.0*f->puu*f->puu*e2
           -(f->puu*f->pvv+2.0*f->puv*f->puv)*e5
           +6.0*pupv*f->puv*(2.0*f->puu-f->pvv));
  f->Lpv.x = (float)(2.0*(f->puuu*f->pv*(e0+e2)-f->puuv*f->pu*e4
                +f->puvv*f->pv*e1)
           -3.0*f->puv*(f->puu*e6+f->pvv*e1)
           +pupv*(f->puu*(4.0*f->pvv-6.0*f->puu)+8.0*f->puv*f->puv));
  f->Lpuu.x = (float)(-f->pu*(f->pvv*e8-6.0*pupv*f->puv)
            -(6.0*f->pu*f->puu +3.0*f->pv*f->puv)*e2);
  f->Lpuv.x = (float)(3.0*f->pu*pupv*(2.0*f->puu-f->pvv) -3.0*f->pv*(f->pvv+f->puu*e2)
            -4.0*f->pu*f->puv*e8);
  f->Lpuv.x = (float)(-3.0*f->pv*(f->puu*e7+f->pvv*e1)-4.0*f->pu*f->puv*e8);
  f->Lpvv.x = (float)(-3.0*f->pv*f->puv*e1-f->pu*f->puu*e8);
  f->Lpuuu.x = (float)(e2*e0);
  f->Lpuuv.x = (float)(-2.0*pupv*e0);
  f->Lpuvv.x = (float)(e1*e0);
  f->Lpvvv.x = 0.0;

  f->Lpu.y = (float)(2.0*(f->pvvv*f->pu*(e0+e1)-f->puvv*f->pv*e3
                +f->puuv*f->pu*e2)
           -3.0*f->puv*(f->pvv*e5+f->puu*e2)
           +pupv*(8.0*f->puv*f->puv+f->pvv*(4.0*f->puu-6.0*f->pvv)));
  f->Lpv.y = (float)(2.0*(f->pvvv*f->pv*e1-f->puvv*f->pu*e4
                +f->puuv*f->pv*(e0+e2))
           -3.0*f->pvv*f->pvv*e1
           -(f->puu*f->pvv+2.0*f->puv*f->puv)*e6
           +6.0*pupv*f->puv*(2.0*f->pvv-f->puu));
  f->Lpuu.y = (float)(-3.0*f->pu*f->puv*e2-f->pv*f->pvv*e7);
  f->Lpuv.y = (float)(-3.0*f->pu*(f->puu*e2+f->pvv*e8) -4.0*f->pv*f->puv*e7);
  f->Lpvv.y = (float)(-f->pv*(f->puu*e7-6.0*pupv*f->puv)
            -(6.0*f->pv*f->pvv+3.0*f->pu*f->puv)*e1);
  f->Lpuuu.y = 0.0;
  f->Lpuuv.y = (float)(e2*e0);
  f->Lpuvv.y = (float)(-2.0*pupv*e0);
  f->Lpvvv.y = (float)(e1*e0);

  f->BLT.x = f->B[0]*f->L.x + f->B[1]*f->L.y;
  f->BLT.y = f->B[1]*f->L.x + f->B[2]*f->L.y;
  f->LBLT = f->L.x*f->BLT.x + f->L.y*f->BLT.y;
/*
  *funct += 0.25*(f->L.x*(f->L.x*f->B[0] + 2.0*f->L.y*f->B[1]) + f->L.y*f->L.y*f->B[2])*d*f->jac;
*/
  *funct += (float)(0.25*f->LBLT*d*f->jac);
} /*_g2h_IntFunc1bf*/

void _g2h_IntFunc1cf ( G2HNLFuncf *f, float *funct )
{
  double a, d, e0, e1, e2, e3, e4, e5, e6, e7, e8, pu2, pv2, pupv;

  f->B[0] = (float)(e2 = 1.0 + (pv2 = f->pv*f->pv));
  f->B[1] = (float)(-(pupv = f->pu*f->pv));
  f->B[2] = (float)(e1 = 1.0 + (pu2 = f->pu*f->pu));

  e3 = 1.0 + 3.0*pu2 + pv2;
  e4 = 1.0 + pu2 + 3.0*pv2;
  e5 = 1.0 + 3.0*pu2 - 2.0*pv2;
  e6 = 1.0 - 2.0*pu2 + 3.0*pv2;
  e7 = 1.0 - 2.0*pu2 + pv2;
  e8 = 1.0 + pu2 - 2.0*pv2;
  e0 = 1.0 + pu2 + pv2;
  a = e0*e0;
  f->D = (float)(d = 1.0/(a*a*e0*sqrt(e0)));  /* e0^(-11/5) */

  f->L.x = (float)((f->puuu*e2 -2.0*f->puuv*pupv +f->puvv*e1)*e0
         +f->pu*f->puu*(6.0*pupv*f->puv-f->pvv*e8-3.0*f->puu*e2)
         -f->pu*f->puv*(2.0*f->puv*e8+3.0*pupv*f->pvv)
         -3.0*f->pv*f->puv*(f->puu*e2+f->pvv));
  f->L.y = (float)((f->puuv*e2 -2.0*f->puvv*pupv +f->pvvv*e1)*e0
         +f->pv*f->pvv*(6.0*pupv*f->puv-f->puu*e7-3.0*f->pvv*e1)
         -f->pv*f->puv*(2.0*f->puv*e7+3.0*pupv*f->puu)
         -3.0*f->pu*f->puv*(f->pvv*e1+f->puu));

  f->Bpu[0] = 0.0;     f->Bpu[1] = -f->pv;  f->Bpu[2] = (float)(2.0*f->pu);
  f->Bpv[0] = (float)(2.0*f->pv);  f->Bpv[1] = -f->pu;  f->Bpv[2] = 0.0;
  f->Dpu = (float)(-11.0*f->pu*d/e0);
  f->Dpv = (float)(-11.0*f->pv*d/e0);
  f->Lpu.x = (float)(2.0*(f->puuu*f->pu*e2-f->puuv*f->pv*e3
                +f->puvv*f->pu*(e0+e1))
           -3.0*f->puu*f->puu*e2
           -(f->puu*f->pvv+2.0*f->puv*f->puv)*e5
           +6.0*pupv*f->puv*(2.0*f->puu-f->pvv));
  f->Lpv.x = (float)(2.0*(f->puuu*f->pv*(e0+e2)-f->puuv*f->pu*e4
                +f->puvv*f->pv*e1)
           -3.0*f->puv*(f->puu*e6+f->pvv*e1)
           +pupv*(f->puu*(4.0*f->pvv-6.0*f->puu)+8.0*f->puv*f->puv));
  f->Lpuu.x = (float)(-f->pu*(f->pvv*e8-6.0*pupv*f->puv)
            -(6.0*f->pu*f->puu +3.0*f->pv*f->puv)*e2);
  f->Lpuv.x = (float)(3.0*f->pu*pupv*(2.0*f->puu-f->pvv) -3.0*f->pv*(f->pvv+f->puu*e2)
            -4.0*f->pu*f->puv*e8);
  f->Lpuv.x = (float)(-3.0*f->pv*(f->puu*e7+f->pvv*e1)-4.0*f->pu*f->puv*e8);
  f->Lpvv.x = (float)(-3.0*f->pv*f->puv*e1-f->pu*f->puu*e8);
  f->Lpuuu.x = (float)(e2*e0);
  f->Lpuuv.x = (float)(-2.0*pupv*e0);
  f->Lpuvv.x = (float)(e1*e0);
  f->Lpvvv.x = 0.0;

  f->Lpu.y = (float)(2.0*(f->pvvv*f->pu*(e0+e1)-f->puvv*f->pv*e3
                +f->puuv*f->pu*e2)
           -3.0*f->puv*(f->pvv*e5+f->puu*e2)
           +pupv*(8.0*f->puv*f->puv+f->pvv*(4.0*f->puu-6.0*f->pvv)));
  f->Lpv.y = (float)(2.0*(f->pvvv*f->pv*e1-f->puvv*f->pu*e4
                +f->puuv*f->pv*(e0+e2))
           -3.0*f->pvv*f->pvv*e1
           -(f->puu*f->pvv+2.0*f->puv*f->puv)*e6
           +6.0*pupv*f->puv*(2.0*f->pvv-f->puu));
  f->Lpuu.y = (float)(-3.0*f->pu*f->puv*e2-f->pv*f->pvv*e7);
  f->Lpuv.y = (float)(-3.0*f->pu*(f->puu*e2+f->pvv*e8) -4.0*f->pv*f->puv*e7);
  f->Lpvv.y = (float)(-f->pv*(f->puu*e7-6.0*pupv*f->puv)
            -(6.0*f->pv*f->pvv+3.0*f->pu*f->puv)*e1);
  f->Lpuuu.y = 0.0;
  f->Lpuuv.y = (float)(e2*e0);
  f->Lpuvv.y = (float)(-2.0*pupv*e0);
  f->Lpvvv.y = (float)(e1*e0);

  f->Bpupu[0] = 0.0;  f->Bpupu[1] =  0.0;  f->Bpupu[2] = 2.0;
  f->Bpupv[0] = 0.0;  f->Bpupv[1] = -1.0;  f->Bpupv[2] = 0.0;
  f->Bpvpv[0] = 2.0;  f->Bpvpv[1] =  0.0;  f->Bpvpv[2] = 0.0;

  f->Dpupu = (float)(-11.0*(1.0-12.0*pu2+pv2)*d/(e0*e0));
  f->Dpupv = (float)(143.0*pupv*d/(e0*e0));
  f->Dpvpv = (float)(-11.0*(1.0+pu2-12.0*pv2)*d/(e0*e0));

  f->Lpupu.x = (float)(2.0*f->puuu*e2 -12.0*pupv*f->puuv + 2.0*f->puvv*(2.0+6.0*pu2+pv2)
             -6.0*f->pu*f->puu*f->pvv +12.0*f->pv*f->puu*f->puv -12.0*f->pu*f->puv*f->puv - 6.0*f->pv*f->puv*f->pvv);
  f->Lpupv.x = (float)(4.0*pupv*(f->puuu+f->puvv) -2.0*f->puuv*(1.0+3.0*pu2+3.0*pv2)
             -f->puv*(6.0*f->pu*f->pvv-8.0*f->puv*f->pv) +f->puu*(12.0*f->pu*f->puv-6.0*f->puu*f->pv+4.0*f->pv*f->pvv));
  f->Lpupuu.x = (float)(12.0*pupv*f->puv-6.0*f->puu*e2-f->pvv*e5);
  f->Lpupuv.x = (float)(6.0*pupv*(2.0*f->puu-f->pvv)-4.0*f->puv*e5);
  f->Lpupvv.x = (float)(-f->puu*e5 - 6.0*pupv*f->puv);
  f->Lpupuuu.x = (float)(2.0*f->pu*e2);
  f->Lpupuuv.x = (float)(-2.0*f->pv*e3);
  f->Lpupuvv.x = (float)(2.0*f->pu*(e0+e1));
  f->Lpupvvv.x = 0.0;
  f->Lpvpv.x = (float)(f->puuu*(4.0+2.0*pu2+12.0*pv2) -12.0*pupv*f->puuv +2.0*f->puvv*e1
             +8.0*f->pu*f->puv*f->puv -18.0*f->pv*f->puu*f->puv -f->pu*f->puu*(6.0*f->puu-4.0*f->pvv));
  f->Lpvpuu.x = (float)(4.0*pupv*(f->pvv-3.0*f->puu) -3.0*f->puv*e6);
  f->Lpvpuv.x = (float)(-3.0*f->puu*e6 -3.0*f->pvv*e1
              +16.0*pupv*f->puv);
  f->Lpvpvv.x = (float)(-3.0*f->puv*e1+4.0*pupv*f->puu);
  f->Lpvpuuu.x = (float)(2.0*f->pv*(e0+e2));
  f->Lpvpuuv.x = (float)(-2.0*f->pu*e4);
  f->Lpvpuvv.x = (float)(2.0*f->pv*e1);
  f->Lpvpvvv.x = 0.0;
  f->Lpuupuu.x = (float)(-3.0*f->Lpupuuu.x);
  f->Lpuupuv.x = (float)(-3.0*f->pv*e7);
  f->Lpuupvv.x = (float)(-f->pu*e8);
  f->Lpuvpuv.x = (float)(4.0*f->Lpuupvv.x);
  f->Lpuvpvv.x = (float)(-1.5*f->Lpvpuvv.x);
  f->Lpvvpvv.x = 0.0;

  f->Lpupu.y = (float)(f->pvvv*(4.0+2.0*pv2+12.0*pu2) -12.0*pupv*f->puvv +2.0*f->puuv*e2
             +8.0*f->pv*f->puv*f->puv -18.0*f->pu*f->pvv*f->puv -f->pv*f->pvv*(6.0*f->pvv-4.0*f->puu));
  f->Lpupv.y = (float)(4.0*pupv*(f->pvvv+f->puuv) -2.0*f->puvv*(1.0+3.0*pv2+3.0*pu2)
             -f->puv*(6.0*f->pv*f->puu-8.0*f->puv*f->pu) +f->pvv*(12.0*f->pv*f->puv-6.0*f->pvv*f->pu+4.0*f->pu*f->puu));
  f->Lpupuu.y = (float)(-3.0*f->puv*e2+4.0*pupv*f->pvv);
  f->Lpupuv.y = (float)(-3.0*f->pvv*e5 -3.0*f->puu*e2
              +16.0*pupv*f->puv);
  f->Lpupvv.y = (float)(4.0*pupv*(f->puu-3.0*f->pvv) -3.0*f->puv*e5);
  f->Lpupuuu.y = 0.0;
  f->Lpupuuv.y = (float)(2.0*f->pu*e2);
  f->Lpupuvv.y = (float)(-2.0*f->pv*e3);
  f->Lpupvvv.y = f->Lpupuvv.x;
  f->Lpvpv.y = (float)(2.0*f->pvvv*e1 -12.0*pupv*f->puvv + 2.0*f->puuv*(2.0+6.0*pv2+pu2)
             -6.0*f->pv*f->pvv*f->puu +12.0*f->pu*f->pvv*f->puv -12.0*f->pv*f->puv*f->puv - 6.0*f->pu*f->puv*f->puu);
  f->Lpvpuu.y = (float)(-f->pvv*e6 - 6.0*pupv*f->puv);
  f->Lpvpuv.y = (float)(6.0*pupv*(2.0*f->pvv-f->puu)-4.0*f->puv*e6);
  f->Lpvpvv.y = (float)(12.0*pupv*f->puv-6.0*f->pvv*e1-f->puu*e6);
  f->Lpvpuuu.y = 0.0;
  f->Lpvpuuv.y = f->Lpvpuuu.x;
  f->Lpvpuvv.y = (float)(-2.0*f->pu*e4);
  f->Lpvpvvv.y = (float)(2.0*f->pv*e1);
  f->Lpuupuu.y = 0.0;
  f->Lpuupuv.y = (float)(-3.0*f->pu*e2);
  f->Lpuupvv.y = (float)(-f->pv*e7);
  f->Lpuvpuv.y = (float)(-4.0*f->pv*e7);
  f->Lpuvpvv.y = (float)(-3.0*f->pu*e8);
  f->Lpvvpvv.y = (float)(-6.0*f->pv*e1);

  f->BLT.x = f->B[0]*f->L.x + f->B[1]*f->L.y;
  f->BLT.y = f->B[1]*f->L.x + f->B[2]*f->L.y;
  f->LBLT = f->L.x*f->BLT.x + f->L.y*f->BLT.y;
/*
  *funct += 0.25*(f->L.x*(f->L.x*f->B[0] + 2.0*f->L.y*f->B[1]) + f->L.y*f->L.y*f->B[2])*d*f->jac;
*/
  *funct += (float)(0.25*f->LBLT*d*f->jac);
} /*_g2h_IntFunc1cf*/

void _g2h_IntFunc2bf ( G2HNLFuncf *f, float *grad )
{
  float    di, Bi[3];
  vector2f Li, BiLT;

  Li.x = f->Lpu.x*f->psiu + f->Lpv.x*f->psiv + f->Lpuu.x*f->psiuu + f->Lpuv.x*f->psiuv + f->Lpvv.x*f->psivv
          + f->Lpuuu.x*f->psiuuu + f->Lpuuv.x*f->psiuuv + f->Lpuvv.x*f->psiuvv;
  Li.y = f->Lpu.y*f->psiu + f->Lpv.y*f->psiv + f->Lpuu.y*f->psiuu + f->Lpuv.y*f->psiuv + f->Lpvv.y*f->psivv
          + f->Lpuuv.y*f->psiuuv + f->Lpuvv.y*f->psiuvv + f->Lpvvv.y*f->psivvv;
  Bi[0] = f->Bpv[0]*f->psiv;
  Bi[1] = f->Bpu[1]*f->psiu + f->Bpv[1]*f->psiv;
  Bi[2] = f->Bpu[2]*f->psiu;
  BiLT.x = Bi[0]*f->L.x + Bi[1]*f->L.y;
  BiLT.y = Bi[1]*f->L.x + Bi[2]*f->L.y;
  di = f->Dpu*f->psiu + f->Dpv*f->psiv;

/*
  *grad += 0.25*(2.0*(Li->x*f->L.x*f->B[0] + (Li->x*f->L.y+Li->y*f->L.x)*f->B[1]
                      + Li->y*f->L.y*f->B[2])*f->D
           +(f->L.x*(f->L.x*Bi[0] + 2.0*f->L.y*Bi[1]) + f->L.y*f->L.y*Bi[2])*f->D
           +(f->L.x*(f->L.x*f->B[0] + 2.0*f->L.y*f->B[1]) + f->L.y*f->L.y*f->B[2])*di)*f->jac;
*/
  *grad += (float)(0.25*((2.0*(Li.x*f->BLT.x+Li.y*f->BLT.y)
                  +(f->L.x*BiLT.x + f->L.y*BiLT.y))*f->D + f->LBLT*di)*f->jac);
} /*_g2h_IntFunc2bf*/

void _g2h_IntFunc2cf ( G2HNLFuncf *f, vector2f *Li, float *Bi, vector2f *BiLT,
                       float *Di, float *grad )
{
  float di;

  Li->x = f->Lpu.x*f->psiu + f->Lpv.x*f->psiv + f->Lpuu.x*f->psiuu + f->Lpuv.x*f->psiuv + f->Lpvv.x*f->psivv
          + f->Lpuuu.x*f->psiuuu + f->Lpuuv.x*f->psiuuv + f->Lpuvv.x*f->psiuvv;
  Li->y = f->Lpu.y*f->psiu + f->Lpv.y*f->psiv + f->Lpuu.y*f->psiuu + f->Lpuv.y*f->psiuv + f->Lpvv.y*f->psivv
          + f->Lpuuv.y*f->psiuuv + f->Lpuvv.y*f->psiuvv + f->Lpvvv.y*f->psivvv;
  Bi[0] = f->Bpv[0]*f->psiv;
  Bi[1] = f->Bpu[1]*f->psiu + f->Bpv[1]*f->psiv;
  Bi[2] = f->Bpu[2]*f->psiu;
  BiLT->x = Bi[0]*f->L.x + Bi[1]*f->L.y;
  BiLT->y = Bi[1]*f->L.x + Bi[2]*f->L.y;
  di = *Di = f->Dpu*f->psiu + f->Dpv*f->psiv;

/*
  *grad += 0.25*(2.0*(Li->x*f->L.x*f->B[0] + (Li->x*f->L.y+Li->y*f->L.x)*f->B[1]
                      + Li->y*f->L.y*f->B[2])*f->D
           +(f->L.x*(f->L.x*Bi[0] + 2.0*f->L.y*Bi[1]) + f->L.y*f->L.y*Bi[2])*f->D
           +(f->L.x*(f->L.x*f->B[0] + 2.0*f->L.y*f->B[1]) + f->L.y*f->L.y*f->B[2])*di)*f->jac;
*/
  *grad += (float)(0.25*((2.0*(Li->x*f->BLT.x+Li->y*f->BLT.y)
                  +(f->L.x*BiLT->x + f->L.y*BiLT->y))*f->D + f->LBLT*di)*f->jac);
} /*_g2h_IntFunc2cf*/

void _g2h_IntFunc3cf ( G2HNLFuncf *f, vector2f *Li, vector2f *Lj,
                       vector2f *BiLT, vector2f *BjLT,
                       float Di, float Dj, float *hessian )
{
  vector2f Lij;
  float    Bij[3], Dij;

  Bij[0] = f->Bpvpv[0]*f->psiv*f->psjv;
  Bij[1] = f->Bpupv[1]*(f->psiu*f->psjv+f->psju*f->psiv);
  Bij[2] = f->Bpupu[2]*f->psiu*f->psju;

  Dij = f->Dpupu*f->psiu*f->psju +
        f->Dpupv*(f->psiu*f->psjv+f->psju*f->psiv) +
        f->Dpvpv*f->psiv*f->psjv;

/*
  Lij.x = f->Lpupu.x*f->psiu*f->psju
         + f->Lpupv.x*(f->psiu*f->psjv+f->psju*f->psiv)
         + f->Lpupuu.x*(f->psiu*f->psjuu+f->psju*f->psiuu)
         + f->Lpupuv.x*(f->psiu*f->psjuv+f->psju*f->psiuv)
         + f->Lpupvv.x*(f->psiu*f->psjvv+f->psju*f->psivv)
         + f->Lpupuuu.x*(f->psiu*f->psjuuu+f->psju*f->psiuuu)
         + f->Lpupuuv.x*(f->psiu*f->psjuuv+f->psju*f->psiuuv)
         + f->Lpupuvv.x*(f->psiu*f->psjuvv+f->psju*f->psiuvv)
         + f->Lpvpv.x*f->psiv*f->psjv
         + f->Lpvpuu.x*(f->psiv*f->psjuu+f->psjv*f->psiuu)
         + f->Lpvpuv.x*(f->psiv*f->psjuv+f->psjv*f->psiuv)
         + f->Lpvpvv.x*(f->psiv*f->psjvv+f->psjv*f->psivv)
         + f->Lpvpuuu.x*(f->psiv*f->psjuuu+f->psjv*f->psiuuu)
         + f->Lpvpuuv.x*(f->psiv*f->psjuuv+f->psjv*f->psiuuv)
         + f->Lpvpuvv.x*(f->psiv*f->psjuvv+f->psjv*f->psiuvv)
         + f->Lpuupuu.x*f->psiuu*f->psjuu
         + f->Lpuupuv.x*(f->psiuu*f->psjuv+f->psjuu*f->psiuv)
         + f->Lpuupvv.x*(f->psiuu*f->psjvv+f->psjuu*f->psivv)
         + f->Lpuvpuv.x*f->psiuv*f->psjuv
         + f->Lpuvpvv.x*(f->psiuv*f->psjvv+f->psjuv*f->psivv);
  Lij.y = f->Lpupu.y*f->psiu*f->psju
         + f->Lpupv.y*(f->psiu*f->psjv+f->psju*f->psiv)
         + f->Lpupuu.y*(f->psiu*f->psjuu+f->psju*f->psiuu)
         + f->Lpupuv.y*(f->psiu*f->psjuv+f->psju*f->psiuv)
         + f->Lpupvv.y*(f->psiu*f->psjvv+f->psju*f->psivv)
         + f->Lpupuuv.y*(f->psiu*f->psjuuv+f->psju*f->psiuuv)
         + f->Lpupuvv.y*(f->psiu*f->psjuvv+f->psju*f->psiuvv)
         + f->Lpupvvv.y*(f->psiu*f->psjvvv+f->psju*f->psivvv)
         + f->Lpvpv.y*f->psiv*f->psjv
         + f->Lpvpuu.y*(f->psiv*f->psjuu+f->psjv*f->psiuu)
         + f->Lpvpuv.y*(f->psiv*f->psjuv+f->psjv*f->psiuv)
         + f->Lpvpvv.y*(f->psiv*f->psjvv+f->psjv*f->psivv)
         + f->Lpvpuuv.y*(f->psiv*f->psjuuv+f->psjv*f->psiuuv)
         + f->Lpvpuvv.y*(f->psiv*f->psjuvv+f->psjv*f->psiuvv)
         + f->Lpvpvvv.y*(f->psiv*f->psjvvv+f->psjv*f->psivvv)
         + f->Lpuupuv.y*(f->psiuu*f->psjuv+f->psjuu*f->psiuv)
         + f->Lpuupvv.y*(f->psiuu*f->psjvv+f->psjuu*f->psivv)
         + f->Lpuvpuv.y*f->psiuv*f->psjuv
         + f->Lpuvpvv.y*(f->psiuv*f->psjvv+f->psjuv*f->psivv)
         + f->Lpvvpvv.y*f->psivv*f->psjvv;
*/
  Lij.x = f->psiu*(f->Lpupu.x*f->psju
                 + f->Lpupv.x*f->psjv
                 + f->Lpupuu.x*f->psjuu
                 + f->Lpupuv.x*f->psjuv
                 + f->Lpupvv.x*f->psjvv
                 + f->Lpupuuu.x*f->psjuuu
                 + f->Lpupuuv.x*f->psjuuv
                 + f->Lpupuvv.x*f->psjuvv)
        + f->psju*(f->Lpupv.x*f->psiv
                 + f->Lpupuu.x*f->psiuu
                 + f->Lpupuv.x*f->psiuv
                 + f->Lpupvv.x*f->psivv
                 + f->Lpupuuu.x*f->psiuuu
                 + f->Lpupuuv.x*f->psiuuv
                 + f->Lpupuvv.x*f->psiuvv)
        + f->psiv*(f->Lpvpv.x*f->psjv
                 + f->Lpvpuu.x*f->psjuu
                 + f->Lpvpuv.x*f->psjuv
                 + f->Lpvpvv.x*f->psjvv
                 + f->Lpvpuuu.x*f->psjuuu
                 + f->Lpvpuuv.x*f->psjuuv
                 + f->Lpvpuvv.x*f->psjuvv)
        + f->psjv*(f->Lpvpuu.x*f->psiuu
                 + f->Lpvpuv.x*f->psiuv
                 + f->Lpvpvv.x*f->psivv
                 + f->Lpvpuuu.x*f->psiuuu
                 + f->Lpvpuuv.x*f->psiuuv
                 + f->Lpvpuvv.x*f->psiuvv)
       + f->psiuu*(f->Lpuupuu.x*f->psjuu
                 + f->Lpuupuv.x*f->psjuv
                 + f->Lpuupvv.x*f->psjvv)
       + f->psiuv*(f->Lpuvpuv.x*f->psjuv
                 + f->Lpuvpvv.x*f->psjvv)
       + f->psjuu*(f->Lpuupuv.x*f->psiuv
                 + f->Lpuupvv.x*f->psivv)
       + f->Lpuvpvv.x*f->psjuv*f->psivv;

  Lij.y = f->psiu*(f->Lpupu.y*f->psju
                 + f->Lpupv.y*f->psjv
                 + f->Lpupuu.y*f->psjuu
                 + f->Lpupuv.y*f->psjuv
                 + f->Lpupvv.y*f->psjvv
                 + f->Lpupuuv.y*f->psjuuv
                 + f->Lpupuvv.y*f->psjuvv
                 + f->Lpupvvv.y*f->psjvvv)
        + f->psju*(f->Lpupv.y*f->psiv
                 + f->Lpupuu.y*f->psiuu
                 + f->Lpupuv.y*f->psiuv
                 + f->Lpupvv.y*f->psivv
                 + f->Lpupuuv.y*f->psiuuv
                 + f->Lpupuvv.y*f->psiuvv
                 + f->Lpupvvv.y*f->psivvv)
        + f->psiv*(f->Lpvpv.y*f->psjv
                 + f->Lpvpuu.y*f->psjuu
                 + f->Lpvpuv.y*f->psjuv
                 + f->Lpvpvv.y*f->psjvv
                 + f->Lpvpuuv.y*f->psjuuv
                 + f->Lpvpuvv.y*f->psjuvv
                 + f->Lpvpvvv.y*f->psjvvv)
        + f->psjv*(f->Lpvpuu.y*f->psiuu
                 + f->Lpvpuv.y*f->psiuv
                 + f->Lpvpvv.y*f->psivv
                 + f->Lpvpuuv.y*f->psiuuv
                 + f->Lpvpuvv.y*f->psiuvv
                 + f->Lpvpvvv.y*f->psivvv)
       + f->psiuu*(f->Lpuupuv.y*f->psjuv
                 + f->Lpuupvv.y*f->psjvv)
       + f->psjuu*(f->Lpuupuv.y*f->psiuv
                 + f->Lpuupvv.y*f->psivv)
       + f->psiuv*(f->Lpuvpuv.y*f->psjuv
                 + f->Lpuvpvv.y*f->psjvv)
       + f->Lpuvpvv.y*f->psjuv*f->psivv
       + f->Lpvvpvv.y*f->psivv*f->psjvv;
/*
  *hessian += 0.25*(
    (2.0*(Lij.x*f->L.x*f->B[0] + (Lij.x*f->L.y+Lij.y*f->L.x)*f->B[1] + Lij.y*f->L.y*f->B[2]
        +Li->x*f->L.x*Bj[0] + (Li->x*f->L.y+Li->y*f->L.x)*Bj[1] + Li->y*f->L.y*Bj[2]
        +Lj->x*f->L.x*Bi[0] + (Lj->x*f->L.y+Lj->y*f->L.x)*Bi[1] + Lj->y*f->L.y*Bi[2]
        +Li->x*Lj->x*f->B[0] + (Li->x*Lj->y+Li->y*Lj->x)*f->B[1] + Li->y*Lj->y*f->B[2])
      +f->L.x*f->L.x*Bij[0] + (f->L.x*f->L.y+f->L.y*f->L.x)*Bij[1] + f->L.y*f->L.y*Bij[2])*f->D
    +(2.0*(Li->x*f->L.x*f->B[0] + (Li->x*f->L.y+Li->y*f->L.x)*f->B[1] + Li->y*f->L.y*f->B[2])
        +f->L.x*f->L.x*Bi[0] + (f->L.x*f->L.y+f->L.y*f->L.x)*Bi[1] + f->L.y*f->L.y*Bi[2])*Dj
    +(2.0*(Lj->x*f->L.x*f->B[0] + (Lj->x*f->L.y+Lj->y*f->L.x)*f->B[1] + Lj->y*f->L.y*f->B[2])
        +f->L.x*f->L.x*Bj[0] + (f->L.x*f->L.y+f->L.y*f->L.x)*Bj[1] + f->L.y*f->L.y*Bj[2])*Di
    +(f->L.x*f->L.x*f->B[0] + (f->L.x*f->L.y+f->L.y*f->L.x)*f->B[1] + f->L.y*f->L.y*f->B[2])*Dij)*f->jac;
*/
  *hessian += (float)(0.25*(
    (2.0*(Lij.x*f->BLT.x+Lij.y*f->BLT.y
        +Li->x*BjLT->x+Li->y*BjLT->y +Lj->x*BiLT->x+Lj->y*BiLT->y
        +Li->x*Lj->x*f->B[0] + (Li->x*Lj->y+Li->y*Lj->x)*f->B[1] + Li->y*Lj->y*f->B[2])
      +f->L.x*f->L.x*Bij[0] + (f->L.x*f->L.y+f->L.y*f->L.x)*Bij[1] + f->L.y*f->L.y*Bij[2])*f->D
    +(2.0*(Li->x*f->BLT.x+Li->y*f->BLT.y) + f->L.x*BiLT->x+f->L.y*BiLT->y)*Dj
    +(2.0*(Lj->x*f->BLT.x+Lj->y*f->BLT.y) + f->L.x*BjLT->x+f->L.y*BjLT->y)*Di
    +f->LBLT*Dij)*f->jac);
} /*_g2h_IntFunc3cf*/

