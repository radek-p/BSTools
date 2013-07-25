
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2005, 2012                            */
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
#include "eg1hprivatef.h"
#include "eg1herror.h"


G1HNLPrivatef *_g1h_nlprivf;

/* ///////////////////////////////////////////////////////////////////////// */
void g1h_ReflectVectorsf ( int n, const vector3f *v, vector3f *w )
{
  vector3f r;
  double   a;
  int      i;

  r = _g1h_nlprivf->reflv;
  for ( i = 0; i < n; i++ ) {
    a = DotProduct3f ( &v[i], &r );
    AddVector3Mf ( &v[i], &r, -2.0*a, &w[i] );
  }
} /*g1h_ReflectVectorsf*/

void g1h_nonlinoutpatchf ( int n, int m, const float *cp, void *usrptr )
{
#define N (n+1)*(m+1)
  memcpy ( &_g1h_nlprivf->nldi[_g1h_nlprivf->auxc*N], cp, N*sizeof(vector3f) );
  _g1h_nlprivf->auxc ++;
#undef N
} /*g1h_nonlinoutpatchf*/

boolean _g1h_StopItf ( int itn, float gn0, float gn,
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
} /*_g1h_StopItf*/

/* ///////////////////////////////////////////////////////////////////////// */
boolean g1h_GetHoleSurrndPatchf ( GHoleDomainf *domain,
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
  if ( !ind || !q )
    goto failure;

  hole_k = domain->hole_k;
  gh_GetBspInd ( hole_k, i, j, ind );
  for ( k = 0; k < 16; k++ )
    q[k] = hole_cp[ind[k]];
  ukn = &domain->hole_knots[11*((i+hole_k-1) % hole_k)+3 ];
  vkn = &domain->hole_knots[11*i+j];
  mbs_BSPatchToBezf ( 3, 3, 7, ukn, 3, 7, vkn, 12, (float*)q,
                      NULL, NULL, NULL, NULL, NULL, NULL, 12, (float*)bcp );

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*g1h_GetHoleSurrndPatchf*/

boolean g1h_ComputeNLNormalf ( GHoleDomainf *domain,
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
  if ( !bcp )
    goto failure;

  SetVector3f ( &nv, 0.0, 0.0, 0.0 );
  for ( i = 0; i < hole_k; i++ )
    for ( j = 1; j < 3; j++ ) {
      g1h_GetHoleSurrndPatchf ( domain, hole_cp, i, j, bcp );
          /* integrate the unit normal vector at the boundary */
      for ( k = 0; k < DENSITY; k++ ) {
        t = (float)(k+k+1)/(float)(2*DENSITY);
        mbs_BCHornerDerP3f ( 3, 3, bcp, 0.0, t, &p, &pu, &pv );
        NormalizeVector3f ( &pv );
        CrossProduct3f ( &pu, &pv, &nu );
        AddVector3f ( &nv, &nu, &nv );
      }
    }
  NormalizeVector3f ( &nv );

        /* now a verification */
  for ( i = 0; i < hole_k; i++ )
    for ( j = 1; j < 3; j++ ) {
      g1h_GetHoleSurrndPatchf ( domain, hole_cp, i, j, bcp );
          /* integrate the unit normal vector at the boundary */
      for ( k = 0; k < DENSITY; k++ ) {
        t = (float)(k+k+1)/(float)(2*DENSITY);
        mbs_BCHornerDerP3f ( 3, 3, bcp, 0.0, t, &p, &pu, &pv );
        CrossProduct3f ( &pu, &pv, &nu );
        if ( DotProduct3f ( &nu, &nv ) <= 0.0 )
          goto failure;
      }
    }
  *anv = nv;

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
#undef DENSITY
} /*g1h_ComputeNLNormalf*/

boolean _g1h_ComputeNLNormalf ( GHoleDomainf *domain,
                                G1HNLPrivatef *nlprivate,
                                const point3f *hole_cp )
{
  vector3f nv;

  if ( g1h_ComputeNLNormalf ( domain, hole_cp, &nv ) ) {
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
} /*_g1h_ComputeNLNormalf*/

boolean _g1h_TabNLDer0f ( int nkn, const float *tkn,
             const float *hfunc, const float *dhfunc, const float *ddhfunc,
             vector2f *diu, vector2f *div,
             vector2f *diuu, vector2f *diuv, vector2f *divv,
             float *fc00, float *fc01, float *fd00, float *fd01,
             float *psiu, float *psiv,
             float *psiuu, float *psiuv, float *psivv )
{
  void  *sp;
  int   i;
  float *hu, *hv, *huu, *huv, *hvv;

  sp = pkv_GetScratchMemTop ();
  hu = pkv_GetScratchMemf ( 5*nkn*nkn );
  if ( !hu )
    goto failure;
  hv = &hu[nkn*nkn];    huu = &hv[nkn*nkn];
  huv = &huu[nkn*nkn];  hvv = &huv[nkn*nkn];

  if ( !mbs_TabBezC1Coons0Der2f ( 1, nkn, tkn, hfunc, dhfunc, ddhfunc,
                     nkn, tkn, hfunc, dhfunc, ddhfunc,
                     G1_CROSS00DEG, fc00, G1_CROSS01DEG, fc01,
                     G1_CROSS00DEG, fd00, G1_CROSS01DEG, fd01,
                     NULL, hu, hv, huu, huv, hvv ) )
    goto failure;

  for ( i = 0; i < nkn*nkn; i++ )
    pkn_Comp2iDerivatives2f ( diu[i].x, diu[i].y, div[i].x, div[i].y,
          diuu[i].x, diuu[i].y, diuv[i].x, diuv[i].y, divv[i].x, divv[i].y,
          1, &hu[i], &hv[i], &huu[i], &huv[i], &hvv[i],
          &psiu[i], &psiv[i], &psiuu[i], &psiuv[i], &psivv[i] );

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*_g1h_TabNLDer0f*/

boolean _g1h_TabNLDerf ( int nkn, float *tkn,
             const float *hfunc, const float *dhfunc, const float *ddhfunc,
             vector2f *diu, vector2f *div,
             vector2f *diuu, vector2f *diuv, vector2f *divv,
             float *fc00, float *fc01, float *fc10, float *fc11,
             float *fd00, float *fd01, float *fd10, float *fd11,
             float *psiu, float *psiv,
             float *psiuu, float *psiuv, float *psivv )
{
  void  *sp;
  int   i;
  float *hu, *hv, *huu, *huv, *hvv;

  sp = pkv_GetScratchMemTop ();
  hu = pkv_GetScratchMemf ( 5*nkn*nkn );
  if ( !hu )
    goto failure;
  hv = &hu[nkn*nkn];    huu = &hv[nkn*nkn];
  huv = &huu[nkn*nkn];  hvv = &huv[nkn*nkn];

  if ( !mbs_TabBezC1CoonsDer2f ( 1, nkn, tkn, hfunc, dhfunc, ddhfunc,
            nkn, tkn, hfunc, dhfunc, ddhfunc,
            G1_CROSS00DEG, fc00, G1_CROSS01DEG, fc01,
            G1_CROSS10DEG, fc10, G1_CROSS11DEG, fc11,
            G1_CROSS00DEG, fd00, G1_CROSS01DEG, fd01,
            G1_CROSS10DEG, fd10, G1_CROSS11DEG, fd11,
            NULL, hu, hv, huu, huv, hvv ) )
    goto failure;

  for ( i = 0; i < nkn*nkn; i++ )
    pkn_Comp2iDerivatives2f ( diu[i].x, diu[i].y, div[i].x, div[i].y,
          diuu[i].x, diuu[i].y, diuv[i].x, diuv[i].y, divv[i].x, divv[i].y,
          1, &hu[i], &hv[i], &huu[i], &huv[i], &hvv[i],
          &psiu[i], &psiv[i], &psiuu[i], &psiuv[i], &psivv[i] );

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*_g1h_TabNLDerf*/

boolean _g1h_TabNLBasisFunctionsf ( GHoleDomainf *domain, int nkn,
                                    G1HNLPrivatef *nlpr )
{
#define N ((G1H_FINALDEG+1)*(G1H_FINALDEG+1))
  void     *sp;
  GHolePrivateRecf  *privateG;
  G1HolePrivateRecf *privateG1;
  float    *tkn;
  float    *fc00, *fc01, *fc10, *fc11, *fd00, *fd01, *fd10, *fd11,
           *bc00, *bc01, *bc10, *bc11, *bd00, *bd01, *bd10, *bd11;
  float    *hfunc, *dhfunc, *ddhfunc;
  int      i, j, k, l, kN, kNQ2, nkn2, hole_k, f, fN, nfunc_a, nfunc_b;
  point2f  *di, ddi;
  float    bvz;
  unsigned char *bfcpn;

  sp      = pkv_GetScratchMemTop ();
  hole_k  = domain->hole_k;
  privateG  = domain->privateG;
  privateG1 = domain->privateG1;
  nfunc_a = privateG1->nfunc_a;
  nfunc_b = privateG1->nfunc_b;
  di      = pkv_GetScratchMem ( N*sizeof(point2f) );
  nkn2    = nkn*nkn;
  tkn     = pkv_GetScratchMemf ( nkn );
  hfunc   = pkv_GetScratchMemf ( 12*nkn );
  if ( !tkn || !di || !hfunc )
    goto failure;

  dhfunc = &hfunc[4*nkn];
  ddhfunc = &dhfunc[4*nkn];

  _gh_PrepareTabKnotsf ( nkn, privateG1->opt_quad, tkn );
  mbs_TabCubicHFuncDer2f ( 0.0, 1.0, nkn, tkn, hfunc, dhfunc, ddhfunc );

  for ( k = kN = kNQ2 = 0;
        k < hole_k;
        k++, kN += N, kNQ2 += nkn2 ) {
    pkv_Selectf ( N, 2, 3, 2, &nlpr->nldi[kN], di );
    for ( i = l = 0;  i < nkn;  i++ )
         /* ***** this is not an optimal algorithm, to be improved */
      for ( j = 0;  j < nkn;  j++, l++ ) {
        mbs_BCHornerDer2Pf ( G1H_FINALDEG, G1H_FINALDEG, 2, (float*)di,
                    tkn[i], tkn[j], (float*)&ddi,
                    (float*)&nlpr->diu[kNQ2+l], (float*)&nlpr->div[kNQ2+l],
                    (float*)&nlpr->diuu[kNQ2+l], (float*)&nlpr->diuv[kNQ2+l],
                    (float*)&nlpr->divv[kNQ2+l] );
        nlpr->jac[kNQ2+l] = (float)det2f ( &nlpr->diu[kNQ2+l], &nlpr->div[kNQ2+l] );
     }
  }

  /* verify the Jacobian */
  if ( nlpr->jac[0] > 0.0 ) {
    for ( i = 1; i < hole_k*nkn2; i++ )
      if ( nlpr->jac[i] <= 0.0 )
        goto failure;
  }
  else {
    for ( i = 0; i < hole_k*nkn2; i++ ) {
      if ( nlpr->jac[i] >= 0.0 )
        goto failure;
      nlpr->jac[i] = -nlpr->jac[i];
    }
  }

  for ( f = fN = 0;  f < nfunc_a;  f++ ) {
    for ( k = kNQ2 = 0;  k < hole_k;  k++, fN += nkn2, kNQ2 += nkn2 ) {
      _g1h_GetBFAPatchCurvesf ( domain, f, k, &fc00, &fc01, &fd00, &fd01 );
      _g1h_TabNLDer0f ( nkn, tkn,
                        hfunc, dhfunc, ddhfunc,
                        &nlpr->diu[kNQ2], &nlpr->div[kNQ2],
                        &nlpr->diuu[kNQ2], &nlpr->diuv[kNQ2], &nlpr->divv[kNQ2],
                        fc00, fc01, fd00, fd01,
                        &nlpr->psiu[fN], &nlpr->psiv[fN],
                        &nlpr->psiuu[fN], &nlpr->psiuv[fN], &nlpr->psivv[fN] );
    }
  }

  bc00 = pkv_GetScratchMemf ( 2*(G1_CROSSDEGSUM+4)*hole_k );
  if ( !bc00 )
    goto failure;
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
    _g1h_TabNLDerf ( nkn, tkn, hfunc, dhfunc, ddhfunc,
             &nlpr->diu[kNQ2], &nlpr->div[kNQ2],
             &nlpr->diuu[kNQ2], &nlpr->diuv[kNQ2], &nlpr->divv[kNQ2],
             &bc00[(G1_CROSS00DEG+1)*k], &bc01[(G1_CROSS01DEG+1)*k],
             &bc10[(G1_CROSS10DEG+1)*k], &bc11[(G1_CROSS11DEG+1)*k],
             &bd00[(G1_CROSS00DEG+1)*k], &bd01[(G1_CROSS01DEG+1)*k],
             &bd10[(G1_CROSS10DEG+1)*k], &bd11[(G1_CROSS11DEG+1)*k],
             &nlpr->psiu[fN], &nlpr->psiv[fN],
             &nlpr->psiuu[fN], &nlpr->psiuv[fN], &nlpr->psivv[fN] );
  }
  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
#undef N
} /*_g1h_TabNLBasisFunctionsf*/

void _g1h_IntFunc1f ( float pu, float pv,
                      float puu, float puv, float pvv, float jac,
                      float *c, float *cs, float *a, float *b,
                      float *funct )
{
  float A, B, C, CS;

  *c = C = (float)(1.0+pu*pu+pv*pv);
  *cs = CS = (float)sqrt(C);
  *a = A = (float)(0.5*(puu*(1.0+pv*pv)+pvv*(1.0+pu*pu)-2.0*puv*pu*pv));
  *b = B = (float)(1.0/(C*C*CS));
  *funct += A*A*B*jac;
} /*_g1h_IntFunc1f*/

void _g1h_IntFunc2f ( float pu, float pv,
                      float puu, float puv, float pvv,
                      float psiu, float psiv,
                      float psiuu, float psiuv, float psivv,
                      float jac, float c, float A, float B,
                      float *ai, float *bi, float *grad )
{
  float Ai, Bi;

  *ai = Ai = (float)((pu*pvv-puv*pv)*psiu - (pu*puv-puu*pv)*psiv
           +0.5*((1.0+pv*pv)*psiuu+(1.0+pu*pu)*psivv) - pu*pv*psiuv);
  *bi = Bi = (float)(-5.0*(pu*psiu+pv*psiv)*B/c);
  *grad += (float)(A*(2.0*Ai*B+A*Bi)*jac);
} /*_g1h_IntFunc2f*/

void _g1h_IntFunc3f ( float pu, float pv,
                      float puu, float puv, float pvv,
                      float psiu, float psiv,
                      float psiuu, float psiuv, float psivv,
                      float psju, float psjv,
                      float psjuu, float psjuv, float psjvv,
                      float jac, float c,
                      float A, float B, float Ai, float Bi,
                      float Aj, float Bj,
                      float *hessian )
{
  float Aij, Bij;

  Aij = pu*(psiu*psjvv+psivv*psju-psiv*psjuv-psiuv*psjv)+
        pv*(psiv*psjuu+psiuu*psjv-psiu*psjuv-psiuv*psju)+
        puu*psiv*psjv-puv*(psiu*psjv+psiv*psju)+pvv*psiu*psju;
/*
  Bij = (float)(-5.0*((1.0-7.0*pu*pv)*(psiu*psju+psiv*psjv)+
              (pv*pv-6.0*pu*pu)*psiu*psju+
              (pu*pu-6.0*pv*pv)*psiv*psjv)*B/(c*c));
*/
  Bij = (float)(5.0*((6.0*pu*pu-pv*pv-1.0)*psiu*psju
                    +7.0*pu*pv*(psiu*psjv+psju*psiv)
                    +(6.0*pv*pv-pu*pu-1.0)*psiv*psjv)*B/(c*c));
  *hessian += (float)((2.0*A*(Aij*B+Ai*Bj+Aj*Bi)+2.0*Ai*Aj*B+A*A*Bij)*jac);
} /*_g1h_IntFunc3f*/

