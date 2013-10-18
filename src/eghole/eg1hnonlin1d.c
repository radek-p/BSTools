
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2005, 2013                            */
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
#include "eg1hprivated.h"
#include "eg1herror.h"


G1HNLPrivated *_g1h_nlprivd;

/* ///////////////////////////////////////////////////////////////////////// */
void g1h_ReflectVectorsd ( int n, const vector3d *v, vector3d *w )
{
  vector3d r;
  double   a;
  int      i;

  r = _g1h_nlprivd->reflv;
  for ( i = 0; i < n; i++ ) {
    a = DotProduct3d ( &v[i], &r );
    AddVector3Md ( &v[i], &r, -2.0*a, &w[i] );
  }
} /*g1h_ReflectVectorsd*/

void g1h_nonlinoutpatchd ( int n, int m, const double *cp, void *usrptr )
{
#define N (n+1)*(m+1)
  memcpy ( &_g1h_nlprivd->nldi[_g1h_nlprivd->auxc*N], cp, N*sizeof(vector3d) );
  _g1h_nlprivd->auxc ++;
#undef N
} /*g1h_nonlinoutpatchd*/

boolean _g1h_StopItd ( int itn, double gn0, double gn,
                       double cn, double dcn, double scf )
{
#define MAXITER 20
#define EPS0    5.0e-10
#define EPS1    1.0e-8

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
} /*_g1h_StopItd*/

/* ///////////////////////////////////////////////////////////////////////// */
boolean g1h_GetHoleSurrndPatchd ( GHoleDomaind *domain,
                                  const point3d *hole_cp,
                                  int i, int j, point3d *bcp )
{
  void    *sp;
  int     *ind;
  double  *ukn, *vkn;
  point3d *q;
  int     hole_k, k;

  sp  = pkv_GetScratchMemTop ();
  ind = pkv_GetScratchMem ( 16*sizeof(int) );
  q   = pkv_GetScratchMem ( 16*sizeof(point3d) );
  if ( !ind || !q )
    goto failure;

  hole_k = domain->hole_k;
  gh_GetBspInd ( hole_k, i, j, ind );
  for ( k = 0; k < 16; k++ )
    q[k] = hole_cp[ind[k]];
  ukn = &domain->hole_knots[11*((i+hole_k-1) % hole_k)+3 ];
  vkn = &domain->hole_knots[11*i+j];
  if ( !mbs_BSPatchToBezd ( 3, 3, 7, ukn, 3, 7, vkn, 12, (double*)q,
                            NULL, NULL, NULL, NULL, NULL, NULL, 12, (double*)bcp ) )
    goto failure;

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*g1h_GetHoleSurrndPatchd*/

boolean g1h_ComputeNLNormald ( GHoleDomaind *domain,
                               const point3d *hole_cp,
                               vector3d *anv )
{
#define DENSITY 16
  void     *sp;
  int      hole_k;
  point3d  *bcp, p;
  int      i, j, k;
  vector3d pu, pv, nu, nv;
  double   t;

  sp = pkv_GetScratchMemTop ();
  hole_k = domain->hole_k;
  bcp = pkv_GetScratchMem ( 16*sizeof(point3d) );
  if ( !bcp )
    goto failure;

  SetVector3d ( &nv, 0.0, 0.0, 0.0 );
  for ( i = 0; i < hole_k; i++ )
    for ( j = 1; j < 3; j++ ) {
      g1h_GetHoleSurrndPatchd ( domain, hole_cp, i, j, bcp );
          /* integrate the unit normal vector at the boundary */
      for ( k = 0; k < DENSITY; k++ ) {
        t = (double)(k+k+1)/(double)(2*DENSITY);
        if ( !mbs_BCHornerDerP3d ( 3, 3, bcp, 0.0, t, &p, &pu, &pv ) )
          goto failure;
        NormalizeVector3d ( &pv );
        CrossProduct3d ( &pu, &pv, &nu );
        AddVector3d ( &nv, &nu, &nv );
      }
    }
  NormalizeVector3d ( &nv );

        /* now a verification */
  for ( i = 0; i < hole_k; i++ )
    for ( j = 1; j < 3; j++ ) {
      g1h_GetHoleSurrndPatchd ( domain, hole_cp, i, j, bcp );
          /* integrate the unit normal vector at the boundary */
      for ( k = 0; k < DENSITY; k++ ) {
        t = (double)(k+k+1)/(double)(2*DENSITY);
        if ( !mbs_BCHornerDerP3d ( 3, 3, bcp, 0.0, t, &p, &pu, &pv ) )
          goto failure;
        CrossProduct3d ( &pu, &pv, &nu );
        if ( DotProduct3d ( &nu, &nv ) <= 0.0 )
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
} /*g1h_ComputeNLNormald*/

boolean _g1h_ComputeNLNormald ( GHoleDomaind *domain,
                                G1HNLPrivated *nlprivate,
                                const point3d *hole_cp )
{
  vector3d nv;

  if ( g1h_ComputeNLNormald ( domain, hole_cp, &nv ) ) {
    nlprivate->nlnv = nv;
    if ( nv.z > 0.0 )
      nv.z += 1.0;
    else
      nv.z -= 1.0;
    NormalizeVector3d ( &nv );
    nlprivate->reflv = nv;
    return true;
  }
  else
    return false;
} /*_g1h_ComputeNLNormald*/

boolean _g1h_TabNLDer0d ( int nkn, const double *tkn,
             const double *hfunc, const double *dhfunc, const double *ddhfunc,
             vector2d *diu, vector2d *div,
             vector2d *diuu, vector2d *diuv, vector2d *divv,
             double *fc00, double *fc01, double *fd00, double *fd01,
             double *psiu, double *psiv,
             double *psiuu, double *psiuv, double *psivv )
{
  void   *sp;
  int    i;
  double *hu, *hv, *huu, *huv, *hvv;

  sp = pkv_GetScratchMemTop ();
  hu = pkv_GetScratchMemd ( 5*nkn*nkn );
  if ( !hu )
    goto failure;
  hv = &hu[nkn*nkn];    huu = &hv[nkn*nkn];
  huv = &huu[nkn*nkn];  hvv = &huv[nkn*nkn];

  if ( !mbs_TabBezC1Coons0Der2d ( 1, nkn, tkn, hfunc, dhfunc, ddhfunc,
                     nkn, tkn, hfunc, dhfunc, ddhfunc,
                     G1_CROSS00DEG, fc00, G1_CROSS01DEG, fc01,
                     G1_CROSS00DEG, fd00, G1_CROSS01DEG, fd01,
                     NULL, hu, hv, huu, huv, hvv ) )
    goto failure;

  for ( i = 0; i < nkn*nkn; i++ )
    if ( !pkn_Comp2iDerivatives2d ( diu[i].x, diu[i].y, div[i].x, div[i].y,
              diuu[i].x, diuu[i].y, diuv[i].x, diuv[i].y, divv[i].x, divv[i].y,
              1, &hu[i], &hv[i], &huu[i], &huv[i], &hvv[i],
              &psiu[i], &psiv[i], &psiuu[i], &psiuv[i], &psivv[i] ) )
      goto failure;

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*_g1h_TabNLDer0d*/

boolean _g1h_TabNLDerd ( int nkn, double *tkn,
             const double *hfunc, const double *dhfunc, const double *ddhfunc,
             vector2d *diu, vector2d *div,
             vector2d *diuu, vector2d *diuv, vector2d *divv,
             double *fc00, double *fc01, double *fc10, double *fc11,
             double *fd00, double *fd01, double *fd10, double *fd11,
             double *psiu, double *psiv,
             double *psiuu, double *psiuv, double *psivv )
{
  void   *sp;
  int    i;
  double *hu, *hv, *huu, *huv, *hvv;

  sp = pkv_GetScratchMemTop ();
  hu = pkv_GetScratchMemd ( 5*nkn*nkn );
  if ( !hu )
    goto failure;
  hv = &hu[nkn*nkn];    huu = &hv[nkn*nkn];
  huv = &huu[nkn*nkn];  hvv = &huv[nkn*nkn];

  if ( !mbs_TabBezC1CoonsDer2d ( 1, nkn, tkn, hfunc, dhfunc, ddhfunc,
            nkn, tkn, hfunc, dhfunc, ddhfunc,
            G1_CROSS00DEG, fc00, G1_CROSS01DEG, fc01,
            G1_CROSS10DEG, fc10, G1_CROSS11DEG, fc11,
            G1_CROSS00DEG, fd00, G1_CROSS01DEG, fd01,
            G1_CROSS10DEG, fd10, G1_CROSS11DEG, fd11,
            NULL, hu, hv, huu, huv, hvv ) )
    goto failure;

  for ( i = 0; i < nkn*nkn; i++ )
    if ( !pkn_Comp2iDerivatives2d ( diu[i].x, diu[i].y, div[i].x, div[i].y,
              diuu[i].x, diuu[i].y, diuv[i].x, diuv[i].y, divv[i].x, divv[i].y,
              1, &hu[i], &hv[i], &huu[i], &huv[i], &hvv[i],
              &psiu[i], &psiv[i], &psiuu[i], &psiuv[i], &psivv[i] ) )
      goto failure;

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*_g1h_TabNLDerd*/

boolean _g1h_TabNLBasisFunctionsd ( GHoleDomaind *domain, int nkn,
                                    G1HNLPrivated *nlpr )
{
#define N ((G1H_FINALDEG+1)*(G1H_FINALDEG+1))
  void     *sp;
  GHolePrivateRecd  *privateG;
  G1HolePrivateRecd *privateG1;
  double   *tkn;
  double   *fc00, *fc01, *fc10, *fc11, *fd00, *fd01, *fd10, *fd11,
           *bc00, *bc01, *bc10, *bc11, *bd00, *bd01, *bd10, *bd11;
  double   *hfunc, *dhfunc, *ddhfunc;
  int      i, j, k, l, kN, kNQ2, nkn2, hole_k, f, fN, nfunc_a, nfunc_b;
  point2d  *di, ddi;
  double   bvz;
  unsigned char *bfcpn;

  sp      = pkv_GetScratchMemTop ();
  hole_k  = domain->hole_k;
  privateG  = domain->privateG;
  privateG1 = domain->privateG1;
  nfunc_a = privateG1->nfunc_a;
  nfunc_b = privateG1->nfunc_b;
  di      = pkv_GetScratchMem ( N*sizeof(point2d) );
  nkn2    = nkn*nkn;
  tkn     = pkv_GetScratchMemd ( nkn );
  hfunc   = pkv_GetScratchMemd ( 12*nkn );
  if ( !tkn || !di || !hfunc )
    goto failure;

  dhfunc = &hfunc[4*nkn];
  ddhfunc = &dhfunc[4*nkn];

  _gh_PrepareTabKnotsd ( nkn, privateG1->opt_quad, tkn );
  if ( !mbs_TabCubicHFuncDer2d ( 0.0, 1.0, nkn, tkn, hfunc, dhfunc, ddhfunc ) )
    goto failure;

  for ( k = kN = kNQ2 = 0;
        k < hole_k;
        k++, kN += N, kNQ2 += nkn2 ) {
    pkv_Selectd ( N, 2, 3, 2, &nlpr->nldi[kN], di );
    for ( i = l = 0;  i < nkn;  i++ )
         /* ***** this is not an optimal algorithm, to be improved */
      for ( j = 0;  j < nkn;  j++, l++ ) {
        if ( !mbs_BCHornerDer2Pd ( G1H_FINALDEG, G1H_FINALDEG, 2, (double*)di,
                    tkn[i], tkn[j], (double*)&ddi,
                    (double*)&nlpr->diu[kNQ2+l], (double*)&nlpr->div[kNQ2+l],
                    (double*)&nlpr->diuu[kNQ2+l], (double*)&nlpr->diuv[kNQ2+l],
                    (double*)&nlpr->divv[kNQ2+l] ) )
          goto failure;
        nlpr->jac[kNQ2+l] = (double)det2d ( &nlpr->diu[kNQ2+l], &nlpr->div[kNQ2+l] );
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
      _g1h_GetBFAPatchCurvesd ( domain, f, k, &fc00, &fc01, &fd00, &fd01 );
      _g1h_TabNLDer0d ( nkn, tkn,
                        hfunc, dhfunc, ddhfunc,
                        &nlpr->diu[kNQ2], &nlpr->div[kNQ2],
                        &nlpr->diuu[kNQ2], &nlpr->diuv[kNQ2], &nlpr->divv[kNQ2],
                        fc00, fc01, fd00, fd01,
                        &nlpr->psiu[fN], &nlpr->psiv[fN],
                        &nlpr->psiuu[fN], &nlpr->psiuv[fN], &nlpr->psivv[fN] );
    }
  }

  bc00 = pkv_GetScratchMemd ( 2*(G1_CROSSDEGSUM+4)*hole_k );
  if ( !bc00 )
    goto failure;
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
    _g1h_TabNLDerd ( nkn, tkn, hfunc, dhfunc, ddhfunc,
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
} /*_g1h_TabNLBasisFunctionsd*/

void _g1h_IntFunc1d ( double pu, double pv,
                      double puu, double puv, double pvv, double jac,
                      double *c, double *cs, double *a, double *b,
                      double *funct )
{
  double A, B, C, CS;

  *c = C = (double)(1.0+pu*pu+pv*pv);
  *cs = CS = (double)sqrt(C);
  *a = A = (double)(0.5*(puu*(1.0+pv*pv)+pvv*(1.0+pu*pu)-2.0*puv*pu*pv));
  *b = B = (double)(1.0/(C*C*CS));
  *funct += A*A*B*jac;
} /*_g1h_IntFunc1d*/

void _g1h_IntFunc2d ( double pu, double pv,
                      double puu, double puv, double pvv,
                      double psiu, double psiv,
                      double psiuu, double psiuv, double psivv,
                      double jac, double c, double A, double B,
                      double *ai, double *bi, double *grad )
{
  double Ai, Bi;

  *ai = Ai = (double)((pu*pvv-puv*pv)*psiu - (pu*puv-puu*pv)*psiv
           +0.5*((1.0+pv*pv)*psiuu+(1.0+pu*pu)*psivv) - pu*pv*psiuv);
  *bi = Bi = (double)(-5.0*(pu*psiu+pv*psiv)*B/c);
  *grad += (double)(A*(2.0*Ai*B+A*Bi)*jac);
} /*_g1h_IntFunc2d*/

void _g1h_IntFunc3d ( double pu, double pv,
                      double puu, double puv, double pvv,
                      double psiu, double psiv,
                      double psiuu, double psiuv, double psivv,
                      double psju, double psjv,
                      double psjuu, double psjuv, double psjvv,
                      double jac, double c,
                      double A, double B, double Ai, double Bi,
                      double Aj, double Bj,
                      double *hessian )
{
  double Aij, Bij;

  Aij = pu*(psiu*psjvv+psivv*psju-psiv*psjuv-psiuv*psjv)+
        pv*(psiv*psjuu+psiuu*psjv-psiu*psjuv-psiuv*psju)+
        puu*psiv*psjv-puv*(psiu*psjv+psiv*psju)+pvv*psiu*psju;
/*
  Bij = (double)(-5.0*((1.0-7.0*pu*pv)*(psiu*psju+psiv*psjv)+
              (pv*pv-6.0*pu*pu)*psiu*psju+
              (pu*pu-6.0*pv*pv)*psiv*psjv)*B/(c*c));
*/
  Bij = (double)(5.0*((6.0*pu*pu-pv*pv-1.0)*psiu*psju
                     +7.0*pu*pv*(psiu*psjv+psju*psiv)
                     +(6.0*pv*pv-pu*pu-1.0)*psiv*psjv)*B/(c*c));
  *hessian += (double)((2.0*A*(Aij*B+Ai*Bj+Aj*Bi)+2.0*Ai*Aj*B+A*A*Bij)*jac);
} /*_g1h_IntFunc3d*/

