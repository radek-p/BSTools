
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

#include "eg1holef.h"
#include "eg1hprivatef.h"
#include "eg1herror.h"


static G1HNLPrivatef *_g1h_InitExtNLprf ( int hole_k, int nfunc_a, int nconstr )
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
    nlpr->diu = pkv_GetScratchMem ( 5*G1_NQUADSQ*hole_k*sizeof(vector2f) );
    nlpr->jac = pkv_GetScratchMemf ( G1_NQUADSQ*hole_k );
    nlpr->psiu = pkv_GetScratchMemf ( 5*((nfunc_a+1)*hole_k+nfunc_c)*G1_NQUADSQ );

    if ( !nlpr->nldi || !nlpr->acoeff || !nlpr->rhole_cp ||
         !nlpr->diu || !nlpr->jac || !nlpr->psiu )
      return NULL;

    nlpr->div = &nlpr->diu[G1_NQUADSQ*hole_k];
    nlpr->diuu = &nlpr->div[G1_NQUADSQ*hole_k];
    nlpr->diuv = &nlpr->diuu[G1_NQUADSQ*hole_k];
    nlpr->divv = &nlpr->diuv[G1_NQUADSQ*hole_k];
    nlpr->psiv = &nlpr->psiu[((nfunc_a+1)*hole_k+nfunc_c)*G1_NQUADSQ];
    nlpr->psiuu = &nlpr->psiv[((nfunc_a+1)*hole_k+nfunc_c)*G1_NQUADSQ];
    nlpr->psiuv = &nlpr->psiuu[((nfunc_a+1)*hole_k+nfunc_c)*G1_NQUADSQ];
    nlpr->psivv = &nlpr->psiuv[((nfunc_a+1)*hole_k+nfunc_c)*G1_NQUADSQ];
  }
  return nlpr;
#undef N
} /*_g1h_InitExtNLprf*/

static boolean _g1h_TabExtNLBasisFunctionsf ( GHoleDomainf *domain,
                                              G1HNLPrivatef *nlpr )
{
#define N ((G1H_FINALDEG+1)*(G1H_FINALDEG+1))
  void     *sp;
  G1HolePrivateRecf *privateG1;
  int      hole_k, nfunc_a, i, j, k, l, f, fN, bN;
  float    *tkn, *tbez, *tbezu, *tbezv, *tbezuu, *tbezuv, *tbezvv,
           *psiu, *psiv, *psiuu, *psiuv, *psivv;
  vector2f *diu, *div, *diuu, *diuv, *divv;

  sp      = pkv_GetScratchMemTop ();
  if ( !_g1h_TabNLBasisFunctionsf ( domain, G1_NQUAD, nlpr ) )
    goto failure;

  hole_k  = domain->hole_k;
  privateG1 = domain->privateG1;
  nfunc_a = privateG1->nfunc_a;
  tkn     = pkv_GetScratchMemf ( G1_NQUAD );
  tbez    = pkv_GetScratchMemf ( 6*G1_NQUADSQ*G1_DBDIM );
  if ( !tkn || !tbez )
    goto failure;
  tbezu = &tbez[G1_NQUADSQ*G1_DBDIM];     tbezv = &tbezu[G1_NQUADSQ*G1_DBDIM];
  tbezuu = &tbezv[G1_NQUADSQ*G1_DBDIM];   tbezuv = &tbezuu[G1_NQUADSQ*G1_DBDIM];
  tbezvv = &tbezuv[G1_NQUADSQ*G1_DBDIM];
  _gh_PrepareTabKnotsf ( G1_NQUAD, privateG1->opt_quad, tkn );

  _g1h_TabTensBezPolyDer2f ( G1_NQUAD, tkn, tbez, tbezu, tbezv,
                             tbezuu, tbezuv, tbezvv );
  for ( k = f = 0; k < hole_k; k++ ) {
    diu = &nlpr->diu[k*G1_NQUADSQ];    div = &nlpr->div[k*G1_NQUADSQ];
    diuu = &nlpr->diuu[k*G1_NQUADSQ];  diuv = &nlpr->diuv[k*G1_NQUADSQ];
    divv = &nlpr->divv[k*G1_NQUADSQ];
    for ( i = bN = 0;  i < G1H_FINALDEG-3;  i++ )
      for ( j = 0;  j < G1H_FINALDEG-3;  j++, f++, bN += G1_NQUADSQ ) {
        fN = ((nfunc_a+1)*hole_k+f)*G1_NQUADSQ;
        psiu = &nlpr->psiu[fN];    psiv = &nlpr->psiv[fN];
        psiuu = &nlpr->psiuu[fN];  psiuv = &nlpr->psiuv[fN];
        psivv = &nlpr->psivv[fN];
        for ( l = 0; l < G1_NQUADSQ; l++ )
          if ( !pkn_Comp2iDerivatives2f ( diu[l].x, diu[l].y, div[l].x, div[l].y,
                     diuu[l].x, diuu[l].y, diuv[l].x, diuv[l].y,
                     divv[l].x, divv[l].y, 1, &tbezu[bN+l], &tbezv[bN+l],
                     &tbezuu[bN+l], &tbezuv[bN+l], &tbezvv[bN+l],
                     &psiu[l], &psiv[l], &psiuu[l], &psiuv[l], &psivv[l] ) )
            goto failure;
      }
  }
  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
#undef N
} /*_g1h_TabExtNLBasisFunctionsf*/

static float _g1h_ComputeExtNLFuncf ( GHoleDomainf *domain,
                                      G1HNLPrivatef *nlprivate,
                                      const float *coeff )
{
#define SCALE (1.0/(float)G1_NQUADSQ)
  G1HolePrivateRecf *privateG1;
  int   i, k, kn, knot, fi, hole_k, nfunc_a, nfunc_c;
  float pu, pv, puu, puv, pvv, A, B, c, cs, funct;
  float *psiu, *psiv, *psiuu, *psiuv, *psivv, *jac;

  privateG1 = domain->privateG1;
  hole_k = domain->hole_k;
  nfunc_a = privateG1->nfunc_a;
  nfunc_c = hole_k*G1_DBDIM;
  psiu = nlprivate->psiu;
  psiv = nlprivate->psiv;
  psiuu = nlprivate->psiuu;
  psiuv = nlprivate->psiuv;
  psivv = nlprivate->psivv;
  jac = nlprivate->jac;

  funct = 0.0;
  for ( k = knot = 0;  k < hole_k;  k++ )
    for ( kn = 0;  kn < G1_NQUADSQ;  kn++, knot++ ) {
        /* evaluate the function and its derivatives */
      fi = nfunc_a*G1_NQUADSQ*hole_k+knot;
      pu = psiu[fi];    pv = psiv[fi];
      puu = psiuu[fi];  puv = psiuv[fi];  pvv = psivv[fi];
      for ( i = 0; i < nfunc_a; i++ ) {
        fi = i*G1_NQUADSQ*hole_k+knot;
        c = coeff[nfunc_c+i];
        pu -= c*psiu[fi];    pv -= c*psiv[fi];
        puu -= c*psiuu[fi];  puv -= c*psiuv[fi];  pvv -= c*psivv[fi];
      }
      for ( i = 0; i < G1_DBDIM; i++ ) {
        fi = ((nfunc_a+1)*hole_k+k*G1_DBDIM+i)*G1_NQUADSQ+kn;
        c = coeff[k*G1_DBDIM+i];
        pu -= c*psiu[fi];    pv -= c*psiv[fi];
        puu -= c*psiuu[fi];  puv -= c*psiuv[fi];  pvv -= c*psivv[fi];
      }

        /* make the integration */
    _g1h_IntFunc1f ( pu, pv, puu, puv, pvv, jac[knot],
                     &c, &cs, &A, &B, &funct );
/*
if ( k == 20 )
  Verify ( 3.0/(2*G1_NQUAD), 9.0/(2*G1_NQUAD), p, pu, pv, puu, puv, pvv, A/(c*cs) );
*/
    }
  return (float)(funct*SCALE);
#undef SCALE
} /*_g1h_ComputeExtNLFuncf*/

static boolean _g1h_ComputeExtNLFuncGradf ( GHoleDomainf *domain,
                                            G1HNLPrivatef *nlprivate,
                                            const float *coeff,
                                            float *func, float *grad )
{
#define SCALE (1.0/(float)G1_NQUADSQ)
  void  *sp;
  G1HolePrivateRecf *privateG1;
  int   i, k, kn, knot, fi, hole_k, nfunc_a, nfunc_c, nfunc;
  float pu, pv, puu, puv, pvv, A, B, c, cs, funct;
  float *psiu, *psiv, *psiuu, *psiuv, *psivv, *jac, *Ai, *Bi;

  sp = pkv_GetScratchMemTop ();
  privateG1 = domain->privateG1;
  hole_k = domain->hole_k;
  nfunc_a = privateG1->nfunc_a;
  nfunc_c = hole_k*G1_DBDIM;
  nfunc = nfunc_a+nfunc_c;
  psiu = nlprivate->psiu;
  psiv = nlprivate->psiv;
  psiuu = nlprivate->psiuu;
  psiuv = nlprivate->psiuv;
  psivv = nlprivate->psivv;
  jac = nlprivate->jac;
  Ai = pkv_GetScratchMemf ( nfunc );
  Bi = pkv_GetScratchMemf ( nfunc );
  if ( !Ai || !Bi )
    goto failure;

  funct = 0.0;
  memset ( grad, 0, nfunc*sizeof(float) );
  for ( k = knot = 0;  k < hole_k;  k++ )
    for ( kn = 0;  kn < G1_NQUADSQ;  kn++, knot++ ) {
        /* evaluate the function and its derivatives */
      fi = nfunc_a*G1_NQUADSQ*hole_k+knot;
      pu = psiu[fi];    pv = psiv[fi];
      puu = psiuu[fi];  puv = psiuv[fi];  pvv = psivv[fi];
      for ( i = 0; i < nfunc_a; i++ ) {
        fi = i*G1_NQUADSQ*hole_k+knot;
        c = coeff[nfunc_c+i];
        pu -= c*psiu[fi];    pv -= c*psiv[fi];
        puu -= c*psiuu[fi];  puv -= c*psiuv[fi];  pvv -= c*psivv[fi];
      }
      for ( i = 0; i < G1_DBDIM; i++ ) {
        fi = ((nfunc_a+1)*hole_k+k*G1_DBDIM+i)*G1_NQUADSQ+kn;
        c = coeff[k*G1_DBDIM+i];
        pu -= c*psiu[fi];    pv -= c*psiv[fi];
        puu -= c*psiuu[fi];  puv -= c*psiuv[fi];  pvv -= c*psivv[fi];
      }

        /* make the integration */
      _g1h_IntFunc1f ( pu, pv, puu, puv, pvv, jac[knot],
                       &c, &cs, &A, &B, &funct );
/*
if ( k == 20 )
  Verify ( 3.0/(2*G1_NQUAD), 9.0/(2*G1_NQUAD), p, pu, pv, puu, puv, pvv, A/(c*cs) );
*/
        /* integrate the gradient */
      for ( i = 0; i < nfunc_a; i++ ) {
        fi = i*G1_NQUADSQ*hole_k+knot;
        _g1h_IntFunc2f ( pu, pv, puu, puv, pvv,
                         psiu[fi], psiv[fi], psiuu[fi], psiuv[fi], psivv[fi],
                         jac[knot], c, A, B,
                         &Ai[nfunc_c+i], &Bi[nfunc_c+i], &grad[nfunc_c+i] );
      }
      for ( i = 0; i < G1_DBDIM; i++ ) {
        fi = ((nfunc_a+1)*hole_k+k*G1_DBDIM+i)*G1_NQUADSQ+kn;
        _g1h_IntFunc2f ( pu, pv, puu, puv, pvv,
                         psiu[fi], psiv[fi], psiuu[fi], psiuv[fi], psivv[fi],
                         jac[knot], c, A, B,
                         &Ai[k*G1_DBDIM+i], &Bi[k*G1_DBDIM+i], &grad[k*G1_DBDIM+i] );
      }
    }

  *func = (float)(funct*SCALE);
  pkn_MultMatrixNumf ( 1, nfunc, 0, grad, SCALE, 0, grad );
  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
#undef SCALE
} /*_g1h_ComputeExtNLFuncGradf*/

static boolean _g1h_ComputeExtNLFuncGradHessianf ( GHoleDomainf *domain,
                          G1HNLPrivatef *nlprivate,
                          const float *coeff,
                          float *func, float *grad, float *hii )
{
#define SCALE (1.0/(float)G1_NQUADSQ)
  void  *sp;
  G1HolePrivateRecf *privateG1;
  int   i, j, k, kn, knot, fi, fj, hole_k, nfunc_a, nfunc_c, nfunc;
  int   diagblsize, sideblsize, akksize;
  float pu, pv, puu, puv, pvv, A, B, c, cs, funct;
  float *psiu, *psiv, *psiuu, *psiuv, *psivv, *jac, *Ai, *Bi;
  float *hki, *hkk;

  sp = pkv_GetScratchMemTop ();
  privateG1 = domain->privateG1;
  hole_k = domain->hole_k;
  nfunc_a = privateG1->nfunc_a;
  nfunc_c = hole_k*G1_DBDIM;
  nfunc = nfunc_a+nfunc_c;
  hki = &hii[pkn_Block1FindBlockPos ( hole_k, G1_DBDIM, nfunc_a, hole_k, 0)];
  hkk = &hii[pkn_Block1FindBlockPos ( hole_k, G1_DBDIM, nfunc_a, hole_k, hole_k)];
  psiu = nlprivate->psiu;
  psiv = nlprivate->psiv;
  psiuu = nlprivate->psiuu;
  psiuv = nlprivate->psiuv;
  psivv = nlprivate->psivv;
  jac = nlprivate->jac;
  Ai = pkv_GetScratchMemf ( nfunc );
  Bi = pkv_GetScratchMemf ( nfunc );
  if ( !Ai || !Bi )
    goto failure;
  diagblsize = G1_DBDIM*(G1_DBDIM+1)/2;
  sideblsize = G1_DBDIM*nfunc_a;
  akksize = nfunc_a*(nfunc_a+1)/2;

  funct = 0.0;
  memset ( grad, 0, nfunc*sizeof(float) );
  memset ( hii, 0, hole_k*diagblsize*sizeof(float) );
  memset ( hki, 0, hole_k*sideblsize*sizeof(float) );
  memset ( hkk, 0, akksize*sizeof(float) );
  for ( k = knot = 0;  k < hole_k;  k++ )
    for ( kn = 0;  kn < G1_NQUADSQ;  kn++, knot++ ) {
        /* evaluate the function and its derivatives */
      fi = nfunc_a*G1_NQUADSQ*hole_k+knot;
      pu = psiu[fi];    pv = psiv[fi];
      puu = psiuu[fi];  puv = psiuv[fi];  pvv = psivv[fi];
      for ( i = 0; i < nfunc_a; i++ ) {
        fi = i*G1_NQUADSQ*hole_k+knot;
        c = coeff[nfunc_c+i];
        pu -= c*psiu[fi];    pv -= c*psiv[fi];
        puu -= c*psiuu[fi];  puv -= c*psiuv[fi];  pvv -= c*psivv[fi];
      }
      for ( i = 0; i < G1_DBDIM; i++ ) {
        fi = ((nfunc_a+1)*hole_k+k*G1_DBDIM+i)*G1_NQUADSQ+kn;
        c = coeff[k*G1_DBDIM+i];
        pu -= c*psiu[fi];    pv -= c*psiv[fi];
        puu -= c*psiuu[fi];  puv -= c*psiuv[fi];  pvv -= c*psivv[fi];
      }

        /* make the integration */
      _g1h_IntFunc1f ( pu, pv, puu, puv, pvv, jac[knot],
                       &c, &cs, &A, &B, &funct );
/*
if ( k == 20 )
  Verify ( 3.0/(2*G1_NQUAD), 9.0/(2*G1_NQUAD), p, pu, pv, puu, puv, pvv, A/(c*cs) );
*/
        /* integrate the gradient */
      for ( i = 0; i < nfunc_a; i++ ) {
        fi = i*G1_NQUADSQ*hole_k+knot;
        _g1h_IntFunc2f ( pu, pv, puu, puv, pvv,
                         psiu[fi], psiv[fi], psiuu[fi], psiuv[fi], psivv[fi],
                         jac[knot], c, A, B,
                         &Ai[nfunc_c+i], &Bi[nfunc_c+i], &grad[nfunc_c+i] );
      }
      for ( i = 0; i < G1_DBDIM; i++ ) {
        fi = ((nfunc_a+1)*hole_k+k*G1_DBDIM+i)*G1_NQUADSQ+kn;
        _g1h_IntFunc2f ( pu, pv, puu, puv, pvv,
                         psiu[fi], psiv[fi], psiuu[fi], psiuv[fi], psivv[fi],
                         jac[knot], c, A, B,
                         &Ai[k*G1_DBDIM+i], &Bi[k*G1_DBDIM+i], &grad[k*G1_DBDIM+i] );
      }

        /* integrate the Hessian */
      for ( i = 0; i < nfunc_a; i++ ) {
        fi = i*G1_NQUADSQ*hole_k+knot;
        for ( j = 0; j <= i; j++ ) {
          fj = j*G1_NQUADSQ*hole_k+knot;
          _g1h_IntFunc3f ( pu, pv, puu, puv, pvv,
                           psiu[fi], psiv[fi], psiuu[fi], psiuv[fi], psivv[fi],
                           psiu[fj], psiv[fj], psiuu[fj], psiuv[fj], psivv[fj],
                           jac[knot], c, A, B,
                           Ai[nfunc_c+i], Bi[nfunc_c+i],
                           Ai[nfunc_c+j], Bi[nfunc_c+j],
                           &hkk[pkn_SymMatIndex(i,j)] );
        }
        for ( j = 0; j < G1_DBDIM; j++ ) {
          fj = ((nfunc_a+1)*hole_k+k*G1_DBDIM+j)*G1_NQUADSQ+kn;
          _g1h_IntFunc3f ( pu, pv, puu, puv, pvv,
                           psiu[fi], psiv[fi], psiuu[fi], psiuv[fi], psivv[fi],
                           psiu[fj], psiv[fj], psiuu[fj], psiuv[fj], psivv[fj],
                           jac[knot], c, A, B,
                           Ai[nfunc_c+i], Bi[nfunc_c+i],
                           Ai[k*G1_DBDIM+j], Bi[k*G1_DBDIM+j],
                           &hki[(k*nfunc_a+i)*G1_DBDIM+j] );
        }
      }
      for ( i = 0; i < G1_DBDIM; i++ ) {
        fi = ((nfunc_a+1)*hole_k+k*G1_DBDIM+i)*G1_NQUADSQ+kn;
        for ( j = 0; j <= i; j++ ) {
          fj = ((nfunc_a+1)*hole_k+k*G1_DBDIM+j)*G1_NQUADSQ+kn;
          _g1h_IntFunc3f ( pu, pv, puu, puv, pvv,
                           psiu[fi], psiv[fi], psiuu[fi], psiuv[fi], psivv[fi],
                           psiu[fj], psiv[fj], psiuu[fj], psiuv[fj], psivv[fj],
                           jac[knot], c, A, B,
                           Ai[k*G1_DBDIM+i], Bi[k*G1_DBDIM+i],
                           Ai[k*G1_DBDIM+j], Bi[k*G1_DBDIM+j],
                           &hii[k*diagblsize+pkn_SymMatIndex(i,j)] );
        }
      }
    }

  *func = (float)(funct*SCALE);
  pkn_MultMatrixNumf ( 1, nfunc, 0, grad, SCALE, 0, grad );
  pkn_MultMatrixNumf ( 1, akksize, 0, hkk, SCALE, 0, hkk );
  pkn_MultMatrixNumf ( 1, hole_k*sideblsize, 0, hki, SCALE, 0, hki );
  pkn_MultMatrixNumf ( 1, hole_k*diagblsize, 0, hii, SCALE, 0, hii );
  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
#undef SCALE
} /*_g1h_ComputeExtNLFuncGradHessianf*/

static boolean g1h_ExtNLNewtonf ( GHoleDomainf *domain, G1HNLPrivatef *nlprivate )
{
#define EPSF 2.0e-4
  void    *sp;
  G1HolePrivateRecf *privateG1;
  int     itn, jtn, ktn, hole_k, nfunc_a, nfunc_c, nfunc, asize;
  int     diagblsize, sideblsize, akksize;
  float   *coeff, *dcoeff, *grad, *hii, *chii;
  float   func, func0, func1, aux, gn, gn0 = 0.0, gn1, dco, dyn;
  boolean positive;

  sp = pkv_GetScratchMemTop ();
  hole_k  = domain->hole_k;
  privateG1 = domain->privateG1;
  nfunc_a = privateG1->nfunc_a;
  nfunc_c = hole_k*G1_DBDIM;
  nfunc   = nfunc_a+nfunc_c;
  diagblsize = G1_DBDIM*(G1_DBDIM+1)/2;
  sideblsize = G1_DBDIM*nfunc_a;
  akksize = nfunc_a*(nfunc_a+1)/2;

  coeff  = pkv_GetScratchMemf ( 3*nfunc );
  asize  = pkn_Block1ArraySize ( hole_k, G1_DBDIM, nfunc_a );
  hii    = pkv_GetScratchMemf ( 2*asize );
  if ( !coeff || !hii )
    goto failure;
  dcoeff = &coeff[nfunc];
  grad   = &dcoeff[nfunc];
  chii   = &hii[asize];

        /* setup the initial point */
  pkv_Selectf ( nfunc, 1, 3, 1, &nlprivate->acoeff[0].z, coeff );

        /* Newton iterations */
  for ( itn = ktn = 0; ; itn++ ) {
    if ( !_g1h_ComputeExtNLFuncGradHessianf ( domain, nlprivate, coeff,
                             &func, grad, hii ) )
      goto failure;
    gn = (float)sqrt ( pkn_ScalarProductf ( nfunc, grad, grad ) );
    if ( itn == 0 ) {
/*
printf ( "func = %f, gn0 = %f\n", func, gn );
*/
      gn0 = gn;
    }
    memcpy ( chii, hii,
            (hole_k*(diagblsize+sideblsize)+akksize)*sizeof(float) );

    if ( (positive = pkn_Block1CholeskyDecompMf ( hole_k, G1_DBDIM, nfunc_a,
                                     hii )) ) {
      pkn_Block1LowerTrMSolvef ( hole_k, G1_DBDIM, nfunc_a, hii, 1, 1, grad );
      pkn_Block1UpperTrMSolvef ( hole_k, G1_DBDIM, nfunc_a, hii, 1, 1, grad );
    }
    else {

printf ( "! " );

      if ( !pkn_Block1SymMatrixMultf ( hole_k, G1_DBDIM, nfunc_a, chii,
                                       1, 1, grad, 1, dcoeff ) )
        goto failure;
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
      func1 = _g1h_ComputeExtNLFuncf ( domain, nlprivate, dcoeff );
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
/*
printf ( "+" );
*/
        _g1h_ComputeExtNLFuncGradf ( domain, nlprivate, coeff, &func0, grad );
        gn1 = (float)sqrt ( pkn_ScalarProductf ( nfunc, grad, grad ) );
        pkn_Block1LowerTrMSolvef ( hole_k, G1_DBDIM, nfunc_a, hii, 1, 1, grad );
        pkn_Block1UpperTrMSolvef ( hole_k, G1_DBDIM, nfunc_a, hii, 1, 1, grad );
        pkn_AddMatrixf ( 1, nfunc, 0, coeff, 0, grad, 0, dcoeff );
        func1 = _g1h_ComputeExtNLFuncf ( domain, nlprivate, dcoeff );
        if ( func1 >= func0 || gn1 >= 0.5*gn )
          break;
/*
printf ( "    func = %f, gn = %f\n", func1, gn );
*/
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
} /*g1h_ExtNLNewtonf*/

/* ///////////////////////////////////////////////////////////////////////// */
boolean g1h_NLExtFillHolef ( GHoleDomainf *domain, const point3f *hole_cp,  
                             float *acoeff, void *usrptr,
                             void (*outpatch) ( int n, int m, const point3f *cp,
                                                void *usrptr ) )
{
  typedef void outscf ( int n, int m, const float *cp, void *usrptr );

  void     *sp;
  G1HNLPrivatef     *nlprivate;
  G1HolePrivateRecf *privateG1;
  float    *fc00, *Bi, *Bk;
  int      hole_k, nfunc_a, nfunc_b, nfunc_c;

  sp = pkv_GetScratchMemTop ();
  privateG1 = domain->privateG1;
  hole_k = domain->hole_k;
  nfunc_a = privateG1->nfunc_a;
  nfunc_b = privateG1->nfunc_b;
  nfunc_c = hole_k*G1_DBDIM;
  _g1h_nlprivf = nlprivate = _g1h_InitExtNLprf ( hole_k, nfunc_a, 0 );
  if ( !nlprivate )
    goto failure;

  if ( !_g1h_ComputeNLNormalf ( domain, nlprivate, hole_cp ) )
    goto failure;
  g1h_ReflectVectorsf ( 12*hole_k+1, hole_cp, nlprivate->rhole_cp );

  nlprivate->auxc = 0;
  if ( !g1h_ExtFillHolef ( domain, 3, (float*)nlprivate->rhole_cp,
                           (float*)nlprivate->acoeff, NULL, g1h_nonlinoutpatchf ) )
    goto failure;

  if ( !_g1h_TabExtNLBasisFunctionsf ( domain, nlprivate ) )
    goto failure;

  if ( !g1h_ExtNLNewtonf ( domain, nlprivate ) )
    goto failure;

  g1h_ReflectVectorsf ( nfunc_a+nfunc_c, nlprivate->acoeff,
                        nlprivate->acoeff );
  if ( acoeff )
    memcpy ( acoeff, nlprivate->acoeff,
            (nfunc_a+nfunc_c)*sizeof(vector3f) );

  fc00 = pkv_GetScratchMem ( (G1_CROSSDEGSUM+4)*2*hole_k*sizeof(vector3f) );
  if ( !fc00 )
    goto failure;
  Bi = privateG1->EBmat;
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
} /*g1h_NLExtFillHolef*/

/* ///////////////////////////////////////////////////////////////////////// */
static boolean g1h_ExtNLConstrNewtonf ( GHoleDomainf *domain,
                                        G1HNLPrivatef *nlprivate, int nconstr,
                                        float *ECmat )
{
#define EPSF 2.0e-4
  void    *sp;
  G1HolePrivateRecf *privateG1;
  int     itn, jtn, ktn, i, j, hole_k, nfunc_a, nfunc_c, nfunc_ac, nfunc;
  int     diagblsize, sideblsize, esideblsize, asize, esize;
  float   *coeff, *grad, *hii, *hkk, *hki;
  float   *cT, *E22ii, *E22kk, *E22ki, *cE22ii;
  float   *aa, *D1, *y, *y1, *M, *f;
  float   func, func0, func1, aux, gn, gn0 = 0.0, gn1, dco, dyn, dyn1;
  boolean positive;

  sp = pkv_GetScratchMemTop ();
  hole_k   = domain->hole_k;
  privateG1  = domain->privateG1;
  nfunc_a  = privateG1->nfunc_a;
  nfunc_c  = hole_k*G1_DBDIM;
  nfunc_ac = nfunc_a+nfunc_c;
  nfunc    = nfunc_ac-nconstr;
  diagblsize  = G1_DBDIM*(G1_DBDIM+1)/2;
  sideblsize  = G1_DBDIM*nfunc_a;
  esideblsize = G1_DBDIM*(nfunc_a-nconstr);

  coeff    = pkv_GetScratchMemf ( nfunc_ac );
  grad     = pkv_GetScratchMemf ( nfunc_ac );
  asize    = pkn_Block1ArraySize ( hole_k, G1_DBDIM, nfunc_a );
  hii      = pkv_GetScratchMemf ( asize );
  cT       = pkv_GetScratchMemf ( nfunc_a*nconstr );
  aa       = pkv_GetScratchMemf ( 2*nconstr );
  D1       = pkv_GetScratchMemf ( (nconstr*(nconstr+1))/2);
  esize    = pkn_Block1ArraySize ( hole_k, G1_DBDIM, nfunc_a-nconstr );
  E22ii    = pkv_GetScratchMemf ( esize );
  cE22ii   = pkv_GetScratchMemf ( esize );
  y        = pkv_GetScratchMemf ( nfunc_ac );
  y1       = pkv_GetScratchMemf ( nfunc_ac );
  M        = pkv_GetScratchMemf ( (nfunc_a*(nfunc_a+1))/2 );
  f        = pkv_GetScratchMemf ( nfunc_ac );
  if ( !coeff || !grad || !hii || !cT || !aa ||
       !D1 || !E22ii || !cE22ii || !y || !y1 || !M || !f )
    goto failure;
  hkk = &hii[pkn_Block1FindBlockPos ( hole_k, G1_DBDIM, nfunc_a, hole_k, hole_k )];
  hki = &hii[pkn_Block1FindBlockPos ( hole_k, G1_DBDIM, nfunc_a, hole_k, 0 )];
  E22kk = &E22ii[pkn_Block1FindBlockPos ( hole_k, G1_DBDIM, nfunc_a-nconstr,
                                          hole_k, hole_k )];
  E22ki = &E22ii[pkn_Block1FindBlockPos ( hole_k, G1_DBDIM, nfunc_a-nconstr,
                                          hole_k, 0 )];

        /* setup the initial point */
  pkv_Selectf ( nfunc_ac, 1, 3, 1, &nlprivate->acoeff[0].z, coeff );

        /* step 1: decompose the constraint equations matrix */
  pkv_TransposeMatrixf ( nconstr, nfunc_a, nfunc_ac, &ECmat[nfunc_c],
                         nconstr, cT );
  pkn_QRDecomposeMatrixf ( nfunc_a, nconstr, cT, aa );
  for ( i = 0; i < nconstr; i++ )
    for ( j = i; j < nconstr; j++ )
      D1[pkn_SymMatIndex(i,j)] = cT[i*nconstr+j];

        /* Newton iterations */
  for ( itn = ktn = 0; ; itn++ ) {
          /* step 2 */
    memcpy ( y, coeff, nfunc_ac*sizeof(float) );
    pkn_multiReflectVectorf ( nfunc_a, nconstr, cT, aa, 1, 1, &y[nfunc_c] );
    pkn_LowerTrMatrixSolvef ( nconstr, D1, 1,
                              3, &nlprivate->rhole_cp[12*hole_k+1].z,
                              1, &y[nfunc_c] );
    for ( i = nfunc_c; i < nfunc_c+nconstr; i++ )
      y[i] = -y[i];

          /* step 3 */
    if ( !_g1h_ComputeExtNLFuncGradHessianf ( domain, nlprivate, coeff,
                                              &func, grad, hii ) )
      goto failure;

    if ( !pkn_ComputeQTSQf ( nfunc_a, hkk, nconstr, cT, aa, M ) )
      goto failure;
    for ( i = 0; i < nfunc_a-nconstr; i++ )
      for ( j = i; j < nfunc_a-nconstr; j++ )
        E22kk[pkn_SymMatIndex(i,j)] = M[pkn_SymMatIndex(nconstr+i,nconstr+j)];
    for ( i = 0; i < hole_k; i++ )
      pkn_multiReflectVectorf ( nfunc_a, nconstr, cT, aa,
                                G1_DBDIM, G1_DBDIM, &hki[i*sideblsize] );
    pkv_Selectf ( hole_k, esideblsize, sideblsize, esideblsize,
                  &hki[nconstr*G1_DBDIM], E22ki );
    memcpy ( E22ii, hii, hole_k*diagblsize*sizeof(float) );
    memcpy ( cE22ii, E22ii, esize*sizeof(float) );

          /* step 4 */
    memcpy ( f, grad, nfunc_ac*sizeof(float) );
    pkn_multiReflectVectorf ( nfunc_a, nconstr, cT, aa, 1, 1, &f[nfunc_c] );
    memmove ( &f[nfunc_c], &f[nfunc_c+nconstr],
              (nfunc_a-nconstr)*sizeof(float) );
    gn = (float)sqrt ( pkn_ScalarProductf ( nfunc, f, f ) );
    if ( itn == 0 ) {
/*
printf ( "func = %f, gn0 = %f\n", func, gn );
*/
      gn0 = gn;
    }
          /* step 5 */
    if ( (positive = pkn_Block1CholeskyDecompMf ( hole_k, G1_DBDIM, nfunc_a-nconstr,
                                     E22ii )) ) {
      pkn_Block1LowerTrMSolvef ( hole_k, G1_DBDIM, nfunc_a-nconstr, E22ii, 1, 1, f );
      pkn_Block1UpperTrMSolvef ( hole_k, G1_DBDIM, nfunc_a-nconstr, E22ii, 1, 1, f );
    }
    else {

printf ( "! " );

      if ( !pkn_Block1SymMatrixMultf ( hole_k, G1_DBDIM, nfunc_a-nconstr, cE22ii,
                                       1, 1, f, 1, y1 ) )
        goto failure;
      aux = (float)pkn_ScalarProductf ( nfunc, f, y1 );
      if ( aux <= 0.0 || aux < EPSF*gn )
        goto failure;
      pkn_MultMatrixNumf ( 1, nfunc, 1, f, gn/aux, 1, f );
    }

          /* step 6 */
    dyn = (float)sqrt ( pkn_ScalarProductf ( nfunc, f, f ) );
    memmove ( &f[nfunc_c+nconstr], &f[nfunc_c],
              (nfunc_a-nconstr)*sizeof(float) );
    dco = (float)sqrt ( pkn_ScalarProductf ( nfunc_ac, coeff, coeff ) );
    for ( aux = 1.0; aux >= EPSF; aux *= 0.5 ) {
      for ( i = 0; i < nfunc_c; i++ )
        y1[i] = y[i]+aux*f[i];
      memcpy ( &y1[nfunc_c], &y[nfunc_c], nconstr*sizeof(float) );
      for ( i = nfunc_c+nconstr; i < nfunc_ac; i++ )
        y1[i] = y[i]+aux*f[i];
      pkn_multiInvReflectVectorf ( nfunc_a, nconstr, cT, aa, 1, 1, &y1[nfunc_c] );
      func1 = _g1h_ComputeExtNLFuncf ( domain, nlprivate, y1 );
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
        pkn_multiReflectVectorf ( nfunc_a, nconstr, cT, aa, 1, 1, &y[nfunc_c] );
        pkn_LowerTrMatrixSolvef ( nconstr, D1, 1,
                                  3, &nlprivate->rhole_cp[12*hole_k+1].z,
                                  1, &y[nfunc_c] );
        for ( i = nfunc_c; i < nfunc_c+nconstr; i++ )
          y[i] = -y[i];
              /* step 3' */
        _g1h_ComputeExtNLFuncGradf ( domain, nlprivate, coeff, &func0, grad );
              /* step 4' */
        memcpy ( f, grad, nfunc_ac*sizeof(float) );
        pkn_multiReflectVectorf ( nfunc_a, nconstr, cT, aa, 1, 1, &f[nfunc_c] );
        memmove ( &f[nfunc_c], &f[nfunc_c+nconstr],
                  (nfunc_a-nconstr)*sizeof(float) );
        gn1 = (float)sqrt ( pkn_ScalarProductf ( nfunc, f, f ) );
              /* step 5' */
        pkn_Block1LowerTrMSolvef ( hole_k, G1_DBDIM, nfunc_a-nconstr,
                                   E22ii, 1, 1, f );
        pkn_Block1UpperTrMSolvef ( hole_k, G1_DBDIM, nfunc_a-nconstr,
                                   E22ii, 1, 1, f );
              /* step 6' */
        dyn1 = (float)sqrt ( pkn_ScalarProductf ( nfunc, f, f ) );
        memmove ( &f[nfunc_c+nconstr], &f[nfunc_c],
                  (nfunc_a-nconstr)*sizeof(float) );
        for ( i = 0; i < nfunc_c; i++ )
          y1[i] = y[i]+f[i];
        memcpy ( &y1[nfunc_c], &y[nfunc_c], nconstr*sizeof(float) );
        for ( i = nfunc_c+nconstr; i < nfunc_ac; i++ )
          y1[i] = y[i]+f[i];
        pkn_multiInvReflectVectorf ( nfunc_a, nconstr, cT, aa, 1, 1, &y1[nfunc_c] );

        func1 = _g1h_ComputeExtNLFuncf ( domain, nlprivate, y1 );
        if ( func1 >= func0 || gn1 >= 0.5*gn )
          break;
/*
printf ( "    func = %f, gn = %f\n", func1, gn );
*/
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
/*
printf ( "func = %f\n", func1 );
*/
  pkv_Selectf ( nfunc_ac, 1, 1, 3, coeff, &nlprivate->acoeff[0].z );

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
#undef EPSF
} /*g1h_ExtNLConstrNewtonf*/

boolean g1h_NLExtFillHoleConstrf ( GHoleDomainf *domain,
                     const point3f *hole_cp,
                     int nconstr, const vector3f *constr,
                     float *acoeff, void *usrptr,
                     void (*outpatch) ( int n, int m, const point3f *cp,
                                        void *usrptr ) )
{
  typedef void outscf ( int n, int m, const float *cp, void *usrptr );

  void     *sp;
  G1HNLPrivatef     *nlprivate;
  G1HolePrivateRecf *privateG1;
  float    *fc00, *Bi, *Bk;
  int      hole_k, nfunc_a, nfunc_b, nfunc_c, nfunc_ac;

  sp = pkv_GetScratchMemTop ();
  privateG1 = domain->privateG1;
  hole_k  = domain->hole_k;
  nfunc_a = privateG1->nfunc_a;
  nfunc_b = privateG1->nfunc_b;
  nfunc_c = hole_k*G1_DBDIM;
  nfunc_ac = nfunc_a+nfunc_c;
  _g1h_nlprivf = nlprivate = _g1h_InitExtNLprf ( hole_k, nfunc_a, nconstr );
  if ( !nlprivate )
    goto failure;

  if ( !_g1h_ComputeNLNormalf ( domain, nlprivate, hole_cp ) )
    goto failure;
  g1h_ReflectVectorsf ( 12*hole_k+1, hole_cp, nlprivate->rhole_cp );
  g1h_ReflectVectorsf ( nconstr, constr, &nlprivate->rhole_cp[12*hole_k+1] );

  nlprivate->auxc = 0;
  if ( !g1h_ExtFillHoleConstrf ( domain, 3, (float*)nlprivate->rhole_cp,
                    nconstr, (float*)&nlprivate->rhole_cp[12*hole_k+1],
                    (float*)nlprivate->acoeff, NULL, g1h_nonlinoutpatchf ) )
    goto failure;

  if ( !_g1h_TabExtNLBasisFunctionsf ( domain, nlprivate ) )
    goto failure;

  if ( !g1h_ExtNLConstrNewtonf ( domain, nlprivate, nconstr, privateG1->ECmat ) )
    goto failure;

  g1h_ReflectVectorsf ( nfunc_ac, nlprivate->acoeff,
                        nlprivate->acoeff );
  if ( acoeff )
    memcpy ( acoeff, nlprivate->acoeff, nfunc_ac*sizeof(vector3f) );

  fc00 = pkv_GetScratchMem ( (G1_CROSSDEGSUM+4)*2*hole_k*sizeof(vector3f) );
  if ( !fc00 )
    goto failure;
  Bi = privateG1->EBmat;
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
} /*g1h_NLExtFillHoleConstrf*/

/* ///////////////////////////////////////////////////////////////////////// */
boolean _g1h_ReflectExtAltConstrMatrixf ( GHoleDomainf *domain,
                                          vector3f *reflv, float *RACmat )
{
  void  *sp;
  G1HolePrivateRecf *privateG1;
  int   i, j, nconstr, hole_k, nfunc_a, nfunc_c, nfunc_ac;
  float *buf;

  sp = pkv_GetScratchMemTop ();
  hole_k  = domain->hole_k;
  privateG1 = domain->privateG1;
  nfunc_a = privateG1->nfunc_a;
  nfunc_c = hole_k*G1_DBDIM;
  nfunc_ac = nfunc_a+nfunc_c;
  nconstr = privateG1->extnaconstr;
  memcpy ( RACmat, privateG1->AECmat, nconstr*nfunc_ac*3*sizeof(float) );
  buf = pkv_GetScratchMemf ( nfunc_ac );
  if ( !buf ) {
    domain->error_code = G1H_ERROR_NO_SCRATCH_MEMORY;
    goto failure;
  }
  for ( i = 0; i < nconstr; i++ ) {
    pkn_MultMatrixNumf ( 1, nfunc_ac, 0, RACmat, reflv->x, 0, buf );
    pkn_AddMatrixMf ( 1, nfunc_ac, 0, buf, 0, &RACmat[nfunc_ac], reflv->y, 0, buf );
    pkn_AddMatrixMf ( 1, nfunc_ac, 0, buf, 0, &RACmat[2*nfunc_ac], reflv->z, 0, buf );
    for ( j = 0; j < nfunc_ac; j++ ) {
      RACmat[j]            -= (float)(2.0*buf[j]*reflv->x);
      RACmat[nfunc_ac+j]   -= (float)(2.0*buf[j]*reflv->y);
      RACmat[2*nfunc_ac+j] -= (float)(2.0*buf[j]*reflv->z);
    }
    RACmat = &RACmat[3*nfunc_ac];
  }
  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*_g1h_ReflectExtAltConstrMatrixf*/

boolean g1h_NLExtFillHoleAltConstrf ( GHoleDomainf *domain, const point3f *hole_cp,
                         int naconstr, const float *constr,
                         float *acoeff, void *usrptr,
                         void (*outpatch) ( int n, int m, const point3f *cp,
                                            void *usrptr ) )
{
  typedef void outscf ( int n, int m, const float *cp, void *usrptr );

  void     *sp;
  G1HNLPrivatef     *nlprivate;
  G1HolePrivateRecf *privateG1;
  float    *fc00, *Bi, *Bk;
  float    *saveconstr = NULL, *newconstr, *ECmat, *rsconstr;
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
  _g1h_nlprivf = nlprivate = _g1h_InitExtNLprf ( hole_k, nfunc_a, naconstr );
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
  if ( !g1h_ExtFillHoleAltConstrf ( domain, 3, (float*)nlprivate->rhole_cp,
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

  if ( !_g1h_TabExtNLBasisFunctionsf ( domain, nlprivate ) )
    goto failure;

  if ( !g1h_ExtNLConstrNewtonf ( domain, nlprivate, naconstr, ECmat ) )
    goto failure;

  g1h_ReflectVectorsf ( nfunc_ac, nlprivate->acoeff, nlprivate->acoeff );
  if ( acoeff )
    memcpy ( acoeff, nlprivate->acoeff, nfunc_ac*sizeof(vector3f) );

  fc00 = pkv_GetScratchMem ( (G1_CROSSDEGSUM+4)*2*hole_k*sizeof(vector3f) );
  if ( !fc00 )
    goto failure;
  Bi = privateG1->EBmat;
  Bk = &Bi[hole_k*G1_DBDIM*nfunc_b];
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
} /*g1h_NLExtFillHoleAltConstrf*/

/* ///////////////////////////////////////////////////////////////////////// */
static float _g1h_ComputeExtLFuncf ( GHoleDomainf *domain,
                                     G1HNLPrivatef *nlprivate, float *coeff )
{
  G1HolePrivateRecf *privateG1;
  int   k, i, kn, knot, fi, hole_k, nfunc_a, nfunc_c;
  float *psiuu, *psivv, *jac;
  float funct, lap;

  privateG1 = domain->privateG1;
  hole_k  = domain->hole_k;
  nfunc_a = privateG1->nfunc_a;
  nfunc_c = hole_k*G1_DBDIM;
  psiuu = nlprivate->psiuu;
  psivv = nlprivate->psivv;
  jac   = nlprivate->jac;
  funct = 0.0;
  for ( k = knot = 0;  k < hole_k;  k++ )
    for ( kn = 0;  kn < G1_NQUADSQ;  kn++, knot++ ) {
      fi = nfunc_a*G1_NQUADSQ*hole_k+knot;
      lap = psiuu[fi]+psivv[fi];
      for ( i = 0; i < nfunc_a; i++ ) {
        fi = i*G1_NQUADSQ*hole_k+knot;
        lap -= coeff[nfunc_c+i]*(psiuu[fi]+psivv[fi]);
      }
      for ( i = 0; i < G1_DBDIM; i++ ) {
        fi = ((nfunc_a+1)*hole_k+k*G1_DBDIM+i)*G1_NQUADSQ+kn;
        lap -= coeff[k*G1_DBDIM+i]*(psiuu[fi]+psivv[fi]);
      }
      funct += lap*lap*jac[knot];
    }

  return (float)(funct/(4.0*(float)G1_NQUADSQ));
} /*_g1h_ComputeExtLFuncf*/

boolean g1h_NLExtFunctionalValuef ( GHoleDomainf *domain,
                                    const point3f *hole_cp, const vector3f *acoeff,
                                    float *funcval )
{
  void  *sp;
  G1HNLPrivatef     *nlprivate;
  G1HolePrivateRecf *privateG1;
  int   hole_k, nfunc_a, nfunc_b, nfunc_c;
  float *fc00, *coeff, *Bi, *Bk;

  sp = pkv_GetScratchMemTop ();
  privateG1 = domain->privateG1;
  hole_k = domain->hole_k;
  nfunc_a = privateG1->nfunc_a;
  nfunc_b = privateG1->nfunc_b;
  nfunc_c = hole_k*G1_DBDIM;
  _g1h_nlprivf = nlprivate = _g1h_InitExtNLprf ( hole_k, nfunc_a, 0 );
  if ( !nlprivate )
    goto failure;

  if ( !_g1h_ComputeNLNormalf ( domain, nlprivate, hole_cp ) )
    goto failure;
  g1h_ReflectVectorsf ( 12*hole_k+1, hole_cp, nlprivate->rhole_cp );
  g1h_ReflectVectorsf ( nfunc_a+nfunc_c, acoeff, nlprivate->acoeff );

  fc00 = pkv_GetScratchMem ( (G1_CROSSDEGSUM+4)*2*hole_k*sizeof(vector3f) );
  coeff = pkv_GetScratchMemf ( nfunc_a+nfunc_c );
  if ( !fc00 || !coeff )
    goto failure;
  Bi = privateG1->EBmat;
  Bk = &Bi[hole_k*G1_DBDIM*nfunc_b];
  if ( !_g1h_SetExtRightSidef ( domain, Bi, Bk, 3, (float*)nlprivate->rhole_cp,
                                fc00, NULL ) )
    goto failure;
  nlprivate->auxc = 0;
  if ( !_g1h_OutputExtPatchesf ( domain, 3,
                (float*)nlprivate->acoeff, fc00, NULL, g1h_nonlinoutpatchf ) )
    goto failure;
  if ( !_g1h_TabExtNLBasisFunctionsf ( domain, nlprivate ) )
    goto failure;
  pkv_Selectf ( nfunc_a+nfunc_c, 1, 3, 1, &nlprivate->acoeff[0].z, coeff );
  funcval[0] = _g1h_ComputeExtLFuncf ( domain, nlprivate, coeff );
  funcval[1] = _g1h_ComputeExtNLFuncf ( domain, nlprivate, coeff );

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*g1h_NLExtFunctionalValuef*/

