
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

#include "eg2holef.h"
#include "eg2hprivatef.h"
#include "eg2herror.h"


/* ///////////////////////////////////////////////////////////////////////// */
static G2HNLPrivatef *_g2h_InitNLprf ( int hole_k, int V0SpaceDim, int nconstr )
{
#define N (G2H_FINALDEG+1)*(G2H_FINALDEG+1)
  G2HNLPrivatef *nlpr;

  if ( (nlpr = pkv_GetScratchMem ( sizeof(G2HNLPrivatef) )) ) {
    nlpr->auxc = 0;
    nlpr->nldi = pkv_GetScratchMem ( N*hole_k*sizeof(point3f) );
    nlpr->acoeff = pkv_GetScratchMem ( V0SpaceDim*sizeof(vector3f) );
    nlpr->rhole_cp = pkv_GetScratchMem ( (12*hole_k+1+nconstr)*sizeof(vector3f) );
    nlpr->diu = pkv_GetScratchMem ( 9*G2_NQUADSQ*hole_k*sizeof(vector2f) );
    nlpr->jac = pkv_GetScratchMemf ( G2_NQUADSQ*hole_k );
    nlpr->psiu = pkv_GetScratchMemf ( 9*(V0SpaceDim+1)*G2_NQUADSQ*hole_k );

    if ( !nlpr->nldi || !nlpr->acoeff || !nlpr->rhole_cp ||
         !nlpr->diu || !nlpr->jac || !nlpr->psiu )
      return NULL;

    nlpr->div = &nlpr->diu[G2_NQUADSQ*hole_k];
    nlpr->diuu = &nlpr->div[G2_NQUADSQ*hole_k];
    nlpr->diuv = &nlpr->diuu[G2_NQUADSQ*hole_k];
    nlpr->divv = &nlpr->diuv[G2_NQUADSQ*hole_k];
    nlpr->diuuu = &nlpr->divv[G2_NQUADSQ*hole_k];
    nlpr->diuuv = &nlpr->diuuu[G2_NQUADSQ*hole_k];
    nlpr->diuvv = &nlpr->diuuv[G2_NQUADSQ*hole_k];
    nlpr->divvv = &nlpr->diuvv[G2_NQUADSQ*hole_k];
    nlpr->psiv = &nlpr->psiu[(V0SpaceDim+1)*G2_NQUADSQ*hole_k];
    nlpr->psiuu = &nlpr->psiv[(V0SpaceDim+1)*G2_NQUADSQ*hole_k];
    nlpr->psiuv = &nlpr->psiuu[(V0SpaceDim+1)*G2_NQUADSQ*hole_k];
    nlpr->psivv = &nlpr->psiuv[(V0SpaceDim+1)*G2_NQUADSQ*hole_k];
    nlpr->psiuuu = &nlpr->psivv[(V0SpaceDim+1)*G2_NQUADSQ*hole_k];
    nlpr->psiuuv = &nlpr->psiuuu[(V0SpaceDim+1)*G2_NQUADSQ*hole_k];
    nlpr->psiuvv = &nlpr->psiuuv[(V0SpaceDim+1)*G2_NQUADSQ*hole_k];
    nlpr->psivvv = &nlpr->psiuvv[(V0SpaceDim+1)*G2_NQUADSQ*hole_k];
  }
  return nlpr;
#undef N
} /*_g2h_InitNLprf*/

static float _g2h_ComputeNLFuncf ( GHoleDomainf *domain, G2HNLPrivatef *nlprivate,
                                   const float *coeff )
{
  G2HolePrivateRecf *privateG2;
  G2HNLFuncf f;
  int   k, i, fi, hole_k, nfunc_a;
  float c, funct;

  privateG2 = domain->privateG2;
  hole_k = domain->hole_k;
  nfunc_a = privateG2->nfunc_a;

  funct = 0.0;
  for ( k = 0; k < hole_k*G2_NQUADSQ; k++ ) {
        /* evaluate the function and its derivatives */
    fi = nfunc_a*G2_NQUADSQ*hole_k+k;
    f.pu = nlprivate->psiu[fi];      f.pv = nlprivate->psiv[fi];
    f.puu = nlprivate->psiuu[fi];    f.puv = nlprivate->psiuv[fi];
    f.pvv = nlprivate->psivv[fi];    f.puuu = nlprivate->psiuuu[fi];
    f.puuv = nlprivate->psiuuv[fi];  f.puvv = nlprivate->psiuvv[fi];
    f.pvvv = nlprivate->psivvv[fi];
    for ( i = 0; i < nfunc_a; i++ ) {
      fi = i*G2_NQUADSQ*hole_k+k;
      c = coeff[i];
      f.pu -= c*nlprivate->psiu[fi];      f.pv -= c*nlprivate->psiv[fi];
      f.puu -= c*nlprivate->psiuu[fi];    f.puv -= c*nlprivate->psiuv[fi];
      f.pvv -= c*nlprivate->psivv[fi];    f.puuu -= c*nlprivate->psiuuu[fi];
      f.puuv -= c*nlprivate->psiuuv[fi];  f.puvv -= c*nlprivate->psiuvv[fi];
      f.pvvv -= c*nlprivate->psivvv[fi];
    }
    f.jac = nlprivate->jac[k];
        /* make the integration */
    _g2h_IntFunc1af ( &f, &funct );
  }
  funct /= (float)G2_NQUADSQ;
  return funct;
} /*_g2h_ComputeNLFuncf*/

static boolean _g2h_ComputeNLFuncGradf ( GHoleDomainf *domain,
                                         G2HNLPrivatef *nlprivate,
                                         const float *coeff,
                                         float *func, float *grad )
{
  void  *sp;
  G2HolePrivateRecf *privateG2;
  G2HNLFuncf f;
  int   i, k, fi, hole_k, nfunc_a;
  float c, funct;

  sp = pkv_GetScratchMemTop ();
  privateG2 = domain->privateG2;
  hole_k = domain->hole_k;
  nfunc_a = privateG2->nfunc_a;

  funct = 0.0;
  memset ( grad, 0, nfunc_a*sizeof(float) );
  for ( k = 0; k < hole_k*G2_NQUADSQ; k++ ) {
        /* evaluate the function and its derivatives */
    fi = nfunc_a*G2_NQUADSQ*hole_k+k;
    f.pu = nlprivate->psiu[fi];      f.pv = nlprivate->psiv[fi];
    f.puu = nlprivate->psiuu[fi];    f.puv = nlprivate->psiuv[fi];
    f.pvv = nlprivate->psivv[fi];    f.puuu = nlprivate->psiuuu[fi];
    f.puuv = nlprivate->psiuuv[fi];  f.puvv = nlprivate->psiuvv[fi];
    f.pvvv = nlprivate->psivvv[fi];
    for ( i = 0; i < nfunc_a; i++ ) {
      fi = i*G2_NQUADSQ*hole_k+k;
      c = coeff[i];
      f.pu -= c*nlprivate->psiu[fi];      f.pv -= c*nlprivate->psiv[fi];
      f.puu -= c*nlprivate->psiuu[fi];    f.puv -= c*nlprivate->psiuv[fi];
      f.pvv -= c*nlprivate->psivv[fi];    f.puuu -= c*nlprivate->psiuuu[fi];
      f.puuv -= c*nlprivate->psiuuv[fi];  f.puvv -= c*nlprivate->psiuvv[fi];
      f.pvvv -= c*nlprivate->psivvv[fi];
    }
    f.jac = nlprivate->jac[k];
        /* make the integration */
    _g2h_IntFunc1bf ( &f, &funct );

        /* integrate the gradient */
    for ( i = 0; i < nfunc_a; i++ ) {
      fi = i*G2_NQUADSQ*hole_k+k;
      f.psiu = nlprivate->psiu[fi];      f.psiv = nlprivate->psiv[fi];
      f.psiuu = nlprivate->psiuu[fi];    f.psiuv = nlprivate->psiuv[fi];
      f.psivv = nlprivate->psivv[fi];    f.psiuuu = nlprivate->psiuuu[fi];
      f.psiuuv = nlprivate->psiuuv[fi];  f.psiuvv = nlprivate->psiuvv[fi];
      f.psivvv = nlprivate->psivvv[fi];
      _g2h_IntFunc2bf ( &f, &grad[i] );
    }
  }
  *func = funct/(float)G2_NQUADSQ;
  pkn_MultMatrixNumf ( 1, nfunc_a, 0, grad, 1.0/(float)G2_NQUADSQ, 0, grad );

  pkv_SetScratchMemTop ( sp );
  return true;

/*
failure:
  pkv_SetScratchMemTop ( sp );
  return false;
*/
} /*_g2h_ComputeNLFuncGradf*/

static boolean _g2h_ComputeNLFuncGradHessianf ( GHoleDomainf *domain,
                          G2HNLPrivatef *nlprivate,
                          const float *coeff,
                          float *func, float *grad, float *hessian )
{
  void  *sp;
  G2HolePrivateRecf *privateG2;
  G2HNLFuncf f;
  int   i, j, k, fi, fj, hole_k, nfunc_a;
  float *Li, *Bi, *BiLT, *Di;
  float c, funct;

  sp = pkv_GetScratchMemTop ();
  privateG2 = domain->privateG2;
  hole_k = domain->hole_k;
  nfunc_a = privateG2->nfunc_a;
  Li = pkv_GetScratchMemf ( 8*nfunc_a );
  if ( !Li ) {
    domain->error_code = G2H_ERROR_NO_SCRATCH_MEMORY;
    goto failure;
  }
  Bi = &Li[2*nfunc_a];  BiLT = &Bi[3*nfunc_a];  Di = &BiLT[2*nfunc_a];

  funct = 0.0;
  memset ( grad, 0, nfunc_a*sizeof(float) );
  memset ( hessian, 0, (nfunc_a*(nfunc_a+1)/2)*sizeof(float) );
  for ( k = 0; k < hole_k*G2_NQUADSQ; k++ ) {
        /* evaluate the function and its derivatives */
    fi = nfunc_a*G2_NQUADSQ*hole_k+k;
    f.pu = nlprivate->psiu[fi];      f.pv = nlprivate->psiv[fi];
    f.puu = nlprivate->psiuu[fi];    f.puv = nlprivate->psiuv[fi];
    f.pvv = nlprivate->psivv[fi];    f.puuu = nlprivate->psiuuu[fi];
    f.puuv = nlprivate->psiuuv[fi];  f.puvv = nlprivate->psiuvv[fi];
    f.pvvv = nlprivate->psivvv[fi];
    for ( i = 0; i < nfunc_a; i++ ) {
      fi = i*G2_NQUADSQ*hole_k+k;
      c = coeff[i];
      f.pu -= c*nlprivate->psiu[fi];      f.pv -= c*nlprivate->psiv[fi];
      f.puu -= c*nlprivate->psiuu[fi];    f.puv -= c*nlprivate->psiuv[fi];
      f.pvv -= c*nlprivate->psivv[fi];    f.puuu -= c*nlprivate->psiuuu[fi];
      f.puuv -= c*nlprivate->psiuuv[fi];  f.puvv -= c*nlprivate->psiuvv[fi];
      f.pvvv -= c*nlprivate->psivvv[fi];
    }
    f.jac = nlprivate->jac[k];
        /* make the integration */
    _g2h_IntFunc1cf ( &f, &funct );

        /* integrate the gradient */
    for ( i = 0; i < nfunc_a; i++ ) {
      fi = i*G2_NQUADSQ*hole_k+k;
      f.psiu = nlprivate->psiu[fi];      f.psiv = nlprivate->psiv[fi];
      f.psiuu = nlprivate->psiuu[fi];    f.psiuv = nlprivate->psiuv[fi];
      f.psivv = nlprivate->psivv[fi];    f.psiuuu = nlprivate->psiuuu[fi];
      f.psiuuv = nlprivate->psiuuv[fi];  f.psiuvv = nlprivate->psiuvv[fi];
      f.psivvv = nlprivate->psivvv[fi];
      _g2h_IntFunc2cf ( &f, (vector2f*)(&Li[2*i]), &Bi[3*i],
                        (vector2f*)(&BiLT[2*i]), &Di[i], &grad[i] );
    }

        /* integrate the Hessian */
    for ( i = 0; i < nfunc_a; i++ ) {
      fi = i*G2_NQUADSQ*hole_k+k;
      f.psiu = nlprivate->psiu[fi];      f.psiv = nlprivate->psiv[fi];
      f.psiuu = nlprivate->psiuu[fi];    f.psiuv = nlprivate->psiuv[fi];
      f.psivv = nlprivate->psivv[fi];    f.psiuuu = nlprivate->psiuuu[fi];
      f.psiuuv = nlprivate->psiuuv[fi];  f.psiuvv = nlprivate->psiuvv[fi];
      f.psivvv = nlprivate->psivvv[fi];
      for ( j = 0; j <= i; j++ ) {
        fj = j*G2_NQUADSQ*hole_k+k;
        f.psju = nlprivate->psiu[fj];      f.psjv = nlprivate->psiv[fj];
        f.psjuu = nlprivate->psiuu[fj];    f.psjuv = nlprivate->psiuv[fj];
        f.psjvv = nlprivate->psivv[fj];    f.psjuuu = nlprivate->psiuuu[fj];
        f.psjuuv = nlprivate->psiuuv[fj];  f.psjuvv = nlprivate->psiuvv[fj];
        f.psjvvv = nlprivate->psivvv[fj];
        _g2h_IntFunc3cf ( &f, (vector2f*)(&Li[2*i]), (vector2f*)(&Li[2*j]),
                          (vector2f*)(&BiLT[2*i]),
                          (vector2f*)(&BiLT[2*j]), Di[i], Di[j],
                          &hessian[pkn_SymMatIndex(i,j)] );
      }
    }
  }
  *func = funct/(float)G2_NQUADSQ;
  pkn_MultMatrixNumf ( 1, nfunc_a, 0, grad, 1.0/(float)G2_NQUADSQ, 0, grad );
  pkn_MultMatrixNumf ( 1, nfunc_a*(nfunc_a+1)/2, 0, hessian,
                       1.0/(float)G2_NQUADSQ, 0, hessian );

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*_g2h_ComputeNLFuncGradHessianf*/

/* ///////////////////////////////////////////////////////////////////////// */
static boolean g2h_NLNewtonf ( GHoleDomainf *domain, G2HNLPrivatef *nlprivate )
{
#define EPSF 2.0e-4
  void    *sp;
  G2HolePrivateRecf *privateG2;
  int     itn, jtn, ktn, nfunc_a;
  float   *coeff, *dcoeff, *grad, *hessian, *chess;
  float   func, func0, func1, aux, gn, gn0 = 0.0, gn1, dco, dyn;
  boolean positive;

  sp = pkv_GetScratchMemTop ();
  privateG2 = domain->privateG2;
  nfunc_a = privateG2->nfunc_a;
  coeff   = pkv_GetScratchMemf ( 3*nfunc_a );
  hessian = pkv_GetScratchMemf ( nfunc_a*(nfunc_a+1) );
  if ( !coeff || !hessian ) {
    domain->error_code = G2H_ERROR_NO_SCRATCH_MEMORY;
    goto failure;
  }
  dcoeff = &coeff[nfunc_a];
  grad   = &dcoeff[nfunc_a];
  chess  = &hessian[(nfunc_a*(nfunc_a+1))/2];

        /* setup the initial point */
  pkv_Selectf ( nfunc_a, 1, 3, 1, &nlprivate->acoeff[0].z, coeff );

        /* Newton iterations */
  for ( itn = ktn = 0; ; itn++ ) {
    if ( !_g2h_ComputeNLFuncGradHessianf ( domain, nlprivate, coeff,
                                           &func, grad, hessian ) )
      goto failure;
    gn = (float)sqrt ( pkn_ScalarProductf ( nfunc_a, grad, grad ) );
    if ( itn == 0 ) {
/*
printf ( "func = %f, gn0 = %f\n", func, gn );
*/
      gn0 = gn;
    }
    memcpy ( chess, hessian, (nfunc_a*(nfunc_a+1))/2*sizeof(float) );

    if ( (positive = pkn_CholeskyDecompf ( nfunc_a, hessian )) ) {
      pkn_LowerTrMatrixSolvef ( nfunc_a, hessian, 1, 1, grad, 1, grad );
      pkn_UpperTrMatrixSolvef ( nfunc_a, hessian, 1, 1, grad, 1, grad );
    }
    else {

printf ( "!" );

      pkn_SymMatrixMultf ( nfunc_a, chess, 1, 1, grad, 1, dcoeff );
      aux = (float)pkn_ScalarProductf ( nfunc_a, grad, dcoeff );
      if ( gn < 0.0 || aux < EPSF*gn ) {
        domain->error_code = G2H_ERROR_NL_MINIMIZATION;
        goto failure;
      }
      pkn_MultMatrixNumf ( 1, nfunc_a, 0, grad, gn/aux, 0, grad );
    }
    dco = (float)sqrt ( pkn_ScalarProductf(nfunc_a, coeff, coeff) );
    dyn = (float)sqrt ( pkn_ScalarProductf(nfunc_a, grad, grad) );

    for ( aux = 1.0; aux >= EPSF; aux *= 0.5 ) {
      pkn_AddMatrixMf ( 1, nfunc_a, 0, coeff, 0, grad, aux, 0, dcoeff );
      func1 = _g2h_ComputeNLFuncf ( domain, nlprivate, dcoeff );
      if ( func1 < func )
        break;
    }
    memcpy ( coeff, dcoeff, nfunc_a*sizeof(float) );
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
        _g2h_ComputeNLFuncGradf ( domain, nlprivate, coeff, &func0, grad );
        gn1 = (float)sqrt ( pkn_ScalarProductf ( nfunc_a, grad, grad ) );
        pkn_LowerTrMatrixSolvef ( nfunc_a, hessian, 1, 1, grad, 1, grad );
        pkn_UpperTrMatrixSolvef ( nfunc_a, hessian, 1, 1, grad, 1, grad );
        pkn_AddMatrixf ( 1, nfunc_a, 0, coeff, 0, grad, 0, dcoeff );
        func1 = _g2h_ComputeNLFuncf ( domain, nlprivate, dcoeff );
        if ( func1 >= func0 || gn1 >= 0.5*gn )
          break;
/*
printf ( "    func = %f, gn = %f\n", func1, gn );
*/
        memcpy ( coeff, dcoeff, nfunc_a*sizeof(float) );
        func = func1;
        ktn ++;
        gn = gn1;
        aux = 1.0;
        dco = (float)sqrt ( pkn_ScalarProductf(nfunc_a, coeff, coeff) );
        dyn = (float)sqrt ( pkn_ScalarProductf(nfunc_a, grad, grad) );
      }
    }

    if ( _g2h_StopItf ( itn, gn0, gn, dco, dyn, aux ) )
      break;
  }
/*
printf ( "itn = %d, ktn = %d\n", itn+1, ktn );
*/
  pkv_Selectf ( nfunc_a, 1, 1, 3, coeff, &nlprivate->acoeff[0].z );

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
#undef EPSF
} /*g2h_NLNewtonf*/

boolean g2h_NLFillHolef ( GHoleDomainf *domain, const point3f *hole_cp,  
                          float *acoeff, void *usrptr,
                          void (*outpatch) ( int n, int m, const point3f *cp,
                                             void *usrptr ) )
{
  typedef void outscf ( int n, int m, const float *cp, void *usrptr );

  void     *sp;
  G2HNLPrivatef     *nlprivate;
  G2HolePrivateRecf *privateG2;
  float    *fc00;
  int      hole_k, nfunc_a;

  sp = pkv_GetScratchMemTop ();
  if ( !g2h_DecomposeMatrixf ( domain ) )
    goto failure;

  privateG2 = domain->privateG2;
  hole_k = domain->hole_k;
  nfunc_a = privateG2->nfunc_a;
  _g2h_nlprivf = nlprivate = _g2h_InitNLprf ( hole_k, nfunc_a, 0 );
  if ( !nlprivate )
    goto failure;

  if ( !_g2h_ComputeNLNormalf ( domain, nlprivate, hole_cp ) )
    goto failure;
  g2h_ReflectVectorsf ( 12*hole_k+1, hole_cp, nlprivate->rhole_cp );

  nlprivate->auxc = 0;
  if ( !g2h_FillHolef ( domain, 3, (float*)nlprivate->rhole_cp,
                        (float*)nlprivate->acoeff, NULL, g2h_nonlinoutpatchf ) )
    goto failure;

  if ( !_g2h_TabNLBasisFunctionsf ( domain, nlprivate ) )
    goto failure;

  if ( !g2h_NLNewtonf ( domain, nlprivate ) )
    goto failure;

  g2h_ReflectVectorsf ( nfunc_a, nlprivate->acoeff, nlprivate->acoeff );
  if ( acoeff )
    memcpy ( acoeff, nlprivate->acoeff, nfunc_a*sizeof(vector3f) );

  fc00 = pkv_GetScratchMem ( (G2_CROSSDEGSUM+6)*2*hole_k*sizeof(vector3f) );
  if ( !fc00 ) {
    domain->error_code = G2H_ERROR_NO_SCRATCH_MEMORY;
    goto failure;
  }
  if ( !_g2h_SetRightSidef ( domain, 3, (float*)hole_cp, fc00, NULL ) )
    goto failure;

  if ( !_g2h_OutputPatchesf ( domain, 3, (float*)nlprivate->acoeff, fc00,
                              usrptr, (outscf*)outpatch ) )
    goto failure;

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*g2h_NLFillHolef*/

/* ///////////////////////////////////////////////////////////////////////// */
static boolean g2h_NLConstrNewtonf ( GHoleDomainf *domain,
                                     G2HNLPrivatef *nlprivate, int nconstr,
                                     float *Cmat )
{
#define EPSF 2.0e-4
  void    *sp;
  G2HolePrivateRecf *privateG2;
  int     itn, jtn, ktn, i, j, hole_k, nfunc_a, nfunc;
  float   *coeff, *grad, *hessian;
  float   *cT, *E22, *cE22, *aa, *D1, *y, *y1, *M, *f, *f2;
  float   func, func0, func1, aux, t, gn, gn0 = 0.0, gn1, dco, dyn;
  boolean positive;

  sp = pkv_GetScratchMemTop ();
  hole_k  = domain->hole_k;
  privateG2 = domain->privateG2;
  nfunc_a = privateG2->nfunc_a;
  nfunc   = nfunc_a-nconstr;
  coeff   = pkv_GetScratchMemf ( nfunc_a );
  grad    = pkv_GetScratchMemf ( nfunc_a );
  hessian = pkv_GetScratchMemf ( (nfunc_a*(nfunc_a+1))/2 );
  cT      = pkv_GetScratchMemf ( nfunc_a*nconstr );
  aa      = pkv_GetScratchMemf ( 2*nconstr );
  D1      = pkv_GetScratchMemf ( (nconstr*(nconstr+1))/2);
  E22     = pkv_GetScratchMemf ( (nfunc*(nfunc+1))/2 );
  cE22    = pkv_GetScratchMemf ( (nfunc*(nfunc+1))/2 );
  y       = pkv_GetScratchMemf ( nfunc_a );
  y1      = pkv_GetScratchMemf ( nfunc_a );
  M       = pkv_GetScratchMemf ( (nfunc_a*(nfunc_a+1))/2 );
  f       = pkv_GetScratchMemf ( nfunc_a );
  if ( !coeff || !grad || !hessian || !cT || !aa ||
       !D1 || !E22 || !cE22 || !y || !y1 || !M || !f ) {
    domain->error_code = G2H_ERROR_NO_SCRATCH_MEMORY;
    goto failure;
  }

        /* setup the initial point */
  pkv_Selectf ( nfunc_a, 1, 3, 1, &nlprivate->acoeff[0].z, coeff );

        /* step 1: decompose the constraint equations matrix */
  pkv_TransposeMatrixf ( nconstr, nfunc_a, nfunc_a, Cmat, nconstr, cT );
  pkn_QRDecomposeMatrixf ( nfunc_a, nconstr, cT, aa );
  for ( i = 0; i < nconstr; i++ )
    for ( j = i; j < nconstr; j++ )
      D1[pkn_SymMatIndex(i,j)] = cT[i*nconstr+j];

        /* Newton iterations */
  for ( itn = ktn = 0; ; itn++ ) {
          /* step 2 */
    memcpy ( y, coeff, nfunc_a*sizeof(float) );
    pkn_multiReflectVectorf ( nfunc_a, nconstr, cT, aa, 1, 1, y );
    pkn_LowerTrMatrixSolvef ( nconstr, D1, 1,
                              3, &nlprivate->rhole_cp[12*hole_k+1].z,
                              1, y );
    for ( i = 0; i < nconstr; i++ )
      y[i] = -y[i];

          /* step 3 */
    if ( !_g2h_ComputeNLFuncGradHessianf ( domain, nlprivate, coeff,
                                           &func, grad, hessian ) )
      goto failure;
    if ( !pkn_ComputeQTSQf ( nfunc_a, hessian, nconstr, cT, aa, M ) )
      goto failure;
    for ( i = 0; i < nfunc; i++ )
      for ( j = i; j < nfunc; j++ )
        E22[pkn_SymMatIndex(i,j)] = M[pkn_SymMatIndex(nconstr+i,nconstr+j)];
    memcpy ( cE22, E22, ((nfunc*(nfunc+1))/2)*sizeof(float) );

          /* step 4 */
    memcpy ( f, grad, nfunc_a*sizeof(float) );
    pkn_multiReflectVectorf ( nfunc_a, nconstr, cT, aa, 1, 1, f );
    f2 = &f[nconstr];
    gn = (float)sqrt ( pkn_ScalarProductf ( nfunc, f2, f2 ) );
    if ( itn == 0 ) {
/*
printf ( "func = %f, gn0 = %f\n", func, gn );
*/
      gn0 = gn;
    }

          /* step 5 */
    if ( (positive = pkn_CholeskyDecompf ( nfunc, E22 )) ) {
      pkn_LowerTrMatrixSolvef ( nfunc, E22, 1, 1, f2, 1, f2 );
      pkn_UpperTrMatrixSolvef ( nfunc, E22, 1, 1, f2, 1, f2 );
    }
    else {
/*
printf ( "! " );
*/
      pkn_SymMatrixMultf ( nfunc, cE22, 1, 1, f2, 1, y1 );
      aux = (float)pkn_ScalarProductf ( nfunc, f2, y1 );
      t = (float)pkn_ScalarProductf ( nfunc, f2, f2 );
      if ( aux <= 0.0 || aux < EPSF*t ) {
        domain->error_code = G2H_ERROR_NL_MINIMIZATION;
        goto failure;
      }
      pkn_MultMatrixNumf ( 1, nfunc, 1, f2, t/aux, 1, f2 );
    }

          /* step 6 */
    dco = (float)sqrt ( pkn_ScalarProductf ( nfunc_a, coeff, coeff ) );
    dyn = (float)sqrt ( pkn_ScalarProductf ( nfunc, f2, f2 ) );
    for ( aux = 1.0; aux >= EPSF; aux *= 0.5 ) {
      memcpy ( y1, y, nconstr*sizeof(float) );
      for ( i = nconstr; i < nfunc_a; i++ )
        y1[i] = y[i]+aux*f[i];
      pkn_multiInvReflectVectorf ( nfunc_a, nconstr, cT, aa, 1, 1, y1 );
      func1 = _g2h_ComputeNLFuncf ( domain, nlprivate, y1 );
      if ( func1 < func )
        break;
    }
    memcpy ( coeff, y1, nfunc_a*sizeof(float) );
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
        memcpy ( y, coeff, nfunc_a*sizeof(float) );
        pkn_multiReflectVectorf ( nfunc_a, nconstr, cT, aa, 1, 1, y );
        pkn_LowerTrMatrixSolvef ( nconstr, D1, 1,
                                  3, &nlprivate->rhole_cp[12*hole_k+1].z,
                                  1, y );
        for ( i = 0; i < nconstr; i++ )
          y[i] = -y[i];
              /* step 3' */
        _g2h_ComputeNLFuncGradf ( domain, nlprivate, coeff, &func0, grad );
              /* step 4' */
        memcpy ( f, grad, nfunc_a*sizeof(float) );
        pkn_multiReflectVectorf ( nfunc_a, nconstr, cT, aa, 1, 1, f );
        f2 = &f[nconstr];
        gn1 = (float)sqrt ( pkn_ScalarProductf ( nfunc, f2, f2 ) );
              /* step 5' */
        pkn_LowerTrMatrixSolvef ( nfunc, E22, 1, 1, f2, 1, f2 );
        pkn_UpperTrMatrixSolvef ( nfunc, E22, 1, 1, f2, 1, f2 );
              /* step 6' */
        memcpy ( y1, y, nconstr*sizeof(float) );
        for ( i = nconstr; i < nfunc_a; i++ )
          y1[i] = y[i]+f[i];
        pkn_multiInvReflectVectorf ( nfunc_a, nconstr, cT, aa, 1, 1, y1 );

        func1 = _g2h_ComputeNLFuncf ( domain, nlprivate, y1 );
        if ( func1 >= func0 || gn1 >= 0.5*gn )
          break;
/*
printf ( "    func = %f, gn = %f\n", func1, gn );
*/
        memcpy ( coeff, y1, nfunc_a*sizeof(float) );
        func = func1;
        ktn ++;
        gn = gn1;
        aux = 1.0;
        dco = (float)sqrt ( pkn_ScalarProductf(nfunc_a, coeff, coeff) );
        dyn = (float)sqrt ( pkn_ScalarProductf(nfunc, f2, f2) );
      }
    }
    if ( _g2h_StopItf ( itn, gn0, gn, dco, dyn, aux ) )
      break;
  }
/*
printf ( "itn = %d, ktn = %d\n", itn+1, ktn );
*/
  pkv_Selectf ( nfunc_a, 1, 1, 3, coeff, &nlprivate->acoeff[0].z );

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
#undef EPSF
} /*g2h_NLConstrNewtonf*/

boolean g2h_NLFillHoleConstrf ( GHoleDomainf *domain, const point3f *hole_cp,
                    int nconstr, const vector3f *constr,
                    float *acoeff, void *usrptr,
                    void (*outpatch) ( int n, int m, const point3f *cp,
                                       void *usrptr ) )
{
  typedef void outscf ( int n, int m, const float *cp, void *usrptr );

  void     *sp;
  G2HNLPrivatef     *nlprivate;
  G2HolePrivateRecf *privateG2;
  float    *fc00;
  int      hole_k, nfunc_a;

  sp = pkv_GetScratchMemTop ();
  if ( !g2h_DecomposeMatrixf ( domain ) )
    goto failure;

  privateG2 = domain->privateG2;
  hole_k  = domain->hole_k;
  nfunc_a = privateG2->nfunc_a;
  _g2h_nlprivf = nlprivate = _g2h_InitNLprf ( hole_k, nfunc_a, nconstr );
  if ( !nlprivate )
    goto failure;

  if ( !_g2h_ComputeNLNormalf ( domain, nlprivate, hole_cp ) )
    goto failure;
  g2h_ReflectVectorsf ( 12*hole_k+1, hole_cp, nlprivate->rhole_cp );
  g2h_ReflectVectorsf ( nconstr, constr, &nlprivate->rhole_cp[12*hole_k+1] );

  nlprivate->auxc = 0;
  if ( !g2h_FillHoleConstrf ( domain, 3, (float*)nlprivate->rhole_cp,
                    nconstr, (float*)&nlprivate->rhole_cp[12*hole_k+1],
                    (float*)nlprivate->acoeff, NULL, g2h_nonlinoutpatchf ) )
    goto failure;

  if ( !_g2h_TabNLBasisFunctionsf ( domain, nlprivate ) )
    goto failure;

  if ( !g2h_NLConstrNewtonf ( domain, nlprivate, nconstr, privateG2->Cmat ) )
    goto failure;

  g2h_ReflectVectorsf ( nfunc_a, nlprivate->acoeff, nlprivate->acoeff );
  if ( acoeff )
    memcpy ( acoeff, nlprivate->acoeff, nfunc_a*sizeof(vector3f) );

  fc00 = pkv_GetScratchMem ( (G2_CROSSDEGSUM+6)*2*hole_k*sizeof(vector3f) );
  if ( !fc00 ) {
    domain->error_code = G2H_ERROR_NO_SCRATCH_MEMORY;
    goto failure;
  }
  if ( !_g2h_SetRightSidef ( domain, 3, (float*)hole_cp, fc00, NULL ) )
    goto failure;

  if ( !_g2h_OutputPatchesf ( domain, 3,
                              (float*)nlprivate->acoeff, fc00,
                              usrptr, (outscf*)outpatch ) )
    goto failure;

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*g2h_NLFillHoleConstrf*/

/* ///////////////////////////////////////////////////////////////////////// */
static boolean _g2h_ReflectAltConstrMatrixf ( GHoleDomainf *domain,
                                              vector3f *reflv,
                                              float *RACmat )
{
  void  *sp;
  G2HolePrivateRecf *privateG2;
  int   i, j, nconstr, nfunc_a;
  float *buf;

  sp = pkv_GetScratchMemTop ();
  privateG2 = domain->privateG2;
  nfunc_a = privateG2->nfunc_a;
  nconstr = privateG2->naconstr;
  memcpy ( RACmat, privateG2->ACmat, nconstr*nfunc_a*3*sizeof(float) );
  buf = pkv_GetScratchMemf ( nfunc_a );
  if ( !buf ) {
    domain->error_code = G2H_ERROR_NO_SCRATCH_MEMORY;
    goto failure;
  }
  for ( i = 0; i < nconstr; i++ ) {
    pkn_MultMatrixNumf ( 1, nfunc_a, 0, RACmat, reflv->x, 0, buf );
    pkn_AddMatrixMf ( 1, nfunc_a, 0, buf, 0, &RACmat[nfunc_a], reflv->y, 0, buf );
    pkn_AddMatrixMf ( 1, nfunc_a, 0, buf, 0, &RACmat[2*nfunc_a], reflv->z, 0, buf );
    for ( j = 0; j < nfunc_a; j++ ) {
      RACmat[j]           -= (float)(2.0*buf[j]*reflv->x);
      RACmat[nfunc_a+j]   -= (float)(2.0*buf[j]*reflv->y);
      RACmat[2*nfunc_a+j] -= (float)(2.0*buf[j]*reflv->z);
    }
    RACmat = &RACmat[3*nfunc_a];
  }
  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*_g2h_ReflectAltConstrMatrixf*/

boolean g2h_NLFillHoleAltConstrf ( GHoleDomainf *domain, const point3f *hole_cp,
                    int nconstr, const float *constr,
                    float *acoeff, void *usrptr,
                    void (*outpatch) ( int n, int m, const point3f *cp,
                                       void *usrptr ) )
{
  typedef void outscf ( int n, int m, const float *cp, void *usrptr );

  void     *sp;
  G2HNLPrivatef     *nlprivate;
  G2HolePrivateRecf *privateG2;
  float    *fc00;
  float    *saveconstr = NULL, *newconstr, *Cmat, *rsconstr;
  int      hole_k, nfunc_a, i, j;
  boolean  restore;

  sp = pkv_GetScratchMemTop ();
  restore = false;
  privateG2 = domain->privateG2;
  if ( nconstr <= 0 || privateG2->acdim != 3 || privateG2->naconstr != nconstr ) {
    domain->error_code = G2H_ERROR_INCONSISTENT_CONSTR;
    goto failure;
  }
  hole_k  = domain->hole_k;
  nfunc_a = privateG2->nfunc_a;
  saveconstr = pkv_GetScratchMemf ( 7*nfunc_a*nconstr );
  _g2h_nlprivf = nlprivate = _g2h_InitNLprf ( hole_k, nfunc_a, nconstr );
  if ( !saveconstr || !nlprivate )
    goto failure;

  newconstr = &saveconstr[3*nfunc_a*nconstr];
  Cmat = &newconstr[3*nfunc_a*nconstr];
  memcpy ( saveconstr, privateG2->ACmat, nfunc_a*nconstr*3*sizeof(float) );
  memcpy ( newconstr, privateG2->ACmat, nfunc_a*nconstr*3*sizeof(float) );
  if ( !_g2h_ComputeNLNormalf ( domain, nlprivate, hole_cp ) )
    goto failure;
  g2h_ReflectVectorsf ( 12*hole_k+1, hole_cp, nlprivate->rhole_cp );
  if ( !_g2h_ReflectAltConstrMatrixf ( domain, &nlprivate->reflv, newconstr ) )
    goto failure;

  restore = true;
  if ( !g2h_SetAltConstraintMatrixf ( domain, 3, nconstr, newconstr ) )
    goto failure;
  nlprivate->auxc = 0;
  if ( !g2h_FillHoleAltConstrf ( domain, 3, (float*)nlprivate->rhole_cp,
                    nconstr, constr,
                    (float*)nlprivate->acoeff, NULL, g2h_nonlinoutpatchf ) )
    goto failure;

  pkv_Selectf ( nconstr, nfunc_a, 3*nfunc_a, nfunc_a,
                &newconstr[2*nfunc_a], Cmat );
  rsconstr = &nlprivate->rhole_cp[12*hole_k+1].z;
  pkv_Selectf ( nconstr, 1, 1, 3, constr, rsconstr );
  for ( j = 0; j < nconstr; j++ )
    for ( i = 0; i < nfunc_a; i++ )
      rsconstr[3*j] += newconstr[3*nfunc_a*j+i]*nlprivate->acoeff[i].x +
                       newconstr[(3*j+1)*nfunc_a+i]*nlprivate->acoeff[i].y;

  if ( !_g2h_TabNLBasisFunctionsf ( domain, nlprivate ) )
    goto failure;

  if ( !g2h_NLConstrNewtonf ( domain, nlprivate, nconstr, Cmat ) )
    goto failure;

  g2h_ReflectVectorsf ( nfunc_a, nlprivate->acoeff, nlprivate->acoeff );
  if ( acoeff )
    memcpy ( acoeff, nlprivate->acoeff, nfunc_a*sizeof(vector3f) );

  fc00 = pkv_GetScratchMem ( (G2_CROSSDEGSUM+6)*2*hole_k*sizeof(vector3f) );
  if ( !fc00 ) {
    domain->error_code = G2H_ERROR_NO_SCRATCH_MEMORY;
    goto failure;
  }
  if ( !_g2h_SetRightSidef ( domain, 3, (float*)hole_cp, fc00, NULL ) )
    goto failure;

  if ( !_g2h_OutputPatchesf ( domain, 3,
                              (float*)nlprivate->acoeff, fc00,
                              usrptr, (outscf*)outpatch ) )
    goto failure;

  g2h_SetAltConstraintMatrixf ( domain, 3, nconstr, saveconstr );
  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  if ( restore )
    g2h_SetAltConstraintMatrixf ( domain, 3, nconstr, saveconstr );
  pkv_SetScratchMemTop ( sp );
  return false;
} /*g2h_NLFillHoleAltConstrf*/

/* ///////////////////////////////////////////////////////////////////////// */
static float _g2h_ComputeLFuncf ( GHoleDomainf *domain,
                                  G2HNLPrivatef *nlprivate, float *coeff )
{
  G2HolePrivateRecf *privateG2;
  int      k, i, fi, hole_k, nfunc_a;
  float    *psiuuu, *psiuuv, *psiuvv, *psivvv, *jac;
  float    funct;
  vector2f lapgrad;

  privateG2 = domain->privateG2;
  hole_k  = domain->hole_k;
  nfunc_a = privateG2->nfunc_a;
  psiuuu = nlprivate->psiuuu;
  psiuuv = nlprivate->psiuuv;
  psiuvv = nlprivate->psiuvv;
  psivvv = nlprivate->psivvv;
  jac    = nlprivate->jac;
  funct = 0.0;
  for ( k = 0; k < hole_k*G2_NQUADSQ; k++ ) {
    fi = nfunc_a*G2_NQUADSQ*hole_k+k;
    lapgrad.x = psiuuu[fi]+psiuvv[fi];
    lapgrad.y = psiuuv[fi]+psivvv[fi];
    for ( i = 0; i < nfunc_a; i++ ) {
      fi = i*G2_NQUADSQ*hole_k+k;
      lapgrad.x -= coeff[i]*(psiuuu[fi]+psiuvv[fi]);
      lapgrad.y -= coeff[i]*(psiuuv[fi]+psivvv[fi]);
    }
    funct += (lapgrad.x*lapgrad.x+lapgrad.y*lapgrad.y)*jac[k];
  }

  return (float)(funct/(4.0*(float)G2_NQUADSQ));
} /*_g2h_ComputeLFuncf*/

boolean g2h_NLFunctionalValuef ( GHoleDomainf *domain,
                                 const point3f *hole_cp, const vector3f *acoeff,
                                 float *funcval )
{
  void  *sp;
  G2HNLPrivatef     *nlprivate;
  G2HolePrivateRecf *privateG2;
  int   hole_k, nfunc_a;
  float *fc00, *coeff;

  sp = pkv_GetScratchMemTop ();
  privateG2 = domain->privateG2;
  hole_k = domain->hole_k;
  nfunc_a = privateG2->nfunc_a;
  _g2h_nlprivf = nlprivate = _g2h_InitNLprf ( hole_k, nfunc_a, 0 );
  if ( !nlprivate )
    goto failure;

  if ( !_g2h_ComputeNLNormalf ( domain, nlprivate, hole_cp ) )
    goto failure;
  g2h_ReflectVectorsf ( 12*hole_k+1, hole_cp, nlprivate->rhole_cp );
  g2h_ReflectVectorsf ( nfunc_a, acoeff, nlprivate->acoeff );

  fc00 = pkv_GetScratchMem ( (G2_CROSSDEGSUM+6)*2*hole_k*sizeof(vector3f) );
  coeff = pkv_GetScratchMemf ( nfunc_a );
  if ( !fc00 || !coeff ) {
    domain->error_code = G2H_ERROR_NO_SCRATCH_MEMORY;
    goto failure;
  }
  if ( !_g2h_SetRightSidef ( domain, 3, (float*)nlprivate->rhole_cp, fc00,
                             NULL ) )
    goto failure;
  nlprivate->auxc = 0;
  if ( !_g2h_OutputPatchesf ( domain, 3,
                (float*)nlprivate->acoeff, fc00, NULL, g2h_nonlinoutpatchf ) )
    goto failure;
  if ( !_g2h_TabNLBasisFunctionsf ( domain, nlprivate ) )
    goto failure;
  pkv_Selectf ( nfunc_a, 1, 3, 1, &nlprivate->acoeff[0].z, coeff );
  funcval[0] = _g2h_ComputeLFuncf ( domain, nlprivate, coeff );
  funcval[1] = _g2h_ComputeNLFuncf ( domain, nlprivate, coeff );

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*g2h_NLFunctionalValuef*/

