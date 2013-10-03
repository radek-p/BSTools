
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


static G1HNLPrivated *_g1h_InitNLprd ( int hole_k, int V0SpaceDim, int nconstr )
{
#define N (G1H_FINALDEG+1)*(G1H_FINALDEG+1)
  G1HNLPrivated *nlpr;

  if ( (nlpr = pkv_GetScratchMem ( sizeof(G1HNLPrivated) )) ) {
    nlpr->auxc = 0;
    nlpr->nldi = pkv_GetScratchMem ( N*hole_k*sizeof(point3d) );
    nlpr->acoeff = pkv_GetScratchMem ( V0SpaceDim*sizeof(vector3d) );
    nlpr->rhole_cp = pkv_GetScratchMem ( (12*hole_k+1+nconstr)*sizeof(vector3d) );
    nlpr->diu = pkv_GetScratchMem ( 5*G1_NQUADSQ*hole_k*sizeof(vector2d) );
    nlpr->jac = pkv_GetScratchMemd ( G1_NQUADSQ*hole_k );
    nlpr->psiu = pkv_GetScratchMemd ( 5*(V0SpaceDim+1)*G1_NQUADSQ*hole_k );

    if ( !nlpr->nldi || !nlpr->acoeff || !nlpr->rhole_cp ||
         !nlpr->diu || !nlpr->jac || !nlpr->psiu )
      return NULL;

    nlpr->div = &nlpr->diu[G1_NQUADSQ*hole_k];
    nlpr->diuu = &nlpr->div[G1_NQUADSQ*hole_k];
    nlpr->diuv = &nlpr->diuu[G1_NQUADSQ*hole_k];
    nlpr->divv = &nlpr->diuv[G1_NQUADSQ*hole_k];
    nlpr->psiv = &nlpr->psiu[(V0SpaceDim+1)*G1_NQUADSQ*hole_k];
    nlpr->psiuu = &nlpr->psiv[(V0SpaceDim+1)*G1_NQUADSQ*hole_k];
    nlpr->psiuv = &nlpr->psiuu[(V0SpaceDim+1)*G1_NQUADSQ*hole_k];
    nlpr->psivv = &nlpr->psiuv[(V0SpaceDim+1)*G1_NQUADSQ*hole_k];
  }
  return nlpr;
#undef N
} /*_g1h_InitNLprd*/

static double _g1h_ComputeNLFuncd ( GHoleDomaind *domain,
                                   G1HNLPrivated *nlprivate,
                                   const double *coeff )
{
  G1HolePrivateRecd *privateG1;
  int    k, i, fi, hole_k, nfunc_a;
  double pu, pv, puu, puv, pvv, A, B, c, cs, funct;
  double *psiu, *psiv, *psiuu, *psiuv, *psivv, *jac;

  privateG1 = domain->privateG1;
  hole_k = domain->hole_k;
  nfunc_a = privateG1->nfunc_a;
  psiu = nlprivate->psiu;
  psiv = nlprivate->psiv;
  psiuu = nlprivate->psiuu;
  psiuv = nlprivate->psiuv;
  psivv = nlprivate->psivv;
  jac = nlprivate->jac;

  funct = 0.0;
  for ( k = 0; k < hole_k*G1_NQUADSQ; k++ ) {
        /* evaluate the function and its derivatives */
    fi = nfunc_a*G1_NQUADSQ*hole_k+k;
    pu = psiu[fi];    pv = psiv[fi];
    puu = psiuu[fi];  puv = psiuv[fi];  pvv = psivv[fi];
    for ( i = 0; i < nfunc_a; i++ ) {
      fi = i*G1_NQUADSQ*hole_k+k;
      c = coeff[i];
      pu -= c*psiu[fi];    pv -= c*psiv[fi];
      puu -= c*psiuu[fi];  puv -= c*psiuv[fi];  pvv -= c*psivv[fi];
    }
        /* make the integration */
    _g1h_IntFunc1d ( pu, pv, puu, puv, pvv, jac[k],
                     &c, &cs, &A, &B, &funct );
  }
  funct /= (double)G1_NQUADSQ;
  return funct;
} /*_g1h_ComputeNLFuncd*/

static boolean _g1h_ComputeNLFuncGradd ( GHoleDomaind *domain,
                                         G1HNLPrivated *nlprivate,
                                         const double *coeff,
                                         double *func, double *grad )
{
  void  *sp;
  G1HolePrivateRecd *privateG1;
  int    i, k, fi, hole_k, nfunc_a;
  double pu, pv, puu, puv, pvv, A, B, c, cs, funct;
  double *psiu, *psiv, *psiuu, *psiuv, *psivv, *jac, *Ai, *Bi;

  sp = pkv_GetScratchMemTop ();
  privateG1 = domain->privateG1;
  hole_k = domain->hole_k;
  nfunc_a = privateG1->nfunc_a;
  psiu = nlprivate->psiu;
  psiv = nlprivate->psiv;
  psiuu = nlprivate->psiuu;
  psiuv = nlprivate->psiuv;
  psivv = nlprivate->psivv;
  jac = nlprivate->jac;
  Ai = pkv_GetScratchMemd ( nfunc_a );
  Bi = pkv_GetScratchMemd ( nfunc_a );
  if ( !Ai || !Bi )
    goto failure;

  funct = 0.0;
  memset ( grad, 0, nfunc_a*sizeof(double) );
  for ( k = 0; k < hole_k*G1_NQUADSQ; k++ ) {
        /* evaluate the function and its derivatives */
    fi = nfunc_a*G1_NQUADSQ*hole_k+k;
    pu = psiu[fi];    pv = psiv[fi];
    puu = psiuu[fi];  puv = psiuv[fi];  pvv = psivv[fi];
    for ( i = 0; i < nfunc_a; i++ ) {
      fi = i*G1_NQUADSQ*hole_k+k;
      c = coeff[i];
      pu -= c*psiu[fi];    pv -= c*psiv[fi];
      puu -= c*psiuu[fi];  puv -= c*psiuv[fi];  pvv -= c*psivv[fi];
    }
        /* make the integration */
    _g1h_IntFunc1d ( pu, pv, puu, puv, pvv, jac[k],
                     &c, &cs, &A, &B, &funct );
/*
if ( k == 20 )
  Verify ( 3.0/(2*G1_NQUAD), 9.0/(2*G1_NQUAD), p, pu, pv, puu, puv, pvv, A/(c*cs) );
*/
        /* integrate the gradient */
    for ( i = 0; i < nfunc_a; i++ ) {
      fi = i*G1_NQUADSQ*hole_k+k;
      _g1h_IntFunc2d ( pu, pv, puu, puv, pvv,
                       psiu[fi], psiv[fi], psiuu[fi], psiuv[fi], psivv[fi],
                       jac[k], c, A, B, &Ai[i], &Bi[i], &grad[i] );
    }
  }
  *func = funct/(double)G1_NQUADSQ;
  pkn_MultMatrixNumd ( 1, nfunc_a, 0, grad, 1.0/(double)G1_NQUADSQ, 0, grad );

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*_g1h_ComputeNLFuncGradd*/

static boolean _g1h_ComputeNLFuncGradHessiand ( GHoleDomaind *domain,
                          G1HNLPrivated *nlprivate,
                          const double *coeff,
                          double *func, double *grad, double *hessian )
{
  void   *sp;
  G1HolePrivateRecd *privateG1;
  int    i, j, k, fi, fj, hole_k, nfunc_a;
  double pu, pv, puu, puv, pvv, A, B, c, cs, funct;
  double *psiu, *psiv, *psiuu, *psiuv, *psivv, *jac, *Ai, *Bi;

  sp = pkv_GetScratchMemTop ();
  privateG1 = domain->privateG1;
  hole_k = domain->hole_k;
  nfunc_a = privateG1->nfunc_a;
  psiu = nlprivate->psiu;
  psiv = nlprivate->psiv;
  psiuu = nlprivate->psiuu;
  psiuv = nlprivate->psiuv;
  psivv = nlprivate->psivv;
  jac = nlprivate->jac;
  Ai = pkv_GetScratchMemd ( nfunc_a );
  Bi = pkv_GetScratchMemd ( nfunc_a );
  if ( !Ai || !Bi )
    goto failure;

  funct = 0.0;
  memset ( grad, 0, nfunc_a*sizeof(double) );
  memset ( hessian, 0, (nfunc_a*(nfunc_a+1)/2)*sizeof(double) );
  for ( k = 0; k < hole_k*G1_NQUADSQ; k++ ) {
        /* evaluate the function and its derivatives */
    fi = nfunc_a*G1_NQUADSQ*hole_k+k;
    pu = psiu[fi];    pv = psiv[fi];
    puu = psiuu[fi];  puv = psiuv[fi];  pvv = psivv[fi];
    for ( i = 0; i < nfunc_a; i++ ) {
      fi = i*G1_NQUADSQ*hole_k+k;
      c = coeff[i];
      pu -= c*psiu[fi];    pv -= c*psiv[fi];
      puu -= c*psiuu[fi];  puv -= c*psiuv[fi];  pvv -= c*psivv[fi];
    }
        /* make the integration */
    _g1h_IntFunc1d ( pu, pv, puu, puv, pvv, jac[k],
                     &c, &cs, &A, &B, &funct );
/*
if ( k == 20 )
  Verify ( 3.0/(2*G1_NQUAD), 9.0/(2*G1_NQUAD), p, pu, pv, puu, puv, pvv, A/(c*cs) );
*/
        /* integrate the gradient */
    for ( i = 0; i < nfunc_a; i++ ) {
      fi = i*G1_NQUADSQ*hole_k+k;
      _g1h_IntFunc2d ( pu, pv, puu, puv, pvv,
                       psiu[fi], psiv[fi], psiuu[fi], psiuv[fi], psivv[fi],
                       jac[k], c, A, B, &Ai[i], &Bi[i], &grad[i] );
    }

        /* integrate the Hessian */
    for ( i = 0; i < nfunc_a; i++ ) {
      fi = i*G1_NQUADSQ*hole_k+k;
      for ( j = 0; j <= i; j++ ) {
        fj = j*G1_NQUADSQ*hole_k+k;
        _g1h_IntFunc3d ( pu, pv, puu, puv, pvv,
                         psiu[fi], psiv[fi], psiuu[fi], psiuv[fi], psivv[fi],
                         psiu[fj], psiv[fj], psiuu[fj], psiuv[fj], psivv[fj],
                         jac[k], c, A, B, Ai[i], Bi[i], Ai[j], Bi[j],
                         &hessian[pkn_SymMatIndex(i,j)] );
      }
    }
  }
  *func = funct/(double)G1_NQUADSQ;
  pkn_MultMatrixNumd ( 1, nfunc_a, 0, grad, 1.0/(double)G1_NQUADSQ, 0, grad );
  pkn_MultMatrixNumd ( 1, nfunc_a*(nfunc_a+1)/2, 0, hessian,
                       1.0/(double)G1_NQUADSQ, 0, hessian );

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*_g1h_ComputeNLFuncGradHessiand*/

/* ///////////////////////////////////////////////////////////////////////// */
static boolean g1h_NLNewtond ( GHoleDomaind *domain, G1HNLPrivated *nlprivate )
{
#define EPSF 2.0e-4
  void    *sp;
  G1HolePrivateRecd *privateG1;
  int     itn, jtn, ktn, nfunc_a;
  double  *coeff, *dcoeff, *grad, *hessian, *chess;
  double  func, func0, func1, aux, gn, gn0 = 0.0, gn1, dco, dyn;
  boolean positive;

  sp = pkv_GetScratchMemTop ();
  privateG1 = domain->privateG1;
  nfunc_a = privateG1->nfunc_a;
  coeff   = pkv_GetScratchMemd ( 3*nfunc_a );
  hessian = pkv_GetScratchMemd ( nfunc_a*(nfunc_a+1) );
  if ( !coeff || !hessian )
    goto failure;
  dcoeff = &coeff[nfunc_a];
  grad   = &dcoeff[nfunc_a];
  chess  = &hessian[(nfunc_a*(nfunc_a+1))/2];

        /* setup the initial point */
  pkv_Selectd ( nfunc_a, 1, 3, 1, &nlprivate->acoeff[0].z, coeff );

        /* Newton iterations */
  for ( itn = ktn = 0; ; itn++ ) {
    if ( !_g1h_ComputeNLFuncGradHessiand ( domain, nlprivate, coeff,
                                           &func, grad, hessian ) )
      goto failure;
    gn = (double)sqrt ( pkn_ScalarProductd ( nfunc_a, grad, grad ) );
    if ( itn == 0 ) {
/*
printf ( "func = %f, gn0 = %f\n", func, gn );
*/
      gn0 = gn;
    }
    memcpy ( chess, hessian, (nfunc_a*(nfunc_a+1))/2*sizeof(double) );

    if ( (positive = pkn_CholeskyDecompd ( nfunc_a, hessian )) ) {
      pkn_LowerTrMatrixSolved ( nfunc_a, hessian, 1, 1, grad, 1, grad );
      pkn_UpperTrMatrixSolved ( nfunc_a, hessian, 1, 1, grad, 1, grad );
    }
    else {

printf ( "!" );

      pkn_SymMatrixMultd ( nfunc_a, chess, 1, 1, grad, 1, dcoeff );
      aux = (double)pkn_ScalarProductd ( nfunc_a, grad, dcoeff );
      if ( gn < 0.0 || aux < EPSF*gn )
        goto failure;
      pkn_MultMatrixNumd ( 1, nfunc_a, 0, grad, gn/aux, 0, grad );
    }
    dco = (double)sqrt ( pkn_ScalarProductd(nfunc_a, coeff, coeff) );
    dyn = (double)sqrt ( pkn_ScalarProductd(nfunc_a, grad, grad) );

    for ( aux = 1.0; aux >= EPSF; aux *= 0.5 ) {
      pkn_AddMatrixMd ( 1, nfunc_a, 0, coeff, 0, grad, aux, 0, dcoeff );
      func1 = _g1h_ComputeNLFuncd ( domain, nlprivate, dcoeff );
      if ( func1 < func )
        break;
    }
    memcpy ( coeff, dcoeff, nfunc_a*sizeof(double) );
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
        _g1h_ComputeNLFuncGradd ( domain, nlprivate, coeff, &func0, grad );
        gn1 = (double)sqrt ( pkn_ScalarProductd ( nfunc_a, grad, grad ) );
        pkn_LowerTrMatrixSolved ( nfunc_a, hessian, 1, 1, grad, 1, grad );
        pkn_UpperTrMatrixSolved ( nfunc_a, hessian, 1, 1, grad, 1, grad );
        pkn_AddMatrixd ( 1, nfunc_a, 0, coeff, 0, grad, 0, dcoeff );
        func1 = _g1h_ComputeNLFuncd ( domain, nlprivate, dcoeff );
        if ( func1 >= func0 || gn1 >= 0.5*gn )
          break;
/*
printf ( "    func = %f, gn = %f\n", func1, gn );
*/
        memcpy ( coeff, dcoeff, nfunc_a*sizeof(double) );
        func = func1;
        ktn ++;
        gn = gn1;
        aux = 1.0;
        dco = (double)sqrt ( pkn_ScalarProductd(nfunc_a, coeff, coeff) );
        dyn = (double)sqrt ( pkn_ScalarProductd(nfunc_a, grad, grad) );
      }
    }

    if ( _g1h_StopItd ( itn, gn0, gn, dco, dyn, aux ) )
      break;
  }

printf ( "itn = %d, ktn = %d\n", itn+1, ktn );
/*
printf ( "func = %f\n", func1 );
*/
  pkv_Selectd ( nfunc_a, 1, 1, 3, coeff, &nlprivate->acoeff[0].z );

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
#undef EPSF
} /*g1h_NLNewtond*/

boolean g1h_NLFillHoled ( GHoleDomaind *domain, const point3d *hole_cp,  
                          double *acoeff, void *usrptr,
                          void (*outpatch) ( int n, int m, const point3d *cp,
                                             void *usrptr ) )
{
  typedef void outscf ( int n, int m, const double *cp, void *usrptr );

  void     *sp;
  G1HNLPrivated     *nlprivate;
  G1HolePrivateRecd *privateG1;
  double   *fc00;
  int      hole_k, nfunc_a;

  sp = pkv_GetScratchMemTop ();
  privateG1 = domain->privateG1;
  hole_k = domain->hole_k;
  nfunc_a = privateG1->nfunc_a;
  _g1h_nlprivd = nlprivate = _g1h_InitNLprd ( hole_k, nfunc_a, 0 );
  if ( !nlprivate )
    goto failure;

  if ( !_g1h_ComputeNLNormald ( domain, nlprivate, hole_cp ) )
    goto failure;
  g1h_ReflectVectorsd ( 12*hole_k+1, hole_cp, nlprivate->rhole_cp );

  nlprivate->auxc = 0;
  if ( !g1h_FillHoled ( domain, 3, (double*)nlprivate->rhole_cp,
                        (double*)nlprivate->acoeff, NULL, g1h_nonlinoutpatchd ) )
    goto failure;

  if ( !_g1h_TabNLBasisFunctionsd ( domain, G1_NQUAD, nlprivate ) )
    goto failure;

  if ( !g1h_NLNewtond ( domain, nlprivate ) )
    goto failure;

  g1h_ReflectVectorsd ( nfunc_a, nlprivate->acoeff, nlprivate->acoeff );
  if ( acoeff )
    memcpy ( acoeff, nlprivate->acoeff, nfunc_a*sizeof(vector3d) );

  fc00 = pkv_GetScratchMem ( (G1_CROSSDEGSUM+4)*2*hole_k*sizeof(vector3d) );
  if ( !fc00 )
    goto failure;
  if ( !_g1h_SetRightSided ( domain, privateG1->Bmat,
                             3, (double*)hole_cp, fc00, NULL ) )
    goto failure;

  if ( !_g1h_OutputPatchesd ( domain, 3, (double*)nlprivate->acoeff, fc00,
                              usrptr, (outscf*)outpatch ) )
    goto failure;

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*g1h_NLFillHoled*/

/* ///////////////////////////////////////////////////////////////////////// */
static boolean g1h_NLConstrNewtond ( GHoleDomaind *domain,
                                     G1HNLPrivated *nlprivate, int nconstr,
                                     double *Cmat )
{
#define EPSF 2.0e-4
  void    *sp;
  G1HolePrivateRecd *privateG1;
  int     itn, jtn, ktn, i, j, hole_k, nfunc_a, nfunc;
  double  *coeff, *grad, *hessian;
  double  *cT, *E22, *cE22, *aa, *D1, *y, *y1, *M, *f, *f2;
  double  func, func0, func1, aux, t, gn, gn0 = 0.0, gn1, dco, dyn;
  boolean positive;

  sp = pkv_GetScratchMemTop ();
  hole_k  = domain->hole_k;
  privateG1 = domain->privateG1;
  nfunc_a = privateG1->nfunc_a;
  nfunc   = nfunc_a-nconstr;
  coeff   = pkv_GetScratchMemd ( nfunc_a );
  grad    = pkv_GetScratchMemd ( nfunc_a );
  hessian = pkv_GetScratchMemd ( (nfunc_a*(nfunc_a+1))/2 );
  cT      = pkv_GetScratchMemd ( nfunc_a*nconstr );
  aa      = pkv_GetScratchMemd ( 2*nconstr );
  D1      = pkv_GetScratchMemd ( (nconstr*(nconstr+1))/2);
  E22     = pkv_GetScratchMemd ( (nfunc*(nfunc+1))/2 );
  cE22    = pkv_GetScratchMemd ( (nfunc*(nfunc+1))/2 );
  y       = pkv_GetScratchMemd ( nfunc_a );
  y1      = pkv_GetScratchMemd ( nfunc_a );
  M       = pkv_GetScratchMemd ( (nfunc_a*(nfunc_a+1))/2 );
  f       = pkv_GetScratchMemd ( nfunc_a );
  if ( !coeff || !grad || !hessian || !cT || !aa ||
       !D1 || !E22 || !cE22 || !y || !y1 || !M || !f )
    goto failure;

        /* setup the initial point */
  pkv_Selectd ( nfunc_a, 1, 3, 1, &nlprivate->acoeff[0].z, coeff );

        /* step 1: decompose the constraint equations matrix */
  pkv_TransposeMatrixd ( nconstr, nfunc_a, nfunc_a, Cmat, nconstr, cT );
  pkn_QRDecomposeMatrixd ( nfunc_a, nconstr, cT, aa );
  for ( i = 0; i < nconstr; i++ )
    for ( j = i; j < nconstr; j++ )
      D1[pkn_SymMatIndex(i,j)] = cT[i*nconstr+j];

        /* Newton iterations */
  for ( itn = ktn = 0; ; itn++ ) {
          /* step 2 */
    memcpy ( y, coeff, nfunc_a*sizeof(double) );
    pkn_multiReflectVectord ( nfunc_a, nconstr, cT, aa, 1, 1, y );
    pkn_LowerTrMatrixSolved ( nconstr, D1, 1,
                              3, &nlprivate->rhole_cp[12*hole_k+1].z,
                              1, y );
    for ( i = 0; i < nconstr; i++ )
      y[i] = -y[i];

          /* step 3 */
    if ( !_g1h_ComputeNLFuncGradHessiand ( domain, nlprivate, coeff,
                                           &func, grad, hessian ) )
      goto failure;
    if ( !pkn_ComputeQTSQd ( nfunc_a, hessian, nconstr, cT, aa, M ) )
      goto failure;
    for ( i = 0; i < nfunc; i++ )
      for ( j = i; j < nfunc; j++ )
        E22[pkn_SymMatIndex(i,j)] = M[pkn_SymMatIndex(nconstr+i,nconstr+j)];
    memcpy ( cE22, E22, ((nfunc*(nfunc+1))/2)*sizeof(double) );

          /* step 4 */
    memcpy ( f, grad, nfunc_a*sizeof(double) );
    pkn_multiReflectVectord ( nfunc_a, nconstr, cT, aa, 1, 1, f );
    f2 = &f[nconstr];
    gn = (double)sqrt ( pkn_ScalarProductd ( nfunc, f2, f2 ) );
    if ( itn == 0 ) {
/*
printf ( "func = %f, gn0 = %f\n", func, gn );
*/
      gn0 = gn;
    }

          /* step 5 */
    if ( (positive = pkn_CholeskyDecompd ( nfunc, E22 )) ) {
      pkn_LowerTrMatrixSolved ( nfunc, E22, 1, 1, f2, 1, f2 );
      pkn_UpperTrMatrixSolved ( nfunc, E22, 1, 1, f2, 1, f2 );
    }
    else {

printf ( "! " );

      pkn_SymMatrixMultd ( nfunc, cE22, 1, 1, f2, 1, y1 );
      aux = (double)pkn_ScalarProductd ( nfunc, f2, y1 );
      t = (double)pkn_ScalarProductd ( nfunc, f2, f2 );
      if ( aux <= 0.0 || aux < EPSF*t )
        goto failure;
      pkn_MultMatrixNumd ( 1, nfunc, 1, f2, t/aux, 1, f2 );
    }

          /* step 6 */
    dco = (double)sqrt ( pkn_ScalarProductd ( nfunc_a, coeff, coeff ) );
    dyn = (double)sqrt ( pkn_ScalarProductd ( nfunc, f2, f2 ) );
    for ( aux = 1.0; aux >= EPSF; aux *= 0.5 ) {
      memcpy ( y1, y, nconstr*sizeof(double) );
      for ( i = nconstr; i < nfunc_a; i++ )
        y1[i] = y[i]+aux*f[i];
      pkn_multiInvReflectVectord ( nfunc_a, nconstr, cT, aa, 1, 1, y1 );
      func1 = _g1h_ComputeNLFuncd ( domain, nlprivate, y1 );
      if ( func1 < func )
        break;
    }
    memcpy ( coeff, y1, nfunc_a*sizeof(double) );
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
        memcpy ( y, coeff, nfunc_a*sizeof(double) );
        pkn_multiReflectVectord ( nfunc_a, nconstr, cT, aa, 1, 1, y );
        pkn_LowerTrMatrixSolved ( nconstr, D1, 1,
                                  3, &nlprivate->rhole_cp[12*hole_k+1].z,
                                  1, y );
        for ( i = 0; i < nconstr; i++ )
          y[i] = -y[i];
              /* step 3' */
        _g1h_ComputeNLFuncGradd ( domain, nlprivate, coeff, &func0, grad );
              /* step 4' */
        memcpy ( f, grad, nfunc_a*sizeof(double) );
        pkn_multiReflectVectord ( nfunc_a, nconstr, cT, aa, 1, 1, f );
        f2 = &f[nconstr];
        gn1 = (double)sqrt ( pkn_ScalarProductd ( nfunc, f2, f2 ) );
              /* step 5' */
        pkn_LowerTrMatrixSolved ( nfunc, E22, 1, 1, f2, 1, f2 );
        pkn_UpperTrMatrixSolved ( nfunc, E22, 1, 1, f2, 1, f2 );
              /* step 6' */
        memcpy ( y1, y, nconstr*sizeof(double) );
        for ( i = nconstr; i < nfunc_a; i++ )
          y1[i] = y[i]+f[i];
        pkn_multiInvReflectVectord ( nfunc_a, nconstr, cT, aa, 1, 1, y1 );

        func1 = _g1h_ComputeNLFuncd ( domain, nlprivate, y1 );
        if ( func1 >= func0 || gn1 >= 0.5*gn )
          break;
/*
printf ( "    func = %f, gn = %f\n", func1, gn );
*/
        memcpy ( coeff, y1, nfunc_a*sizeof(double) );
        func = func1;
        ktn ++;
        gn = gn1;
        aux = 1.0;
        dco = (double)sqrt ( pkn_ScalarProductd(nfunc_a, coeff, coeff) );
        dyn = (double)sqrt ( pkn_ScalarProductd(nfunc, f2, f2) );
      }
    }
    if ( _g1h_StopItd ( itn, gn0, gn, dco, dyn, aux ) )
      break;
  }

printf ( "itn = %d, ktn = %d\n", itn+1, ktn );
/*
printf ( "func = %f\n", func1 );
*/
  pkv_Selectd ( nfunc_a, 1, 1, 3, coeff, &nlprivate->acoeff[0].z );

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
#undef EPSF
} /*g1h_NLConstrNewtond*/

boolean g1h_NLFillHoleConstrd ( GHoleDomaind *domain, const point3d *hole_cp,
                    int nconstr, const vector3d *constr,
                    double *acoeff, void *usrptr,
                    void (*outpatch) ( int n, int m, const point3d *cp,
                                       void *usrptr ) )
{
  typedef void outscf ( int n, int m, const double *cp, void *usrptr );

  void     *sp;
  G1HNLPrivated     *nlprivate;
  G1HolePrivateRecd *privateG1;
  double   *fc00;
  int      hole_k, nfunc_a;

  sp = pkv_GetScratchMemTop ();
  privateG1 = domain->privateG1;
  hole_k  = domain->hole_k;
  nfunc_a = privateG1->nfunc_a;
  _g1h_nlprivd = nlprivate = _g1h_InitNLprd ( hole_k, nfunc_a, nconstr );
  if ( !nlprivate )
    goto failure;

  if ( !_g1h_ComputeNLNormald ( domain, nlprivate, hole_cp ) )
    goto failure;
  g1h_ReflectVectorsd ( 12*hole_k+1, hole_cp, nlprivate->rhole_cp );
  g1h_ReflectVectorsd ( nconstr, constr, &nlprivate->rhole_cp[12*hole_k+1] );

  nlprivate->auxc = 0;
  if ( !g1h_FillHoleConstrd ( domain, 3, (double*)nlprivate->rhole_cp,
                    nconstr, (double*)&nlprivate->rhole_cp[12*hole_k+1],
                    (double*)nlprivate->acoeff, NULL, g1h_nonlinoutpatchd ) )
    goto failure;

  if ( !_g1h_TabNLBasisFunctionsd ( domain, G1_NQUAD, nlprivate ) )
    goto failure;

  if ( !g1h_NLConstrNewtond ( domain, nlprivate, nconstr, privateG1->Cmat ) )
    goto failure;

  g1h_ReflectVectorsd ( nfunc_a, nlprivate->acoeff, nlprivate->acoeff );
  if ( acoeff )
    memcpy ( acoeff, nlprivate->acoeff, nfunc_a*sizeof(vector3d) );

  fc00 = pkv_GetScratchMem ( (G1_CROSSDEGSUM+4)*2*hole_k*sizeof(vector3d) );
  if ( !fc00 )
    goto failure;
  if ( !_g1h_SetRightSided ( domain, privateG1->Bmat,
                             3, (double*)hole_cp, fc00, NULL ) )
    goto failure;

  if ( !_g1h_OutputPatchesd ( domain, 3, (double*)nlprivate->acoeff, fc00,
                              usrptr, (outscf*)outpatch ) )
    goto failure;

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*g1h_NLFillHoleConstrd*/

/* ///////////////////////////////////////////////////////////////////////// */
boolean _g1h_ReflectAltConstrMatrixd ( GHoleDomaind *domain,
                                       vector3d *reflv, double *RACmat )
{
  void   *sp;
  G1HolePrivateRecd *privateG1;
  int    i, j, nconstr, nfunc_a;
  double *buf;

  sp = pkv_GetScratchMemTop ();
  privateG1 = domain->privateG1;
  nfunc_a = privateG1->nfunc_a;
  nconstr = privateG1->naconstr;
  memcpy ( RACmat, privateG1->ACmat, nconstr*nfunc_a*3*sizeof(double) );
  buf = pkv_GetScratchMemd ( nfunc_a );
  if ( !buf ) {
    domain->error_code = G1H_ERROR_NO_SCRATCH_MEMORY;
    goto failure;
  }
  for ( i = 0; i < nconstr; i++ ) {
    pkn_MultMatrixNumd ( 1, nfunc_a, 0, RACmat, reflv->x, 0, buf );
    pkn_AddMatrixMd ( 1, nfunc_a, 0, buf, 0, &RACmat[nfunc_a], reflv->y, 0, buf );
    pkn_AddMatrixMd ( 1, nfunc_a, 0, buf, 0, &RACmat[2*nfunc_a], reflv->z, 0, buf );
    for ( j = 0; j < nfunc_a; j++ ) {
      RACmat[j]           -= (double)(2.0*buf[j]*reflv->x);
      RACmat[nfunc_a+j]   -= (double)(2.0*buf[j]*reflv->y);
      RACmat[2*nfunc_a+j] -= (double)(2.0*buf[j]*reflv->z);
    }
    RACmat = &RACmat[3*nfunc_a];
  }
  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*_g1h_ReflectAltConstrMatrixd*/

boolean g1h_NLFillHoleAltConstrd ( GHoleDomaind *domain, const point3d *hole_cp,
                    int nconstr, const double *constr,
                    double *acoeff, void *usrptr,
                    void (*outpatch) ( int n, int m, const point3d *cp,
                                       void *usrptr ) )
{
  typedef void outscf ( int n, int m, const double *cp, void *usrptr );

  void     *sp;
  G1HNLPrivated     *nlprivate;
  G1HolePrivateRecd *privateG1;
  double   *fc00;
  double   *saveconstr = NULL, *newconstr, *Cmat, *rsconstr;
  int      hole_k, nfunc_a, i, j;
  boolean  restore;

  sp = pkv_GetScratchMemTop ();
  restore = false;
  privateG1 = domain->privateG1;
  if ( nconstr <= 0 || privateG1->acdim != 3 || privateG1->naconstr != nconstr ) {
    domain->error_code = G1H_ERROR_INCONSISTENT_CONSTR;
    goto failure;
  }
  hole_k  = domain->hole_k;
  nfunc_a = privateG1->nfunc_a;
  saveconstr = pkv_GetScratchMemd ( 7*nfunc_a*nconstr );
  _g1h_nlprivd = nlprivate = _g1h_InitNLprd ( hole_k, nfunc_a, nconstr );
  if ( !saveconstr || !nlprivate )
    goto failure;

  newconstr = &saveconstr[3*nfunc_a*nconstr];
  Cmat = &newconstr[3*nfunc_a*nconstr];
  memcpy ( saveconstr, privateG1->ACmat, nfunc_a*nconstr*3*sizeof(double) );
  memcpy ( newconstr, privateG1->ACmat, nfunc_a*nconstr*3*sizeof(double) );
  if ( !_g1h_ComputeNLNormald ( domain, nlprivate, hole_cp ) )
    goto failure;
  g1h_ReflectVectorsd ( 12*hole_k+1, hole_cp, nlprivate->rhole_cp );
  if ( !_g1h_ReflectAltConstrMatrixd ( domain, &nlprivate->reflv, newconstr ) )
    goto failure;

  restore = true;
  if ( !g1h_SetAltConstraintMatrixd ( domain, 3, nconstr, newconstr ) )
    goto failure;
  nlprivate->auxc = 0;
  if ( !g1h_FillHoleAltConstrd ( domain, 3, (double*)nlprivate->rhole_cp,
                    nconstr, constr,
                    (double*)nlprivate->acoeff, NULL, g1h_nonlinoutpatchd ) )
    goto failure;

  pkv_Selectd ( nconstr, nfunc_a, 3*nfunc_a, nfunc_a,
                &newconstr[2*nfunc_a], Cmat );
  rsconstr = &nlprivate->rhole_cp[12*hole_k+1].z;
  pkv_Selectd ( nconstr, 1, 1, 3, constr, rsconstr );
  for ( j = 0; j < nconstr; j++ )
    for ( i = 0; i < nfunc_a; i++ )
      rsconstr[3*j] += newconstr[3*nfunc_a*j+i]*nlprivate->acoeff[i].x +
                       newconstr[(3*j+1)*nfunc_a+i]*nlprivate->acoeff[i].y;

  if ( !_g1h_TabNLBasisFunctionsd ( domain, G1_NQUAD, nlprivate ) )
    goto failure;

  if ( !g1h_NLConstrNewtond ( domain, nlprivate, nconstr, Cmat ) )
    goto failure;

  g1h_ReflectVectorsd ( nfunc_a, nlprivate->acoeff, nlprivate->acoeff );
  if ( acoeff )
    memcpy ( acoeff, nlprivate->acoeff, nfunc_a*sizeof(vector3d) );

  fc00 = pkv_GetScratchMem ( (G1_CROSSDEGSUM+4)*2*hole_k*sizeof(vector3d) );
  if ( !fc00 )
    goto failure;
  if ( !_g1h_SetRightSided ( domain, privateG1->Bmat,
                             3, (double*)hole_cp, fc00, NULL ) )
    goto failure;

  if ( !_g1h_OutputPatchesd ( domain, 3, (double*)nlprivate->acoeff, fc00,
                              usrptr, (outscf*)outpatch ) )
    goto failure;

  g1h_SetAltConstraintMatrixd ( domain, 3, nconstr, saveconstr );
  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  if ( restore )
    g1h_SetAltConstraintMatrixd ( domain, 3, nconstr, saveconstr );
  pkv_SetScratchMemTop ( sp );
  return false;
} /*g1h_NLFillHoleAltConstrd*/

/* ///////////////////////////////////////////////////////////////////////// */
static double _g1h_ComputeLFuncd ( GHoleDomaind *domain,
                                  G1HNLPrivated *nlprivate, double *coeff )
{
  G1HolePrivateRecd *privateG1;
  int    k, i, fi, hole_k, nfunc_a;
  double *psiuu, *psivv, *jac;
  double funct, lap;

  privateG1 = domain->privateG1;
  hole_k  = domain->hole_k;
  nfunc_a = privateG1->nfunc_a;
  psiuu = nlprivate->psiuu;
  psivv = nlprivate->psivv;
  jac   = nlprivate->jac;
  funct = 0.0;
  for ( k = 0; k < hole_k*G1_NQUADSQ; k++ ) {
    fi = nfunc_a*G1_NQUADSQ*hole_k+k;
    lap = psiuu[fi]+psivv[fi];
    for ( i = 0; i < nfunc_a; i++ ) {
      fi = i*G1_NQUADSQ*hole_k+k;
      lap -= coeff[i]*(psiuu[fi]+psivv[fi]);
    }
    funct += lap*lap*jac[k];
  }

  return (double)(funct/(4.0*(double)G1_NQUADSQ));
} /*_g1h_ComputeLFuncd*/

boolean g1h_NLFunctionalValued ( GHoleDomaind *domain,
                                 const point3d *hole_cp, const vector3d *acoeff,
                                 double *funcval )
{
  void   *sp;
  G1HNLPrivated     *nlprivate;
  G1HolePrivateRecd *privateG1;
  int    hole_k, nfunc_a;
  double *fc00, *coeff;

  sp = pkv_GetScratchMemTop ();
  privateG1 = domain->privateG1;
  hole_k = domain->hole_k;
  nfunc_a = privateG1->nfunc_a;
  _g1h_nlprivd = nlprivate = _g1h_InitNLprd ( hole_k, nfunc_a, 0 );
  if ( !nlprivate )
    goto failure;

  if ( !_g1h_ComputeNLNormald ( domain, nlprivate, hole_cp ) )
    goto failure;
  g1h_ReflectVectorsd ( 12*hole_k+1, hole_cp, nlprivate->rhole_cp );
  g1h_ReflectVectorsd ( nfunc_a, acoeff, nlprivate->acoeff );

  fc00 = pkv_GetScratchMem ( (G1_CROSSDEGSUM+4)*2*hole_k*sizeof(vector3d) );
  coeff = pkv_GetScratchMemd ( nfunc_a );
  if ( !fc00 || !coeff )
    goto failure;
  if ( !_g1h_SetRightSided ( domain, privateG1->Bmat,
                             3, (double*)nlprivate->rhole_cp, fc00, NULL ) )
    goto failure;
  nlprivate->auxc = 0;
  if ( !_g1h_OutputPatchesd ( domain, 3,
                (double*)nlprivate->acoeff, fc00, NULL, g1h_nonlinoutpatchd ) )
    goto failure;
  if ( !_g1h_TabNLBasisFunctionsd ( domain, G1_NQUAD, nlprivate ) )
    goto failure;
  pkv_Selectd ( nfunc_a, 1, 3, 1, &nlprivate->acoeff[0].z, coeff );
  funcval[0] = _g1h_ComputeLFuncd ( domain, nlprivate, coeff );
  funcval[1] = _g1h_ComputeNLFuncd ( domain, nlprivate, coeff );

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*g1h_NLFunctionalValued*/

