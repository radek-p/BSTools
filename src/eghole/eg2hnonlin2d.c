
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2005, 2009                            */
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

#include "eg2holed.h"
#include "eg2hprivated.h"
#include "eg2herror.h"


/* ///////////////////////////////////////////////////////////////////////// */
static G2HNLPrivated *_g2h_InitNLprd ( int hole_k, int V0SpaceDim, int nconstr )
{
#define N (G2H_FINALDEG+1)*(G2H_FINALDEG+1)
  G2HNLPrivated *nlpr;

  if ( (nlpr = pkv_GetScratchMem ( sizeof(G2HNLPrivated) )) ) {
    nlpr->auxc = 0;
    nlpr->nldi = pkv_GetScratchMem ( N*hole_k*sizeof(point3d) );
    nlpr->acoeff = pkv_GetScratchMem ( V0SpaceDim*sizeof(vector3d) );
    nlpr->rhole_cp = pkv_GetScratchMem ( (12*hole_k+1+nconstr)*sizeof(vector3d) );
    nlpr->diu = pkv_GetScratchMem ( 9*G2_NQUADSQ*hole_k*sizeof(vector2d) );
    nlpr->jac = pkv_GetScratchMemd ( G2_NQUADSQ*hole_k );
    nlpr->psiu = pkv_GetScratchMemd ( 9*(V0SpaceDim+1)*G2_NQUADSQ*hole_k );

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
} /*_g2h_InitNLprd*/

static double _g2h_ComputeNLFuncd ( GHoleDomaind *domain, G2HNLPrivated *nlprivate,
                                    const double *coeff )
{
  G2HolePrivateRecd *privateG2;
  G2HNLFuncd f;
  int    k, i, fi, hole_k, nfunc_a;
  double c, funct;

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
    _g2h_IntFunc1ad ( &f, &funct );
  }
  funct /= (double)G2_NQUADSQ;
  return funct;
} /*_g2h_ComputeNLFuncd*/

static boolean _g2h_ComputeNLFuncGradd ( GHoleDomaind *domain,
                                         G2HNLPrivated *nlprivate,
                                         const double *coeff,
                                         double *func, double *grad )
{
  void   *sp;
  G2HolePrivateRecd *privateG2;
  G2HNLFuncd f;
  int    i, k, fi, hole_k, nfunc_a;
  double c, funct;

  sp = pkv_GetScratchMemTop ();
  privateG2 = domain->privateG2;
  hole_k = domain->hole_k;
  nfunc_a = privateG2->nfunc_a;

  funct = 0.0;
  memset ( grad, 0, nfunc_a*sizeof(double) );
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
    _g2h_IntFunc1bd ( &f, &funct );

        /* integrate the gradient */
    for ( i = 0; i < nfunc_a; i++ ) {
      fi = i*G2_NQUADSQ*hole_k+k;
      f.psiu = nlprivate->psiu[fi];      f.psiv = nlprivate->psiv[fi];
      f.psiuu = nlprivate->psiuu[fi];    f.psiuv = nlprivate->psiuv[fi];
      f.psivv = nlprivate->psivv[fi];    f.psiuuu = nlprivate->psiuuu[fi];
      f.psiuuv = nlprivate->psiuuv[fi];  f.psiuvv = nlprivate->psiuvv[fi];
      f.psivvv = nlprivate->psivvv[fi];
      _g2h_IntFunc2bd ( &f, &grad[i] );
    }
  }
  *func = funct/(double)G2_NQUADSQ;
  pkn_MultMatrixNumd ( 1, nfunc_a, 0, grad, 1.0/(double)G2_NQUADSQ, 0, grad );

  pkv_SetScratchMemTop ( sp );
  return true;

/*
failure:
  pkv_SetScratchMemTop ( sp );
  return false;
*/
} /*_g2h_ComputeNLFuncGradd*/

static boolean _g2h_ComputeNLFuncGradHessiand ( GHoleDomaind *domain,
                          G2HNLPrivated *nlprivate,
                          const double *coeff,
                          double *func, double *grad, double *hessian )
{
  void   *sp;
  G2HolePrivateRecd *privateG2;
  G2HNLFuncd f;
  int    i, j, k, fi, fj, hole_k, nfunc_a;
  double *Li, *Bi, *BiLT, *Di;
  double c, funct;

  sp = pkv_GetScratchMemTop ();
  privateG2 = domain->privateG2;
  hole_k = domain->hole_k;
  nfunc_a = privateG2->nfunc_a;
  Li = pkv_GetScratchMemd ( 8*nfunc_a );
  if ( !Li ) {
    domain->error_code = G2H_ERROR_NO_SCRATCH_MEMORY;
    goto failure;
  }
  Bi = &Li[2*nfunc_a];  BiLT = &Bi[3*nfunc_a];  Di = &BiLT[2*nfunc_a];

  funct = 0.0;
  memset ( grad, 0, nfunc_a*sizeof(double) );
  memset ( hessian, 0, (nfunc_a*(nfunc_a+1)/2)*sizeof(double) );
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
    _g2h_IntFunc1cd ( &f, &funct );

        /* integrate the gradient */
    for ( i = 0; i < nfunc_a; i++ ) {
      fi = i*G2_NQUADSQ*hole_k+k;
      f.psiu = nlprivate->psiu[fi];      f.psiv = nlprivate->psiv[fi];
      f.psiuu = nlprivate->psiuu[fi];    f.psiuv = nlprivate->psiuv[fi];
      f.psivv = nlprivate->psivv[fi];    f.psiuuu = nlprivate->psiuuu[fi];
      f.psiuuv = nlprivate->psiuuv[fi];  f.psiuvv = nlprivate->psiuvv[fi];
      f.psivvv = nlprivate->psivvv[fi];
      _g2h_IntFunc2cd ( &f, (vector2d*)(&Li[2*i]), &Bi[3*i],
                        (vector2d*)(&BiLT[2*i]), &Di[i], &grad[i] );
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
        _g2h_IntFunc3cd ( &f, (vector2d*)(&Li[2*i]), (vector2d*)(&Li[2*j]),
                          (vector2d*)(&BiLT[2*i]),
                          (vector2d*)(&BiLT[2*j]), Di[i], Di[j],
                          &hessian[pkn_SymMatIndex(i,j)] );
      }
    }
  }
  *func = funct/(double)G2_NQUADSQ;
  pkn_MultMatrixNumd ( 1, nfunc_a, 0, grad, 1.0/(double)G2_NQUADSQ, 0, grad );
  pkn_MultMatrixNumd ( 1, nfunc_a*(nfunc_a+1)/2, 0, hessian,
                       1.0/(double)G2_NQUADSQ, 0, hessian );

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*_g2h_ComputeNLFuncGradHessiand*/

/* ///////////////////////////////////////////////////////////////////////// */
static boolean g2h_NLNewtond ( GHoleDomaind *domain, G2HNLPrivated *nlprivate )
{
#define EPSF 2.0e-4
  void    *sp;
  G2HolePrivateRecd *privateG2;
  int     itn, jtn, ktn, nfunc_a;
  double  *coeff, *dcoeff, *grad, *hessian, *chess;
  double  func, func0, func1, aux, gn, gn0 = 0.0, gn1, dco, dyn;
  boolean positive;

  sp = pkv_GetScratchMemTop ();
  privateG2 = domain->privateG2;
  nfunc_a = privateG2->nfunc_a;
  coeff   = pkv_GetScratchMemd ( 3*nfunc_a );
  hessian = pkv_GetScratchMemd ( nfunc_a*(nfunc_a+1) );
  if ( !coeff || !hessian ) {
    domain->error_code = G2H_ERROR_NO_SCRATCH_MEMORY;
    goto failure;
  }
  dcoeff = &coeff[nfunc_a];
  grad   = &dcoeff[nfunc_a];
  chess  = &hessian[(nfunc_a*(nfunc_a+1))/2];

        /* setup the initial point */
  pkv_Selectd ( nfunc_a, 1, 3, 1, &nlprivate->acoeff[0].z, coeff );

        /* Newton iterations */
  for ( itn = ktn = 0; ; itn++ ) {
    if ( !_g2h_ComputeNLFuncGradHessiand ( domain, nlprivate, coeff,
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
      if ( gn < 0.0 || aux < EPSF*gn ) {
        domain->error_code = G2H_ERROR_NL_MINIMIZATION;
        goto failure;
      }
      pkn_MultMatrixNumd ( 1, nfunc_a, 0, grad, gn/aux, 0, grad );
    }
    dco = (double)sqrt ( pkn_ScalarProductd(nfunc_a, coeff, coeff) );
    dyn = (double)sqrt ( pkn_ScalarProductd(nfunc_a, grad, grad) );

    for ( aux = 1.0; aux >= EPSF; aux *= 0.5 ) {
      pkn_AddMatrixMd ( 1, nfunc_a, 0, coeff, 0, grad, aux, 0, dcoeff );
      func1 = _g2h_ComputeNLFuncd ( domain, nlprivate, dcoeff );
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
        _g2h_ComputeNLFuncGradd ( domain, nlprivate, coeff, &func0, grad );
        gn1 = (double)sqrt ( pkn_ScalarProductd ( nfunc_a, grad, grad ) );
        pkn_LowerTrMatrixSolved ( nfunc_a, hessian, 1, 1, grad, 1, grad );
        pkn_UpperTrMatrixSolved ( nfunc_a, hessian, 1, 1, grad, 1, grad );
        pkn_AddMatrixd ( 1, nfunc_a, 0, coeff, 0, grad, 0, dcoeff );
        func1 = _g2h_ComputeNLFuncd ( domain, nlprivate, dcoeff );
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

    if ( _g2h_StopItd ( itn, gn0, gn, dco, dyn, aux ) )
      break;
  }
/*
printf ( "itn = %d, ktn = %d\n", itn+1, ktn );
*/
  pkv_Selectd ( nfunc_a, 1, 1, 3, coeff, &nlprivate->acoeff[0].z );

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
#undef EPSF
} /*g2h_NLNewtond*/

boolean g2h_NLFillHoled ( GHoleDomaind *domain, const point3d *hole_cp,  
                          double *acoeff, void *usrptr,
                          void (*outpatch) ( int n, int m, const point3d *cp,
                                             void *usrptr ) )
{
  typedef void outscf ( int n, int m, const double *cp, void *usrptr );

  void     *sp;
  G2HNLPrivated     *nlprivate;
  G2HolePrivateRecd *privateG2;
  double   *fc00;
  int      hole_k, nfunc_a;

  sp = pkv_GetScratchMemTop ();
  if ( !g2h_DecomposeMatrixd ( domain ) )
    goto failure;

  privateG2 = domain->privateG2;
  hole_k = domain->hole_k;
  nfunc_a = privateG2->nfunc_a;
  _g2h_nlprivd = nlprivate = _g2h_InitNLprd ( hole_k, nfunc_a, 0 );
  if ( !nlprivate )
    goto failure;

  if ( !_g2h_ComputeNLNormald ( domain, nlprivate, hole_cp ) )
    goto failure;
  g2h_ReflectVectorsd ( 12*hole_k+1, hole_cp, nlprivate->rhole_cp );

  nlprivate->auxc = 0;
  if ( !g2h_FillHoled ( domain, 3, (double*)nlprivate->rhole_cp,
                        (double*)nlprivate->acoeff, NULL, g2h_nonlinoutpatchd ) )
    goto failure;

  if ( !_g2h_TabNLBasisFunctionsd ( domain, nlprivate ) )
    goto failure;

  if ( !g2h_NLNewtond ( domain, nlprivate ) )
    goto failure;

  g2h_ReflectVectorsd ( nfunc_a, nlprivate->acoeff, nlprivate->acoeff );
  if ( acoeff )
    memcpy ( acoeff, nlprivate->acoeff, nfunc_a*sizeof(vector3d) );

  fc00 = pkv_GetScratchMem ( (G2_CROSSDEGSUM+6)*2*hole_k*sizeof(vector3d) );
  if ( !fc00 ) {
    domain->error_code = G2H_ERROR_NO_SCRATCH_MEMORY;
    goto failure;
  }
  if ( !_g2h_SetRightSided ( domain, 3, (double*)hole_cp, fc00, NULL ) )
    goto failure;

  if ( !_g2h_OutputPatchesd ( domain, 3, (double*)nlprivate->acoeff, fc00,
                              usrptr, (outscf*)outpatch ) )
    goto failure;

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*g2h_NLFillHoled*/

/* ///////////////////////////////////////////////////////////////////////// */
static boolean g2h_NLConstrNewtond ( GHoleDomaind *domain,
                                     G2HNLPrivated *nlprivate, int nconstr,
                                     double *Cmat )
{
#define EPSF 2.0e-4
  void    *sp;
  G2HolePrivateRecd *privateG2;
  int     itn, jtn, ktn, i, j, hole_k, nfunc_a, nfunc;
  double  *coeff, *grad, *hessian;
  double  *cT, *E22, *cE22, *aa, *D1, *y, *y1, *M, *f, *f2;
  double  func, func0, func1, aux, t, gn, gn0 = 0.0, gn1, dco, dyn;
  boolean positive;

  sp = pkv_GetScratchMemTop ();
  hole_k  = domain->hole_k;
  privateG2 = domain->privateG2;
  nfunc_a = privateG2->nfunc_a;
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
       !D1 || !E22 || !cE22 || !y || !y1 || !M || !f ) {
    domain->error_code = G2H_ERROR_NO_SCRATCH_MEMORY;
    goto failure;
  }

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
    if ( !_g2h_ComputeNLFuncGradHessiand ( domain, nlprivate, coeff,
                                           &func, grad, hessian ) )
      goto failure;
    pkn_ComputeQTSQd ( nfunc_a, hessian, nconstr, cT, aa, M );
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
/*
printf ( "! " );
*/
      pkn_SymMatrixMultd ( nfunc, cE22, 1, 1, f2, 1, y1 );
      aux = (double)pkn_ScalarProductd ( nfunc, f2, y1 );
      t = (double)pkn_ScalarProductd ( nfunc, f2, f2 );
      if ( aux <= 0.0 || aux < EPSF*t ) {
        domain->error_code = G2H_ERROR_NL_MINIMIZATION;
        goto failure;
      }
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
      func1 = _g2h_ComputeNLFuncd ( domain, nlprivate, y1 );
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
        _g2h_ComputeNLFuncGradd ( domain, nlprivate, coeff, &func0, grad );
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

        func1 = _g2h_ComputeNLFuncd ( domain, nlprivate, y1 );
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
    if ( _g2h_StopItd ( itn, gn0, gn, dco, dyn, aux ) )
      break;
  }
/*
printf ( "itn = %d, ktn = %d\n", itn+1, ktn );
*/
  pkv_Selectd ( nfunc_a, 1, 1, 3, coeff, &nlprivate->acoeff[0].z );

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
#undef EPSF
} /*g2h_NLConstrNewtond*/

boolean g2h_NLFillHoleConstrd ( GHoleDomaind *domain, const point3d *hole_cp,
                    int nconstr, const vector3d *constr,
                    double *acoeff, void *usrptr,
                    void (*outpatch) ( int n, int m, const point3d *cp,
                                       void *usrptr ) )
{
  typedef void outscf ( int n, int m, const double *cp, void *usrptr );

  void     *sp;
  G2HNLPrivated     *nlprivate;
  G2HolePrivateRecd *privateG2;
  double   *fc00;
  int      hole_k, nfunc_a;

  sp = pkv_GetScratchMemTop ();
  if ( !g2h_DecomposeMatrixd ( domain ) )
    goto failure;

  privateG2 = domain->privateG2;
  hole_k  = domain->hole_k;
  nfunc_a = privateG2->nfunc_a;
  _g2h_nlprivd = nlprivate = _g2h_InitNLprd ( hole_k, nfunc_a, nconstr );
  if ( !nlprivate )
    goto failure;

  if ( !_g2h_ComputeNLNormald ( domain, nlprivate, hole_cp ) )
    goto failure;
  g2h_ReflectVectorsd ( 12*hole_k+1, hole_cp, nlprivate->rhole_cp );
  g2h_ReflectVectorsd ( nconstr, constr, &nlprivate->rhole_cp[12*hole_k+1] );

  nlprivate->auxc = 0;
  if ( !g2h_FillHoleConstrd ( domain, 3, (double*)nlprivate->rhole_cp,
                    nconstr, (double*)&nlprivate->rhole_cp[12*hole_k+1],
                    (double*)nlprivate->acoeff, NULL, g2h_nonlinoutpatchd ) )
    goto failure;

  if ( !_g2h_TabNLBasisFunctionsd ( domain, nlprivate ) )
    goto failure;

  if ( !g2h_NLConstrNewtond ( domain, nlprivate, nconstr, privateG2->Cmat ) )
    goto failure;

  g2h_ReflectVectorsd ( nfunc_a, nlprivate->acoeff, nlprivate->acoeff );
  if ( acoeff )
    memcpy ( acoeff, nlprivate->acoeff, nfunc_a*sizeof(vector3d) );

  fc00 = pkv_GetScratchMem ( (G2_CROSSDEGSUM+6)*2*hole_k*sizeof(vector3d) );
  if ( !fc00 ) {
    domain->error_code = G2H_ERROR_NO_SCRATCH_MEMORY;
    goto failure;
  }
  if ( !_g2h_SetRightSided ( domain, 3, (double*)hole_cp, fc00, NULL ) )
    goto failure;

  if ( !_g2h_OutputPatchesd ( domain, 3,
                              (double*)nlprivate->acoeff, fc00,
                              usrptr, (outscf*)outpatch ) )
    goto failure;

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*g2h_NLFillHoleConstrd*/

/* ///////////////////////////////////////////////////////////////////////// */
static boolean _g2h_ReflectAltConstrMatrixd ( GHoleDomaind *domain,
                                              vector3d *reflv,
                                              double *RACmat )
{
  void   *sp;
  G2HolePrivateRecd *privateG2;
  int    i, j, nconstr, nfunc_a;
  double *buf;

  sp = pkv_GetScratchMemTop ();
  privateG2 = domain->privateG2;
  nfunc_a = privateG2->nfunc_a;
  nconstr = privateG2->naconstr;
  memcpy ( RACmat, privateG2->ACmat, nconstr*nfunc_a*3*sizeof(double) );
  buf = pkv_GetScratchMemd ( nfunc_a );
  if ( !buf ) {
    domain->error_code = G2H_ERROR_NO_SCRATCH_MEMORY;
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
} /*_g2h_ReflectAltConstrMatrixd*/

boolean g2h_NLFillHoleAltConstrd ( GHoleDomaind *domain, const point3d *hole_cp,
                    int nconstr, const double *constr,
                    double *acoeff, void *usrptr,
                    void (*outpatch) ( int n, int m, const point3d *cp,
                                       void *usrptr ) )
{
  typedef void outscf ( int n, int m, const double *cp, void *usrptr );

  void     *sp;
  G2HNLPrivated     *nlprivate;
  G2HolePrivateRecd *privateG2;
  double   *fc00;
  double   *saveconstr = NULL, *newconstr, *Cmat, *rsconstr;
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
  saveconstr = pkv_GetScratchMemd ( 7*nfunc_a*nconstr );
  _g2h_nlprivd = nlprivate = _g2h_InitNLprd ( hole_k, nfunc_a, nconstr );
  if ( !saveconstr || !nlprivate )
    goto failure;

  newconstr = &saveconstr[3*nfunc_a*nconstr];
  Cmat = &newconstr[3*nfunc_a*nconstr];
  memcpy ( saveconstr, privateG2->ACmat, nfunc_a*nconstr*3*sizeof(double) );
  memcpy ( newconstr, privateG2->ACmat, nfunc_a*nconstr*3*sizeof(double) );
  if ( !_g2h_ComputeNLNormald ( domain, nlprivate, hole_cp ) )
    goto failure;
  g2h_ReflectVectorsd ( 12*hole_k+1, hole_cp, nlprivate->rhole_cp );
  if ( !_g2h_ReflectAltConstrMatrixd ( domain, &nlprivate->reflv, newconstr ) )
    goto failure;

  restore = true;
  if ( !g2h_SetAltConstraintMatrixd ( domain, 3, nconstr, newconstr ) )
    goto failure;
  nlprivate->auxc = 0;
  if ( !g2h_FillHoleAltConstrd ( domain, 3, (double*)nlprivate->rhole_cp,
                    nconstr, constr,
                    (double*)nlprivate->acoeff, NULL, g2h_nonlinoutpatchd ) )
    goto failure;

  pkv_Selectd ( nconstr, nfunc_a, 3*nfunc_a, nfunc_a,
                &newconstr[2*nfunc_a], Cmat );
  rsconstr = &nlprivate->rhole_cp[12*hole_k+1].z;
  pkv_Selectd ( nconstr, 1, 1, 3, constr, rsconstr );
  for ( j = 0; j < nconstr; j++ )
    for ( i = 0; i < nfunc_a; i++ )
      rsconstr[3*j] += newconstr[3*nfunc_a*j+i]*nlprivate->acoeff[i].x +
                       newconstr[(3*j+1)*nfunc_a+i]*nlprivate->acoeff[i].y;

  if ( !_g2h_TabNLBasisFunctionsd ( domain, nlprivate ) )
    goto failure;

  if ( !g2h_NLConstrNewtond ( domain, nlprivate, nconstr, Cmat ) )
    goto failure;

  g2h_ReflectVectorsd ( nfunc_a, nlprivate->acoeff, nlprivate->acoeff );
  if ( acoeff )
    memcpy ( acoeff, nlprivate->acoeff, nfunc_a*sizeof(vector3d) );

  fc00 = pkv_GetScratchMem ( (G2_CROSSDEGSUM+6)*2*hole_k*sizeof(vector3d) );
  if ( !fc00 ) {
    domain->error_code = G2H_ERROR_NO_SCRATCH_MEMORY;
    goto failure;
  }
  if ( !_g2h_SetRightSided ( domain, 3, (double*)hole_cp, fc00, NULL ) )
    goto failure;

  if ( !_g2h_OutputPatchesd ( domain, 3,
                              (double*)nlprivate->acoeff, fc00,
                              usrptr, (outscf*)outpatch ) )
    goto failure;

  g2h_SetAltConstraintMatrixd ( domain, 3, nconstr, saveconstr );
  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  if ( restore )
    g2h_SetAltConstraintMatrixd ( domain, 3, nconstr, saveconstr );
  pkv_SetScratchMemTop ( sp );
  return false;
} /*g2h_NLFillHoleAltConstrd*/

/* ///////////////////////////////////////////////////////////////////////// */
static double _g2h_ComputeLFuncd ( GHoleDomaind *domain,
                                  G2HNLPrivated *nlprivate, double *coeff )
{
  G2HolePrivateRecd *privateG2;
  int      k, i, fi, hole_k, nfunc_a;
  double   *psiuuu, *psiuuv, *psiuvv, *psivvv, *jac;
  double   funct;
  vector2d lapgrad;

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

  return (double)(funct/(4.0*(double)G2_NQUADSQ));
} /*_g2h_ComputeLFuncd*/

boolean g2h_NLFunctionalValued ( GHoleDomaind *domain,
                                 const point3d *hole_cp, const vector3d *acoeff,
                                 double *funcval )
{
  void   *sp;
  G2HNLPrivated     *nlprivate;
  G2HolePrivateRecd *privateG2;
  int    hole_k, nfunc_a;
  double *fc00, *coeff;

  sp = pkv_GetScratchMemTop ();
  privateG2 = domain->privateG2;
  hole_k = domain->hole_k;
  nfunc_a = privateG2->nfunc_a;
  _g2h_nlprivd = nlprivate = _g2h_InitNLprd ( hole_k, nfunc_a, 0 );
  if ( !nlprivate )
    goto failure;

  if ( !_g2h_ComputeNLNormald ( domain, nlprivate, hole_cp ) )
    goto failure;
  g2h_ReflectVectorsd ( 12*hole_k+1, hole_cp, nlprivate->rhole_cp );
  g2h_ReflectVectorsd ( nfunc_a, acoeff, nlprivate->acoeff );

  fc00 = pkv_GetScratchMem ( (G2_CROSSDEGSUM+6)*2*hole_k*sizeof(vector3d) );
  coeff = pkv_GetScratchMemd ( nfunc_a );
  if ( !fc00 || !coeff ) {
    domain->error_code = G2H_ERROR_NO_SCRATCH_MEMORY;
    goto failure;
  }
  if ( !_g2h_SetRightSided ( domain, 3, (double*)nlprivate->rhole_cp, fc00,
                             NULL ) )
    goto failure;
  nlprivate->auxc = 0;
  if ( !_g2h_OutputPatchesd ( domain, 3,
                (double*)nlprivate->acoeff, fc00, NULL, g2h_nonlinoutpatchd ) )
    goto failure;
  if ( !_g2h_TabNLBasisFunctionsd ( domain, nlprivate ) )
    goto failure;
  pkv_Selectd ( nfunc_a, 1, 3, 1, &nlprivate->acoeff[0].z, coeff );
  funcval[0] = _g2h_ComputeLFuncd ( domain, nlprivate, coeff );
  funcval[1] = _g2h_ComputeNLFuncd ( domain, nlprivate, coeff );

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*g2h_NLFunctionalValued*/

