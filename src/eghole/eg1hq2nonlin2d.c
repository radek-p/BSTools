
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

#include "eg1holed.h"
#include "eg2holed.h"
#include "eg1hprivated.h"
#include "eg2hprivated.h"
#include "eg1herror.h"


static G1HNLPrivated *_g1hq2_InitNLprd ( GHoleDomaind *domain,
                                         int hole_k, int V0SpaceDim,
                                         int nconstr )
{
#define N (G1H_FINALDEG+1)*(G1H_FINALDEG+1)
  G1HNLPrivated *nlpr;

  if ( (nlpr = pkv_GetScratchMem ( sizeof(G1HNLPrivated) )) ) {
    nlpr->auxc = 0;
    nlpr->nldi = pkv_GetScratchMem ( N*hole_k*sizeof(point3d) );
    nlpr->acoeff = pkv_GetScratchMem ( V0SpaceDim*sizeof(vector3d) );
    nlpr->rhole_cp = pkv_GetScratchMem ( (12*hole_k+1+nconstr)*sizeof(vector3d) );
    nlpr->diu = pkv_GetScratchMem ( 9*G1_NQUADSQ*hole_k*sizeof(vector2d) );
    nlpr->jac = pkv_GetScratchMemd ( G1_NQUADSQ*hole_k );
    nlpr->psiu = pkv_GetScratchMemd ( 9*(V0SpaceDim+1)*G1_NQUADSQ*hole_k );
    nlpr->ctang = pkv_GetScratchMem ( 3*G1_NQUAD*hole_k*sizeof(vector2d) );
    nlpr->cpsiu = pkv_GetScratchMemd ( 5*3*G1_NQUAD*(V0SpaceDim+1)*hole_k );

    if ( !nlpr->nldi || !nlpr->acoeff || !nlpr->rhole_cp ||
         !nlpr->diu || !nlpr->jac || !nlpr->psiu ||
         !nlpr->ctang || !nlpr->cpsiu ) {
      domain->error_code = G1H_ERROR_NO_SCRATCH_MEMORY;
      return NULL;
    }

    nlpr->div = &nlpr->diu[G1_NQUADSQ*hole_k];
    nlpr->diuu = &nlpr->div[G1_NQUADSQ*hole_k];
    nlpr->diuv = &nlpr->diuu[G1_NQUADSQ*hole_k];
    nlpr->divv = &nlpr->diuv[G1_NQUADSQ*hole_k];
    nlpr->diuuu = &nlpr->divv[G1_NQUADSQ*hole_k];
    nlpr->diuuv = &nlpr->diuuu[G1_NQUADSQ*hole_k];
    nlpr->diuvv = &nlpr->diuuv[G1_NQUADSQ*hole_k];
    nlpr->divvv = &nlpr->diuvv[G1_NQUADSQ*hole_k];
    nlpr->psiv = &nlpr->psiu[(V0SpaceDim+1)*G1_NQUADSQ*hole_k];
    nlpr->psiuu = &nlpr->psiv[(V0SpaceDim+1)*G1_NQUADSQ*hole_k];
    nlpr->psiuv = &nlpr->psiuu[(V0SpaceDim+1)*G1_NQUADSQ*hole_k];
    nlpr->psivv = &nlpr->psiuv[(V0SpaceDim+1)*G1_NQUADSQ*hole_k];
    nlpr->psiuuu = &nlpr->psivv[(V0SpaceDim+1)*G1_NQUADSQ*hole_k];
    nlpr->psiuuv = &nlpr->psiuuu[(V0SpaceDim+1)*G1_NQUADSQ*hole_k];
    nlpr->psiuvv = &nlpr->psiuuv[(V0SpaceDim+1)*G1_NQUADSQ*hole_k];
    nlpr->psivvv = &nlpr->psiuvv[(V0SpaceDim+1)*G1_NQUADSQ*hole_k];
    nlpr->cpsiv = &nlpr->cpsiu[3*G1_NQUAD*(V0SpaceDim+1)*hole_k];
    nlpr->cpsiuu = &nlpr->cpsiv[3*G1_NQUAD*(V0SpaceDim+1)*hole_k];
    nlpr->cpsiuv = &nlpr->cpsiuu[3*G1_NQUAD*(V0SpaceDim+1)*hole_k];
    nlpr->cpsivv = &nlpr->cpsiuv[3*G1_NQUAD*(V0SpaceDim+1)*hole_k];
  }
  return nlpr;
#undef N
} /*_g1hq2_InitNLprd*/

static double _g1hq2_ComputeNLFuncd ( GHoleDomaind *domain,
                                      G1HNLPrivated *nlprivate,
                                      const double *coeff,
                                      double C )
{
  G1HolePrivateRecd *privateG1;
  G2HNLFuncd   f;
  G1Q2HNLFuncd cf;
  int    k, i, fi, hole_k, nfunc_a;
  double c, funct, cfunct;

  privateG1 = domain->privateG1;
  hole_k = domain->hole_k;
  nfunc_a = privateG1->nfunc_a;

  funct = cfunct = 0.0;
        /* integration in the area Qmega */
  for ( k = 0; k < hole_k*G1_NQUADSQ; k++ ) {
        /* evaluate the function and its derivatives */
    fi = nfunc_a*G1_NQUADSQ*hole_k+k;
    f.pu = nlprivate->psiu[fi];      f.pv = nlprivate->psiv[fi];
    f.puu = nlprivate->psiuu[fi];    f.puv = nlprivate->psiuv[fi];
    f.pvv = nlprivate->psivv[fi];    f.puuu = nlprivate->psiuuu[fi];
    f.puuv = nlprivate->psiuuv[fi];  f.puvv = nlprivate->psiuvv[fi];
    f.pvvv = nlprivate->psivvv[fi];
    for ( i = 0; i < nfunc_a; i++ ) {
      fi = i*G1_NQUADSQ*hole_k+k;
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
  funct /= (double)G1_NQUADSQ;

        /* integration along the curves */
  for ( k = 0; k < hole_k*G1_NQUAD*3; k++ ) {
        /* evaluate the function and its derivatives */
    fi = nfunc_a*3*G1_NQUAD*hole_k+k;
    cf.pu = nlprivate->cpsiu[fi];     cf.pv = nlprivate->cpsiv[fi];
    cf.jpuu = nlprivate->cpsiuu[fi];  cf.jpuv = nlprivate->cpsiuv[fi];
    cf.jpvv = nlprivate->cpsivv[fi];
    for ( i = 0; i < nfunc_a; i++ ) {
      fi = i*3*G1_NQUAD*hole_k+k;
      c = coeff[i];
      cf.pu -= c*nlprivate->cpsiu[fi];     cf.pv -= c*nlprivate->cpsiv[fi];
      cf.jpuu -= c*nlprivate->cpsiuu[fi];  cf.jpuv -= c*nlprivate->cpsiuv[fi];
      cf.jpvv -= c*nlprivate->cpsivv[fi];
    }
    cf.tang = nlprivate->ctang[k];
        /* make the integration */
    _g1hq2_IntFunc1ad ( &cf, &cfunct );
  }
  cfunct /= (double)G1_NQUAD;

  return funct + C*cfunct;
} /*_g1hq2_ComputeNLFuncd*/

static boolean _g1hq2_ComputeNLFuncGradd ( GHoleDomaind *domain,
                                           G1HNLPrivated *nlprivate,
                                           const double *coeff,
                                           double C,
                                           double *func, double *grad )
{
  void   *sp;
  G1HolePrivateRecd *privateG1;
  G2HNLFuncd   f;
  G1Q2HNLFuncd cf;
  int    i, k, fi, hole_k, nfunc_a;
  double c, funct, cfunct, *cgrad;

  sp = pkv_GetScratchMemTop ();
  privateG1 = domain->privateG1;
  hole_k = domain->hole_k;
  nfunc_a = privateG1->nfunc_a;
  cgrad = pkv_GetScratchMemd ( nfunc_a );
  if ( !cgrad ) {
    domain->error_code = G1H_ERROR_NO_SCRATCH_MEMORY;
    goto failure;
  }

  funct = cfunct = 0.0;
  memset ( grad, 0, nfunc_a*sizeof(double) );
  memset ( cgrad, 0, nfunc_a*sizeof(double) );
        /* integration in the area Omega */
  for ( k = 0; k < hole_k*G1_NQUADSQ; k++ ) {
        /* evaluate the function and its derivatives */
    fi = nfunc_a*G1_NQUADSQ*hole_k+k;
    f.pu = nlprivate->psiu[fi];      f.pv = nlprivate->psiv[fi];
    f.puu = nlprivate->psiuu[fi];    f.puv = nlprivate->psiuv[fi];
    f.pvv = nlprivate->psivv[fi];    f.puuu = nlprivate->psiuuu[fi];
    f.puuv = nlprivate->psiuuv[fi];  f.puvv = nlprivate->psiuvv[fi];
    f.pvvv = nlprivate->psivvv[fi];
    for ( i = 0; i < nfunc_a; i++ ) {
      fi = i*G1_NQUADSQ*hole_k+k;
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
      fi = i*G1_NQUADSQ*hole_k+k;
      f.psiu = nlprivate->psiu[fi];      f.psiv = nlprivate->psiv[fi];
      f.psiuu = nlprivate->psiuu[fi];    f.psiuv = nlprivate->psiuv[fi];
      f.psivv = nlprivate->psivv[fi];    f.psiuuu = nlprivate->psiuuu[fi];
      f.psiuuv = nlprivate->psiuuv[fi];  f.psiuvv = nlprivate->psiuvv[fi];
      f.psivvv = nlprivate->psivvv[fi];
      _g2h_IntFunc2bd ( &f, &grad[i] );
    }
  }
  funct /= (double)G1_NQUADSQ;
  pkn_MultMatrixNumd ( 1, nfunc_a, 0, grad, 1.0/(double)G1_NQUADSQ, 0, grad );

        /* integration along the curves */
  for ( k = 0; k < hole_k*G1_NQUAD*3; k++ ) {
        /* evaluate the function and its derivatives */
    fi = nfunc_a*3*G1_NQUAD*hole_k+k;
    cf.pu = nlprivate->cpsiu[fi];     cf.pv = nlprivate->cpsiv[fi];
    cf.jpuu = nlprivate->cpsiuu[fi];  cf.jpuv = nlprivate->cpsiuv[fi];
    cf.jpvv = nlprivate->cpsivv[fi];
    for ( i = 0; i < nfunc_a; i++ ) {
      fi = i*3*G1_NQUAD*hole_k+k;
      c = coeff[i];
      cf.pu -= c*nlprivate->cpsiu[fi];     cf.pv -= c*nlprivate->cpsiv[fi];
      cf.jpuu -= c*nlprivate->cpsiuu[fi];  cf.jpuv -= c*nlprivate->cpsiuv[fi];
      cf.jpvv -= c*nlprivate->cpsivv[fi];
    }
    cf.tang = nlprivate->ctang[k];
        /* make the integration */
    _g1hq2_IntFunc1bd ( &cf, &cfunct );
        /* integrate the gradient */
    for ( i = 0; i < nfunc_a; i++ ) {
      fi = i*3*G1_NQUAD*hole_k+k;
      cf.psiu = nlprivate->cpsiu[fi];     cf.psiv = nlprivate->cpsiv[fi];
      cf.jpsiuu = nlprivate->cpsiuu[fi];  cf.jpsiuv = nlprivate->cpsiuv[fi];
      cf.jpsivv = nlprivate->cpsivv[fi];
      _g1hq2_IntFunc2bd ( &cf, &cgrad[i] );
    }
  }
  cfunct /= (double)G1_NQUAD;
  *func = funct + C*cfunct;
  pkn_AddMatrixMd ( 1, nfunc_a, 0, grad, 0, cgrad, C/(double)G1_NQUAD, 0, grad );

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*_g1hq2_ComputeNLFuncGradd*/

static boolean _g1hq2_ComputeNLFuncGradHessiand ( GHoleDomaind *domain,
                          G1HNLPrivated *nlprivate,
                          const double *coeff,
                          double C,
                          double *func, double *grad, double *hessian,
                          double *chessian )
{
  void   *sp;
  G1HolePrivateRecd *privateG1;
  G2HNLFuncd   f;
  G1Q2HNLFuncd cf;
  int    i, j, k, fi, fj, hole_k, nfunc_a;
  double *Li, *Bi, *BiLT, *Di;
  double c, funct, cfunct, *cgrad;

  sp = pkv_GetScratchMemTop ();
  privateG1 = domain->privateG1;
  hole_k = domain->hole_k;
  nfunc_a = privateG1->nfunc_a;
  cgrad = pkv_GetScratchMemd ( nfunc_a );
  Li = pkv_GetScratchMemd ( 8*nfunc_a );
  if ( !cgrad || !Li ) {
    domain->error_code = G1H_ERROR_NO_SCRATCH_MEMORY;
    goto failure;
  }
  Bi = &Li[2*nfunc_a];  BiLT = &Bi[3*nfunc_a];  Di = &BiLT[2*nfunc_a];

  funct = cfunct = 0.0;
  memset ( grad, 0, nfunc_a*sizeof(double) );
  memset ( cgrad, 0, nfunc_a*sizeof(double) );
  memset ( hessian, 0, (nfunc_a*(nfunc_a+1)/2)*sizeof(double) );
  memset ( chessian, 0, (nfunc_a*(nfunc_a+1)/2)*sizeof(double) );
        /* integration in the area Omega */
  for ( k = 0; k < hole_k*G1_NQUADSQ; k++ ) {
        /* evaluate the function and its derivatives */
    fi = nfunc_a*G1_NQUADSQ*hole_k+k;
    f.pu = nlprivate->psiu[fi];      f.pv = nlprivate->psiv[fi];
    f.puu = nlprivate->psiuu[fi];    f.puv = nlprivate->psiuv[fi];
    f.pvv = nlprivate->psivv[fi];    f.puuu = nlprivate->psiuuu[fi];
    f.puuv = nlprivate->psiuuv[fi];  f.puvv = nlprivate->psiuvv[fi];
    f.pvvv = nlprivate->psivvv[fi];
    for ( i = 0; i < nfunc_a; i++ ) {
      fi = i*G1_NQUADSQ*hole_k+k;
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
      fi = i*G1_NQUADSQ*hole_k+k;
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
      fi = i*G1_NQUADSQ*hole_k+k;
      f.psiu = nlprivate->psiu[fi];      f.psiv = nlprivate->psiv[fi];
      f.psiuu = nlprivate->psiuu[fi];    f.psiuv = nlprivate->psiuv[fi];
      f.psivv = nlprivate->psivv[fi];    f.psiuuu = nlprivate->psiuuu[fi];
      f.psiuuv = nlprivate->psiuuv[fi];  f.psiuvv = nlprivate->psiuvv[fi];
      f.psivvv = nlprivate->psivvv[fi];
      for ( j = 0; j <= i; j++ ) {
        fj = j*G1_NQUADSQ*hole_k+k;
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
  funct /= (double)G1_NQUADSQ;
  pkn_MultMatrixNumd ( 1, nfunc_a, 0, grad, 1.0/(double)G1_NQUADSQ, 0, grad );
  pkn_MultMatrixNumd ( 1, nfunc_a*(nfunc_a+1)/2, 0, hessian,
                       1.0/(double)G1_NQUADSQ, 0, hessian );

        /* integration along the curves */
  for ( k = 0; k < hole_k*G1_NQUAD*3; k++ ) {
        /* evaluate the function and its derivatives */
    fi = nfunc_a*3*G1_NQUAD*hole_k+k;
    cf.pu = nlprivate->cpsiu[fi];     cf.pv = nlprivate->cpsiv[fi];
    cf.jpuu = nlprivate->cpsiuu[fi];  cf.jpuv = nlprivate->cpsiuv[fi];
    cf.jpvv = nlprivate->cpsivv[fi];
    for ( i = 0; i < nfunc_a; i++ ) {
      fi = i*3*G1_NQUAD*hole_k+k;
      c = coeff[i];
      cf.pu -= c*nlprivate->cpsiu[fi];     cf.pv -= c*nlprivate->cpsiv[fi];
      cf.jpuu -= c*nlprivate->cpsiuu[fi];  cf.jpuv -= c*nlprivate->cpsiuv[fi];
      cf.jpvv -= c*nlprivate->cpsivv[fi];
    }
    cf.tang = nlprivate->ctang[k];
        /* make the integration */
    _g1hq2_IntFunc1cd ( &cf, &cfunct );
        /* integrate the gradient */
    for ( i = 0; i < nfunc_a; i++ ) {
      fi = i*3*G1_NQUAD*hole_k+k;
      cf.psiu = nlprivate->cpsiu[fi];     cf.psiv = nlprivate->cpsiv[fi];
      cf.jpsiuu = nlprivate->cpsiuu[fi];  cf.jpsiuv = nlprivate->cpsiuv[fi];
      cf.jpsivv = nlprivate->cpsivv[fi];
      _g1hq2_IntFunc2cd ( &cf, &Di[i], &Bi[i], &cgrad[i] );
    }
        /* integrate the Hessian */
    for ( i = 0; i < nfunc_a; i++ ) {
      fi = i*3*G1_NQUAD*hole_k+k;
      cf.psiu = nlprivate->cpsiu[fi];     cf.psiv = nlprivate->cpsiv[fi];
      cf.jpsiuu = nlprivate->cpsiuu[fi];  cf.jpsiuv = nlprivate->cpsiuv[fi];
      cf.jpsivv = nlprivate->cpsivv[fi];
      for ( j = 0; j <= i; j++ ) {
        fj = j*3*G1_NQUAD*hole_k+k;
        cf.psju = nlprivate->cpsiu[fj];     cf.psjv = nlprivate->cpsiv[fj];
        cf.jpsjuu = nlprivate->cpsiuu[fj];  cf.jpsjuv = nlprivate->cpsiuv[fj];
        cf.jpsjvv = nlprivate->cpsivv[fj];
        _g1hq2_IntFunc3cd ( &cf, Di[i], Bi[i], Di[j], Bi[j],
                            &chessian[pkn_SymMatIndex(i,j)] );
      }
    }
  }
  cfunct /= (double)G1_NQUAD;
/*
printf ( "f1 = %10.6f, f2 = %10.6f, C = %10.6f, f = %10.6f\n",
         funct, cfunct, C, funct+C*cfunct );
*/
  *func = funct + C*cfunct;
  pkn_AddMatrixMd ( 1, nfunc_a, 0, grad, 0, cgrad, C/(double)G1_NQUAD, 0, grad );
  pkn_AddMatrixMd ( 1, nfunc_a*(nfunc_a+1)/2, 0, hessian, 0, chessian,
                    C/(double)G1_NQUAD, 0, hessian );

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*_g1hq2_ComputeNLFuncGradHessiand*/

/* ///////////////////////////////////////////////////////////////////////// */
static boolean g1hq2_NLNewtond ( GHoleDomaind *domain,
                                 G1HNLPrivated *nlprivate )
{
#define EPSF 2.0e-4
  void    *sp;
  G1HolePrivateRecd *privateG1;
  int     itn, jtn, ktn, nfunc_a;
  double  *coeff, *dcoeff, *grad, *hessian, *chess;
  double  func, func0, func1, aux, gn, gn0 = 0.0, gn1, dco, dyn;
  double  C;
  boolean positive;

  sp = pkv_GetScratchMemTop ();
  privateG1 = domain->privateG1;
        /* the penalty constant C1b is set by the procedure */
        /* computing the initial approximation of the solution */
  C = 5.0*privateG1->C1b/nlprivate->ddiam;
  nfunc_a = privateG1->nfunc_a;
  coeff   = pkv_GetScratchMemd ( 3*nfunc_a );
  hessian = pkv_GetScratchMemd ( nfunc_a*(nfunc_a+1) );
  if ( !coeff || !hessian ) {
    domain->error_code = G1H_ERROR_NO_SCRATCH_MEMORY;
    goto failure;
  }
  dcoeff = &coeff[nfunc_a];
  grad   = &dcoeff[nfunc_a];
  chess  = &hessian[(nfunc_a*(nfunc_a+1))/2];

        /* setup the initial point */
  pkv_Selectd ( nfunc_a, 1, 3, 1, &nlprivate->acoeff[0].z, coeff );

        /* Newton iterations */
  func0 = 1.0e+8;
  for ( itn = ktn = 0; ; itn++ ) {
    if ( !_g1hq2_ComputeNLFuncGradHessiand ( domain, nlprivate, coeff, C,
                                             &func, grad, hessian, chess ) ) 
      goto failure;
    gn = (double)sqrt ( pkn_ScalarProductd ( nfunc_a, grad, grad ) );
    if ( itn == 0 ) {

printf ( "func = %f, gn0 = %f\n", func, gn );

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
      func1 = _g1hq2_ComputeNLFuncd ( domain, nlprivate, dcoeff, C );
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

printf ( "+" );

        _g1hq2_ComputeNLFuncGradd ( domain, nlprivate, coeff,
                                    C, &func0, grad );
        gn1 = (double)sqrt ( pkn_ScalarProductd ( nfunc_a, grad, grad ) );
        pkn_LowerTrMatrixSolved ( nfunc_a, hessian, 1, 1, grad, 1, grad );
        pkn_UpperTrMatrixSolved ( nfunc_a, hessian, 1, 1, grad, 1, grad );
        pkn_AddMatrixd ( 1, nfunc_a, 0, coeff, 0, grad, 0, dcoeff );
        func1 = _g1hq2_ComputeNLFuncd ( domain, nlprivate, dcoeff, C );
        if ( func1 >= func0 || gn1 >= 0.5*gn )
          break;

printf ( "    func = %f, gn = %f\n", func1, gn );

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

printf ( "func = %f\n", func1 );

  pkv_Selectd ( nfunc_a, 1, 1, 3, coeff, &nlprivate->acoeff[0].z );

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
#undef EPSF
} /*g1hq2_NLNewtond*/

boolean g1h_Q2NLFillHoled ( GHoleDomaind *domain,
                    const point3d *hole_cp,
                    double *acoeff, void *usrptr,
                    void (*outpatch) ( int n, int m, const point3d *cp,
                                       void *usrptr ) )
{
  typedef void outscf ( int n, int m, const double *cp, void *usrptr );

  void   *sp, *sp1;
  G1HNLPrivated     *nlprivate;
  G1HolePrivateRecd *privateG1;
  double *fc00, *ctrd, *bc00;
  int    hole_k, nfunc_a;

  sp = pkv_GetScratchMemTop ();
  privateG1 = domain->privateG1;
  hole_k = domain->hole_k;
  nfunc_a = privateG1->nfunc_a;
  _g1h_nlprivd = nlprivate = _g1hq2_InitNLprd ( domain, hole_k, nfunc_a, 0 );
  if ( !nlprivate )
    goto failure;

  if ( !_g1h_ComputeNLNormald ( domain, nlprivate, hole_cp ) )
    goto failure;
  g1h_ReflectVectorsd ( 12*hole_k+1, hole_cp, nlprivate->rhole_cp );

  nlprivate->auxc = 0;
  if ( !g1h_Q2FillHoled ( domain, 3, (double*)nlprivate->rhole_cp,
                          (double*)nlprivate->acoeff, NULL, g1h_nonlinoutpatchd ) )
    goto failure;

  sp1 = pkv_GetScratchMemTop ();
  if ( !(bc00 = pkv_GetScratchMemd ( 2*(G1_CROSSDEGSUM+4)*hole_k )) ) {
    domain->error_code = G1H_ERROR_NO_SCRATCH_MEMORY;
    goto failure;
  }
  if ( !_g1hq2_TabNLBasisFunctionsOmegad ( domain, G1_NQUAD, nlprivate, bc00 ) )
    goto failure;
  if ( !(ctrd = pkv_GetScratchMemd ( 3*G1_NQUAD*38*hole_k )) ) {
    domain->error_code = G1H_ERROR_NO_SCRATCH_MEMORY;
    goto failure;
  }
  if ( !_g1hq2_TabNLBasisFunctionsGammad ( domain, G1_NQUAD, nlprivate, ctrd, bc00 ) )
    goto failure;
  pkv_SetScratchMemTop ( sp1 );

  if ( !g1hq2_NLNewtond ( domain, nlprivate ) )
    goto failure;

  g1h_ReflectVectorsd ( nfunc_a, nlprivate->acoeff, nlprivate->acoeff );
  if ( acoeff )
    memcpy ( acoeff, nlprivate->acoeff, nfunc_a*sizeof(vector3d) );

  fc00 = pkv_GetScratchMem ( (G1_CROSSDEGSUM+4)*2*hole_k*sizeof(vector3d) );
  if ( !fc00 ) {
    domain->error_code = G1H_ERROR_NO_SCRATCH_MEMORY;
    goto failure;
  }
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
} /*g1h_Q2NLFillHoled*/

/* ///////////////////////////////////////////////////////////////////////// */
static boolean g1h_Q2NLConstrNewtond ( GHoleDomaind *domain,
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
  double  C;
  boolean positive;

  sp = pkv_GetScratchMemTop ();
  hole_k = domain->hole_k;
  privateG1 = domain->privateG1;
        /* the penalty constant C1b is set by the procedure */
        /* computing the initial approximation of the solution */
  C = 5.0*privateG1->C1b/nlprivate->ddiam;
  nfunc_a = privateG1->nfunc_a;
  nfunc   = nfunc_a-nconstr;
  coeff   = pkv_GetScratchMemd ( nfunc_a );
  grad    = pkv_GetScratchMemd ( nfunc_a );
  hessian = pkv_GetScratchMemd ( (nfunc_a*(nfunc_a+1))/2 );
  cT      = pkv_GetScratchMemd ( nfunc_a*nconstr );
  aa      = pkv_GetScratchMemd ( 2*nconstr );
  D1      = pkv_GetScratchMemd ( (nconstr*(nconstr+1))/2 );
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
    if ( !_g1hq2_ComputeNLFuncGradHessiand ( domain, nlprivate, coeff, C,
                                             &func, grad, hessian, M ) )
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
      func1 = _g1hq2_ComputeNLFuncd ( domain, nlprivate, y1, C );
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
        _g1hq2_ComputeNLFuncGradd ( domain, nlprivate, coeff, C, &func0, grad );
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

        func1 = _g1hq2_ComputeNLFuncd ( domain, nlprivate, y1, C );
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
} /*g1h_Q2NLConstrNewtond*/

boolean g1h_Q2NLFillHoleConstrd ( GHoleDomaind *domain, CONST_ point3d *hole_cp,
                                  int nconstr, CONST_ vector3d *constr,
                                  double *acoeff, void *usrptr,
                                  void (*outpatch) ( int n, int m,
                                                     const point3d *cp,
                                                     void *usrptr ) )
{
  typedef void outscf ( int n, int m, const double *cp, void *usrptr );

  void   *sp, *sp1;
  G1HNLPrivated     *nlprivate;
  G1HolePrivateRecd *privateG1;
  double *fc00, *ctrd, *bc00;
  int    hole_k, nfunc_a;

  sp = pkv_GetScratchMemTop ();
  privateG1 = domain->privateG1;
  hole_k = domain->hole_k;
  nfunc_a = privateG1->nfunc_a;
  _g1h_nlprivd = nlprivate = _g1hq2_InitNLprd ( domain, hole_k, nfunc_a, nconstr );
  if ( !nlprivate )
    goto failure;

  if ( !_g1h_ComputeNLNormald ( domain, nlprivate, hole_cp ) )
    goto failure;
  g1h_ReflectVectorsd ( 12*hole_k+1, hole_cp, nlprivate->rhole_cp );
  g1h_ReflectVectorsd ( nconstr, constr, &nlprivate->rhole_cp[12*hole_k+1] );

  nlprivate->auxc = 0;
  if ( !g1h_Q2FillHoleConstrd ( domain, 3, (double*)nlprivate->rhole_cp,
                    nconstr, (double*)&nlprivate->rhole_cp[12*hole_k+1],
                    (double*)nlprivate->acoeff, NULL, g1h_nonlinoutpatchd ) )
    goto failure;

  sp1 = pkv_GetScratchMemTop ();
  if ( !(bc00 = pkv_GetScratchMemd ( 2*(G1_CROSSDEGSUM+4)*hole_k )) ) {
    domain->error_code = G1H_ERROR_NO_SCRATCH_MEMORY;
    goto failure;
  }
  if ( !_g1hq2_TabNLBasisFunctionsOmegad ( domain, G1_NQUAD, nlprivate, bc00 ) )
    goto failure;
  if ( !(ctrd = pkv_GetScratchMemd ( 3*G1_NQUAD*38*hole_k )) ) {
    domain->error_code = G1H_ERROR_NO_SCRATCH_MEMORY;
    goto failure;
  }
  if ( !_g1hq2_TabNLBasisFunctionsGammad ( domain, G1_NQUAD, nlprivate, ctrd, bc00 ) )
    goto failure;
  pkv_SetScratchMemTop ( sp1 );

  if ( !g1h_Q2NLConstrNewtond ( domain, nlprivate, nconstr, privateG1->Cmat ) )
    goto failure;

  g1h_ReflectVectorsd ( nfunc_a, nlprivate->acoeff, nlprivate->acoeff );
  if ( acoeff )
    memcpy ( acoeff, nlprivate->acoeff, nfunc_a*sizeof(vector3d) );

  fc00 = pkv_GetScratchMem ( (G1_CROSSDEGSUM+4)*2*hole_k*sizeof(vector3d) );
  if ( !fc00 )
    goto failure;
  if ( !_g1h_SetRightSided ( domain, privateG1->Q2BMat,
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
} /*g1h_Q2NLFillHoleConstrd*/

/* ///////////////////////////////////////////////////////////////////////// */
boolean g1h_Q2NLFillHoleAltConstrd ( GHoleDomaind *domain, CONST_ point3d *hole_cp,
                                     int naconstr, CONST_ double *constr,
                                     double *acoeff, void *usrptr,
                                     void (*outpatch) ( int n, int m,
                                                        const point3d *cp,
                                                        void *usrptr ) )
{
  typedef void outscf ( int n, int m, const double *cp, void *usrptr );

  void     *sp, *sp1;
  G1HNLPrivated     *nlprivate;
  G1HolePrivateRecd *privateG1;
  double   *fc00, *ctrd, *bc00;
  double   *saveconstr = NULL, *newconstr, *Cmat, *rsconstr;
  int      hole_k, nfunc_a, i, j;
  boolean  restore;

  sp = pkv_GetScratchMemTop ();
  restore = false;
  privateG1 = domain->privateG1;
  if ( naconstr <= 0 || privateG1->acdim != 3 || privateG1->naconstr != naconstr ) {
    domain->error_code = G1H_ERROR_INCONSISTENT_CONSTR;
    goto failure;
  }
  hole_k  = domain->hole_k;
  nfunc_a = privateG1->nfunc_a;
  saveconstr = pkv_GetScratchMemd ( 7*nfunc_a*naconstr );
  _g1h_nlprivd = nlprivate = _g1hq2_InitNLprd ( domain, hole_k, nfunc_a, naconstr );
  if ( !saveconstr || !nlprivate )
    goto failure;

  newconstr = &saveconstr[3*nfunc_a*naconstr];
  Cmat = &newconstr[3*nfunc_a*naconstr];
  memcpy ( saveconstr, privateG1->ACmat, nfunc_a*naconstr*3*sizeof(double) );
  memcpy ( newconstr, privateG1->ACmat, nfunc_a*naconstr*3*sizeof(double) );
  if ( !_g1h_ComputeNLNormald ( domain, nlprivate, hole_cp ) )
    goto failure;
  g1h_ReflectVectorsd ( 12*hole_k+1, hole_cp, nlprivate->rhole_cp );
  if ( !_g1h_ReflectAltConstrMatrixd ( domain, &nlprivate->reflv, newconstr ) )
    goto failure;

  restore = true;
  if ( !g1h_SetAltConstraintMatrixd ( domain, 3, naconstr, newconstr ) )
    goto failure;
  nlprivate->auxc = 0;
  if ( !g1h_Q2FillHoleAltConstrd ( domain, 3, (double*)nlprivate->rhole_cp,
                    naconstr, constr,
                    (double*)nlprivate->acoeff, NULL, g1h_nonlinoutpatchd ) )
    goto failure;

  pkv_Selectd ( naconstr, nfunc_a, 3*nfunc_a, nfunc_a,
                &newconstr[2*nfunc_a], Cmat );
  rsconstr = &nlprivate->rhole_cp[12*hole_k+1].z;
  pkv_Selectd ( naconstr, 1, 1, 3, constr, rsconstr );
  for ( j = 0; j < naconstr; j++ )
    for ( i = 0; i < nfunc_a; i++ )
      rsconstr[3*j] += newconstr[3*nfunc_a*j+i]*nlprivate->acoeff[i].x +
                       newconstr[(3*j+1)*nfunc_a+i]*nlprivate->acoeff[i].y;

  sp1 = pkv_GetScratchMemTop ();
  if ( !(bc00 = pkv_GetScratchMemd ( 2*(G1_CROSSDEGSUM+4)*hole_k)) ) {
    domain->error_code = G1H_ERROR_NO_SCRATCH_MEMORY;
    goto failure;
  }
  if ( !_g1hq2_TabNLBasisFunctionsOmegad ( domain, G1_NQUAD, nlprivate, bc00 ) )
    goto failure;
  if ( !(ctrd = pkv_GetScratchMemd ( 3*G1_NQUAD*38*hole_k)) ) {
    domain->error_code = G1H_ERROR_NO_SCRATCH_MEMORY;
    goto failure;
  }
  if ( !_g1hq2_TabNLBasisFunctionsGammad ( domain, G1_NQUAD, nlprivate, ctrd, bc00 ) )
    goto failure;
  pkv_SetScratchMemTop ( sp1 );

  if ( !g1h_Q2NLConstrNewtond ( domain, nlprivate, naconstr, Cmat ) )
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

  g1h_SetAltConstraintMatrixd ( domain, 3, naconstr, saveconstr );
  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  if ( restore )
    g1h_SetAltConstraintMatrixd ( domain, 3, naconstr, saveconstr );
  pkv_SetScratchMemTop ( sp );
  return false;
} /*g1h_Q2NLFillHoleAltConstrd*/

