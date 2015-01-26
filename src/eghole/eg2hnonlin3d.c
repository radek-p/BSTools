
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

#include "eg2holed.h"
#include "eg2hprivated.h"
#include "eg2herror.h"


static G2HNLPrivated *_g2h_InitExtNLprd ( int hole_k, int nfunc_a, int nconstr )
{
#define N (G2H_FINALDEG+1)*(G2H_FINALDEG+1)
  G2HNLPrivated *nlpr;
  int  nfunc_c;

  if ( (nlpr = pkv_GetScratchMem ( sizeof(G2HNLPrivated) )) ) {
    nfunc_c = hole_k*G2_DBDIM;
    nlpr->auxc = 0;
    nlpr->nldi = pkv_GetScratchMem ( N*hole_k*sizeof(point3d) );
    nlpr->acoeff = pkv_GetScratchMem ( (nfunc_a+nfunc_c)*sizeof(vector3d) );
    nlpr->rhole_cp = pkv_GetScratchMem ( (12*hole_k+1+nconstr)*sizeof(vector3d) );
    nlpr->diu  = pkv_GetScratchMem ( 9*G2_NQUADSQ*hole_k*sizeof(vector2d) );
    nlpr->jac = pkv_GetScratchMemd ( G2_NQUADSQ*hole_k );
    nlpr->psiu = pkv_GetScratchMemd ( 9*((nfunc_a+1)*hole_k+nfunc_c)*G2_NQUADSQ );

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
    nlpr->psiv = &nlpr->psiu[((nfunc_a+1)*hole_k+nfunc_c)*G2_NQUADSQ];
    nlpr->psiuu = &nlpr->psiv[((nfunc_a+1)*hole_k+nfunc_c)*G2_NQUADSQ];
    nlpr->psiuv = &nlpr->psiuu[((nfunc_a+1)*hole_k+nfunc_c)*G2_NQUADSQ];
    nlpr->psivv = &nlpr->psiuv[((nfunc_a+1)*hole_k+nfunc_c)*G2_NQUADSQ];
    nlpr->psiuuu = &nlpr->psivv[((nfunc_a+1)*hole_k+nfunc_c)*G2_NQUADSQ];
    nlpr->psiuuv = &nlpr->psiuuu[((nfunc_a+1)*hole_k+nfunc_c)*G2_NQUADSQ];
    nlpr->psiuvv = &nlpr->psiuuv[((nfunc_a+1)*hole_k+nfunc_c)*G2_NQUADSQ];
    nlpr->psivvv = &nlpr->psiuvv[((nfunc_a+1)*hole_k+nfunc_c)*G2_NQUADSQ];
  }
  return nlpr;
#undef N
} /*_g2h_InitExtNLprd*/

static boolean _g2h_TabExtNLBasisFunctionsd ( GHoleDomaind *domain,
                                              G2HNLPrivated *nlpr )
{
#define N ((G2H_FINALDEG+1)*(G2H_FINALDEG+1))
  void     *sp;
  G2HolePrivateRecd *privateG2;
  int      hole_k, nfunc_a, i, j, k, l, f, fN, bN;
  double   *tkn, *tbez, *tbezu, *tbezv, *tbezuu, *tbezuv, *tbezvv,
           *tbezuuu, *tbezuuv, *tbezuvv, *tbezvvv,
           *psiu, *psiv, *psiuu, *psiuv, *psivv,
           *psiuuu, *psiuuv, *psiuvv, *psivvv;
  vector2d *diu, *div, *diuu, *diuv, *divv, *diuuu, *diuuv, *diuvv, *divvv;

  sp      = pkv_GetScratchMemTop ();
  if ( !_g2h_TabNLBasisFunctionsd ( domain, nlpr ) )
    goto failure;

  hole_k  = domain->hole_k;
  privateG2 = domain->privateG2;
  nfunc_a = privateG2->nfunc_a;
  tkn     = pkv_GetScratchMemd ( G2_NQUAD );
  tbez    = pkv_GetScratchMemd ( 10*G2_NQUADSQ*G2_DBDIM );
  if ( !tkn || !tbez ) {
    domain->error_code = G2H_ERROR_NO_SCRATCH_MEMORY;
    goto failure;
  }
  tbezu = &tbez[G2_NQUADSQ*G2_DBDIM];       tbezv = &tbezu[G2_NQUADSQ*G2_DBDIM];
  tbezuu = &tbezv[G2_NQUADSQ*G2_DBDIM];     tbezuv = &tbezuu[G2_NQUADSQ*G2_DBDIM];
  tbezvv = &tbezuv[G2_NQUADSQ*G2_DBDIM];    tbezuuu = &tbezvv[G2_NQUADSQ*G2_DBDIM];
  tbezuuv = &tbezuuu[G2_NQUADSQ*G2_DBDIM];  tbezuvv = &tbezuuv[G2_NQUADSQ*G2_DBDIM];
  tbezvvv = &tbezuvv[G2_NQUADSQ*G2_DBDIM];
  _gh_PrepareTabKnotsd ( G2_NQUAD, privateG2->opt_quad, tkn );

  _g2h_TabTensBezPolyDer3d ( G2_NQUAD, tkn, tbez, tbezu, tbezv,
                             tbezuu, tbezuv, tbezvv,
                             tbezuuu, tbezuuv, tbezuvv, tbezvvv );
  for ( k = f = 0; k < hole_k; k++ ) {
    diu = &nlpr->diu[k*G2_NQUADSQ];      div = &nlpr->div[k*G2_NQUADSQ];
    diuu = &nlpr->diuu[k*G2_NQUADSQ];    diuv = &nlpr->diuv[k*G2_NQUADSQ];
    divv = &nlpr->divv[k*G2_NQUADSQ];    diuuu = &nlpr->diuuu[k*G2_NQUADSQ];
    diuuv = &nlpr->diuuv[k*G2_NQUADSQ];  diuvv = &nlpr->diuvv[k*G2_NQUADSQ];
    divvv = &nlpr->divvv[k*G2_NQUADSQ];
    for ( i = bN = 0;  i < G2H_FINALDEG-5;  i++ )
      for ( j = 0;  j < G2H_FINALDEG-5;  j++, f++, bN += G2_NQUADSQ ) {
        fN = ((nfunc_a+1)*hole_k+f)*G2_NQUADSQ;
        psiu = &nlpr->psiu[fN];      psiv = &nlpr->psiv[fN];
        psiuu = &nlpr->psiuu[fN];    psiuv = &nlpr->psiuv[fN];
        psivv = &nlpr->psivv[fN];
        psiuuu = &nlpr->psiuuu[fN];  psiuuv = &nlpr->psiuuv[fN];
        psiuvv = &nlpr->psiuvv[fN];  psivvv = &nlpr->psivvv[fN];
        for ( l = 0; l < G2_NQUADSQ; l++ )
          if ( !pkn_Comp2iDerivatives3d ( diu[l].x, diu[l].y, div[l].x, div[l].y,
                  diuu[l].x, diuu[l].y, diuv[l].x, diuv[l].y,
                  divv[l].x, divv[l].y, diuuu[l].x, diuuu[l].y,
                  diuuv[l].x, diuuv[l].y, diuvv[l].x, diuvv[l].y,
                  divvv[l].x, divvv[l].y, 1, &tbezu[bN+l], &tbezv[bN+l],
                  &tbezuu[bN+l], &tbezuv[bN+l], &tbezvv[bN+l],
                  &tbezuuu[bN+l], &tbezuuv[bN+l], &tbezuvv[bN+l], &tbezvvv[bN+l],
                  &psiu[l], &psiv[l], &psiuu[l], &psiuv[l], &psivv[l],
                  &psiuuu[l], &psiuuv[l], &psiuvv[l], &psivvv[l] ) )
            goto failure;
      }
  }
  pkv_SetScratchMemTop ( sp );
/*
DrawFGraph ( 1, G2_NQUAD, &nlpr->di[G2_NQUADSQ], &nlpr->psiuu[((nfunc_a+1)*hole_k+4)*G2_NQUADSQ] );
*/
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
#undef N
} /*_g2h_TabExtNLBasisFunctionsd*/

static double _g2h_ComputeExtNLFuncd ( GHoleDomaind *domain,
                                       G2HNLPrivated *nlprivate,
                                       const double *coeff )
{
#define SCALE (1.0/(double)G2_NQUADSQ)
  G2HolePrivateRecd *privateG2;
  G2HNLFuncd f;
  int    i, k, kn, knot, fi, hole_k, nfunc_a, nfunc_c;
  double c, funct;

  privateG2 = domain->privateG2;
  hole_k = domain->hole_k;
  nfunc_a = privateG2->nfunc_a;
  nfunc_c = hole_k*G2_DBDIM;

  funct = 0.0;
  for ( k = knot = 0;  k < hole_k;  k++ )
    for ( kn = 0;  kn < G2_NQUADSQ;  kn++, knot++ ) {
        /* evaluate the function and its derivatives */
      fi = nfunc_a*G2_NQUADSQ*hole_k+knot;
      f.pu = nlprivate->psiu[fi];      f.pv = nlprivate->psiv[fi];
      f.puu = nlprivate->psiuu[fi];    f.puv = nlprivate->psiuv[fi];
      f.pvv = nlprivate->psivv[fi];    f.puuu = nlprivate->psiuuu[fi];
      f.puuv = nlprivate->psiuuv[fi];  f.puvv = nlprivate->psiuvv[fi];
      f.pvvv = nlprivate->psivvv[fi];
      for ( i = 0; i < nfunc_a; i++ ) {
        fi = i*G2_NQUADSQ*hole_k+knot;
        c = coeff[nfunc_c+i];
        f.pu -= c*nlprivate->psiu[fi];      f.pv -= c*nlprivate->psiv[fi];
        f.puu -= c*nlprivate->psiuu[fi];    f.puv -= c*nlprivate->psiuv[fi];
        f.pvv -= c*nlprivate->psivv[fi];    f.puuu -= c*nlprivate->psiuuu[fi];
        f.puuv -= c*nlprivate->psiuuv[fi];  f.puvv -= c*nlprivate->psiuvv[fi];
        f.pvvv -= c*nlprivate->psivvv[fi];
      }
      for ( i = 0; i < G2_DBDIM; i++ ) {
        fi = ((nfunc_a+1)*hole_k+k*G2_DBDIM+i)*G2_NQUADSQ+kn;
        c = coeff[k*G2_DBDIM+i];
        f.pu -= c*nlprivate->psiu[fi];      f.pv -= c*nlprivate->psiv[fi];
        f.puu -= c*nlprivate->psiuu[fi];    f.puv -= c*nlprivate->psiuv[fi];
        f.pvv -= c*nlprivate->psivv[fi];    f.puuu -= c*nlprivate->psiuuu[fi];
        f.puuv -= c*nlprivate->psiuuv[fi];  f.puvv -= c*nlprivate->psiuvv[fi];
        f.pvvv -= c*nlprivate->psivvv[fi];
      }
      f.jac = nlprivate->jac[knot];
        /* make the integration */
      _g2h_IntFunc1ad ( &f, &funct );
    }
  return (double)(funct*SCALE);
#undef SCALE
} /*_g2h_ComputeExtNLFuncd*/

static boolean _g2h_ComputeExtNLFuncGradd ( GHoleDomaind *domain,
                                            G2HNLPrivated *nlprivate,
                                            const double *coeff,
                                            double *func, double *grad )
{
#define SCALE (1.0/(double)G2_NQUADSQ)
  void   *sp;
  G2HolePrivateRecd *privateG2;
  G2HNLFuncd f;
  int    i, k, kn, knot, fi, hole_k, nfunc_a, nfunc_c, nfunc;
  double c, funct;

  sp = pkv_GetScratchMemTop ();
  privateG2 = domain->privateG2;
  hole_k = domain->hole_k;
  nfunc_a = privateG2->nfunc_a;
  nfunc_c = hole_k*G2_DBDIM;
  nfunc = nfunc_a+nfunc_c;

  funct = 0.0;
  memset ( grad, 0, nfunc*sizeof(double) );
  for ( k = knot = 0;  k < hole_k;  k++ )
    for ( kn = 0;  kn < G2_NQUADSQ;  kn++, knot++ ) {
        /* evaluate the function and its derivatives */
      fi = nfunc_a*G2_NQUADSQ*hole_k+knot;
      f.pu = nlprivate->psiu[fi];      f.pv = nlprivate->psiv[fi];
      f.puu = nlprivate->psiuu[fi];    f.puv = nlprivate->psiuv[fi];
      f.pvv = nlprivate->psivv[fi];    f.puuu = nlprivate->psiuuu[fi];
      f.puuv = nlprivate->psiuuv[fi];  f.puvv = nlprivate->psiuvv[fi];
      f.pvvv = nlprivate->psivvv[fi];
      for ( i = 0; i < nfunc_a; i++ ) {
        fi = i*G2_NQUADSQ*hole_k+knot;
        c = coeff[nfunc_c+i];
        f.pu -= c*nlprivate->psiu[fi];      f.pv -= c*nlprivate->psiv[fi];
        f.puu -= c*nlprivate->psiuu[fi];    f.puv -= c*nlprivate->psiuv[fi];
        f.pvv -= c*nlprivate->psivv[fi];    f.puuu -= c*nlprivate->psiuuu[fi];
        f.puuv -= c*nlprivate->psiuuv[fi];  f.puvv -= c*nlprivate->psiuvv[fi];
        f.pvvv -= c*nlprivate->psivvv[fi];
      }
      for ( i = 0; i < G2_DBDIM; i++ ) {
        fi = ((nfunc_a+1)*hole_k+k*G2_DBDIM+i)*G2_NQUADSQ+kn;
        c = coeff[k*G2_DBDIM+i];
        f.pu -= c*nlprivate->psiu[fi];      f.pv -= c*nlprivate->psiv[fi];
        f.puu -= c*nlprivate->psiuu[fi];    f.puv -= c*nlprivate->psiuv[fi];
        f.pvv -= c*nlprivate->psivv[fi];    f.puuu -= c*nlprivate->psiuuu[fi];
        f.puuv -= c*nlprivate->psiuuv[fi];  f.puvv -= c*nlprivate->psiuvv[fi];
        f.pvvv -= c*nlprivate->psivvv[fi];
      }
      f.jac = nlprivate->jac[knot];

        /* make the integration */
      _g2h_IntFunc1bd ( &f, &funct );

        /* integrate the gradient */
      for ( i = 0; i < nfunc_a; i++ ) {
        fi = i*G2_NQUADSQ*hole_k+knot;
        f.psiu = nlprivate->psiu[fi];      f.psiv = nlprivate->psiv[fi];
        f.psiuu = nlprivate->psiuu[fi];    f.psiuv = nlprivate->psiuv[fi];
        f.psivv = nlprivate->psivv[fi];    f.psiuuu = nlprivate->psiuuu[fi];
        f.psiuuv = nlprivate->psiuuv[fi];  f.psiuvv = nlprivate->psiuvv[fi];
        f.psivvv = nlprivate->psivvv[fi];
        _g2h_IntFunc2bd ( &f, &grad[nfunc_c+i] );
      }
      for ( i = 0; i < G2_DBDIM; i++ ) {
        fi = ((nfunc_a+1)*hole_k+k*G2_DBDIM+i)*G2_NQUADSQ+kn;
        f.psiu = nlprivate->psiu[fi];      f.psiv = nlprivate->psiv[fi];
        f.psiuu = nlprivate->psiuu[fi];    f.psiuv = nlprivate->psiuv[fi];
        f.psivv = nlprivate->psivv[fi];    f.psiuuu = nlprivate->psiuuu[fi];
        f.psiuuv = nlprivate->psiuuv[fi];  f.psiuvv = nlprivate->psiuvv[fi];
        f.psivvv = nlprivate->psivvv[fi];
        _g2h_IntFunc2bd ( &f, &grad[k*G2_DBDIM+i] );
      }
    }

  *func = (double)(funct*SCALE);
  pkn_MultMatrixNumd ( 1, nfunc, 0, grad, SCALE, 0, grad );
  pkv_SetScratchMemTop ( sp );
  return true;

/*
failure:
  pkv_SetScratchMemTop ( sp );
  return false;
*/
#undef SCALE
} /*_g2h_ComputeExtNLFuncGradd*/

static boolean _g2h_ComputeExtNLFuncGradHessiand ( GHoleDomaind *domain,
                          G2HNLPrivated *nlprivate,
                          const double *coeff,
                          double *func, double *grad, double *hii )
{
#define SCALE (1.0/(double)G2_NQUADSQ)
  void   *sp;
  G2HolePrivateRecd *privateG2;
  G2HNLFuncd f;
  int    i, j, k, kn, knot, fi, fj, hole_k, nfunc_a, nfunc_c, nfunc;
  int    diagblsize, sideblsize, akksize;
  double c, funct;
  double *Li, *Bi, *BiLT, *Di, *hki, *hkk;

  sp = pkv_GetScratchMemTop ();
  privateG2 = domain->privateG2;
  hole_k = domain->hole_k;
  nfunc_a = privateG2->nfunc_a;
  nfunc_c = hole_k*G2_DBDIM;
  nfunc = nfunc_a+nfunc_c;
  hki = &hii[pkn_Block1FindBlockPos (hole_k, G2_DBDIM, nfunc_a, hole_k, 0)];
  hkk = &hii[pkn_Block1FindBlockPos (hole_k, G2_DBDIM, nfunc_a, hole_k, hole_k)];
  Li = pkv_GetScratchMemd ( 8*nfunc );
  if ( !Li ) {
    domain->error_code = G2H_ERROR_NO_SCRATCH_MEMORY;
    goto failure;
  }
  Bi = &Li[2*nfunc];  BiLT = &Bi[3*nfunc];  Di = &BiLT[2*nfunc];

  diagblsize = G2_DBDIM*(G2_DBDIM+1)/2;
  sideblsize = G2_DBDIM*nfunc_a;
  akksize = nfunc_a*(nfunc_a+1)/2;

  funct = 0.0;
  memset ( grad, 0, nfunc*sizeof(double) );
  memset ( hii, 0, hole_k*diagblsize*sizeof(double) );
  memset ( hki, 0, hole_k*sideblsize*sizeof(double) );
  memset ( hkk, 0, akksize*sizeof(double) );
  for ( k = knot = 0;  k < hole_k;  k++ )
    for ( kn = 0;  kn < G2_NQUADSQ;  kn++, knot++ ) {
        /* evaluate the function and its derivatives */
      fi = nfunc_a*G2_NQUADSQ*hole_k+knot;
      f.pu = nlprivate->psiu[fi];      f.pv = nlprivate->psiv[fi];
      f.puu = nlprivate->psiuu[fi];    f.puv = nlprivate->psiuv[fi];
      f.pvv = nlprivate->psivv[fi];    f.puuu = nlprivate->psiuuu[fi];
      f.puuv = nlprivate->psiuuv[fi];  f.puvv = nlprivate->psiuvv[fi];
      f.pvvv = nlprivate->psivvv[fi];
      for ( i = 0; i < nfunc_a; i++ ) {
        fi = i*G2_NQUADSQ*hole_k+knot;
        c = coeff[nfunc_c+i];
        f.pu -= c*nlprivate->psiu[fi];      f.pv -= c*nlprivate->psiv[fi];
        f.puu -= c*nlprivate->psiuu[fi];    f.puv -= c*nlprivate->psiuv[fi];
        f.pvv -= c*nlprivate->psivv[fi];    f.puuu -= c*nlprivate->psiuuu[fi];
        f.puuv -= c*nlprivate->psiuuv[fi];  f.puvv -= c*nlprivate->psiuvv[fi];
        f.pvvv -= c*nlprivate->psivvv[fi];
      }
      for ( i = 0; i < G2_DBDIM; i++ ) {
        fi = ((nfunc_a+1)*hole_k+k*G2_DBDIM+i)*G2_NQUADSQ+kn;
        c = coeff[k*G2_DBDIM+i];
        f.pu -= c*nlprivate->psiu[fi];      f.pv -= c*nlprivate->psiv[fi];
        f.puu -= c*nlprivate->psiuu[fi];    f.puv -= c*nlprivate->psiuv[fi];
        f.pvv -= c*nlprivate->psivv[fi];    f.puuu -= c*nlprivate->psiuuu[fi];
        f.puuv -= c*nlprivate->psiuuv[fi];  f.puvv -= c*nlprivate->psiuvv[fi];
        f.pvvv -= c*nlprivate->psivvv[fi];
      }
      f.jac = nlprivate->jac[knot];

        /* make the integration */
      _g2h_IntFunc1cd ( &f, &funct );

        /* integrate the gradient */
      for ( i = 0; i < nfunc_a; i++ ) {
        fi = i*G2_NQUADSQ*hole_k+knot;
        f.psiu = nlprivate->psiu[fi];      f.psiv = nlprivate->psiv[fi];
        f.psiuu = nlprivate->psiuu[fi];    f.psiuv = nlprivate->psiuv[fi];
        f.psivv = nlprivate->psivv[fi];    f.psiuuu = nlprivate->psiuuu[fi];
        f.psiuuv = nlprivate->psiuuv[fi];  f.psiuvv = nlprivate->psiuvv[fi];
        f.psivvv = nlprivate->psivvv[fi];
        _g2h_IntFunc2cd ( &f, (vector2d*)(&Li[2*(nfunc_c+i)]),
                          &Bi[3*(nfunc_c+i)], (vector2d*)(&BiLT[2*(nfunc_c+i)]),
                          &Di[nfunc_c+i], &grad[nfunc_c+i] );
      }
      for ( i = 0; i < G2_DBDIM; i++ ) {
        fi = ((nfunc_a+1)*hole_k+k*G2_DBDIM+i)*G2_NQUADSQ+kn;
        f.psiu = nlprivate->psiu[fi];      f.psiv = nlprivate->psiv[fi];
        f.psiuu = nlprivate->psiuu[fi];    f.psiuv = nlprivate->psiuv[fi];
        f.psivv = nlprivate->psivv[fi];    f.psiuuu = nlprivate->psiuuu[fi];
        f.psiuuv = nlprivate->psiuuv[fi];  f.psiuvv = nlprivate->psiuvv[fi];
        f.psivvv = nlprivate->psivvv[fi];
        _g2h_IntFunc2cd ( &f, (vector2d*)(&Li[2*(k*G2_DBDIM+i)]),
                          &Bi[3*(k*G2_DBDIM+i)], (vector2d*)(&BiLT[2*(k*G2_DBDIM+i)]),
                          &Di[k*G2_DBDIM+i], &grad[k*G2_DBDIM+i] );
      }

        /* integrate the Hessian */
      for ( i = 0; i < nfunc_a; i++ ) {
        fi = i*G2_NQUADSQ*hole_k+knot;
        f.psiu = nlprivate->psiu[fi];      f.psiv = nlprivate->psiv[fi];
        f.psiuu = nlprivate->psiuu[fi];    f.psiuv = nlprivate->psiuv[fi];
        f.psivv = nlprivate->psivv[fi];    f.psiuuu = nlprivate->psiuuu[fi];
        f.psiuuv = nlprivate->psiuuv[fi];  f.psiuvv = nlprivate->psiuvv[fi];
        f.psivvv = nlprivate->psivvv[fi];
        for ( j = 0; j <= i; j++ ) {
          fj = j*G2_NQUADSQ*hole_k+knot;
          f.psju = nlprivate->psiu[fj];      f.psjv = nlprivate->psiv[fj];
          f.psjuu = nlprivate->psiuu[fj];    f.psjuv = nlprivate->psiuv[fj];
          f.psjvv = nlprivate->psivv[fj];    f.psjuuu = nlprivate->psiuuu[fj];
          f.psjuuv = nlprivate->psiuuv[fj];  f.psjuvv = nlprivate->psiuvv[fj];
          f.psjvvv = nlprivate->psivvv[fj];
          _g2h_IntFunc3cd ( &f, (vector2d*)(&Li[2*(nfunc_c+i)]),
                                (vector2d*)(&Li[2*(nfunc_c+j)]),
                                (vector2d*)(&BiLT[2*(nfunc_c+i)]),
                                (vector2d*)(&BiLT[2*(nfunc_c+j)]),
                                Di[nfunc_c+i], Di[nfunc_c+j],
                                &hkk[pkn_SymMatIndex(i,j)] );
        }
        for ( j = 0; j < G2_DBDIM; j++ ) {
          fj = ((nfunc_a+1)*hole_k+k*G2_DBDIM+j)*G2_NQUADSQ+kn;
          f.psju = nlprivate->psiu[fj];      f.psjv = nlprivate->psiv[fj];
          f.psjuu = nlprivate->psiuu[fj];    f.psjuv = nlprivate->psiuv[fj];
          f.psjvv = nlprivate->psivv[fj];    f.psjuuu = nlprivate->psiuuu[fj];
          f.psjuuv = nlprivate->psiuuv[fj];  f.psjuvv = nlprivate->psiuvv[fj];
          f.psjvvv = nlprivate->psivvv[fj];
          _g2h_IntFunc3cd ( &f, (vector2d*)(&Li[2*(nfunc_c+i)]),
                                (vector2d*)(&Li[2*(k*G2_DBDIM+j)]),
                                (vector2d*)(&BiLT[2*(nfunc_c+i)]),
                                (vector2d*)(&BiLT[2*(k*G2_DBDIM+j)]),
                                Di[nfunc_c+i], Di[k*G2_DBDIM+j],
                                &hki[(k*nfunc_a+i)*G2_DBDIM+j] );
        }
      }
      for ( i = 0; i < G2_DBDIM; i++ ) {
        fi = ((nfunc_a+1)*hole_k+k*G2_DBDIM+i)*G2_NQUADSQ+kn;
        f.psiu = nlprivate->psiu[fi];      f.psiv = nlprivate->psiv[fi];
        f.psiuu = nlprivate->psiuu[fi];    f.psiuv = nlprivate->psiuv[fi];
        f.psivv = nlprivate->psivv[fi];    f.psiuuu = nlprivate->psiuuu[fi];
        f.psiuuv = nlprivate->psiuuv[fi];  f.psiuvv = nlprivate->psiuvv[fi];
        f.psivvv = nlprivate->psivvv[fi];
        for ( j = 0; j <= i; j++ ) {
          fj = ((nfunc_a+1)*hole_k+k*G2_DBDIM+j)*G2_NQUADSQ+kn;
          f.psju = nlprivate->psiu[fj];      f.psjv = nlprivate->psiv[fj];
          f.psjuu = nlprivate->psiuu[fj];    f.psjuv = nlprivate->psiuv[fj];
          f.psjvv = nlprivate->psivv[fj];    f.psjuuu = nlprivate->psiuuu[fj];
          f.psjuuv = nlprivate->psiuuv[fj];  f.psjuvv = nlprivate->psiuvv[fj];
          f.psjvvv = nlprivate->psivvv[fj];
          _g2h_IntFunc3cd ( &f, (vector2d*)(&Li[2*(k*G2_DBDIM+i)]),
                                (vector2d*)(&Li[2*(k*G2_DBDIM+j)]),
                                (vector2d*)(&BiLT[2*(k*G2_DBDIM+i)]),
                                (vector2d*)(&BiLT[2*(k*G2_DBDIM+j)]),
                                Di[k*G2_DBDIM+i], Di[k*G2_DBDIM+j],
                                &hii[k*diagblsize+pkn_SymMatIndex(i,j)] );
        }
      }
    }

  *func = (double)(funct*SCALE);
  pkn_MultMatrixNumd ( 1, nfunc, 0, grad, SCALE, 0, grad );
  pkn_MultMatrixNumd ( 1, akksize, 0, hkk, SCALE, 0, hkk );
  pkn_MultMatrixNumd ( 1, hole_k*sideblsize, 0, hki, SCALE, 0, hki );
  pkn_MultMatrixNumd ( 1, hole_k*diagblsize, 0, hii, SCALE, 0, hii );
  pkv_SetScratchMemTop ( sp );
/*
DrawExtMatrix ( hole_k, G2_DBDIM, nfunc_a, hii, hki, hkk );
*/
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
#undef SCALE
} /*_g2h_ComputeExtNLFuncGradHessiand*/

static boolean g2h_ExtNLNewtond ( GHoleDomaind *domain, G2HNLPrivated *nlprivate )
{
#define EPSF 2.0e-4
  void    *sp;
  G2HolePrivateRecd *privateG2;
  int     itn, jtn, ktn, hole_k, nfunc_a, nfunc_c, nfunc, asize;
  double  *coeff, *dcoeff, *grad, *hii, *chii;
  double  func, func0, func1, aux, gn, gn0 = 0.0, gn1, dco, dyn;
  boolean positive;

  sp = pkv_GetScratchMemTop ();
  hole_k  = domain->hole_k;
  privateG2 = domain->privateG2;
  nfunc_a = privateG2->nfunc_a;
  nfunc_c = hole_k*G2_DBDIM;
  nfunc   = nfunc_a+nfunc_c;

  coeff  = pkv_GetScratchMemd ( 3*nfunc );
  asize  = pkn_Block1ArraySize ( hole_k, G2_DBDIM, nfunc_a );
  hii    = pkv_GetScratchMemd ( 2*asize );
  if ( !coeff || !hii ) {
    domain->error_code = G2H_ERROR_NO_SCRATCH_MEMORY;
    goto failure;
  }
  dcoeff = &coeff[nfunc];
  grad   = &dcoeff[nfunc];
  chii   = &hii[asize];

        /* setup the initial point */
  pkv_Selectd ( nfunc, 1, 3, 1, &nlprivate->acoeff[0].z, coeff );

        /* Newton iterations */
  for ( itn = ktn = 0; ; itn++ ) {
    if ( !_g2h_ComputeExtNLFuncGradHessiand ( domain, nlprivate, coeff,
                             &func, grad, hii ) )
      goto failure;
    gn = (double)sqrt ( pkn_ScalarProductd ( nfunc, grad, grad ) );
    if ( itn == 0 ) {
/*
printf ( "func = %f, gn0 = %f\n", func, gn );
*/
      gn0 = gn;
    }
    memcpy ( chii, hii, asize );

    if ( (positive = pkn_Block1CholeskyDecompMd ( hole_k, G2_DBDIM, nfunc_a,
                                                  hii )) ) {
      pkn_Block1LowerTrMSolved ( hole_k, G2_DBDIM, nfunc_a, hii, 1, 1, grad );
      pkn_Block1UpperTrMSolved ( hole_k, G2_DBDIM, nfunc_a, hii, 1, 1, grad );
    }
    else {
/*
printf ( "! " );
*/
      if ( !pkn_Block1SymMatrixMultd ( hole_k, G2_DBDIM, nfunc_a, chii,
                                       1, 1, grad, 1, dcoeff ) )
        goto failure;
      aux = (double)pkn_ScalarProductd ( nfunc, grad, dcoeff );
      if ( gn < 0.0 || aux < EPSF*gn ) {
        domain->error_code = G2H_ERROR_NL_MINIMIZATION;
        goto failure;
      }
      pkn_MultMatrixNumd ( 1, nfunc, 0, grad, gn/aux, 0, grad );
    }
    dco = (double)sqrt ( pkn_ScalarProductd(nfunc, coeff, coeff) );
    dyn = (double)sqrt ( pkn_ScalarProductd(nfunc, grad, grad) );

    for ( aux = 1.0; aux > EPSF; aux *= 0.5 ) {
      pkn_AddMatrixMd ( 1, nfunc, 0, coeff, 0, grad, aux, 0, dcoeff );
      func1 = _g2h_ComputeExtNLFuncd ( domain, nlprivate, dcoeff );
      if ( func1 < func )
        break;
    }

    memcpy ( coeff, dcoeff, nfunc*sizeof(double) );
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
        _g2h_ComputeExtNLFuncGradd ( domain, nlprivate, coeff, &func0, grad );
        gn1 = (double)sqrt ( pkn_ScalarProductd ( nfunc, grad, grad ) );
        pkn_Block1LowerTrMSolved ( hole_k, G2_DBDIM, nfunc_a, hii, 1, 1, grad );
        pkn_Block1UpperTrMSolved ( hole_k, G2_DBDIM, nfunc_a, hii, 1, 1, grad );
        pkn_AddMatrixd ( 1, nfunc, 0, coeff, 0, grad, 0, dcoeff );
        func1 = _g2h_ComputeExtNLFuncd ( domain, nlprivate, dcoeff );
        if ( func1 >= func0 || gn1 >= 0.5*gn )
          break;
/*
printf ( "    func = %f, gn = %f\n", func1, gn );
*/
        memcpy ( coeff, dcoeff, nfunc*sizeof(double) );
        func = func1;
        ktn ++;
        gn = gn1;
        aux = 1.0;
        dco = (double)sqrt ( pkn_ScalarProductd(nfunc, coeff, coeff) );
        dyn = (double)sqrt ( pkn_ScalarProductd(nfunc, grad, grad) );
      }
    }

    if ( _g2h_StopItd ( itn, gn0, gn, dco, dyn, aux ) )
      break;
  }

printf ( "itn = %d, ktn = %d\n", itn+1, ktn );
printf ( "func = %f\n", func1 );

  pkv_Selectd ( nfunc, 1, 1, 3, coeff, &nlprivate->acoeff[0].z );

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
#undef EPSF
} /*g2h_ExtNLNewtond*/

/* ///////////////////////////////////////////////////////////////////////// */
boolean g2h_NLExtFillHoled ( GHoleDomaind *domain, const point3d *hole_cp,  
                             double *acoeff, void *usrptr,
                             void (*outpatch) ( int n, int m, const point3d *cp,
                                                void *usrptr ) )
{
  typedef void outscf ( int n, int m, const double *cp, void *usrptr );

  void     *sp;
  G2HNLPrivated     *nlprivate;
  G2HolePrivateRecd *privateG2;
  double   *fc00, *Bi, *Bk;
  int      hole_k, nfunc_a, nfunc_b, nfunc_c;

  sp = pkv_GetScratchMemTop ();
  if ( !g2h_DecomposeExtMatrixd ( domain ) )
    goto failure;

  privateG2 = domain->privateG2;
  hole_k = domain->hole_k;
  nfunc_a = privateG2->nfunc_a;
  nfunc_b = privateG2->nfunc_b;
  nfunc_c = hole_k*G2_DBDIM;
  _g2h_nlprivd = nlprivate = _g2h_InitExtNLprd ( hole_k, nfunc_a, 0 );
  if ( !nlprivate )
    goto failure;

  if ( !_g2h_ComputeNLNormald ( domain, nlprivate, hole_cp ) )
    goto failure;
  g2h_ReflectVectorsd ( 12*hole_k+1, hole_cp, nlprivate->rhole_cp );

  nlprivate->auxc = 0;
  if ( !g2h_ExtFillHoled ( domain, 3, (double*)nlprivate->rhole_cp,
                           (double*)nlprivate->acoeff, NULL, g2h_nonlinoutpatchd ) )
    goto failure;

  if ( !_g2h_TabExtNLBasisFunctionsd ( domain, nlprivate ) )
    goto failure;

  if ( !g2h_ExtNLNewtond ( domain, nlprivate ) )
    goto failure;

  g2h_ReflectVectorsd ( nfunc_a+nfunc_c, nlprivate->acoeff,
                        nlprivate->acoeff );
  if ( acoeff )
    memcpy ( acoeff, nlprivate->acoeff,
            (nfunc_a+nfunc_c)*sizeof(vector3d) );

  fc00 = pkv_GetScratchMem ( (G2_CROSSDEGSUM+6)*2*hole_k*sizeof(vector3d) );
  if ( !fc00 ) {
    domain->error_code = G2H_ERROR_NO_SCRATCH_MEMORY;
    goto failure;
  }
  Bi = privateG2->EBmat;
  Bk = &Bi[hole_k*G2_DBDIM*nfunc_b];
  if ( !_g2h_SetExtRightSided ( domain, Bi, Bk, 3, (double*)hole_cp, fc00, NULL ) )
    goto failure;

  if ( !_g2h_OutputExtPatchesd ( domain, 3,
                                (double*)nlprivate->acoeff, fc00,
                                usrptr, (outscf*)outpatch ) )
    goto failure;

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*g2h_NLExtFillHoled*/

/* ///////////////////////////////////////////////////////////////////////// */
static boolean g2h_ExtNLConstrNewtond ( GHoleDomaind *domain,
                                        G2HNLPrivated *nlprivate, int nconstr,
                                        double *ECmat )
{
#define EPSF 2.0e-4
  void    *sp;
  G2HolePrivateRecd *privateG2;
  int     itn, jtn, ktn, i, j, hole_k, nfunc_a, nfunc_c, nfunc_ac, nfunc;
  int     diagblsize, sideblsize, esideblsize, asize, esize;
  double  *coeff, *grad, *hii, *hkk, *hki;
  double  *cT, *E22kk, *E22ki, *E22ii, *cE22ii;
  double  *aa, *D1, *y, *y1, *M, *f;
  double  func, func0, func1, aux, gn, gn0 = 0.0, gn1, dco, dyn, dyn1;
  boolean positive;

  sp = pkv_GetScratchMemTop ();
  hole_k   = domain->hole_k;
  privateG2  = domain->privateG2;
  nfunc_a  = privateG2->nfunc_a;
  nfunc_c  = hole_k*G2_DBDIM;
  nfunc_ac = nfunc_a+nfunc_c;
  nfunc    = nfunc_ac-nconstr;
  diagblsize  = G2_DBDIM*(G2_DBDIM+1)/2;
  sideblsize  = G2_DBDIM*nfunc_a;
  esideblsize = G2_DBDIM*(nfunc_a-nconstr);

  coeff    = pkv_GetScratchMemd ( nfunc_ac );
  grad     = pkv_GetScratchMemd ( nfunc_ac );
  asize    = pkn_Block1ArraySize ( hole_k, G2_DBDIM, nfunc_a );
  hii      = pkv_GetScratchMemd ( asize );
  cT       = pkv_GetScratchMemd ( nfunc_a*nconstr );
  aa       = pkv_GetScratchMemd ( 2*nconstr );
  D1       = pkv_GetScratchMemd ( (nconstr*(nconstr+1))/2);
  esize    = pkn_Block1ArraySize ( hole_k, G2_DBDIM, nfunc_a-nconstr );
  E22ii    = pkv_GetScratchMemd ( esize );
  cE22ii   = pkv_GetScratchMemd ( esize );
  y        = pkv_GetScratchMemd ( nfunc_ac );
  y1       = pkv_GetScratchMemd ( nfunc_ac );
  M        = pkv_GetScratchMemd ( (nfunc_a*(nfunc_a+1))/2 );
  f        = pkv_GetScratchMemd ( nfunc_ac );
  if ( !coeff || !grad || !hii || !cT || !aa ||
       !D1 || !E22ii || !cE22ii || !y || !y1 || !M || !f ) {
    domain->error_code = G2H_ERROR_NO_SCRATCH_MEMORY;
    goto failure;
  }
  hkk = &hii[pkn_Block1FindBlockPos ( hole_k, G2_DBDIM, nfunc_a, hole_k, hole_k )];
  hki = &hii[pkn_Block1FindBlockPos ( hole_k, G2_DBDIM, nfunc_a, hole_k, 0 )];
  E22kk = &E22ii[pkn_Block1FindBlockPos ( hole_k, G2_DBDIM, nfunc_a-nconstr,
                                          hole_k, hole_k )];
  E22ki = &E22ii[pkn_Block1FindBlockPos ( hole_k, G2_DBDIM, nfunc_a-nconstr,
                                          hole_k, 0 )];

        /* setup the initial point */
  pkv_Selectd ( nfunc_ac, 1, 3, 1, &nlprivate->acoeff[0].z, coeff );

        /* step 1: decompose the constraint equations matrix */
  pkv_TransposeMatrixd ( nconstr, nfunc_a, nfunc_ac, &ECmat[nfunc_c],
                         nconstr, cT );
  pkn_QRDecomposeMatrixd ( nfunc_a, nconstr, cT, aa );
  for ( i = 0; i < nconstr; i++ )
    for ( j = i; j < nconstr; j++ )
      D1[pkn_SymMatIndex(i,j)] = cT[i*nconstr+j];

        /* Newton iterations */
  for ( itn = ktn = 0; ; itn++ ) {
printf ( "*" );
          /* step 2 */
    memcpy ( y, coeff, nfunc_ac*sizeof(double) );
    pkn_multiReflectVectord ( nfunc_a, nconstr, cT, aa, 1, 1, &y[nfunc_c] );
    pkn_LowerTrMatrixSolved ( nconstr, D1, 1,
                              3, &nlprivate->rhole_cp[12*hole_k+1].z,
                              1, &y[nfunc_c] );
    for ( i = nfunc_c; i < nfunc_c+nconstr; i++ )
      y[i] = -y[i];

          /* step 3 */
    if ( !_g2h_ComputeExtNLFuncGradHessiand ( domain, nlprivate, coeff,
                                              &func, grad, hii ) )
      goto failure;

    if ( !pkn_ComputeQTSQd ( nfunc_a, hkk, nconstr, cT, aa, M ) )
      goto failure;
    for ( i = 0; i < nfunc_a-nconstr; i++ )
      for ( j = i; j < nfunc_a-nconstr; j++ )
        E22kk[pkn_SymMatIndex(i,j)] = M[pkn_SymMatIndex(nconstr+i,nconstr+j)];
    for ( i = 0; i < hole_k; i++ )
      pkn_multiReflectVectord ( nfunc_a, nconstr, cT, aa,
                                G2_DBDIM, G2_DBDIM, &hki[i*sideblsize] );
    pkv_Selectd ( hole_k, esideblsize, sideblsize, esideblsize,
                  &hki[nconstr*G2_DBDIM], E22ki );
    memcpy ( E22ii, hii, hole_k*diagblsize*sizeof(double) );
    memcpy ( cE22ii, E22ii, esize*sizeof(double) );

          /* step 4 */
    memcpy ( f, grad, nfunc_ac*sizeof(double) );
    pkn_multiReflectVectord ( nfunc_a, nconstr, cT, aa, 1, 1, &f[nfunc_c] );
    memmove ( &f[nfunc_c], &f[nfunc_c+nconstr],
              (nfunc_a-nconstr)*sizeof(double) );

    gn = (double)sqrt ( pkn_ScalarProductd ( nfunc, f, f ) );
    if ( itn == 0 ) {
/*
printf ( "func = %f, gn0 = %f\n", func, gn );
*/
      gn0 = gn;
    }
          /* step 5 */
    if ( (positive = pkn_Block1CholeskyDecompMd ( hole_k, G2_DBDIM, nfunc_a-nconstr,
                                     E22ii )) ) {
      pkn_Block1LowerTrMSolved ( hole_k, G2_DBDIM, nfunc_a-nconstr, E22ii, 1, 1, f );
      pkn_Block1UpperTrMSolved ( hole_k, G2_DBDIM, nfunc_a-nconstr, E22ii, 1, 1, f );
    }
    else {
/*
printf ( "! " );
*/
      if ( !pkn_Block1SymMatrixMultd ( hole_k, G2_DBDIM, nfunc_a-nconstr, cE22ii,
                                       1, 1, f, 1, y1 ) )
        goto failure;
      aux = (double)pkn_ScalarProductd ( nfunc, f, y1 );
      if ( aux <= 0.0 || aux < EPSF*gn ) {
        domain->error_code = G2H_ERROR_NL_MINIMIZATION;
        goto failure;
      }
      pkn_MultMatrixNumd ( 1, nfunc, 1, f, gn/aux, 1, f );
    }

          /* step 6 */
    dyn = (double)sqrt ( pkn_ScalarProductd ( nfunc, f, f ) );
    memmove ( &f[nfunc_c+nconstr], &f[nfunc_c],
              (nfunc_a-nconstr)*sizeof(double) );
    dco = (double)sqrt ( pkn_ScalarProductd ( nfunc_ac, coeff, coeff ) );
    for ( aux = 1.0; aux >= EPSF; aux *= 0.5 ) {
      for ( i = 0; i < nfunc_c; i++ )
        y1[i] = y[i]+aux*f[i];
      memcpy ( &y1[nfunc_c], &y[nfunc_c], nconstr*sizeof(double) );
      for ( i = nfunc_c+nconstr; i < nfunc_ac; i++ )
        y1[i] = y[i]+aux*f[i];
      pkn_multiInvReflectVectord ( nfunc_a, nconstr, cT, aa, 1, 1, &y1[nfunc_c] );
      func1 = _g2h_ComputeExtNLFuncd ( domain, nlprivate, y1 );
      if ( func1 < func )
        break;
    }
    memcpy ( coeff, y1, nfunc_ac*sizeof(double) );
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
        memcpy ( y, coeff, nfunc_ac*sizeof(double) );
        pkn_multiReflectVectord ( nfunc_a, nconstr, cT, aa, 1, 1, &y[nfunc_c] );
        pkn_LowerTrMatrixSolved ( nconstr, D1, 1,
                                  3, &nlprivate->rhole_cp[12*hole_k+1].z,
                                  1, &y[nfunc_c] );
        for ( i = nfunc_c; i < nfunc_c+nconstr; i++ )
          y[i] = -y[i];
              /* step 3' */
        _g2h_ComputeExtNLFuncGradd ( domain, nlprivate, coeff, &func0, grad );
              /* step 4' */
        memcpy ( f, grad, nfunc_ac*sizeof(double) );
        pkn_multiReflectVectord ( nfunc_a, nconstr, cT, aa, 1, 1, &f[nfunc_c] );
        memmove ( &f[nfunc_c], &f[nfunc_c+nconstr],
                  (nfunc_a-nconstr)*sizeof(double) );
        gn1 = (double)sqrt ( pkn_ScalarProductd ( nfunc, f, f ) );
              /* step 5' */
        pkn_Block1LowerTrMSolved ( hole_k, G2_DBDIM, nfunc_a-nconstr, E22ii, 1, 1, f );
        pkn_Block1UpperTrMSolved ( hole_k, G2_DBDIM, nfunc_a-nconstr, E22ii, 1, 1, f );
              /* step 6' */
        dyn1 = (double)sqrt ( pkn_ScalarProductd(nfunc, f, f) );
        memmove ( &f[nfunc_c+nconstr], &f[nfunc_c],
                  (nfunc_a-nconstr)*sizeof(double) );
        for ( i = 0; i < nfunc_c; i++ )
          y1[i] = y[i]+f[i];
        memcpy ( &y1[nfunc_c], &y[nfunc_c], nconstr*sizeof(double) );
        for ( i = nfunc_c+nconstr; i < nfunc_ac; i++ )
          y1[i] = y[i]+f[i];
        pkn_multiInvReflectVectord ( nfunc_a, nconstr, cT, aa, 1, 1, &y1[nfunc_c] );

        func1 = _g2h_ComputeExtNLFuncd ( domain, nlprivate, y1 );
        if ( func1 >= func0 || gn1 >= 0.5*gn )
          break;
/*
printf ( "    func = %f, gn = %f\n", func1, gn );
*/
        memcpy ( coeff, y1, nfunc_ac*sizeof(double) );
        func = func1;
        ktn ++;
        gn = gn1;
        aux = 1.0;
        dco = (double)sqrt ( pkn_ScalarProductd(nfunc_ac, coeff, coeff) );
        dyn = dyn1;
      }
    }
    if ( _g2h_StopItd ( itn, gn0, gn, dco, dyn, aux ) )
      break;
  }
/*
printf ( "itn = %d, ktn = %d\n", itn+1, ktn );
*/
  pkv_Selectd ( nfunc_ac, 1, 1, 3, coeff, &nlprivate->acoeff[0].z );

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
#undef EPSF
} /*g2h_ExtNLConstrNewtond*/

boolean g2h_NLExtFillHoleConstrd ( GHoleDomaind *domain,
                     const point3d *hole_cp,
                     int nconstr, const vector3d *constr,
                     double *acoeff, void *usrptr,
                     void (*outpatch) ( int n, int m, const point3d *cp,
                                        void *usrptr ) )
{
  typedef void outscf ( int n, int m, const double *cp, void *usrptr );

  void     *sp;
  G2HNLPrivated     *nlprivate;
  G2HolePrivateRecd *privateG2;
  double   *fc00, *Bi, *Bk;
  int      hole_k, nfunc_a, nfunc_b, nfunc_c, nfunc_ac;

  sp = pkv_GetScratchMemTop ();
  if ( !g2h_DecomposeExtMatrixd ( domain ) )
    goto failure;

  privateG2 = domain->privateG2;
  hole_k  = domain->hole_k;
  nfunc_a = privateG2->nfunc_a;
  nfunc_b = privateG2->nfunc_b;
  nfunc_c = hole_k*G2_DBDIM;
  nfunc_ac = nfunc_a+nfunc_c;
  _g2h_nlprivd = nlprivate = _g2h_InitExtNLprd ( hole_k, nfunc_a, nconstr );
  if ( !nlprivate )
    goto failure;

  if ( !_g2h_ComputeNLNormald ( domain, nlprivate, hole_cp ) )
    goto failure;
  g2h_ReflectVectorsd ( 12*hole_k+1, hole_cp, nlprivate->rhole_cp );
  g2h_ReflectVectorsd ( nconstr, constr, &nlprivate->rhole_cp[12*hole_k+1] );

  nlprivate->auxc = 0;
  if ( !g2h_ExtFillHoleConstrd ( domain, 3, (double*)nlprivate->rhole_cp,
                    nconstr, (double*)&nlprivate->rhole_cp[12*hole_k+1],
                    (double*)nlprivate->acoeff, NULL, g2h_nonlinoutpatchd ) )
    goto failure;

  if ( !_g2h_TabExtNLBasisFunctionsd ( domain, nlprivate ) )
    goto failure;

  if ( !g2h_ExtNLConstrNewtond ( domain, nlprivate, nconstr, privateG2->ECmat ) )
    goto failure;

  g2h_ReflectVectorsd ( nfunc_ac, nlprivate->acoeff,
                        nlprivate->acoeff );
  if ( acoeff )
    memcpy ( acoeff, nlprivate->acoeff, nfunc_ac*sizeof(vector3d) );

  fc00 = pkv_GetScratchMem ( (G2_CROSSDEGSUM+6)*2*hole_k*sizeof(vector3d) );
  if ( !fc00 ) {
    domain->error_code = G2H_ERROR_NO_SCRATCH_MEMORY;
    goto failure;
  }
  Bi = privateG2->EBmat;
  Bk = &Bi[hole_k*G2_DBDIM*nfunc_b];
  if ( !_g2h_SetExtRightSided ( domain, Bi, Bk, 3, (double*)hole_cp, fc00, NULL ) )
    goto failure;

  if ( !_g2h_OutputExtPatchesd ( domain, 3,
                                 (double*)nlprivate->acoeff, fc00,
                                 usrptr, (outscf*)outpatch ) )
    goto failure;

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*g2h_NLExtFillHoleConstrd*/

/* ///////////////////////////////////////////////////////////////////////// */
static boolean _g2h_ReflectExtAltConstrMatrixd ( GHoleDomaind *domain,
                                                 vector3d *reflv,
                                                 double *RACmat )
{
  void   *sp;
  G2HolePrivateRecd *privateG2;
  int    i, j, nconstr, hole_k, nfunc_a, nfunc_c, nfunc_ac;
  double *buf;

  sp = pkv_GetScratchMemTop ();
  hole_k  = domain->hole_k;
  privateG2 = domain->privateG2;
  nfunc_a = privateG2->nfunc_a;
  nfunc_c = hole_k*G2_DBDIM;
  nfunc_ac = nfunc_a+nfunc_c;
  nconstr = privateG2->extnaconstr;
  memcpy ( RACmat, privateG2->AECmat, nconstr*nfunc_ac*3*sizeof(double) );
  buf = pkv_GetScratchMemd ( nfunc_ac );
  if ( !buf ) {
    domain->error_code = G2H_ERROR_NO_SCRATCH_MEMORY;
    goto failure;
  }
  for ( i = 0; i < nconstr; i++ ) {
    pkn_MultMatrixNumd ( 1, nfunc_ac, 0, RACmat, reflv->x, 0, buf );
    pkn_AddMatrixMd ( 1, nfunc_ac, 0, buf, 0, &RACmat[nfunc_ac], reflv->y, 0, buf );
    pkn_AddMatrixMd ( 1, nfunc_ac, 0, buf, 0, &RACmat[2*nfunc_ac], reflv->z, 0, buf );
    for ( j = 0; j < nfunc_ac; j++ ) {
      RACmat[j]            -= (double)(2.0*buf[j]*reflv->x);
      RACmat[nfunc_ac+j]   -= (double)(2.0*buf[j]*reflv->y);
      RACmat[2*nfunc_ac+j] -= (double)(2.0*buf[j]*reflv->z);
    }
    RACmat = &RACmat[3*nfunc_ac];
  }
  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*_g2h_ReflectExtAltConstrMatrixd*/

boolean g2h_NLExtFillHoleAltConstrd ( GHoleDomaind *domain, const point3d *hole_cp,
                         int naconstr, const double *constr,
                         double *acoeff, void *usrptr,
                         void (*outpatch) ( int n, int m, const point3d *cp,
                                            void *usrptr ) )
{
  typedef void outscf ( int n, int m, const double *cp, void *usrptr );

  void     *sp;
  G2HNLPrivated     *nlprivate;
  G2HolePrivateRecd *privateG2;
  double   *fc00, *Bi, *Bk;
  double   *saveconstr = NULL, *newconstr, *ECmat, *rsconstr;
  int      hole_k, nfunc_a, nfunc_b, nfunc_c, nfunc_ac, i, j;
  boolean  restore;

  sp = pkv_GetScratchMemTop ();
  restore = false;
  privateG2 = domain->privateG2;
  if ( naconstr <= 0 || privateG2->extacdim != 3 ||
       privateG2->extnaconstr != naconstr ) {
    domain->error_code = G2H_ERROR_INCONSISTENT_CONSTR;
    goto failure;
  }
  hole_k  = domain->hole_k;
  nfunc_a = privateG2->nfunc_a;
  nfunc_b = privateG2->nfunc_b;
  nfunc_c = hole_k*G2_DBDIM;
  nfunc_ac = nfunc_a+nfunc_c;
  saveconstr = pkv_GetScratchMemd ( 7*nfunc_ac*naconstr );
  _g2h_nlprivd = nlprivate = _g2h_InitExtNLprd ( hole_k, nfunc_a, naconstr );
  if ( !saveconstr || !nlprivate )
    goto failure;

  newconstr = &saveconstr[3*nfunc_ac*naconstr];
  ECmat = &newconstr[3*nfunc_ac*naconstr];
  memcpy ( saveconstr, privateG2->AECmat, nfunc_ac*naconstr*3*sizeof(double) );
  memcpy ( newconstr, privateG2->AECmat, nfunc_ac*naconstr*3*sizeof(double) );
  if ( !_g2h_ComputeNLNormald ( domain, nlprivate, hole_cp ) )
    goto failure;
  g2h_ReflectVectorsd ( 12*hole_k+1, hole_cp, nlprivate->rhole_cp );
  if ( !_g2h_ReflectExtAltConstrMatrixd ( domain, &nlprivate->reflv, newconstr ) )
    goto failure;

  restore = true;
  if ( !g2h_SetExtAltConstraintMatrixd ( domain, 3, naconstr, newconstr ) )
    goto failure;
  nlprivate->auxc = 0;
  if ( !g2h_ExtFillHoleAltConstrd ( domain, 3, (double*)nlprivate->rhole_cp,
                    naconstr, constr,
                    (double*)nlprivate->acoeff, NULL, g2h_nonlinoutpatchd ) )
    goto failure;

  pkv_Selectd ( naconstr, nfunc_ac, 3*nfunc_ac, nfunc_ac,
                &newconstr[2*nfunc_ac], ECmat );
  rsconstr = &nlprivate->rhole_cp[12*hole_k+1].z;
  pkv_Selectd ( naconstr, 1, 1, 3, constr, rsconstr );
  for ( j = 0; j < naconstr; j++ )
    for ( i = 0; i < nfunc_ac; i++ )
      rsconstr[3*j] += newconstr[3*nfunc_ac*j+i]*nlprivate->acoeff[i].x +
                       newconstr[(3*j+1)*nfunc_ac+i]*nlprivate->acoeff[i].y;

  if ( !_g2h_TabExtNLBasisFunctionsd ( domain, nlprivate ) )
    goto failure;

  if ( !g2h_ExtNLConstrNewtond ( domain, nlprivate, naconstr, ECmat ) )
    goto failure;

  g2h_ReflectVectorsd ( nfunc_ac, nlprivate->acoeff, nlprivate->acoeff );
  if ( acoeff )
    memcpy ( acoeff, nlprivate->acoeff, nfunc_ac*sizeof(vector3d) );

  fc00 = pkv_GetScratchMem ( (G2_CROSSDEGSUM+6)*2*hole_k*sizeof(vector3d) );
  if ( !fc00 ) {
    domain->error_code = G2H_ERROR_NO_SCRATCH_MEMORY;
    goto failure;
  }
  Bi = privateG2->EBmat;
  Bk = &Bi[hole_k*G2_DBDIM*nfunc_b];
  if ( !_g2h_SetExtRightSided ( domain, Bi, Bk, 3, (double*)hole_cp, fc00, NULL ) )
    goto failure;

  if ( !_g2h_OutputExtPatchesd ( domain, 3,
                                 (double*)nlprivate->acoeff, fc00,
                                 usrptr, (outscf*)outpatch ) )
    goto failure;

  g2h_SetExtAltConstraintMatrixd ( domain, 3, naconstr, saveconstr );
  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  if ( restore )
    g2h_SetExtAltConstraintMatrixd ( domain, 3, naconstr, saveconstr );
  pkv_SetScratchMemTop ( sp );
  return false;
} /*g2h_NLExtFillHoleAltConstrd*/

/* ///////////////////////////////////////////////////////////////////////// */
static double _g2h_ComputeExtLFuncd ( GHoleDomaind *domain,
                                      G2HNLPrivated *nlprivate, double *coeff )
{
  G2HolePrivateRecd *privateG2;
  int      k, i, kn, knot, fi, hole_k, nfunc_a, nfunc_c;
  double   *psiuuu, *psiuuv, *psiuvv, *psivvv, *jac;
  double   funct;
  vector2d lapgrad;

  privateG2 = domain->privateG2;
  hole_k  = domain->hole_k;
  nfunc_a = privateG2->nfunc_a;
  nfunc_c = hole_k*G2_DBDIM;
  psiuuu = nlprivate->psiuuu;
  psiuuv = nlprivate->psiuuv;
  psiuvv = nlprivate->psiuvv;
  psivvv = nlprivate->psivvv;
  jac   = nlprivate->jac;
  funct = 0.0;
  for ( k = knot = 0;  k < hole_k;  k++ )
    for ( kn = 0;  kn < G2_NQUADSQ;  kn++, knot++ ) {
      fi = nfunc_a*G2_NQUADSQ*hole_k+knot;
      lapgrad.x = psiuuu[fi]+psiuvv[fi];
      lapgrad.y = psiuuv[fi]+psivvv[fi];
      for ( i = 0; i < nfunc_a; i++ ) {
        fi = i*G2_NQUADSQ*hole_k+knot;
        lapgrad.x -= coeff[nfunc_c+i]*(psiuuu[fi]+psiuvv[fi]);
        lapgrad.y -= coeff[nfunc_c+i]*(psiuuv[fi]+psivvv[fi]);
      }
      for ( i = 0; i < G2_DBDIM; i++ ) {
        fi = ((nfunc_a+1)*hole_k+k*G2_DBDIM+i)*G2_NQUADSQ+kn;
        lapgrad.x -= coeff[k*G2_DBDIM+i]*(psiuuu[fi]+psiuvv[fi]);
        lapgrad.y -= coeff[k*G2_DBDIM+i]*(psiuuv[fi]+psivvv[fi]);
      }
      funct += (lapgrad.x*lapgrad.x+lapgrad.y*lapgrad.y)*jac[knot];
    }

  return (double)(funct/(4.0*(double)G2_NQUADSQ));
} /*_g2h_ComputeExtLFuncd*/

boolean g2h_NLExtFunctionalValued ( GHoleDomaind *domain,
                                    const point3d *hole_cp, const vector3d *acoeff,
                                    double *funcval )
{
  void   *sp;
  G2HNLPrivated     *nlprivate;
  G2HolePrivateRecd *privateG2;
  int    hole_k, nfunc_a, nfunc_b, nfunc_c;
  double *fc00, *coeff, *Bi, *Bk;

  sp = pkv_GetScratchMemTop ();
  privateG2 = domain->privateG2;
  hole_k = domain->hole_k;
  nfunc_a = privateG2->nfunc_a;
  nfunc_b = privateG2->nfunc_b;
  nfunc_c = hole_k*G2_DBDIM;
  _g2h_nlprivd = nlprivate = _g2h_InitExtNLprd ( hole_k, nfunc_a, 0 );
  if ( !nlprivate )
    goto failure;

  if ( !_g2h_ComputeNLNormald ( domain, nlprivate, hole_cp ) )
    goto failure;
  g2h_ReflectVectorsd ( 12*hole_k+1, hole_cp, nlprivate->rhole_cp );
  g2h_ReflectVectorsd ( nfunc_a+nfunc_c, acoeff, nlprivate->acoeff );

  fc00 = pkv_GetScratchMem ( (G2_CROSSDEGSUM+6)*2*hole_k*sizeof(vector3d) );
  coeff = pkv_GetScratchMemd ( nfunc_a+nfunc_c );
  if ( !fc00 || !coeff ) {
    domain->error_code = G2H_ERROR_NO_SCRATCH_MEMORY;
    goto failure;
  }
  Bi = privateG2->EBmat;
  Bk = &Bi[hole_k*G2_DBDIM*nfunc_b];
  if ( !_g2h_SetExtRightSided ( domain, Bi, Bk, 3, (double*)nlprivate->rhole_cp,
                                fc00, NULL ) )
    goto failure;
  nlprivate->auxc = 0;
  if ( !_g2h_OutputExtPatchesd ( domain, 3,
                (double*)nlprivate->acoeff, fc00, NULL, g2h_nonlinoutpatchd ) )
    goto failure;
  if ( !_g2h_TabExtNLBasisFunctionsd ( domain, nlprivate ) )
    goto failure;
  pkv_Selectd ( nfunc_a+nfunc_c, 1, 3, 1, &nlprivate->acoeff[0].z, coeff );
  funcval[0] = _g2h_ComputeExtLFuncd ( domain, nlprivate, coeff );
  funcval[1] = _g2h_ComputeExtNLFuncd ( domain, nlprivate, coeff );

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*g2h_NLExtFunctionalValued*/

