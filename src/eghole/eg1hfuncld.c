
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

#undef CONST_
#define CONST_

#include "eg1holed.h"
#include "eg1hprivated.h"
#include "eg1herror.h"

static boolean g1h_ComputeBBMatrixd ( GHoleDomaind *domain )
{
  void     *sp;
  GHolePrivateRecd  *privateG;
  G1HolePrivateRecd *privateG1;
  int      hole_k, nfunc_b;
  int      i, j, l, n;
  double   *trd, *jac, *BBmat;
  double   *tkn, *hfunc, *dhfunc, *ddhfunc;
  vector2d *c00, *c01, *c10, *c11, *d00, *d01, *d10, *d11;
  double   *fc00, *fc01, *fc10, *fc11, *fd00, *fd01, *fd10, *fd11;
  double   *lap;
  unsigned short *support_b;

  sp = pkv_GetScratchMemTop ();
  hole_k = domain->hole_k;
  privateG  = domain->privateG;
  privateG1 = domain->privateG1;
  if ( !privateG1 )
    goto failure;
  if ( !privateG1->Amat && !privateG1->EAmat )
    goto failure;

  if ( !privateG1->BBmat ) {
    nfunc_b   = privateG1->nfunc_b;
    support_b = privateG->support_b;
    BBmat = privateG1->BBmat = malloc ( (nfunc_b*(nfunc_b+1)/2)*sizeof(double) );
    if ( !BBmat )
      goto failure;

    n = hole_k*G1_NQUADSQ;

    tkn = pkv_GetScratchMemd ( 19*G1_NQUAD );
    jac = pkv_GetScratchMemd ( n );
    trd = pkv_GetScratchMemd ( 5*n );
    lap = pkv_GetScratchMemd ( (nfunc_b+1)*n );
    if ( !tkn || !jac || !trd || !lap )
      goto failure;
    hfunc = &tkn[G1_NQUAD];         dhfunc = &hfunc[6*G1_NQUAD];
    ddhfunc = &dhfunc[6*G1_NQUAD];
    _gh_PrepareTabKnotsd ( G1_NQUAD, privateG1->opt_quad, tkn );
    if ( !mbs_TabCubicHFuncDer2d ( 0.0, 1.0, G1_NQUAD, tkn, hfunc, dhfunc, ddhfunc ) )
      goto failure;

    for ( i = 0; i < hole_k; i++ ) {
      _g1h_GetDiPatchCurvesd ( domain, i,
                               &c00, &c01, &c10, &c11,
                               &d00, &d01, &d10, &d11 );
      _g1h_TabDiPatchJac2d ( G1_NQUAD, tkn, hfunc, dhfunc, ddhfunc,
                             c00, c01, c10, c11, d00, d01, d10, d11,
                             &jac[i*G1_NQUADSQ], &trd[i*G1_NQUADSQ*5] );
    }
    for ( j = 0; j < nfunc_b; j++ )
      for ( i = 0; i < hole_k; i++ )
        if ( support_b[j] & (0x0001 << i) ) {
          _g1h_GetBFBPatchCurvesd ( domain, j, i,
                                    &fc00, &fc01, &fc10, &fc11,
                                    &fd00, &fd01, &fd10, &fd11 );
          _g1h_TabLaplaciand ( G1_NQUAD, tkn, hfunc, dhfunc, ddhfunc,
                               fc00, fc01, fc10, fc11, fd00, fd01, fd10, fd11,
                               &trd[i*G1_NQUADSQ*5], &lap[(j*hole_k+i)*G1_NQUADSQ] );
        }
    for ( i = l = 0;  i < nfunc_b;  i++ ) {
      for ( j = 0;  j <= i;  j++, l++ )
        BBmat[l] = _g1h_Integrald ( hole_k, G1_NQUADSQ, jac,
                                    support_b[i], &lap[i*n],
                                    support_b[j], &lap[j*n] );
    }
  }

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*g1h_ComputeBBMatrixd*/

double g1h_FunctionalValued ( GHoleDomaind *domain, int spdimen,
                             const double *hole_cp, const double *acoeff )
{
  void   *sp;
  GHolePrivateRecd  *privateG;
  G1HolePrivateRecd *privateG1;
  unsigned char *bfcpn;
  int    nfunc_a, nfunc_b, j;
  double fval;
  double *Amat, *Bmat, *BBmat, *a, *b, *c;

  sp = pkv_GetScratchMemTop ();
  if ( !g1h_ComputeBBMatrixd ( domain ) )
    goto failure;

  privateG  = domain->privateG;
  privateG1 = domain->privateG1;
  nfunc_a = privateG1->nfunc_a;
  nfunc_b = privateG1->nfunc_b;
  bfcpn   = privateG->bfcpn;
  Amat    = privateG1->Amat;
  Bmat    = privateG1->Bmat;
  BBmat   = privateG1->BBmat;
  if ( !Amat || !Bmat )
    goto failure;

  a = pkv_GetScratchMemd ( spdimen*max(nfunc_a,nfunc_b) );
  b = pkv_GetScratchMemd ( spdimen*nfunc_a );
  c = pkv_GetScratchMemd ( spdimen*max(nfunc_a,nfunc_b) );
  if ( !a || !b || !c )
    goto failure;

  for ( j = 0; j < nfunc_b; j++ )
    memcpy ( &a[j*spdimen], &hole_cp[bfcpn[j]*spdimen], spdimen*sizeof(double) );
  pkn_SymMatrixMultd ( nfunc_b, BBmat, spdimen, spdimen, a, spdimen, c );
  fval = pkn_ScalarProductd ( spdimen*nfunc_b, a, c );

  pkn_SymMatrixMultd ( nfunc_a, Amat, spdimen, spdimen, acoeff, spdimen, c );
  pkn_MultMatrixd ( nfunc_a, nfunc_b, nfunc_b, Bmat,
                    spdimen, spdimen, a, spdimen, b );
  fval += pkn_ScalarProductd ( spdimen*nfunc_a, acoeff, c )
           -2.0*pkn_ScalarProductd ( spdimen*nfunc_a, acoeff, b );

  pkv_SetScratchMemTop ( sp );
  return (double)max ( 0.0, 0.25*fval );

failure:
  pkv_SetScratchMemTop ( sp );
  return -1.0;
} /*g1h_FunctionalValued*/

/*
double _g1h_ExtFunctionalValued ( GHoleDomaind *domain, int spdimen,
                                CONST_ double *hole_cp, CONST_ double *acoeff )
{
  void   *sp;
  G1HolePrivateRecd *privateG1;
  unsigned char *bfcpn;
  int    hole_k, nfunc_a, nfunc_b, nfunc_c, nbf, j;
  double fval;
  double *Aii, *Aki, *Akk, *Bi, *Bk, *BBmat, *Lii, *Lki, *Lkk, *a, *b, *c;

  sp = pkv_GetScratchMemTop ();
  if ( !g1h_ComputeBBMatrixd ( domain ) )
    goto failure;

  privateG1 = domain->privateG1;
  if ( !_g1h_GetExtBlockAddressesd ( domain, &Aii, &Aki, &Akk,
                                     &Bi, &Bk, &Lii, &Lki, &Lkk ) )
    goto failure;

  hole_k  = domain->hole_k;
  nfunc_a = privateG1->nfunc_a;
  nfunc_b = privateG1->nfunc_b;
  nfunc_c = hole_k*G1_DBDIM;
  nbf     = nfunc_a+nfunc_c;
  bfcpn   = privateG1->bfcpn;
  BBmat   = privateG1->BBmat;

  a = pkv_GetScratchMemd ( spdimen*nbf );
  b = pkv_GetScratchMemd ( spdimen*nbf );
  c = pkv_GetScratchMemd ( spdimen*nbf );
  if ( !a || !b || !c )
    goto failure;

  for ( j = 0; j < nfunc_b; j++ )
    memcpy ( &a[j*spdimen], &hole_cp[bfcpn[j]*spdimen], spdimen*sizeof(double) );
  pkn_SymMatrixMultd ( nfunc_b, BBmat, spdimen, spdimen, a, spdimen, c );
  fval = pkn_ScalarProductd ( spdimen*nfunc_b, a, c );

  pkn_BSymMatrixMultd ( hole_k, G1_DBDIM, nfunc_a, Aii, Aki, Akk,
                        spdimen, spdimen, acoeff, spdimen, c );
  pkn_MultMatrixd ( nfunc_c, nfunc_b, nfunc_b, Bi,
                    spdimen, spdimen, a, spdimen, b );
  pkn_MultMatrixd ( nfunc_a, nfunc_b, nfunc_b, Bk,
                    spdimen, spdimen, a, spdimen, &b[nfunc_c*spdimen] );
  fval += pkn_ScalarProductd ( spdimen*nbf, acoeff, c )
            -2.0*pkn_ScalarProductd ( spdimen*nbf, acoeff, b );

  pkv_SetScratchMemTop ( sp );
  return max ( 0.0, fval );

failure:
  pkv_SetScratchMemTop ( sp );
  return -1.0;
}*/ /*g1h_ExtFunctionalValued*/

/* The function above returns incorrect results. Below is a temporary */
/* replacement, to be used before the function above is fixed.  */

static int _np, _spdimen;
static double *patches;
static void myoutpatch ( int n, int m, const double *cp, void *usrptr )
{
  memcpy ( &patches[_np*_spdimen*(n+1)*(m+1)], cp,
           (n+1)*(m+1)*_spdimen*sizeof(double) );
  _np ++;
} /*myoutpatch*/

double g1h_ExtFunctionalValued ( GHoleDomaind *domain, int spdimen,
                                 CONST_ double *hole_cp, CONST_ double *acoeff )
{
  void     *sp;
  G1HolePrivateRecd *privateG1;
  double   *fc00, *x, *Bi, *Bk, fct;
  double   *tkn, jac, lap, lp;
  vector2d *c00, *c01, *c10, *c11, *d00, *d01, *d10, *d11;
  point2d  *bdi;
  vector2d di, diu, div, diuu, diuv, divv;
  double   *pi, *piu, *piv, *piuu, *piuv, *pivv;
  double   *fi, *fiu, *fiv, *fiuu, *fiuv, *fivv;
  int      hole_k, nfunc_a, nfunc_b, nfunc_c, nbf;
  int      i, j, k, l, dn, dm;


  sp = pkv_GetScratchMemTop ();

  hole_k = domain->hole_k;
  privateG1 = domain->privateG1;
  nfunc_a = privateG1->nfunc_a;
  nfunc_b = privateG1->nfunc_b;
  nfunc_c = hole_k*G1_DBDIM;
  nbf     = nfunc_a+nfunc_c;

  if ( !privateG1->EAmat )
    if ( !g1h_DecomposeExtMatrixd ( domain ) )
      goto failure;

  x = pkv_GetScratchMemd ( spdimen*nbf );
  fc00 = pkv_GetScratchMemd ( (G1_CROSSDEGSUM+6)*2*hole_k*spdimen );
  patches = pkv_GetScratchMemd ( hole_k*spdimen*(G1H_FINALDEG+1)*(G1H_FINALDEG+1) );
  tkn = pkv_GetScratchMemd ( G1_NQUAD );
  bdi = pkv_GetScratchMem ( (G1H_FINALDEG+1)*(G1H_FINALDEG+1)*sizeof(point2d) );
  pi  = pkv_GetScratchMemd ( 12*spdimen );
  if ( !x || !fc00 || !patches|| !tkn || !bdi || !pi )
    goto failure;
  piu = &pi[spdimen];    piv = &piu[spdimen];
  piuu = &piv[spdimen];  piuv = &piuu[spdimen];  pivv = &piuv[spdimen];
  fi = &pivv[spdimen];   fiu = &fi[spdimen];     fiv = &fiu[spdimen];
  fiuu = &fiv[spdimen];  fiuv = &fiuu[spdimen];  fivv = &fiuv[spdimen];

  Bi = privateG1->EBmat;
  Bk = &Bi[hole_k*G1_DBDIM*nfunc_b];
  if ( !_g1h_SetExtRightSided ( domain, Bi, Bk, spdimen, hole_cp, fc00, x ) )
    goto failure;

  _spdimen = spdimen;
  _np = 0;
  if ( !_g1h_OutputExtPatchesd ( domain, spdimen, acoeff, fc00, NULL, myoutpatch ) )
    goto failure;

  _gh_PrepareTabKnotsd ( G1_NQUAD, privateG1->opt_quad, tkn );
  fct = 0.0;
  for ( k = 0; k < hole_k; k++ ) {
    _g1h_GetDiPatchCurvesd ( domain, k, &c00, &c01, &c10, &c11,
                             &d00, &d01, &d10, &d11 );
    if ( !mbs_BezC1CoonsToBezd ( 2,
               G1_CROSS00DEG, (double*)c00, G1_CROSS01DEG, (double*)c01,
               G1_CROSS10DEG, (double*)c10, G1_CROSS11DEG, (double*)c11,
               G1_CROSS00DEG, (double*)d00, G1_CROSS01DEG, (double*)d01,
               G1_CROSS10DEG, (double*)d10, G1_CROSS11DEG, (double*)d11,
               &dn, &dm, (double*)bdi ) )
      goto failure;
    for ( i = 0; i < G1_NQUAD; i++ )
      for ( j = 0; j < G1_NQUAD; j++ ) {
        if ( !mbs_BCHornerDer2Pd ( G1H_FINALDEG, G1H_FINALDEG, 2, (double*)bdi,
                             tkn[i], tkn[j],
                             &di.x, &diu.x, &div.x, &diuu.x, &diuv.x, &divv.x ) )
          goto failure;
        if ( !mbs_BCHornerDer2Pd ( G1H_FINALDEG, G1H_FINALDEG, spdimen,
                             &patches[k*(G1H_FINALDEG+1)*(G1H_FINALDEG+1)*spdimen],
                             tkn[i], tkn[j],
                             pi, piu, piv, piuu, piuv, pivv ) )
          goto failure;
        if ( !pkn_Comp2iDerivatives2d ( diu.x, diu.y, div.x, div.y, diuu.x, diuu.y,
                                        diuv.x, diuv.y, divv.x, divv.y,
                                        spdimen, piu, piv, piuu, piuv, pivv,
                                        fiu, fiv, fiuu, fiuv, fivv ) )
          goto failure;
        jac = diu.x*div.y - diu.y*div.x;
        for ( l = 0, lap = 0.0;  l < spdimen;  l++ ) {
          lp = fiuu[l]+fivv[l];
          lap += lp*lp;
        }
        fct += jac*lap;
      }
  }

  pkv_SetScratchMemTop ( sp );
  return (double)(fct/(4.0*(double)G1_NQUADSQ));

failure:
  pkv_SetScratchMemTop ( sp );
  return -1.0;
} /*g1h_ExtFunctionalValued*/

