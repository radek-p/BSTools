
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

#include "eg2holef.h"
#include "eg2hprivatef.h"
#include "eg2herror.h"

static boolean g2h_ComputeBBMatrixf ( GHoleDomainf *domain )
{
  void     *sp;
  GHolePrivateRecf  *privateG;
  G2HolePrivateRecf *privateG2;
  int      hole_k, nfunc_b;
  int      i, j, l, n;
  float    *trd, *jac, *BBmat;
  float    *tkn, *hfunc, *dhfunc, *ddhfunc, *dddhfunc;
  vector2f *c00, *c01, *c02, *c10, *c11, *c12,
           *d00, *d01, *d02, *d10, *d11, *d12;
  float    *fc00, *fc01, *fc02, *fc10, *fc11, *fc12,
           *fd00, *fd01, *fd02, *fd10, *fd11, *fd12;
  vector2f *lgr;
  unsigned short *support_b;

  sp = pkv_GetScratchMemTop ();
  hole_k = domain->hole_k;
  privateG  = domain->privateG;
  privateG2 = domain->privateG2;
  if ( !privateG2 )
    goto failure;
  if ( !privateG2->Amat && !privateG2->EAmat )
    goto failure;

  if ( !privateG2->BBmat ) {
    nfunc_b   = privateG2->nfunc_b;
    support_b = privateG->support_b;
    BBmat = privateG2->BBmat = malloc ( (nfunc_b*(nfunc_b+1)/2)*sizeof(float) );
    if ( !BBmat )
      goto failure;

    n = hole_k*G2_NQUADSQ;

    tkn = pkv_GetScratchMemf ( 25*G2_NQUAD );
    jac = pkv_GetScratchMemf ( n );
    trd = pkv_GetScratchMemf ( 18*n );
    lgr = (vector2f*)pkv_GetScratchMem ( (nfunc_b+1)*n*sizeof(vector2f) );
    if ( !tkn || !jac || !trd || !lgr )
      goto failure;
    hfunc = &tkn[G2_NQUAD];         dhfunc = &hfunc[6*G2_NQUAD];
    ddhfunc = &dhfunc[6*G2_NQUAD];  dddhfunc = &ddhfunc[6*G2_NQUAD];
    _gh_PrepareTabKnotsf ( G2_NQUAD, privateG2->opt_quad, tkn );
    mbs_TabQuinticHFuncDer3f ( 0.0, 1.0, G2_NQUAD, tkn,
                               hfunc, dhfunc, ddhfunc, dddhfunc );

    for ( i = 0; i < hole_k; i++ ) {
      _g2h_GetDiPatchCurvesf ( domain, i,
                               &c00, &c01, &c02, &c10, &c11, &c12,
                               &d00, &d01, &d02, &d10, &d11, &d12 );
      _g2h_TabDiPatchJac3f ( G2_NQUAD, tkn, hfunc, dhfunc, ddhfunc, dddhfunc,
                             c00, c01, c02, c10, c11, c12,
                             d00, d01, d02, d10, d11, d12,
                             &jac[i*G2_NQUADSQ], &trd[i*G2_NQUADSQ*18] );
    }
    for ( j = 0; j < nfunc_b; j++ )
      for ( i = 0; i < hole_k; i++ )
        if ( support_b[j] & (0x0001 << i) ) {
          _g2h_GetBFBPatchCurvesf ( domain, j, i,
                                    &fc00, &fc01, &fc02, &fc10, &fc11, &fc12,
                                    &fd00, &fd01, &fd02, &fd10, &fd11, &fd12 );
          _g2h_TabLaplacianGradf ( G2_NQUAD, tkn, hfunc, dhfunc, ddhfunc, dddhfunc,
                                   fc00, fc01, fc02, fc10, fc11, fc12,
                                   fd00, fd01, fd02, fd10, fd11, fd12,
                                   &trd[i*G2_NQUADSQ*18],
                                   &lgr[(j*hole_k+i)*G2_NQUADSQ] );
        }
    for ( i = l = 0;  i < nfunc_b;  i++ ) {
      for ( j = 0;  j <= i;  j++, l++ )
        BBmat[l] = _g2h_Integralf ( hole_k, G2_NQUADSQ, jac,
                                    support_b[i], &lgr[i*n],
                                    support_b[j], &lgr[j*n] );
    }
  }

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*g2h_ComputeBBMatrixf*/

float g2h_FunctionalValuef ( GHoleDomainf *domain, int spdimen,
                             CONST_ float *hole_cp, CONST_ float *acoeff )
{
  void   *sp;
  GHolePrivateRecf  *privateG;
  G2HolePrivateRecf *privateG2;
  unsigned char *bfcpn;
  int    nfunc_a, nfunc_b, j;
  double fval;
  float  *Amat, *Bmat, *BBmat, *a, *b, *c;

  sp = pkv_GetScratchMemTop ();
  if ( !g2h_ComputeBBMatrixf ( domain ) )
    goto failure;

  privateG  = domain->privateG;
  privateG2 = domain->privateG2;
  nfunc_a = privateG2->nfunc_a;
  nfunc_b = privateG2->nfunc_b;
  bfcpn   = privateG->bfcpn;
  Amat    = privateG2->Amat;
  Bmat    = privateG2->Bmat;
  BBmat   = privateG2->BBmat;
  if ( !Amat || !Bmat )
    goto failure;

  a = pkv_GetScratchMemf ( spdimen*max(nfunc_a,nfunc_b) );
  b = pkv_GetScratchMemf ( spdimen*nfunc_a );
  c = pkv_GetScratchMemf ( spdimen*max(nfunc_a,nfunc_b) );
  if ( !a || !b || !c )
    goto failure;

  for ( j = 0; j < nfunc_b; j++ )
    memcpy ( &a[j*spdimen], &hole_cp[bfcpn[j]*spdimen], spdimen*sizeof(float) );
  pkn_SymMatrixMultf ( nfunc_b, BBmat, spdimen, spdimen, a, spdimen, c );
  fval = pkn_ScalarProductf ( spdimen*nfunc_b, a, c );

  pkn_SymMatrixMultf ( nfunc_a, Amat, spdimen, spdimen, acoeff, spdimen, c );
  pkn_MultMatrixf ( nfunc_a, nfunc_b, nfunc_b, Bmat,
                    spdimen, spdimen, a, spdimen, b );
  fval += pkn_ScalarProductf ( spdimen*nfunc_a, acoeff, c )
          -2.0*pkn_ScalarProductf ( spdimen*nfunc_a, acoeff, b );

  pkv_SetScratchMemTop ( sp );
  return (float)max ( 0.0, 0.25*fval );

failure:
  pkv_SetScratchMemTop ( sp );
  return -1.0;
} /*g2h_FunctionalValuef*/

/*
float _g2h_ExtFunctionalValuef ( GHoleDomainf *domain, int spdimen,
                                CONST_ float *hole_cp, CONST_ float *acoeff )
{
  void   *sp;
  G2HolePrivateRecf *privateG2;
  unsigned char *bfcpn;
  int    hole_k, nfunc_a, nfunc_b, nfunc_c, nbf, j;
  double fval;
  float  *Aii, *Aki, *Akk, *Bi, *Bk, *BBmat, *Lii, *Lki, *Lkk, *a, *b, *c;

  sp = pkv_GetScratchMemTop ();
  if ( !g2h_ComputeBBMatrixf ( domain ) )
    goto failure;

  privateG2 = domain->privateG2;
  if ( !_g2h_GetExtBlockAddressesf ( domain, &Aii, &Aki, &Akk,
                                     &Bi, &Bk, &Lii, &Lki, &Lkk ) )
    goto failure;

  hole_k  = domain->hole_k;
  nfunc_a = privateG2->nfunc_a;
  nfunc_b = privateG2->nfunc_b;
  nfunc_c = hole_k*G2_DBDIM;
  nbf     = nfunc_a+nfunc_c;
  bfcpn   = privateG2->bfcpn;
  BBmat   = privateG2->BBmat;

  a = pkv_GetScratchMemf ( spdimen*nbf );
  b = pkv_GetScratchMemf ( spdimen*nbf );
  c = pkv_GetScratchMemf ( spdimen*nbf );
  if ( !a || !b || !c )
    goto failure;

  for ( j = 0; j < nfunc_b; j++ )
    memcpy ( &a[j*spdimen], &hole_cp[bfcpn[j]*spdimen], spdimen*sizeof(float) );
  pkn_SymMatrixMultf ( nfunc_b, BBmat, spdimen, spdimen, a, spdimen, c );
  fval = pkn_ScalarProductf ( spdimen*nfunc_b, a, c );

  pkn_BSymMatrixMultf ( hole_k, G2_DBDIM, nfunc_a, Aii, Aki, Akk,
                        spdimen, spdimen, acoeff, spdimen, c );
  pkn_MultMatrixf ( nfunc_c, nfunc_b, nfunc_b, Bi,
                    spdimen, spdimen, a, spdimen, b );
  pkn_MultMatrixf ( nfunc_a, nfunc_b, nfunc_b, Bk,
                    spdimen, spdimen, a, spdimen, &b[nfunc_c*spdimen] );
  fval += pkn_ScalarProductf ( spdimen*nbf, acoeff, c )
            -2.0*pkn_ScalarProductf ( spdimen*nbf, acoeff, b );

  pkv_SetScratchMemTop ( sp );
  return max ( 0.0, fval );

failure:
  pkv_SetScratchMemTop ( sp );
  return -1.0;
}*/ /*g2h_ExtFunctionalValuef*/

/* The function above returns incorrect results. Below is a temporary */
/* replacement, to be used before the function above is fixed.  */

static int _np, _spdimen;
static float *patches;
static void myoutpatch ( int n, int m, const float *cp, void *usrptr )
{
  memcpy ( &patches[_np*_spdimen*(n+1)*(m+1)], cp,
           (n+1)*(m+1)*_spdimen*sizeof(float) );
  _np ++;
} /*myoutpatch*/

float g2h_ExtFunctionalValuef ( GHoleDomainf *domain, int spdimen,
                                 CONST_ float *hole_cp, CONST_ float *acoeff )
{
  void     *sp;
  G2HolePrivateRecf *privateG2;
  float   *fc00, *x, *Bi, *Bk, fct;
  float   *tkn, jac, lapg, lpgx, lpgy;
  vector2f *c00, *c01, *c02, *c10, *c11, *c12, *d00, *d01, *d02, *d10, *d11, *d12;
  point2f  *bdi;
  vector2f di, diu, div, diuu, diuv, divv, diuuu, diuuv, diuvv, divvv;
  float   *pi, *piu, *piv, *piuu, *piuv, *pivv, *piuuu, *piuuv, *piuvv, *pivvv;
  float   *fi, *fiu, *fiv, *fiuu, *fiuv, *fivv, *fiuuu, *fiuuv, *fiuvv, *fivvv;
  int      hole_k, nfunc_a, nfunc_b, nfunc_c, nbf;
  int      i, j, k, l, dn, dm;


  sp = pkv_GetScratchMemTop ();

  hole_k = domain->hole_k;
  privateG2 = domain->privateG2;
  nfunc_a = privateG2->nfunc_a;
  nfunc_b = privateG2->nfunc_b;
  nfunc_c = hole_k*G2_DBDIM;
  nbf     = nfunc_a+nfunc_c;

  if ( !privateG2->EAmat )
    if ( !g2h_DecomposeExtMatrixf ( domain ) )
      goto failure;

  x = pkv_GetScratchMemf ( spdimen*nbf );
  fc00 = pkv_GetScratchMemf ( (G2_CROSSDEGSUM+6)*2*hole_k*spdimen );
  patches = pkv_GetScratchMemf ( hole_k*spdimen*(G2H_FINALDEG+1)*(G2H_FINALDEG+1) );
  tkn = pkv_GetScratchMemf ( G2_NQUAD );
  bdi = pkv_GetScratchMem ( (G2H_FINALDEG+1)*(G2H_FINALDEG+1)*sizeof(point2f) );
  pi  = pkv_GetScratchMemf ( 20*spdimen );
  if ( !x || !fc00 || !patches|| !tkn || !bdi || !pi )
    goto failure;
  piu = &pi[spdimen];    piv = &piu[spdimen];
  piuu = &piv[spdimen];  piuv = &piuu[spdimen];  pivv = &piuv[spdimen];
  piuuu = &pivv[spdimen];  piuuv = &piuuu[spdimen];
  piuvv = &piuuv[spdimen];  pivvv = &piuvv[spdimen];
  fi = &pivvv[spdimen];   fiu = &fi[spdimen];     fiv = &fiu[spdimen];
  fiuu = &fiv[spdimen];  fiuv = &fiuu[spdimen];  fivv = &fiuv[spdimen];
  fiuuu = &fivv[spdimen];  fiuuv = &fiuuu[spdimen];
  fiuvv = &fiuuv[spdimen];  fivvv = &fiuvv[spdimen];

  Bi = privateG2->EBmat;
  Bk = &Bi[hole_k*G2_DBDIM*nfunc_b];
  if ( !_g2h_SetExtRightSidef ( domain, Bi, Bk, spdimen, hole_cp, fc00, x ) )
    goto failure;

  _spdimen = spdimen;
  _np = 0;
  if ( !_g2h_OutputExtPatchesf ( domain, spdimen, acoeff, fc00, NULL, myoutpatch ) )
    goto failure;

  _gh_PrepareTabKnotsf ( G2_NQUAD, privateG2->opt_quad, tkn );
  fct = 0.0;
  for ( k = 0; k < hole_k; k++ ) {
    _g2h_GetDiPatchCurvesf ( domain, k, &c00, &c01, &c02, &c10, &c11, &c12,
                             &d00, &d01, &d02, &d10, &d11, &d12 );
    if ( !mbs_BezC2CoonsToBezf ( 2,
               G2_CROSS00DEG, (float*)c00, G2_CROSS01DEG, (float*)c01, G2_CROSS02DEG, (float*)c02,
               G2_CROSS10DEG, (float*)c10, G2_CROSS11DEG, (float*)c11, G2_CROSS12DEG, (float*)c12,
               G2_CROSS00DEG, (float*)d00, G2_CROSS01DEG, (float*)d01, G2_CROSS02DEG, (float*)d02,
               G2_CROSS10DEG, (float*)d10, G2_CROSS11DEG, (float*)d11, G2_CROSS12DEG, (float*)d12,
               &dn, &dm, (float*)bdi ) )
      goto failure;
    for ( i = 0; i < G2_NQUAD; i++ )
      for ( j = 0; j < G2_NQUAD; j++ ) {
        mbs_BCHornerDer3Pf ( G2H_FINALDEG, G2H_FINALDEG, 2, (float*)bdi,
                             tkn[i], tkn[j],
                             &di.x, &diu.x, &div.x, &diuu.x, &diuv.x, &divv.x,
                             &diuuu.x, &diuuv.x, &diuvv.x, &divvv.x );
        mbs_BCHornerDer3Pf ( G2H_FINALDEG, G2H_FINALDEG, spdimen,
                             &patches[k*(G2H_FINALDEG+1)*(G2H_FINALDEG+1)*spdimen],
                             tkn[i], tkn[j],
                             pi, piu, piv, piuu, piuv, pivv,
                             piuuu, piuuv, piuvv, pivvv );
        if ( !pkn_Comp2iDerivatives3f ( diu.x, diu.y, div.x, div.y, diuu.x, diuu.y,
                                  diuv.x, diuv.y, divv.x, divv.y,
                                  diuuu.x, diuuu.y, diuuv.x, diuuv.y,
                                  diuvv.x, diuvv.y, divvv.x, divvv.y,
                                  spdimen, piu, piv, piuu, piuv, pivv,
                                  piuuu, piuuv, piuvv, pivvv,
                                  fiu, fiv, fiuu, fiuv, fivv,
                                  fiuuu, fiuuv, fiuvv, fivvv ) )
          goto failure;
        jac = diu.x*div.y - diu.y*div.x;
        for ( l = 0, lapg = 0.0;  l < spdimen;  l++ ) {
          lpgx = fiuuu[l]+fiuvv[l];
          lpgy = fiuuv[l]+fivvv[l];
          lapg += lpgx*lpgx + lpgy*lpgy;
        }
        fct += jac*lapg;
      }
  }

  pkv_SetScratchMemTop ( sp );
  return (float)(fct/(4.0*(float)G2_NQUADSQ));

failure:
  pkv_SetScratchMemTop ( sp );
  return -1.0;
} /*g2h_ExtFunctionalValuef*/

