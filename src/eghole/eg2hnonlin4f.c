
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

#include "eg2holef.h"
#include "eg2hprivatef.h"
#include "eg2herror.h"

#define _DEBUG_HESSIAN

/* ///////////////////////////////////////////////////////////////////////// */
static void g2h_splnloutpatchf ( int n, int lknu, const float *knu,
                                 int m, int lknv, const float *knv,
                                 const float *cp, void *usrptr )
{
  int size;

  size = (lknu-n)*(lknv-m);
  memcpy ( &_g2h_nlprivf->nldi[_g2h_nlprivf->auxc*size], cp,
           size*sizeof(vector3f) );
  _g2h_nlprivf->auxc ++;
} /*g2h_splnloutpatchf*/

/* ///////////////////////////////////////////////////////////////////////// */
static G2HNLSPrivatef *_g2h_AllocSplNLPrf ( GHoleDomainf *domain )
{
  G2HolePrivateRecf  *privateG2;
  G2HoleSPrivateRecf *sprivateG2;
  G2HNLSPrivatef     *nlspr;
  int   hole_k, nk, m1, m2, rr, nkn, psize;
  int   nfunc_a, nfunc_c, nfunc_d;
  int   *fkn, *lkn, *cfuncvi, *dfuncvi;
  float *tkn;
  int   i, j, k, fn, nr, nc;
  int   i0, j0, i1, j1, nzc;

  if ( (nlspr = pkv_GetScratchMem ( sizeof(G2HNLSPrivatef) )) ) {
    privateG2 = domain->privateG2;
    sprivateG2 = domain->SprivateG2;
    hole_k = domain->hole_k;
    nk  = sprivateG2->nk;
    m1  = sprivateG2->m1;
    m2  = sprivateG2->m2;
    nfunc_a = privateG2->nfunc_a;
    nfunc_c = sprivateG2->nsfunc_c;
    nfunc_d = sprivateG2->nsfunc_d;
    rr  = G2H_FINALDEG-5+nk*m2;
    nkn = nlspr->nkn = (nk+1)*G2_QUAD_FACTOR;
    tkn = nlspr->tkn = pkv_GetScratchMemf ( nkn );
    fkn = nlspr->fkn = pkv_GetScratchMem ( 2*rr*sizeof(int) );
    nlspr->cb = pkv_GetScratchMemf ( 4*rr*nkn );
    cfuncvi = nlspr->cfuncvi = pkv_GetScratchMem ( nfunc_c*sizeof(int) );
    dfuncvi = nlspr->dfuncvi = pkv_GetScratchMem ( nfunc_d*sizeof(int) );
    if ( !tkn || !fkn || !nlspr->cb || !cfuncvi || !dfuncvi )
      return false;
    lkn = nlspr->lkn = &fkn[rr];
    nlspr->cbt = &nlspr->cb[rr*nkn];
    nlspr->cbtt = &nlspr->cbt[rr*nkn];
    nlspr->cbttt = &nlspr->cbtt[rr*nkn];
    _gh_PrepareTabKnotsf ( nkn, privateG2->opt_quad, tkn );
    psize = sprivateG2->lastfpknot-G2H_FINALDEG;
    nlspr->psize = psize*psize;

        /* compute the numbers of quadrature knots for each C block function */
    _g2h_TabBSFuncDer3f ( G2H_FINALDEG,
                  sprivateG2->lastcknot, sprivateG2->cknots,
                  3, G2H_FINALDEG+nk*m2-3, nkn, tkn, nlspr->fkn, nlspr->lkn,
                  nlspr->cb, nlspr->cbt, nlspr->cbtt, nlspr->cbttt );
    for ( j = fn = 0;  j < rr;  j++ ) {
      nr = lkn[j]-fkn[j]+1;
      for ( k = 0;  k < rr;  k++, fn++ ) {
        nc = lkn[k]-fkn[k]+1;
        cfuncvi[fn] = nr*nc;
      }
    }
    rr *= rr;
    for ( i = 1; i < hole_k; i++ )
      memcpy ( &cfuncvi[i*rr], cfuncvi, rr*sizeof(int) );

        /* compute the numbers of quadrature knots for each D block function */
    rr = 3*nk*m1;  /* == nfunc_d/hole_k */
    for ( fn = 0;  fn < rr;  fn++ ) {
      _g2h_FuncDSuppf ( hole_k, nk, m1, fn, 0, &nzc, &i0, &i1, &j0, &j1 );
      dfuncvi[fn] = 2*(i1-i0+1)*(j1-j0+1);
    }
    for ( i = 1; i < hole_k; i++ )
      memcpy ( &dfuncvi[i*rr], dfuncvi, rr*sizeof(int) );

        /* now compute the prefix sums and the array size */
    k = (nfunc_a+1)*hole_k*nkn*nkn;
    for ( i = 0; i < nfunc_c; i++ ) {
      j = cfuncvi[i];
      cfuncvi[i] = k;
      k += j;
    }
    for ( i = 0; i < nfunc_d; i++ ) {
      j = dfuncvi[i];
      dfuncvi[i] = k;
      k += j;
    }
    nlspr->ftabsize = k;
printf ( "ftabsize = %d\n", k );
  }
  return nlspr;
} /*_g2h_AllocSplNLPrf*/

static boolean _g2h_InitSplNLprf ( GHoleDomainf *domain,
                                   G2HNLSPrivatef *nlspr,
                                   int nconstr )
{
  G2HolePrivateRecf  *privateG2;
  G2HoleSPrivateRecf *sprivateG2;
  G2HNLPrivatef *nlpr;
  int hole_k, nfunc_a, nfunc_c, nfunc_d, nk;
  int nfunc, nkn, nkn2, ftabsize;

  privateG2 = domain->privateG2;
  sprivateG2 = domain->SprivateG2;
  hole_k = domain->hole_k;
  nlpr = &nlspr->nlpr;
  nfunc_a = privateG2->nfunc_a;
  nfunc_c = sprivateG2->nsfunc_c;
  nfunc_d = sprivateG2->nsfunc_d;
  nfunc = nfunc_a+nfunc_c+nfunc_d;
  ftabsize = nlspr->ftabsize;
  nk    = sprivateG2->nk;
  nkn   = (nk+1)*G2_QUAD_FACTOR;
  nkn2  = nkn*nkn;
  nlpr->auxc = 0;
  nlpr->nldi = pkv_GetScratchMem ( nlspr->psize*hole_k*sizeof(point3f) );
  nlpr->acoeff = pkv_GetScratchMem ( nfunc*sizeof(vector3f) );
  nlpr->rhole_cp = pkv_GetScratchMem ( (12*hole_k+1+nconstr)*sizeof(vector3f) );
  nlpr->diu = pkv_GetScratchMem ( 9*nkn2*hole_k*sizeof(vector2f) );
  nlpr->jac = pkv_GetScratchMemf ( nkn2*hole_k );
  nlpr->psiu = pkv_GetScratchMemf ( 9*ftabsize );

  if ( !nlpr->nldi || !nlpr->acoeff || !nlpr->rhole_cp ||
       !nlpr->diu || !nlpr->jac || !nlpr->psiu )
    return false;

  nlpr->div = &nlpr->diu[nkn2*hole_k];
  nlpr->diuu = &nlpr->div[nkn2*hole_k];
  nlpr->diuv = &nlpr->diuu[nkn2*hole_k];
  nlpr->divv = &nlpr->diuv[nkn2*hole_k];
  nlpr->diuuu = &nlpr->divv[nkn2*hole_k];
  nlpr->diuuv = &nlpr->diuuu[nkn2*hole_k];
  nlpr->diuvv = &nlpr->diuuv[nkn2*hole_k];
  nlpr->divvv = &nlpr->diuvv[nkn2*hole_k];
  nlpr->psiv = &nlpr->psiu[ftabsize];
  nlpr->psiuu = &nlpr->psiv[ftabsize];
  nlpr->psiuv = &nlpr->psiuu[ftabsize];
  nlpr->psivv = &nlpr->psiuv[ftabsize];
  nlpr->psiuuu = &nlpr->psivv[ftabsize];
  nlpr->psiuuv = &nlpr->psiuuu[ftabsize];
  nlpr->psiuvv = &nlpr->psiuuv[ftabsize];
  nlpr->psivvv = &nlpr->psiuvv[ftabsize];
  return true;
} /*_g2h_InitSplNLprf*/

static boolean _g2h_TabSplNLBasisFunctionsf ( GHoleDomainf *domain,
                                              G2HNLSPrivatef *nlspr )
{
  void     *sp;
  G2HNLPrivatef      *nlpr;
  GHolePrivateRecf   *privateG;
  G2HolePrivateRecf  *privateG2;
  G2HoleSPrivateRecf *sprivate;
  int      hole_k, nk, m1, m2;
  int      nfunc_a, nfunc_b, nfunc_c, nfabc;
  int      i, j, k, l, f, rr, fn, nkn, nkn2, kN, kNQ2,
           i0, i1, j0, j1, ik, jk, nzc, sn, dn, ii, jj, fN;
  int      *fkn, *lkn;
  point2f  ddi;
  vector2f *diu, *div, *diuu, *diuv, *divv, *diuuu, *diuuv, *diuvv, *divvv;
  float    *fc00, *fc01, *fc02, *fc10, *fc11, *fc12,
           *fd00, *fd01, *fd02, *fd10, *fd11, *fd12,
           *bc00, *bc01, *bc02, *bc10, *bc11, *bc12,
           *bd00, *bd01, *bd02, *bd10, *bd11, *bd12;
  float    *psiu, *psiv, *psiuu, *psiuv, *psivv,
           *psiuuu, *psiuuv, *psiuvv, *psivvv;
  float    *cb, *cbt, *cbtt, *cbttt;
  int      lastomcknot, lastpvknot, lastpvvknot, lastfpknot, disize;
  float    *omcknots, *pvknots, *pvvknots,
           *fcomc, *fpv, *fpvv, *fpu, *fpuu, *tkn;
  float    *hfunc, *dhfunc, *ddhfunc, *dddhfunc, *di;
  float    p, pu, puu, puuu, q, qv, qvv, qvvv, pq[9], bvz;
  float    p0, p0u, p0uu, p0uuu, p1, p1u, p1uu, p1uuu, p2, p2u, p2uu, p2uuu;
  unsigned char *bfcpn;

  sp = pkv_GetScratchMemTop ();
  hole_k = domain->hole_k;
  nlpr = &nlspr->nlpr;
  privateG  = domain->privateG;
  privateG2 = domain->privateG2;
  sprivate  = domain->SprivateG2;
  nk = sprivate->nk;
  m1 = sprivate->m1;
  m2 = sprivate->m2;
  nkn = nlspr->nkn;
  nkn2 = nkn*nkn;
  tkn = nlspr->tkn;
  nfunc_a = privateG2->nfunc_a;
  nfunc_b = privateG2->nfunc_b;
  nfunc_c = sprivate->nsfunc_c;
  nfabc = nfunc_a+nfunc_b+nfunc_c;
  hfunc = pkv_GetScratchMemf ( 24*nkn );
  lastfpknot = sprivate->lastfpknot;
  disize = lastfpknot-G2H_FINALDEG;
  disize *= disize;
  di = pkv_GetScratchMem ( disize*sizeof(point2f) );
  if ( !hfunc || !di ) {
    domain->error_code = G2H_ERROR_NO_SCRATCH_MEMORY;
    goto failure;
  }
  dhfunc = &hfunc[6*nkn];
  ddhfunc = &dhfunc[6*nkn];
  dddhfunc = &ddhfunc[6*nkn];
  if ( !mbs_TabQuinticHFuncDer3f ( 0.0, 1.0, nkn, tkn,
                                   hfunc, dhfunc, ddhfunc, dddhfunc ) )
    goto failure;

        /* compute the Jacobian */
  for ( k = kN = kNQ2 = 0;
        k < hole_k; 
        k++, kN += disize, kNQ2 += nkn2 ) {
    pkv_Selectf ( disize, 2, 3, 2, &nlpr->nldi[kN], di );
    for ( i = l = 0; i < nkn; i++ )
          /* this is not an optimal algorithm, to be improved */
      for ( j = 0;  j < nkn;  j++, l++ ) {
        mbs_deBoorDer3Pf ( G2H_FINALDEG, lastfpknot, sprivate->fpknots,
                       G2H_FINALDEG, lastfpknot, sprivate->fpknots, 2,
                       2*(lastfpknot-G2H_FINALDEG), di,
                       tkn[i], tkn[j], (float*)&ddi,
                       (float*)&nlpr->diu[kNQ2+l], (float*)&nlpr->div[kNQ2+l],
                       (float*)&nlpr->diuu[kNQ2+l], (float*)&nlpr->diuv[kNQ2+l],
                       (float*)&nlpr->divv[kNQ2+l],
                       (float*)&nlpr->diuuu[kNQ2+l], (float*)&nlpr->diuuv[kNQ2+l],
                       (float*)&nlpr->diuvv[kNQ2+l], (float*)&nlpr->divvv[kNQ2+l] );
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

        /* deal with the block A functions */
  for ( f = fN = 0;  f < nfunc_a;  f++ ) {
    for ( k = kNQ2 = 0; k < hole_k; k++, fN += nkn2, kNQ2 += nkn2 ) {
      _g2h_GetBFAPatchCurvesf ( domain, f, k, &fc00, &fc01, &fc02,
                                              &fd00, &fd01, &fd02 );
      _g2h_TabNLDer0f ( domain, nkn, tkn, hfunc, dhfunc, ddhfunc, dddhfunc,
                        &nlpr->diu[kNQ2], &nlpr->div[kNQ2],
                        &nlpr->diuu[kNQ2], &nlpr->diuv[kNQ2], &nlpr->divv[kNQ2],
                        &nlpr->diuuu[kNQ2], &nlpr->diuuv[kNQ2],
                        &nlpr->diuvv[kNQ2], &nlpr->divvv[kNQ2],
                        fc00, fc01, fc02, fd00, fd01, fd02,
                        &nlpr->psiu[fN], &nlpr->psiv[fN],
                        &nlpr->psiuu[fN], &nlpr->psiuv[fN], &nlpr->psivv[fN],
                        &nlpr->psiuuu[fN], &nlpr->psiuuv[fN],
                        &nlpr->psiuvv[fN], &nlpr->psivvv[fN] );
    }
  }
        /* the fixed combination of the block B functions */
  bc00 = pkv_GetScratchMemf ( 2*(G2_CROSSDEGSUM+6)*hole_k );
  if ( !bc00 ) {
    domain->error_code = G2H_ERROR_NO_SCRATCH_MEMORY;
    goto failure;
  }
  bc01 = &bc00[(G2_CROSS00DEG+1)*hole_k];  bc02 = &bc01[(G2_CROSS01DEG+1)*hole_k];
  bc10 = &bc02[(G2_CROSS02DEG+1)*hole_k];  bc11 = &bc10[(G2_CROSS10DEG+1)*hole_k];
  bc12 = &bc11[(G2_CROSS11DEG+1)*hole_k];  bd00 = &bc12[(G2_CROSS12DEG+1)*hole_k];
  bd01 = &bd00[(G2_CROSS00DEG+1)*hole_k];  bd02 = &bd01[(G2_CROSS01DEG+1)*hole_k];
  bd10 = &bd02[(G2_CROSS02DEG+1)*hole_k];  bd11 = &bd10[(G2_CROSS10DEG+1)*hole_k];
  bd12 = &bd11[(G2_CROSS11DEG+1)*hole_k];

  memset ( bc00, 0, 2*(G2_CROSSDEGSUM+6)*hole_k*sizeof(float) );
  bfcpn = privateG->bfcpn;
  for ( f = 0; f < nfunc_b; f++ ) {
        /* find the Coons representation of the constant part of the solution */
    bvz = nlpr->rhole_cp[bfcpn[f]].z;
#ifdef DEBUG_HESSIAN
  /* the test is done for the zero boundary conditions */
bvz = 0.0;
#endif
    for ( k = 0; k < hole_k; k++ ) {
      _g2h_GetBFBPatchCurvesf ( domain, f, k,
                                &fc00, &fc01, &fc02, &fc10, &fc11, &fc12,
                                &fd00, &fd01, &fd02, &fd10, &fd11, &fd12 );
      pkn_AddMatrixMf ( 1, G2_CROSS00DEG+1, 0, &bc00[(G2_CROSS00DEG+1)*k], 0, fc00,
                        bvz, 0, &bc00[(G2_CROSS00DEG+1)*k] );
      pkn_AddMatrixMf ( 1, G2_CROSS01DEG+1, 0, &bc01[(G2_CROSS01DEG+1)*k], 0, fc01,
                        bvz, 0, &bc01[(G2_CROSS01DEG+1)*k] );
      pkn_AddMatrixMf ( 1, G2_CROSS02DEG+1, 0, &bc02[(G2_CROSS02DEG+1)*k], 0, fc02,
                        bvz, 0, &bc02[(G2_CROSS02DEG+1)*k] );
      pkn_AddMatrixMf ( 1, G2_CROSS10DEG+1, 0, &bc10[(G2_CROSS10DEG+1)*k], 0, fc10,
                        bvz, 0, &bc10[(G2_CROSS10DEG+1)*k] );
      pkn_AddMatrixMf ( 1, G2_CROSS11DEG+1, 0, &bc11[(G2_CROSS11DEG+1)*k], 0, fc11,
                        bvz, 0, &bc11[(G2_CROSS11DEG+1)*k] );
      pkn_AddMatrixMf ( 1, G2_CROSS12DEG+1, 0, &bc12[(G2_CROSS12DEG+1)*k], 0, fc12,
                        bvz, 0, &bc12[(G2_CROSS12DEG+1)*k] );
      pkn_AddMatrixMf ( 1, G2_CROSS00DEG+1, 0, &bd00[(G2_CROSS00DEG+1)*k], 0, fd00,
                        bvz, 0, &bd00[(G2_CROSS00DEG+1)*k] );
      pkn_AddMatrixMf ( 1, G2_CROSS01DEG+1, 0, &bd01[(G2_CROSS01DEG+1)*k], 0, fd01,
                        bvz, 0, &bd01[(G2_CROSS01DEG+1)*k] );
      pkn_AddMatrixMf ( 1, G2_CROSS02DEG+1, 0, &bd02[(G2_CROSS02DEG+1)*k], 0, fd02,
                        bvz, 0, &bd02[(G2_CROSS02DEG+1)*k] );
      pkn_AddMatrixMf ( 1, G2_CROSS10DEG+1, 0, &bd10[(G2_CROSS10DEG+1)*k], 0, fd10,
                        bvz, 0, &bd10[(G2_CROSS10DEG+1)*k] );
      pkn_AddMatrixMf ( 1, G2_CROSS11DEG+1, 0, &bd11[(G2_CROSS11DEG+1)*k], 0, fd11,
                        bvz, 0, &bd11[(G2_CROSS11DEG+1)*k] );
      pkn_AddMatrixMf ( 1, G2_CROSS12DEG+1, 0, &bd12[(G2_CROSS12DEG+1)*k], 0, fd12,
                        bvz, 0, &bd12[(G2_CROSS12DEG+1)*k] );
    }
  }
  for ( k = kNQ2 = 0, fN = nfunc_a*nkn2*hole_k;
        k < hole_k;
        k++, kNQ2 += nkn2, fN += nkn2 ) {
    _g2h_TabNLDerf ( domain, nkn, tkn, hfunc, dhfunc, ddhfunc, dddhfunc,
             &nlpr->diu[kNQ2], &nlpr->div[kNQ2],
             &nlpr->diuu[kNQ2], &nlpr->diuv[kNQ2], &nlpr->divv[kNQ2],
             &nlpr->diuuu[kNQ2], &nlpr->diuuv[kNQ2],
             &nlpr->diuvv[kNQ2], &nlpr->divvv[kNQ2],
             &bc00[(G2_CROSS00DEG+1)*k], &bc01[(G2_CROSS01DEG+1)*k],
             &bc02[(G2_CROSS02DEG+1)*k],
             &bc10[(G2_CROSS10DEG+1)*k], &bc11[(G2_CROSS11DEG+1)*k],
             &bc12[(G2_CROSS12DEG+1)*k],
             &bd00[(G2_CROSS00DEG+1)*k], &bd01[(G2_CROSS01DEG+1)*k],
             &bd02[(G2_CROSS02DEG+1)*k],
             &bd10[(G2_CROSS10DEG+1)*k], &bd11[(G2_CROSS11DEG+1)*k],
             &bd12[(G2_CROSS12DEG+1)*k],
             &nlpr->psiu[fN], &nlpr->psiv[fN],
             &nlpr->psiuu[fN], &nlpr->psiuv[fN], &nlpr->psivv[fN],
             &nlpr->psiuuu[fN], &nlpr->psiuuv[fN],
             &nlpr->psiuvv[fN], &nlpr->psivvv[fN] );
  }

        /* now the block C functions */
  rr    = G2H_FINALDEG-5+nk*m2;
  fkn   = nlspr->fkn;
  lkn   = nlspr->lkn;
  cb    = nlspr->cb;
  cbt   = nlspr->cbt;
  cbtt  = nlspr->cbtt;
  cbttt = nlspr->cbttt;
  for ( i = fn = 0;  i < hole_k;  i++ ) {
    diu = &nlpr->diu[i*nkn2];
    div = &nlpr->div[i*nkn2];
    diuu = &nlpr->diuu[i*nkn2];
    diuv = &nlpr->diuv[i*nkn2];
    divv = &nlpr->divv[i*nkn2];
    diuuu = &nlpr->diuuu[i*nkn2];
    diuuv = &nlpr->diuuv[i*nkn2];
    diuvv = &nlpr->diuvv[i*nkn2];
    divvv = &nlpr->divvv[i*nkn2];
    for ( j = 0; j < rr; j++ ) {
      i0 = fkn[j];  i1 = lkn[j];
      for ( k = 0;  k < rr;  k++, fn++ ) {
        jk = nlspr->cfuncvi[fn];
        psiu = &nlpr->psiu[jk];
        psiv = &nlpr->psiv[jk];
        psiuu = &nlpr->psiuu[jk];
        psiuv = &nlpr->psiuv[jk];
        psivv = &nlpr->psivv[jk];
        psiuuu = &nlpr->psiuuu[jk];
        psiuuv = &nlpr->psiuuv[jk];
        psiuvv = &nlpr->psiuvv[jk];
        psivvv = &nlpr->psivvv[jk];
        j0 = fkn[k];  j1 = lkn[k];
        for ( ik = i0, sn = 0;  ik <= i1;  ik++ ) {
          p = cb[j*nkn+ik];      pu = cbt[j*nkn+ik];
          puu = cbtt[j*nkn+ik];  puuu = cbttt[j*nkn+ik];
          for ( jk = j0, dn = ik*nkn+j0;
                jk <= j1;
                jk++, sn++, dn++ ) {
            q = cb[k*nkn+jk];      qv = cbt[k*nkn+jk];
            qvv = cbtt[k*nkn+jk];  qvvv = cbttt[k*nkn+jk];
            _g2h_TensDer3f ( p, pu, puu, puuu, q, qv, qvv, qvvv, pq );
            if ( !pkn_Comp2iDerivatives3f ( diu[dn].x, diu[dn].y, div[dn].x, div[dn].y,
                        diuu[dn].x, diuu[dn].y, diuv[dn].x, diuv[dn].y,
                        divv[dn].x, divv[dn].y,
                        diuuu[dn].x, diuuu[dn].y, diuuv[dn].x, diuuv[dn].y,
                        diuvv[dn].x, diuvv[dn].y, divvv[dn].x, divvv[dn].y,
                        1, &pq[0], &pq[1], &pq[2], &pq[3], &pq[4],
                        &pq[5], &pq[6], &pq[7], &pq[8],
                        &psiu[sn], &psiv[sn], &psiuu[sn], &psiuv[sn], &psivv[sn],
                        &psiuuu[sn], &psiuuv[sn], &psiuvv[sn], &psivvv[sn] ) )
              goto failure;
          }
        }
      }
    }
  }
        /* and the block D functions */
  lastomcknot = sprivate->lastomcknot;
  omcknots    = sprivate->omcknots;
  lastpvknot  = sprivate->lastpvknot;
  pvknots     = sprivate->pvknots;
  lastpvvknot = sprivate->lastpvvknot;
  pvvknots    = sprivate->pvvknots;
  fcomc = pkv_GetScratchMemf ( lastomcknot+2*(lastpvknot+lastpvvknot)-
                        (G2_CROSS00DEG+2*(G2_CROSS01DEG+G2_CROSS02DEG)) );
  if ( !fcomc ) {
    domain->error_code = G2H_ERROR_NO_SCRATCH_MEMORY;
    goto failure;
  }
  fpv = &fcomc[lastomcknot-G2_CROSS00DEG];
  fpvv = &fpv[lastpvknot-G2_CROSS01DEG];
  fpu = &fpvv[lastpvvknot-G2_CROSS02DEG];
  fpuu = &fpu[lastpvknot-G2_CROSS01DEG];
  rr = 3*nk*m1;  /* == nfunc_d/hole_k */
          /* this is done in two stages; first the functions related to */
          /* the curve Gamma_i are evaluated at the quadrature knots */
          /* located in Omega_i */
  for ( i = fn = 0;  i < hole_k;  i++ ) {
    diu = &nlpr->diu[i*nkn2];
    div = &nlpr->div[i*nkn2];
    diuu = &nlpr->diuu[i*nkn2];
    diuv = &nlpr->diuv[i*nkn2];
    divv = &nlpr->divv[i*nkn2];
    diuuu = &nlpr->diuuu[i*nkn2];
    diuuv = &nlpr->diuuv[i*nkn2];
    diuvv = &nlpr->diuvv[i*nkn2];
    divvv = &nlpr->divvv[i*nkn2];
    for ( j = 0;  j < rr;  j++, fn++ ) {
      _g2h_FuncDSuppf ( hole_k, nk, m1, fn, i, &nzc, &i0, &i1, &j0, &j1 );
      _g2h_GetSplDBasisCrossDerf ( domain, nfabc+fn, i,
                                   fcomc, fpv, fpvv, fpu, fpuu );
      jk = nlspr->dfuncvi[fn];
      psiu = &nlpr->psiu[jk];
      psiv = &nlpr->psiv[jk];
      psiuu = &nlpr->psiuu[jk];
      psiuv = &nlpr->psiuv[jk];
      psivv = &nlpr->psivv[jk];
      psiuuu = &nlpr->psiuuu[jk];
      psiuuv = &nlpr->psiuuv[jk];
      psiuvv = &nlpr->psiuvv[jk];
      psivvv = &nlpr->psivvv[jk];
      for ( ik = i0, sn = 0, dn = i0*nkn;  ik <= i1;  ik++ ) {
        if ( nzc == 0 )
          mbs_deBoorDer3C1f ( G2_CROSS00DEG, lastomcknot, omcknots, fcomc,
                              tkn[ik], &p0, &p0u, &p0uu, &p0uuu );
        if ( nzc <= 1 )
          mbs_deBoorDer3C1f ( G2_CROSS01DEG, lastpvknot, pvknots, fpv,
                              tkn[ik], &p1, &p1u, &p1uu, &p1uuu );
        mbs_deBoorDer3C1f ( G2_CROSS02DEG, lastpvvknot, pvvknots, fpvv,
                            tkn[ik], &p2, &p2u, &p2uu, &p2uuu );
        for ( jk = j0, jj = 6*j0;  jk <= j1;  jk++, jj += 6, sn++, dn++ ) {
          switch ( nzc ) {
        case 0:
            pq[0] = p0u*hfunc[jj]   + p1u*hfunc[jj+2]   + p2u*hfunc[jj+4];
            pq[1] = p0*dhfunc[jj]   + p1*dhfunc[jj+2]   + p2*dhfunc[jj+4];
            pq[2] = p0uu*hfunc[jj]  + p1uu*hfunc[jj+2]  + p2uu*hfunc[jj+4];
            pq[3] = p0u*dhfunc[jj]  + p1u*dhfunc[jj+2]  + p2u*dhfunc[jj+4];
            pq[4] = p0*ddhfunc[jj]  + p1*ddhfunc[jj+2]  + p2*ddhfunc[jj+4];
            pq[5] = p0uuu*hfunc[jj] + p1uuu*hfunc[jj+2] + p2uuu*hfunc[jj+4];
            pq[6] = p0uu*dhfunc[jj] + p1uu*dhfunc[jj+2] + p2uu*dhfunc[jj+4];
            pq[7] = p0u*ddhfunc[jj] + p1u*ddhfunc[jj+2] + p2u*ddhfunc[jj+4];
            pq[8] = p0*dddhfunc[jj] + p1*dddhfunc[jj+2] + p2*dddhfunc[jj+4];
            break;
        case 1:
            pq[0] = p1u*hfunc[jj+2]   + p2u*hfunc[jj+4];
            pq[1] = p1*dhfunc[jj+2]   + p2*dhfunc[jj+4];
            pq[2] = p1uu*hfunc[jj+2]  + p2uu*hfunc[jj+4];
            pq[3] = p1u*dhfunc[jj+2]  + p2u*dhfunc[jj+4];
            pq[4] = p1*ddhfunc[jj+2]  + p2*ddhfunc[jj+4];
            pq[5] = p1uuu*hfunc[jj+2] + p2uuu*hfunc[jj+4];
            pq[6] = p1uu*dhfunc[jj+2] + p2uu*dhfunc[jj+4];
            pq[7] = p1u*ddhfunc[jj+2] + p2u*ddhfunc[jj+4];
            pq[8] = p1*dddhfunc[jj+2] + p2*dddhfunc[jj+4];
            break;
        case 2:
            pq[0] = p2u*hfunc[jj+4];
            pq[1] = p2*dhfunc[jj+4];
            pq[2] = p2uu*hfunc[jj+4];
            pq[3] = p2u*dhfunc[jj+4];
            pq[4] = p2*ddhfunc[jj+4];
            pq[5] = p2uuu*hfunc[jj+4];
            pq[6] = p2uu*dhfunc[jj+4];
            pq[7] = p2u*ddhfunc[jj+4];
            pq[8] = p2*dddhfunc[jj+4];
            break;
          }
          if ( !pkn_Comp2iDerivatives3f ( diu[dn].x, diu[dn].y, div[dn].x, div[dn].y,
                        diuu[dn].x, diuu[dn].y, diuv[dn].x, diuv[dn].y,
                        divv[dn].x, divv[dn].y,
                        diuuu[dn].x, diuuu[dn].y, diuuv[dn].x, diuuv[dn].y,
                        diuvv[dn].x, diuvv[dn].y, divvv[dn].x, divvv[dn].y,
                        1, &pq[0], &pq[1], &pq[2], &pq[3], &pq[4],
                        &pq[5], &pq[6], &pq[7], &pq[8],
                        &psiu[sn], &psiv[sn], &psiuu[sn], &psiuv[sn], &psivv[sn],
                        &psiuuu[sn], &psiuuv[sn], &psiuvv[sn], &psivvv[sn] ) )
            goto failure;
        }
      }
    }
  }
          /* and then in Omega_{i-1} */
  for ( ii = hole_k-1, i = fn = 0;  i < hole_k; ii = i++ ) {
    diu = &nlpr->diu[ii*nkn2];
    div = &nlpr->div[ii*nkn2];
    diuu = &nlpr->diuu[ii*nkn2];
    diuv = &nlpr->diuv[ii*nkn2];
    divv = &nlpr->divv[ii*nkn2];
    diuuu = &nlpr->diuuu[ii*nkn2];
    diuuv = &nlpr->diuuv[ii*nkn2];
    diuvv = &nlpr->diuvv[ii*nkn2];
    divvv = &nlpr->divvv[ii*nkn2];
    for ( j = 0;  j < rr;  j++, fn++ ) {
      _g2h_FuncDSuppf ( hole_k, nk, m1, fn, ii, &nzc, &i0, &i1, &j0, &j1 );
      _g2h_GetSplDBasisCrossDerf ( domain, nfabc+fn, i,
                                   fcomc, fpv, fpvv, fpu, fpuu );
      jk = nlspr->dfuncvi[fn] + (i1-i0+1)*(j1-j0+1);
      psiu = &nlpr->psiu[jk];
      psiv = &nlpr->psiv[jk];
      psiuu = &nlpr->psiuu[jk];
      psiuv = &nlpr->psiuv[jk];
      psivv = &nlpr->psivv[jk];
      psiuuu = &nlpr->psiuuu[jk];
      psiuuv = &nlpr->psiuuv[jk];
      psiuvv = &nlpr->psiuvv[jk];
      psivvv = &nlpr->psivvv[jk];
      for ( jk = j0; jk <= j1; jk++ ) {
        if ( nzc == 0 )
          mbs_deBoorDer3C1f ( G2_CROSS00DEG, lastomcknot, omcknots, fcomc,
                              tkn[jk], &p0, &p0u, &p0uu, &p0uuu );
        if ( nzc <= 1 )
          mbs_deBoorDer3C1f ( G2_CROSS01DEG, lastpvknot, pvknots, fpu,
                              tkn[jk], &p1, &p1u, &p1uu, &p1uuu );
        mbs_deBoorDer3C1f ( G2_CROSS02DEG, lastpvvknot, pvvknots, fpuu,
                            tkn[jk], &p2, &p2u, &p2uu, &p2uuu );
        for ( ik = i0, jj = 6*i0, sn = jk-j0, dn = i0*nkn+jk;
              ik <= i1;
              ik++, jj += 6, sn += (j1-j0+1), dn += nkn ) {
          switch ( nzc ) {
        case 0:
            pq[0] = p0*dhfunc[jj]   + p1*dhfunc[jj+2]   + p2*dhfunc[jj+4];
            pq[1] = p0u*hfunc[jj]   + p1u*hfunc[jj+2]   + p2u*hfunc[jj+4];
            pq[2] = p0*ddhfunc[jj]  + p1*ddhfunc[jj+2]  + p2*ddhfunc[jj+4];
            pq[3] = p0u*dhfunc[jj]  + p1u*dhfunc[jj+2]  + p2u*dhfunc[jj+4];
            pq[4] = p0uu*hfunc[jj]  + p1uu*hfunc[jj+2]  + p2uu*hfunc[jj+4];
            pq[5] = p0*dddhfunc[jj] + p1*dddhfunc[jj+2] + p2*dddhfunc[jj+4];
            pq[6] = p0u*ddhfunc[jj] + p1u*ddhfunc[jj+2] + p2u*ddhfunc[jj+4];
            pq[7] = p0uu*dhfunc[jj] + p1uu*dhfunc[jj+2] + p2uu*dhfunc[jj+4];
            pq[8] = p0uuu*hfunc[jj] + p1uuu*hfunc[jj+2] + p2uuu*hfunc[jj+4];
            break;
        case 1:
            pq[0] = p1*dhfunc[jj+2]   + p2*dhfunc[jj+4];
            pq[1] = p1u*hfunc[jj+2]   + p2u*hfunc[jj+4];
            pq[2] = p1*ddhfunc[jj+2]  + p2*ddhfunc[jj+4];
            pq[3] = p1u*dhfunc[jj+2]  + p2u*dhfunc[jj+4];
            pq[4] = p1uu*hfunc[jj+2]  + p2uu*hfunc[jj+4];
            pq[5] = p1*dddhfunc[jj+2] + p2*dddhfunc[jj+4];
            pq[6] = p1u*ddhfunc[jj+2] + p2u*ddhfunc[jj+4];
            pq[7] = p1uu*dhfunc[jj+2] + p2uu*dhfunc[jj+4];
            pq[8] = p1uuu*hfunc[jj+2] + p2uuu*hfunc[jj+4];
            break;
        case 2:
            pq[0] = p2*dhfunc[jj+4];
            pq[1] = p2u*hfunc[jj+4];
            pq[2] = p2*ddhfunc[jj+4];
            pq[3] = p2u*dhfunc[jj+4];
            pq[4] = p2uu*hfunc[jj+4];
            pq[5] = p2*dddhfunc[jj+4];
            pq[6] = p2u*ddhfunc[jj+4];
            pq[7] = p2uu*dhfunc[jj+4];
            pq[8] = p2uuu*hfunc[jj+4];
            break;
          }
          if ( !pkn_Comp2iDerivatives3f ( diu[dn].x, diu[dn].y, div[dn].x, div[dn].y,
                     diuu[dn].x, diuu[dn].y, diuv[dn].x, diuv[dn].y,
                     divv[dn].x, divv[dn].y,
                     diuuu[dn].x, diuuu[dn].y, diuuv[dn].x, diuuv[dn].y,
                     diuvv[dn].x, diuvv[dn].y, divvv[dn].x, divvv[dn].y,
                     1, &pq[0], &pq[1], &pq[2], &pq[3], &pq[4],
                     &pq[5], &pq[6], &pq[7], &pq[8],
                     &psiu[sn], &psiv[sn], &psiuu[sn], &psiuv[sn], &psivv[sn],
                     &psiuuu[sn], &psiuuv[sn], &psiuvv[sn], &psivvv[sn] ) )
            goto failure;
        }
      }
    }
  }

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*_g2h_TabSplNLBasisFunctionsf*/

/* ///////////////////////////////////////////////////////////////////////// */
static float _g2h_ComputeSplNLFuncf ( GHoleDomainf *domain,
                                      G2HNLSPrivatef *nlspr,
                                      const float *coeff )
{
  G2HolePrivateRecf  *privateG2;
  G2HoleSPrivateRecf *sprivate;
  G2HNLPrivatef      *nlpr;
  G2HNLFuncf         f;
  int   hole_k, nfunc_a, nfunc_c, nfunc_d, nfcd;
  float c, funct;
  int   *fkn, *lkn, *cfuncvi, *dfuncvi;
  int   nk, m1, m2, nkn, nkn2, fi, rrc, rrd, bs1,
        fic, fjc, fn, i0, i1, j0, j1, nzc;
  int   k, kk, knot, i, j;

  privateG2 = domain->privateG2;
  sprivate = domain->SprivateG2;
  nlpr = &nlspr->nlpr;
  hole_k = domain->hole_k;
  nfunc_a = privateG2->nfunc_a;
  nfunc_c = sprivate->nsfunc_c;
  nfunc_d = sprivate->nsfunc_d;
  nfcd  = nfunc_c+nfunc_d;
  nk    = sprivate->nk;
  m1    = sprivate->m1;
  m2    = sprivate->m2;
  nkn   = nlspr->nkn;
  nkn2  = nkn*nkn;
  fkn   = nlspr->fkn;
  lkn   = nlspr->lkn;
  cfuncvi = nlspr->cfuncvi;
  dfuncvi = nlspr->dfuncvi;

  rrc = G2H_FINALDEG-5+nk*m2;
  bs1 = rrc*rrc;
  rrd = 3*nk*m1;
  funct = 0.0;
  for ( k = knot = 0, kk = 1;
        k < hole_k;
        k++, kk = (kk+1) % hole_k ) {

        /* evaluate the function and its derivatives */
    for ( i = 0; i < nkn; i++ )
      for ( j = 0;  j < nkn;  j++, knot++ ) {
          /* get the fixed linear combination of the block B functions */
        fi = nfunc_a*nkn2*hole_k + knot;
        f.pu = nlpr->psiu[fi];      f.pv = nlpr->psiv[fi];
        f.puu = nlpr->psiuu[fi];    f.puv = nlpr->psiuv[fi];
        f.pvv = nlpr->psivv[fi];    f.puuu = nlpr->psiuuu[fi];
        f.puuv = nlpr->psiuuv[fi];  f.puvv = nlpr->psiuvv[fi];
        f.pvvv = nlpr->psivvv[fi];
          /* add the linear combination of the block A functions */
        for ( fn = 0; fn < nfunc_a; fn++ ) {
          fi = fn*nkn2*hole_k + knot;
          c = coeff[nfcd+fn];
          f.pu -= c*nlpr->psiu[fi];      f.pv -= c*nlpr->psiv[fi];
          f.puu -= c*nlpr->psiuu[fi];    f.puv -= c*nlpr->psiuv[fi];
          f.pvv -= c*nlpr->psivv[fi];    f.puuu -= c*nlpr->psiuuu[fi];
          f.puuv -= c*nlpr->psiuuv[fi];  f.puvv -= c*nlpr->psiuvv[fi];
          f.pvvv -= c*nlpr->psivvv[fi];
        }
          /* add the linear combination of the block C functions */
        for ( fic = 0; fic < rrc; fic++ )
          if ( i >= fkn[fic] && i <= lkn[fic] )
            for ( fjc = 0; fjc < rrc; fjc++ )
              if ( j >= fkn[fjc] && j <= lkn[fjc] ) {
                fn = bs1*k + fic*rrc + fjc;  /* function number */
                fi = cfuncvi[fn] + (i-fkn[fic])*(lkn[fjc]-fkn[fjc]+1) +
                                   (j-fkn[fjc]);
                c = coeff[fn];
                f.pu -= c*nlpr->psiu[fi];      f.pv -= c*nlpr->psiv[fi];
                f.puu -= c*nlpr->psiuu[fi];    f.puv -= c*nlpr->psiuv[fi];
                f.pvv -= c*nlpr->psivv[fi];    f.puuu -= c*nlpr->psiuuu[fi];
                f.puuv -= c*nlpr->psiuuv[fi];  f.puvv -= c*nlpr->psiuvv[fi];
                f.pvvv -= c*nlpr->psivvv[fi];
              }

          /* add the linear combination of the block D functions */
        for ( fic = 0; fic < rrd; fic++ ) {
          fn = k*rrd + fic;   /* function number */
          _g2h_FuncDSuppf ( hole_k, nk, m1, fn, k,
                            &nzc, &i0, &i1, &j0, &j1 );
          if ( i >= i0 && i <= i1 && j >= j0 && j <= j1 ) {
            fi = dfuncvi[fn] + (i-i0)*(j1-j0+1) + (j-j0);
            c = coeff[nfunc_c+fn];
            f.pu -= c*nlpr->psiu[fi];      f.pv -= c*nlpr->psiv[fi];
            f.puu -= c*nlpr->psiuu[fi];    f.puv -= c*nlpr->psiuv[fi];
            f.pvv -= c*nlpr->psivv[fi];    f.puuu -= c*nlpr->psiuuu[fi];
            f.puuv -= c*nlpr->psiuuv[fi];  f.puvv -= c*nlpr->psiuvv[fi];
            f.pvvv -= c*nlpr->psivvv[fi];
          }
          fn = kk*rrd + fic;  /* function number */
          _g2h_FuncDSuppf ( hole_k, nk, m1, fn, k,
                            &nzc, &i0, &i1, &j0, &j1 );
          if ( i >= i0 && i <= i1 && j >= j0 && j <= j1 ) {
            fi = dfuncvi[fn] + (i1-i0+1)*(j1-j0+1) + (i-i0)*(j1-j0+1) + (j-j0);
            c = coeff[nfunc_c+fn];
            f.pu -= c*nlpr->psiu[fi];      f.pv -= c*nlpr->psiv[fi];
            f.puu -= c*nlpr->psiuu[fi];    f.puv -= c*nlpr->psiuv[fi];
            f.pvv -= c*nlpr->psivv[fi];    f.puuu -= c*nlpr->psiuuu[fi];
            f.puuv -= c*nlpr->psiuuv[fi];  f.puvv -= c*nlpr->psiuvv[fi];
            f.pvvv -= c*nlpr->psivvv[fi];
          }
        }
        f.jac = nlpr->jac[knot];
        /* integrate */
        _g2h_IntFunc1af ( &f, &funct );
      }
  }

  return funct/(float)nkn2;
} /*_g2h_ComputeSplNLFuncf*/

static void _g2h_ComputeSplNLFuncGradf ( GHoleDomainf *domain,
                                            G2HNLSPrivatef *nlspr,
                                            const float *coeff,
                                            float *func, float *grad )
{
  G2HolePrivateRecf  *privateG2;
  G2HoleSPrivateRecf *sprivate;
  G2HNLPrivatef      *nlpr;
  G2HNLFuncf         f;
  int   hole_k, nfunc_a, nfunc_c, nfunc_d, nfcd, nfacd;
  float c, funct;
  int   *fkn, *lkn, *cfuncvi, *dfuncvi;
  int   nk, m1, m2, nkn, nkn2, fi, rrc, rrd, bs1,
        fic, fjc, fn, i0, i1, j0, j1, nzc;
  int   k, kk, knot, i, j;

  privateG2 = domain->privateG2;
  sprivate = domain->SprivateG2;
  nlpr = &nlspr->nlpr;
  hole_k = domain->hole_k;
  nfunc_a = privateG2->nfunc_a;
  nfunc_c = sprivate->nsfunc_c;
  nfunc_d = sprivate->nsfunc_d;
  nfcd  = nfunc_c+nfunc_d;
  nfacd = nfcd+nfunc_a;
  nk    = sprivate->nk;
  m1    = sprivate->m1;
  m2    = sprivate->m2;
  nkn   = nlspr->nkn;
  nkn2  = nkn*nkn;
  fkn   = nlspr->fkn;
  lkn   = nlspr->lkn;
  cfuncvi = nlspr->cfuncvi;
  dfuncvi = nlspr->dfuncvi;

  rrc = G2H_FINALDEG-5+nk*m2;
  bs1 = rrc*rrc;
  rrd = 3*nk*m1;

  funct = 0.0;
  memset ( grad, 0, nfacd*sizeof(float) );

  for ( k = knot = 0, kk = 1;
        k < hole_k;
        k++, kk = (kk+1) % hole_k ) {

        /* evaluate the function and its derivatives */
    for ( i = 0; i < nkn; i++ )
      for ( j = 0;  j < nkn;  j++, knot++ ) {
          /* get the fixed linear combination of the block B functions */
        fi = nfunc_a*nkn2*hole_k + knot;
        f.pu = nlpr->psiu[fi];      f.pv = nlpr->psiv[fi];
        f.puu = nlpr->psiuu[fi];    f.puv = nlpr->psiuv[fi];
        f.pvv = nlpr->psivv[fi];    f.puuu = nlpr->psiuuu[fi];
        f.puuv = nlpr->psiuuv[fi];  f.puvv = nlpr->psiuvv[fi];
        f.pvvv = nlpr->psivvv[fi];
          /* add the linear combination of the block A functions */
        for ( fn = 0; fn < nfunc_a; fn++ ) {
          fi = fn*nkn2*hole_k + knot;
          c = coeff[nfcd+fn];
          f.pu -= c*nlpr->psiu[fi];      f.pv -= c*nlpr->psiv[fi];
          f.puu -= c*nlpr->psiuu[fi];    f.puv -= c*nlpr->psiuv[fi];
          f.pvv -= c*nlpr->psivv[fi];    f.puuu -= c*nlpr->psiuuu[fi];
          f.puuv -= c*nlpr->psiuuv[fi];  f.puvv -= c*nlpr->psiuvv[fi];
          f.pvvv -= c*nlpr->psivvv[fi];
        }
          /* add the linear combination of the block C functions */
        for ( fic = 0; fic < rrc; fic++ )
          if ( i >= fkn[fic] && i <= lkn[fic] )
            for ( fjc = 0; fjc < rrc; fjc++ )
              if ( j >= fkn[fjc] && j <= lkn[fjc] ) {
                fn = bs1*k + fic*rrc + fjc;  /* function number */
                fi = cfuncvi[fn] + (i-fkn[fic])*(lkn[fjc]-fkn[fjc]+1) +
                                   (j-fkn[fjc]);
                c = coeff[fn];
                f.pu -= c*nlpr->psiu[fi];      f.pv -= c*nlpr->psiv[fi];
                f.puu -= c*nlpr->psiuu[fi];    f.puv -= c*nlpr->psiuv[fi];
                f.pvv -= c*nlpr->psivv[fi];    f.puuu -= c*nlpr->psiuuu[fi];
                f.puuv -= c*nlpr->psiuuv[fi];  f.puvv -= c*nlpr->psiuvv[fi];
                f.pvvv -= c*nlpr->psivvv[fi];
              }

          /* add the linear combination of the block D functions */
        for ( fic = 0; fic < rrd; fic++ ) {
          fn = k*rrd + fic;   /* function number */
          _g2h_FuncDSuppf ( hole_k, nk, m1, fn, k,
                            &nzc, &i0, &i1, &j0, &j1 );
          if ( i >= i0 && i <= i1 && j >= j0 && j <= j1 ) {
            fi = dfuncvi[fn] + (i-i0)*(j1-j0+1) + (j-j0);
            c = coeff[nfunc_c+fn];
            f.pu -= c*nlpr->psiu[fi];      f.pv -= c*nlpr->psiv[fi];
            f.puu -= c*nlpr->psiuu[fi];    f.puv -= c*nlpr->psiuv[fi];
            f.pvv -= c*nlpr->psivv[fi];    f.puuu -= c*nlpr->psiuuu[fi];
            f.puuv -= c*nlpr->psiuuv[fi];  f.puvv -= c*nlpr->psiuvv[fi];
            f.pvvv -= c*nlpr->psivvv[fi];
          }
          fn = kk*rrd + fic;  /* function number */
          _g2h_FuncDSuppf ( hole_k, nk, m1, fn, k,
                            &nzc, &i0, &i1, &j0, &j1 );
          if ( i >= i0 && i <= i1 && j >= j0 && j <= j1 ) {
            fi = dfuncvi[fn] + (i1-i0+1)*(j1-j0+1) + (i-i0)*(j1-j0+1) + (j-j0);
            c = coeff[nfunc_c+fn];
            f.pu -= c*nlpr->psiu[fi];      f.pv -= c*nlpr->psiv[fi];
            f.puu -= c*nlpr->psiuu[fi];    f.puv -= c*nlpr->psiuv[fi];
            f.pvv -= c*nlpr->psivv[fi];    f.puuu -= c*nlpr->psiuuu[fi];
            f.puuv -= c*nlpr->psiuuv[fi];  f.puvv -= c*nlpr->psiuvv[fi];
            f.pvvv -= c*nlpr->psivvv[fi];
          }
        }
        f.jac = nlpr->jac[knot];
/* 1. Integration of the functional value */
        _g2h_IntFunc1bf ( &f, &funct );

/* 2. Integration of the functional gradient */
          /* derivatives with respect to the block A coefficients */
        for ( fn = 0; fn < nfunc_a; fn++ ) {
          fi = fn*nkn2*hole_k + knot;
          f.psiu = nlpr->psiu[fi];      f.psiv = nlpr->psiv[fi];
          f.psiuu = nlpr->psiuu[fi];    f.psiuv = nlpr->psiuv[fi];
          f.psivv = nlpr->psivv[fi];    f.psiuuu = nlpr->psiuuu[fi];
          f.psiuuv = nlpr->psiuuv[fi];  f.psiuvv = nlpr->psiuvv[fi];
          f.psivvv = nlpr->psivvv[fi];
          _g2h_IntFunc2bf ( &f, &grad[nfcd+fn] );
        }
          /* derivatives with respect to the block C coefficients */
        for ( fic = 0; fic < rrc; fic++ )
          if ( i >= fkn[fic] && i <= lkn[fic] )
            for ( fjc = 0; fjc < rrc; fjc++ )
              if ( j >= fkn[fjc] && j <= lkn[fjc] ) {
                fn = bs1*k + fic*rrc + fjc;  /* function number */
                fi = cfuncvi[fn] + (i-fkn[fic])*(lkn[fjc]-fkn[fjc]+1) +
                                   (j-fkn[fjc]);
                f.psiu = nlpr->psiu[fi];      f.psiv = nlpr->psiv[fi];
                f.psiuu = nlpr->psiuu[fi];    f.psiuv = nlpr->psiuv[fi];
                f.psivv = nlpr->psivv[fi];    f.psiuuu = nlpr->psiuuu[fi];  
                f.psiuuv = nlpr->psiuuv[fi];  f.psiuvv = nlpr->psiuvv[fi];
                f.psivvv = nlpr->psivvv[fi];
                _g2h_IntFunc2bf ( &f, &grad[fn] );
              }
          /* derivatives with respect to the block D coefficients */
        for ( fic = 0; fic < rrd; fic++ ) {
          fn = k*rrd + fic;   /* function number */
          _g2h_FuncDSuppf ( hole_k, nk, m1, fn, k,
                            &nzc, &i0, &i1, &j0, &j1 );
          if ( i >= i0 && i <= i1 && j >= j0 && j <= j1 ) {
            fi = dfuncvi[fn] + (i-i0)*(j1-j0+1) + (j-j0);
            f.psiu = nlpr->psiu[fi];      f.psiv = nlpr->psiv[fi];
            f.psiuu = nlpr->psiuu[fi];    f.psiuv = nlpr->psiuv[fi];
            f.psivv = nlpr->psivv[fi];    f.psiuuu = nlpr->psiuuu[fi];  
            f.psiuuv = nlpr->psiuuv[fi];  f.psiuvv = nlpr->psiuvv[fi];
            f.psivvv = nlpr->psivvv[fi];
            _g2h_IntFunc2bf ( &f, &grad[nfunc_c+fn] );
          }
          fn = kk*rrd + fic;  /* function number */
          _g2h_FuncDSuppf ( hole_k, nk, m1, fn, k,
                            &nzc, &i0, &i1, &j0, &j1 );
          if ( i >= i0 && i <= i1 && j >= j0 && j <= j1 ) {
            fi = dfuncvi[fn] + (i1-i0+1)*(j1-j0+1) + (i-i0)*(j1-j0+1) + (j-j0);
            f.psiu = nlpr->psiu[fi];      f.psiv = nlpr->psiv[fi];
            f.psiuu = nlpr->psiuu[fi];    f.psiuv = nlpr->psiuv[fi];
            f.psivv = nlpr->psivv[fi];    f.psiuuu = nlpr->psiuuu[fi];  
            f.psiuuv = nlpr->psiuuv[fi];  f.psiuvv = nlpr->psiuvv[fi];
            f.psivvv = nlpr->psivvv[fi];
            _g2h_IntFunc2bf ( &f, &grad[nfunc_c+fn] );
          }
        }
      }
  }
  *func = funct/(float)nkn2;
  pkn_MultMatrixNumf ( 1, nfacd, 0, grad, 1.0/(float)nkn2, 0, grad );
} /*_g2h_ComputeSplNLFuncGradf*/

/* The Hessian is represented as the Block1 structured symmetric matrix, */
/* instead of the more compact Block2. This is to allow imposing         */
/* constraints, which affect directly the block D basis functions        */
static boolean _g2h_ComputeSplNLFuncGradHessianf ( GHoleDomainf *domain,
                                   G2HNLSPrivatef *nlspr,
                                   const float *coeff,
                                   float *func, float *grad, float *hessian )
{
  void *sp;
  G2HolePrivateRecf  *privateG2;
  G2HoleSPrivateRecf *sprivate;
  G2HNLPrivatef      *nlpr;
  G2HNLFuncf         f;
  int   hole_k, nfunc_a, nfunc_c, nfunc_d, nfad, nfcd, nfacd;
  int   hsize;
  float *Li, *Bi, *BiLT, *Di, *hii, *hki, *hkk;
  float c, funct;
  int   *fkn, *lkn, *cfuncvi, *dfuncvi;
  int   nk, m1, m2, nkn, nkn2, fi, fj, rrc, rrd, bs1,
        fic, fjc, fn, fim, fjm, fm, pos,
        i00, i01, j00, j01, i10, i11, j10, j11, nzc;
  int   k, kk, knot, i, j;

  sp = pkv_GetScratchMemTop ();
  privateG2 = domain->privateG2;
  sprivate = domain->SprivateG2;
  nlpr = &nlspr->nlpr;
  hole_k = domain->hole_k;
  nfunc_a = privateG2->nfunc_a;
  nfunc_c = sprivate->nsfunc_c;
  nfunc_d = sprivate->nsfunc_d;
  nfad  = nfunc_a+nfunc_d;
  nfcd  = nfunc_c+nfunc_d;
  nfacd = nfcd+nfunc_a;
  Li = pkv_GetScratchMemf ( 8*nfacd );
  if ( !Li ) {
    domain->error_code = G2H_ERROR_NO_SCRATCH_MEMORY;
    goto failure;
  }
  Bi = &Li[2*nfacd]; BiLT = &Bi[3*nfacd];  Di = &BiLT[2*nfacd];

  nk    = sprivate->nk;
  m1    = sprivate->m1;
  m2    = sprivate->m2;
  nkn   = nlspr->nkn;
  nkn2  = nkn*nkn;
  fkn   = nlspr->fkn;
  lkn   = nlspr->lkn;
  cfuncvi = nlspr->cfuncvi;
  dfuncvi = nlspr->dfuncvi;

  rrc = G2H_FINALDEG-5+nk*m2;
  bs1 = rrc*rrc;
  rrd = 3*nk*m1;
  hsize = pkn_Block1ArraySize ( hole_k, bs1, nfad );
  hkk = &hessian[pkn_Block1FindBlockPos(hole_k,bs1,nfad,hole_k,hole_k)];

  funct = 0.0;
  memset ( grad, 0, nfacd*sizeof(float) );
  memset ( hessian, 0, hsize*sizeof(float) );

  for ( k = knot = 0, kk = 1;
        k < hole_k;
        k++, kk = (kk+1) % hole_k ) {
    hii = &hessian[pkn_Block1FindBlockPos(hole_k,bs1,nfad,k,k)];
    hki = &hessian[pkn_Block1FindBlockPos(hole_k,bs1,nfad,hole_k,k)];

        /* evaluate the function and its derivatives */
    for ( i = 0; i < nkn; i++ )
      for ( j = 0;  j < nkn;  j++, knot++ ) {
          /* get the fixed linear combination of the block B functions */
        fi = nfunc_a*nkn2*hole_k + knot;
        f.pu = nlpr->psiu[fi];      f.pv = nlpr->psiv[fi];
        f.puu = nlpr->psiuu[fi];    f.puv = nlpr->psiuv[fi];
        f.pvv = nlpr->psivv[fi];    f.puuu = nlpr->psiuuu[fi];
        f.puuv = nlpr->psiuuv[fi];  f.puvv = nlpr->psiuvv[fi];
        f.pvvv = nlpr->psivvv[fi];
          /* add the linear combination of the block A functions */
        for ( fn = 0; fn < nfunc_a; fn++ ) {
          fi = fn*nkn2*hole_k + knot;
          c = coeff[nfcd+fn];
          f.pu -= c*nlpr->psiu[fi];      f.pv -= c*nlpr->psiv[fi];
          f.puu -= c*nlpr->psiuu[fi];    f.puv -= c*nlpr->psiuv[fi];
          f.pvv -= c*nlpr->psivv[fi];    f.puuu -= c*nlpr->psiuuu[fi];
          f.puuv -= c*nlpr->psiuuv[fi];  f.puvv -= c*nlpr->psiuvv[fi];
          f.pvvv -= c*nlpr->psivvv[fi];
        }
          /* add the linear combination of the block C functions */
        for ( fic = 0; fic < rrc; fic++ )
          if ( i >= fkn[fic] && i <= lkn[fic] )
            for ( fjc = 0; fjc < rrc; fjc++ )
              if ( j >= fkn[fjc] && j <= lkn[fjc] ) {
                fn = bs1*k + fic*rrc + fjc;  /* function number */
                fi = cfuncvi[fn] + (i-fkn[fic])*(lkn[fjc]-fkn[fjc]+1) +
                                   (j-fkn[fjc]);
                c = coeff[fn];
                f.pu -= c*nlpr->psiu[fi];      f.pv -= c*nlpr->psiv[fi];
                f.puu -= c*nlpr->psiuu[fi];    f.puv -= c*nlpr->psiuv[fi];
                f.pvv -= c*nlpr->psivv[fi];    f.puuu -= c*nlpr->psiuuu[fi];
                f.puuv -= c*nlpr->psiuuv[fi];  f.puvv -= c*nlpr->psiuvv[fi];
                f.pvvv -= c*nlpr->psivvv[fi];
              }

          /* add the linear combination of the block D functions */
        for ( fic = 0; fic < rrd; fic++ ) {
          fn = k*rrd + fic;   /* function number */
          _g2h_FuncDSuppf ( hole_k, nk, m1, fn, k,
                            &nzc, &i00, &i01, &j00, &j01 );
          if ( i >= i00 && i <= i01 && j >= j00 && j <= j01 ) {
            fi = dfuncvi[fn] + (i-i00)*(j01-j00+1) + (j-j00);
            c = coeff[nfunc_c+fn];
            f.pu -= c*nlpr->psiu[fi];      f.pv -= c*nlpr->psiv[fi];
            f.puu -= c*nlpr->psiuu[fi];    f.puv -= c*nlpr->psiuv[fi];
            f.pvv -= c*nlpr->psivv[fi];    f.puuu -= c*nlpr->psiuuu[fi];
            f.puuv -= c*nlpr->psiuuv[fi];  f.puvv -= c*nlpr->psiuvv[fi];
            f.pvvv -= c*nlpr->psivvv[fi];
          }
          fn = kk*rrd + fic;  /* function number */
          _g2h_FuncDSuppf ( hole_k, nk, m1, fn, k,
                            &nzc, &i00, &i01, &j00, &j01 );
          if ( i >= i00 && i <= i01 && j >= j00 && j <= j01 ) {
            fi = dfuncvi[fn] + (i01-i00+1)*(j01-j00+1) + (i-i00)*(j01-j00+1) + (j-j00);
            c = coeff[nfunc_c+fn];
            f.pu -= c*nlpr->psiu[fi];      f.pv -= c*nlpr->psiv[fi];
            f.puu -= c*nlpr->psiuu[fi];    f.puv -= c*nlpr->psiuv[fi];
            f.pvv -= c*nlpr->psivv[fi];    f.puuu -= c*nlpr->psiuuu[fi];
            f.puuv -= c*nlpr->psiuuv[fi];  f.puvv -= c*nlpr->psiuvv[fi];
            f.pvvv -= c*nlpr->psivvv[fi];
          }
        }
        f.jac = nlpr->jac[knot];
/* 1. Integration of the functional value */
        _g2h_IntFunc1cf ( &f, &funct );

/* 2. Integration of the functional gradient */
          /* derivatives with respect to the block A coefficients */
        for ( fn = 0; fn < nfunc_a; fn++ ) {
          fi = fn*nkn2*hole_k + knot;
          f.psiu = nlpr->psiu[fi];      f.psiv = nlpr->psiv[fi];
          f.psiuu = nlpr->psiuu[fi];    f.psiuv = nlpr->psiuv[fi];
          f.psivv = nlpr->psivv[fi];    f.psiuuu = nlpr->psiuuu[fi];
          f.psiuuv = nlpr->psiuuv[fi];  f.psiuvv = nlpr->psiuvv[fi];
          f.psivvv = nlpr->psivvv[fi];
          _g2h_IntFunc2cf ( &f, (vector2f*)&Li[2*(nfcd+fn)],
                            &Bi[3*(nfcd+fn)], (vector2f*)&BiLT[2*(nfcd+fn)],
                            &Di[nfcd+fn], &grad[nfcd+fn] );
        }
          /* derivatives with respect to the block C coefficients */
        for ( fic = 0; fic < rrc; fic++ )
          if ( i >= fkn[fic] && i <= lkn[fic] )
            for ( fjc = 0; fjc < rrc; fjc++ )
              if ( j >= fkn[fjc] && j <= lkn[fjc] ) {
                fn = bs1*k + fic*rrc + fjc;  /* function number */
                fi = cfuncvi[fn] + (i-fkn[fic])*(lkn[fjc]-fkn[fjc]+1) +
                                   (j-fkn[fjc]);
                f.psiu = nlpr->psiu[fi];      f.psiv = nlpr->psiv[fi];
                f.psiuu = nlpr->psiuu[fi];    f.psiuv = nlpr->psiuv[fi];
                f.psivv = nlpr->psivv[fi];    f.psiuuu = nlpr->psiuuu[fi];
                f.psiuuv = nlpr->psiuuv[fi];  f.psiuvv = nlpr->psiuvv[fi];
                f.psivvv = nlpr->psivvv[fi];
                _g2h_IntFunc2cf ( &f, (vector2f*)&Li[2*fn],
                                  &Bi[3*fn], (vector2f*)&BiLT[2*fn],
                                  &Di[fn], &grad[fn] );
              }
          /* derivatives with respect to the block D coefficients */
        for ( fic = 0; fic < rrd; fic++ ) {
          fn = k*rrd + fic;   /* function number */
          _g2h_FuncDSuppf ( hole_k, nk, m1, fn, k,
                            &nzc, &i00, &i01, &j00, &j01 );
          if ( i >= i00 && i <= i01 && j >= j00 && j <= j01 ) {
            fi = dfuncvi[fn] + (i-i00)*(j01-j00+1) + (j-j00);
            f.psiu = nlpr->psiu[fi];      f.psiv = nlpr->psiv[fi];
            f.psiuu = nlpr->psiuu[fi];    f.psiuv = nlpr->psiuv[fi];
            f.psivv = nlpr->psivv[fi];    f.psiuuu = nlpr->psiuuu[fi];
            f.psiuuv = nlpr->psiuuv[fi];  f.psiuvv = nlpr->psiuvv[fi];
            f.psivvv = nlpr->psivvv[fi];
            _g2h_IntFunc2cf ( &f, (vector2f*)&Li[2*(nfunc_c+fn)],
                              &Bi[3*(nfunc_c+fn)], (vector2f*)&BiLT[2*(nfunc_c+fn)],
                              &Di[nfunc_c+fn], &grad[nfunc_c+fn] );
          }
          fn = kk*rrd + fic;  /* function number */
          _g2h_FuncDSuppf ( hole_k, nk, m1, fn, k,
                            &nzc, &i00, &i01, &j00, &j01 );
          if ( i >= i00 && i <= i01 && j >= j00 && j <= j01 ) {
            fi = dfuncvi[fn] + (i01-i00+1)*(j01-j00+1) +
                               (i-i00)*(j01-j00+1) + (j-j00);
            f.psiu = nlpr->psiu[fi];      f.psiv = nlpr->psiv[fi];
            f.psiuu = nlpr->psiuu[fi];    f.psiuv = nlpr->psiuv[fi];
            f.psivv = nlpr->psivv[fi];    f.psiuuu = nlpr->psiuuu[fi];
            f.psiuuv = nlpr->psiuuv[fi];  f.psiuvv = nlpr->psiuvv[fi];
            f.psivvv = nlpr->psivvv[fi];
            _g2h_IntFunc2cf ( &f, (vector2f*)&Li[2*(nfunc_c+fn)],
                              &Bi[3*(nfunc_c+fn)], (vector2f*)&BiLT[2*(nfunc_c+fn)],
                              &Di[nfunc_c+fn], &grad[nfunc_c+fn] );
          }
        }

/* 3. Integration of the functional Hessian */
        for ( fn = 0; fn < nfunc_a; fn++ ) {
          fi = fn*nkn2*hole_k + knot;
          f.psiu = nlpr->psiu[fi];      f.psiv = nlpr->psiv[fi];
          f.psiuu = nlpr->psiuu[fi];    f.psiuv = nlpr->psiuv[fi];
          f.psivv = nlpr->psivv[fi];    f.psiuuu = nlpr->psiuuu[fi];
          f.psiuuv = nlpr->psiuuv[fi];  f.psiuvv = nlpr->psiuvv[fi];
          f.psivvv = nlpr->psivvv[fi];
          /* A x A */
          for ( fm = 0; fm <= fn; fm++ ) {
            fj = fm*nkn2*hole_k + knot;
            f.psju = nlpr->psiu[fj];      f.psjv = nlpr->psiv[fj];
            f.psjuu = nlpr->psiuu[fj];    f.psjuv = nlpr->psiuv[fj];
            f.psjvv = nlpr->psivv[fj];    f.psjuuu = nlpr->psiuuu[fj];
            f.psjuuv = nlpr->psiuuv[fj];  f.psjuvv = nlpr->psiuvv[fj];
            f.psjvvv = nlpr->psivvv[fj];
            pos = pkn_SymMatIndex ( nfunc_d+fn, nfunc_d+fm );
            _g2h_IntFunc3cf ( &f, (vector2f*)&Li[2*(nfcd+fn)],
                              (vector2f*)&Li[2*(nfcd+fm)],
                              (vector2f*)&BiLT[2*(nfcd+fn)],
                              (vector2f*)&BiLT[2*(nfcd+fm)],
                              Di[nfcd+fn], Di[nfcd+fm], &hkk[pos] );
          }
          /* A x C */
          for ( fic = 0; fic < rrc; fic++ )
            if ( i >= fkn[fic] && i <= lkn[fic] )
              for ( fjc = 0; fjc < rrc; fjc++ )
                if ( j >= fkn[fjc] && j <= lkn[fjc] ) {
                  fm = (k*rrc + fic)*rrc + fjc;  /* function number */
                  fj = cfuncvi[fm] + (i-fkn[fic])*(lkn[fjc]-fkn[fjc]+1) + j-fkn[fjc];
                  f.psju = nlpr->psiu[fj];      f.psjv = nlpr->psiv[fj];
                  f.psjuu = nlpr->psiuu[fj];    f.psjuv = nlpr->psiuv[fj];
                  f.psjvv = nlpr->psivv[fj];    f.psjuuu = nlpr->psiuuu[fj];
                  f.psjuuv = nlpr->psiuuv[fj];  f.psjuvv = nlpr->psiuvv[fj];
                  f.psjvvv = nlpr->psivvv[fj];
                  pos = (nfunc_d+fn)*bs1 + fic*rrc + fjc;
                  _g2h_IntFunc3cf ( &f, (vector2f*)&Li[2*(nfcd+fn)],
                                    (vector2f*)&Li[2*fm],
                                    (vector2f*)&BiLT[2*(nfcd+fn)],
                                    (vector2f*)&BiLT[2*fm],
                                    Di[nfcd+fn], Di[fm], &hki[pos] );
                }
          /* A x D */
            /* block D functions related with Gamma_i */
          for ( fic = 0, fm = k*rrd;  fic < rrd;  fic++, fm++ ) {
            _g2h_FuncDSuppf ( hole_k, nk, m1, fm, k,
                              &nzc, &i10, &i11, &j10, &j11 );
            if ( i >= i10 && i <= i11 && j >= j10 && j <= j11 ) {
              fj = dfuncvi[fm] + (i-i10)*(j11-j10+1) + j-j10;
              f.psju = nlpr->psiu[fj];      f.psjv = nlpr->psiv[fj];
              f.psjuu = nlpr->psiuu[fj];    f.psjuv = nlpr->psiuv[fj];
              f.psjvv = nlpr->psivv[fj];    f.psjuuu = nlpr->psiuuu[fj];
              f.psjuuv = nlpr->psiuuv[fj];  f.psjuvv = nlpr->psiuvv[fj];
              f.psjvvv = nlpr->psivvv[fj];
              pos = pkn_SymMatIndex ( nfunc_d+fn, fm );
              _g2h_IntFunc3cf ( &f, (vector2f*)&Li[2*(nfcd+fn)],
                              (vector2f*)&Li[2*(nfunc_c+fm)],
                              (vector2f*)&BiLT[2*(nfcd+fn)],
                              (vector2f*)&BiLT[2*(nfunc_c+fm)],
                              Di[nfcd+fn], Di[nfunc_c+fm], &hkk[pos] );
            }
          }
            /* block D functions related with Gamma_{i+1} */
          for ( fic = 0, fm = kk*rrd;  fic < rrd;  fic++, fm++ ) {
            _g2h_FuncDSuppf ( hole_k, nk, m1, fm, k,
                              &nzc, &i10, &i11, &j10, &j11 );
            if ( i >= i10 && i <= i11 && j >= j10 && j <= j11 ) {
              fj = dfuncvi[fm] + (i11-i10+1)*(j11-j10+1) + (i-i10)*(j11-j10+1) + j-j10;
              f.psju = nlpr->psiu[fj];      f.psjv = nlpr->psiv[fj];
              f.psjuu = nlpr->psiuu[fj];    f.psjuv = nlpr->psiuv[fj];
              f.psjvv = nlpr->psivv[fj];    f.psjuuu = nlpr->psiuuu[fj];
              f.psjuuv = nlpr->psiuuv[fj];  f.psjuvv = nlpr->psiuvv[fj];
              f.psjvvv = nlpr->psivvv[fj];
              pos = pkn_SymMatIndex ( nfunc_d+fn, fm );
              _g2h_IntFunc3cf ( &f, (vector2f*)&Li[2*(nfcd+fn)],
                              (vector2f*)&Li[2*(nfunc_c+fm)],
                              (vector2f*)&BiLT[2*(nfcd+fn)],
                              (vector2f*)&BiLT[2*(nfunc_c+fm)],
                              Di[nfcd+fn], Di[nfunc_c+fm], &hkk[pos] );
            }
          }
        }

        for ( fic = 0; fic < rrc; fic++ )
          if ( i >= fkn[fic] && i <= lkn[fic] )
            for ( fjc = 0; fjc < rrc; fjc++ )
              if ( j >= fkn[fjc] && j <= lkn[fjc] ) {
                fn = (k*rrc + fic)*rrc + fjc;  /* function number */
                fi = cfuncvi[fn] + (i-fkn[fic])*(lkn[fjc]-fkn[fjc]+1) + j-fkn[fjc];
                f.psiu = nlpr->psiu[fi];      f.psiv = nlpr->psiv[fi];
                f.psiuu = nlpr->psiuu[fi];    f.psiuv = nlpr->psiuv[fi];
                f.psivv = nlpr->psivv[fi];    f.psiuuu = nlpr->psiuuu[fi];
                f.psiuuv = nlpr->psiuuv[fi];  f.psiuvv = nlpr->psiuvv[fi];
                f.psivvv = nlpr->psivvv[fi];
          /* C x C */
                for ( fim = 0; fim < rrc; fim++ )
                  if ( i >= fkn[fim] && i <= lkn[fim] )
                    for ( fjm = 0; fjm < rrc; fjm++ )
                      if ( j >= fkn[fjm] && j <= lkn[fjm] ) {
                        fm = (k*rrc + fim)*rrc + fjm;
                        if ( fm <= fn ) {
                          fj = cfuncvi[fm] + (i-fkn[fim])*(lkn[fjm]-fkn[fjm]+1) + j-fkn[fjm];
                          f.psju = nlpr->psiu[fj];      f.psjv = nlpr->psiv[fj];
                          f.psjuu = nlpr->psiuu[fj];    f.psjuv = nlpr->psiuv[fj];
                          f.psjvv = nlpr->psivv[fj];    f.psjuuu = nlpr->psiuuu[fj];
                          f.psjuuv = nlpr->psiuuv[fj];  f.psjuvv = nlpr->psiuvv[fj];
                          f.psjvvv = nlpr->psivvv[fj];
                          pos = pkn_SymMatIndex ( fic*rrc+fjc, fim*rrc+fjm );
                          _g2h_IntFunc3cf ( &f, (vector2f*)&Li[2*fn],
                              (vector2f*)&Li[2*fm],
                              (vector2f*)&BiLT[2*fn],
                              (vector2f*)&BiLT[2*fm],
                              Di[fn], Di[fm], &hii[pos] );
                        }
                        else goto cont_c;
                      }
cont_c: ;
          /* C x D */
            /* block D functions related with Gamma_i */
                for ( fim = 0, fm = k*rrd;  fim < rrd;  fim++, fm++ ) {
                  _g2h_FuncDSuppf ( hole_k, nk, m1, fm, k,
                                    &nzc, &i10, &i11, &j10, &j11 );
                  if ( i >= i10 && i <= i11 && j >= j10 && j <= j11 ) {
                    fj = dfuncvi[fm] + (i-i10)*(j11-j10+1) + j-j10;
                    f.psju = nlpr->psiu[fj];      f.psjv = nlpr->psiv[fj];
                    f.psjuu = nlpr->psiuu[fj];    f.psjuv = nlpr->psiuv[fj];
                    f.psjvv = nlpr->psivv[fj];    f.psjuuu = nlpr->psiuuu[fj];
                    f.psjuuv = nlpr->psiuuv[fj];  f.psjuvv = nlpr->psiuvv[fj];
                    f.psjvvv = nlpr->psivvv[fj];
                    pos = bs1*fm + fic*rrc+fjc;
                    _g2h_IntFunc3cf ( &f, (vector2f*)&Li[2*fn],
                              (vector2f*)&Li[2*(nfunc_c+fm)],
                              (vector2f*)&BiLT[2*fn],
                              (vector2f*)&BiLT[2*(nfunc_c+fm)],
                              Di[fn], Di[nfunc_c+fm], &hki[pos] );
                  }
                }
            /* block D functions related with Gamma_{i+1} */
                for ( fim = 0, fm = kk*rrd;  fim < rrd;  fim++, fm++ ) {
                  _g2h_FuncDSuppf ( hole_k, nk, m1, fm, k,
                                    &nzc, &i10, &i11, &j10, &j11 );
                  if ( i >= i10 && i <= i11 && j >= j10 && j <= j11 ) {
                    fj = dfuncvi[fm] + (i11-i10+1)*(j11-j10+1) + (i-i10)*(j11-j10+1) + j-j10;
                    f.psju = nlpr->psiu[fj];      f.psjv = nlpr->psiv[fj];
                    f.psjuu = nlpr->psiuu[fj];    f.psjuv = nlpr->psiuv[fj];
                    f.psjvv = nlpr->psivv[fj];    f.psjuuu = nlpr->psiuuu[fj];
                    f.psjuuv = nlpr->psiuuv[fj];  f.psjuvv = nlpr->psiuvv[fj];
                    f.psjvvv = nlpr->psivvv[fj];
                    pos = bs1*fm + fic*rrc + fjc;
                    _g2h_IntFunc3cf ( &f, (vector2f*)&Li[2*fn],
                              (vector2f*)&Li[2*(nfunc_c+fm)],
                              (vector2f*)&BiLT[2*fn],
                              (vector2f*)&BiLT[2*(nfunc_c+fm)],
                              Di[fn], Di[nfunc_c+fm], &hki[pos] );
                  }
                }
              }
          /* D x D */
        for ( fic = 0, fn = k*rrd;  fic < rrd;  fic++, fn++ ) {
          _g2h_FuncDSuppf ( hole_k, nk, m1, fn, k,
                            &nzc, &i00, &i01, &j00, &j01 );
          if ( i >= i00 && i <= i01 && j >= j00 && j <= j01 ) {
            fi = dfuncvi[fn] + (i-i00)*(j01-j00+1) + j-j00;
            f.psiu = nlpr->psiu[fi];      f.psiv = nlpr->psiv[fi];
            f.psiuu = nlpr->psiuu[fi];    f.psiuv = nlpr->psiuv[fi];
            f.psivv = nlpr->psivv[fi];    f.psiuuu = nlpr->psiuuu[fi];
            f.psiuuv = nlpr->psiuuv[fi];  f.psiuvv = nlpr->psiuvv[fi];
            f.psivvv = nlpr->psivvv[fi];
            for ( fim = 0, fm = k*rrd;  fim < rrd;  fim++, fm++ )
              if ( fm <= fn ) {
                _g2h_FuncDSuppf ( hole_k, nk, m1, fm, k,
                                  &nzc, &i10, &i11, &j10, &j11 );
                if ( i >= i10 && i <= i11 && j >= j10 && j <= j11 ) {
                  fj = dfuncvi[fm] + (i-i10)*(j11-j10+1) + j-j10;
                  f.psju = nlpr->psiu[fj];      f.psjv = nlpr->psiv[fj];
                  f.psjuu = nlpr->psiuu[fj];    f.psjuv = nlpr->psiuv[fj];
                  f.psjvv = nlpr->psivv[fj];    f.psjuuu = nlpr->psiuuu[fj];
                  f.psjuuv = nlpr->psiuuv[fj];  f.psjuvv = nlpr->psiuvv[fj];
                  f.psjvvv = nlpr->psivvv[fj];
                  pos = pkn_SymMatIndex ( fn, fm );
                  _g2h_IntFunc3cf ( &f, (vector2f*)&Li[2*(nfunc_c+fn)],
                              (vector2f*)&Li[2*(nfunc_c+fm)],
                              (vector2f*)&BiLT[2*(nfunc_c+fn)],
                              (vector2f*)&BiLT[2*(nfunc_c+fm)],
                              Di[nfunc_c+fn], Di[nfunc_c+fm], &hkk[pos] );
                }
              }
            for ( fim = 0, fm = kk*rrd;  fim < rrd;  fim++, fm++ )
              if ( fm <= fn ) {
                _g2h_FuncDSuppf ( hole_k, nk, m1, fm, k,
                                  &nzc, &i10, &i11, &j10, &j11 );
                if ( i >= i10 && i <= i11 && j >= j10 && j <= j11 ) {
                  fj = dfuncvi[fm] + (i11-i10+1)*(j11-j10+1) + (i-i10)*(j11-j10+1) + j-j10;
                  f.psju = nlpr->psiu[fj];      f.psjv = nlpr->psiv[fj];
                  f.psjuu = nlpr->psiuu[fj];    f.psjuv = nlpr->psiuv[fj];
                  f.psjvv = nlpr->psivv[fj];    f.psjuuu = nlpr->psiuuu[fj];
                  f.psjuuv = nlpr->psiuuv[fj];  f.psjuvv = nlpr->psiuvv[fj];
                  f.psjvvv = nlpr->psivvv[fj];
                  pos = pkn_SymMatIndex ( fn, fm );
                  _g2h_IntFunc3cf ( &f, (vector2f*)&Li[2*(nfunc_c+fn)],
                              (vector2f*)&Li[2*(nfunc_c+fm)],
                              (vector2f*)&BiLT[2*(nfunc_c+fn)],
                              (vector2f*)&BiLT[2*(nfunc_c+fm)],
                              Di[nfunc_c+fn], Di[nfunc_c+fm], &hkk[pos] );
                }
              }
          }
        }
        for ( fic = 0, fn = kk*rrd;  fic < rrd;  fic++, fn++ ) {
          _g2h_FuncDSuppf ( hole_k, nk, m1, fn, k,
                            &nzc, &i00, &i01, &j00, &j01 );
          if ( i >= i00 && i <= i01 && j >= j00 && j <= j01 ) {
            fi = dfuncvi[fn] + (i01-i00+1)*(j01-j00+1) + (i-i00)*(j01-j00+1) + j-j00;
            f.psiu = nlpr->psiu[fi];      f.psiv = nlpr->psiv[fi];
            f.psiuu = nlpr->psiuu[fi];    f.psiuv = nlpr->psiuv[fi];
            f.psivv = nlpr->psivv[fi];    f.psiuuu = nlpr->psiuuu[fi];
            f.psiuuv = nlpr->psiuuv[fi];  f.psiuvv = nlpr->psiuvv[fi];
            f.psivvv = nlpr->psivvv[fi];
            for ( fim = 0, fm = k*rrd;  fim < rrd;  fim++, fm++ )
              if ( fm <= fn ) {
                _g2h_FuncDSuppf ( hole_k, nk, m1, fm, k,
                                  &nzc, &i10, &i11, &j10, &j11 );
                if ( i >= i10 && i <= i11 && j >= j10 && j <= j11 ) {
                  fj = dfuncvi[fm] + (i-i10)*(j11-j10+1) + j-j10;
                  f.psju = nlpr->psiu[fj];      f.psjv = nlpr->psiv[fj];
                  f.psjuu = nlpr->psiuu[fj];    f.psjuv = nlpr->psiuv[fj];
                  f.psjvv = nlpr->psivv[fj];    f.psjuuu = nlpr->psiuuu[fj];
                  f.psjuuv = nlpr->psiuuv[fj];  f.psjuvv = nlpr->psiuvv[fj];
                  f.psjvvv = nlpr->psivvv[fj];
                  pos = pkn_SymMatIndex ( fn, fm );
                  _g2h_IntFunc3cf ( &f, (vector2f*)&Li[2*(nfunc_c+fn)],
                              (vector2f*)&Li[2*(nfunc_c+fm)],
                              (vector2f*)&BiLT[2*(nfunc_c+fn)],
                              (vector2f*)&BiLT[2*(nfunc_c+fm)],
                              Di[nfunc_c+fn], Di[nfunc_c+fm], &hkk[pos] );
                }
              }
            for ( fim = 0, fm = kk*rrd;  fim < rrd;  fim++, fm++ )
              if ( fm <= fn ) {
                _g2h_FuncDSuppf ( hole_k, nk, m1, fm, k,
                                  &nzc, &i10, &i11, &j10, &j11 );
                if ( i >= i10 && i <= i11 && j >= j10 && j <= j11 ) {
                  fj = dfuncvi[fm] + (i11-i10+1)*(j11-j10+1) + (i-i10)*(j11-j10+1) + j-j10;
                  f.psju = nlpr->psiu[fj];      f.psjv = nlpr->psiv[fj];
                  f.psjuu = nlpr->psiuu[fj];    f.psjuv = nlpr->psiuv[fj];
                  f.psjvv = nlpr->psivv[fj];    f.psjuuu = nlpr->psiuuu[fj];
                  f.psjuuv = nlpr->psiuuv[fj];  f.psjuvv = nlpr->psiuvv[fj];
                  f.psjvvv = nlpr->psivvv[fj];
                  pos = pkn_SymMatIndex ( fn, fm );
                  _g2h_IntFunc3cf ( &f, (vector2f*)&Li[2*(nfunc_c+fn)],
                              (vector2f*)&Li[2*(nfunc_c+fm)],
                              (vector2f*)&BiLT[2*(nfunc_c+fn)],
                              (vector2f*)&BiLT[2*(nfunc_c+fm)],
                              Di[nfunc_c+fn], Di[nfunc_c+fm], &hkk[pos] );
                }
              }
          }
        }

      }
  }
  *func = funct/(float)nkn2;
  pkn_MultMatrixNumf ( 1, nfacd, 0, grad, 1.0/(float)nkn2, 0, grad );
  pkn_MultMatrixNumf ( 1, hsize, 0, hessian, 1.0/(float)nkn2, 0, hessian );

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*_g2h_ComputeSplNLFuncGradHessianf*/

static boolean g2h_SplNLNewtonf ( GHoleDomainf *domain, G2HNLSPrivatef *nlsprivate )
{
#define EPSF 2.0e-4
  void    *sp;
  G2HolePrivateRecf  *privateG2;
  G2HoleSPrivateRecf *sprivate;
  G2HNLPrivatef      *nlpr;
  int     itn, jtn, ktn, hole_k, nfunc_a, nfunc_c, nfunc_d, nfunc;
  int     bs1, bs2, asize;
  float   func, *grad, *coeff, *dcoeff, *hessian, *chess;
  float   func0, func1, aux, gn, gn0 = 0.0, gn1, dco, dyn;
  boolean positive;

  sp = pkv_GetScratchMemTop ();
  hole_k = domain->hole_k;
  privateG2 = domain->privateG2;
  sprivate  = domain->SprivateG2;
  nlpr      = &nlsprivate->nlpr;
  nfunc_a = privateG2->nfunc_a;
  nfunc_c = sprivate->nsfunc_c;
  nfunc_d = sprivate->nsfunc_d;
  nfunc = nfunc_a+nfunc_c+nfunc_d;

  coeff = pkv_GetScratchMemf ( 3*nfunc );
  bs1 = nfunc_c/hole_k;
  bs2 = nfunc_d+nfunc_a;
  asize = pkn_Block1ArraySize ( hole_k, bs1, bs2 );
  hessian = pkv_GetScratchMemf ( 2*asize );
  if ( !coeff || !hessian ) {
    domain->error_code = G2H_ERROR_NO_SCRATCH_MEMORY;
    goto failure;
  }
  dcoeff = &coeff[nfunc];
  grad   = &dcoeff[nfunc];
  chess  = &hessian[asize];

        /* setup the initial point */
  pkv_Selectf ( nfunc, 1, 3, 1, &nlpr->acoeff[0].z, coeff );

        /* Newton iterations */
  for ( itn = ktn = 0; ; itn++ ) {
    if ( !_g2h_ComputeSplNLFuncGradHessianf ( domain, nlsprivate, coeff,
                              &func, grad, hessian ) )
      goto failure;
    gn = (float)sqrt ( pkn_ScalarProductf ( nfunc, grad, grad ) );
    if ( itn == 0 ) {
printf ( "func = %f, gn0 = %f\n", func, gn );
      gn0 = gn;
    }
    memcpy ( chess, hessian, asize*sizeof(float) );
    if ( (positive = pkn_Block1CholeskyDecompMf ( hole_k, bs1, bs2,
                                                  hessian ) ) ) {
      pkn_Block1LowerTrMSolvef ( hole_k, bs1, bs2, hessian, 1, 1, grad );
      pkn_Block1UpperTrMSolvef ( hole_k, bs1, bs2, hessian, 1, 1, grad );
    }
    else {
printf ( "! " );

      if ( !pkn_Block1SymMatrixMultf ( hole_k, bs1, bs2, chess,
                                       1, 1, grad, 1, dcoeff ) )
        goto failure;
      aux = (float)pkn_ScalarProductf ( nfunc, grad, dcoeff );
      if ( gn < 0.0 || aux < EPSF*gn ) {
        domain->error_code = G2H_ERROR_NL_MINIMIZATION;
        goto failure;
      }
      pkn_MultMatrixNumf ( 1, nfunc, 0, grad, gn/aux, 0, grad );
    }
    dco = (float)sqrt ( pkn_ScalarProductf ( nfunc, coeff, coeff ) );
    dyn = (float)sqrt ( pkn_ScalarProductf ( nfunc, grad, grad ) );

    for ( aux = 1.0; aux > EPSF; aux *= 0.5 ) {
      pkn_AddMatrixMf ( 1, nfunc, 0, coeff, 0, grad, aux, 0, dcoeff );
      func1 = _g2h_ComputeSplNLFuncf ( domain, nlsprivate, dcoeff );
      if ( func1 < func )
        break;
    }

    memcpy ( coeff, dcoeff, nfunc*sizeof(float) );
    func = func1;
    ktn ++;

    if ( positive && aux > 0.1 ) {
        /* With the positive-definite Hessian we try to make some */
        /* extra iterations */

      for ( jtn = 0; jtn < 10; jtn ++ ) {

printf ( "+" );

        _g2h_ComputeSplNLFuncGradf ( domain, nlsprivate, coeff, &func0, grad );
        gn1 = (float)sqrt ( pkn_ScalarProductf ( nfunc, grad, grad ) );
        pkn_Block1LowerTrMSolvef ( hole_k, bs1, bs2, hessian, 1, 1, grad );
        pkn_Block1UpperTrMSolvef ( hole_k, bs1, bs2, hessian, 1, 1, grad );
        pkn_AddMatrixf ( 1, nfunc, 0, coeff, 0, grad, 0, dcoeff );
        func1 = _g2h_ComputeSplNLFuncf ( domain, nlsprivate, dcoeff );
        if ( func1 >= func0 || gn1 >= 0.5*gn )
          break;

printf ( "    func = %f, gn = %f\n", func1, gn );

        memcpy ( coeff, dcoeff, nfunc*sizeof(float) );
        func = func1;
        ktn ++;
        gn = gn1;
        aux = 1.0;
        dco = (float)sqrt ( pkn_ScalarProductf ( nfunc, coeff, coeff ) );
        dyn = (float)sqrt ( pkn_ScalarProductf ( nfunc, grad, grad ) );
      }
    }

    if ( _g2h_StopItf ( itn, gn0, gn, dco, dyn, aux ) )
      break;
  }

printf ( "itn = %d, ktn = %d, func = %f\n", itn+1, ktn, func1 );

  pkv_Selectf ( nfunc, 1, 1, 3, coeff, &nlpr->acoeff[0].z );

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
#undef EPSF
} /*g2h_SplNLNewtonf*/

/* ////////////////////////////////////////////////////////////////////////// */
#ifdef DEBUG_HESSIAN
static boolean _TestHessian ( GHoleDomainf *domain, G2HNLSPrivatef *nlsprivate )
{
  void *sp;
  G2HNLPrivatef      *nlprivate;
  G2HolePrivateRecf  *privateG2;
  G2HoleSPrivateRecf *sprivate;
  int   hole_k, nfunc_a, nfunc_b, nfunc_c, nfunc_d, nfunc, bs1, bs2, asize;
  float func, *grad, *hessian, *coeff;
  float *amat;
  int   i, j, pos1, pos2;
  FILE  *f;

  sp = pkv_GetScratchMemTop ();
  privateG2 = domain->privateG2;
  sprivate  = domain->SprivateG2;
  nlprivate = &nlsprivate->nlpr;
  hole_k  = domain->hole_k;
  nfunc_a = privateG2->nfunc_a;
  nfunc_b = privateG2->nfunc_b;
  nfunc_c = sprivate->nsfunc_c;
  nfunc_d = sprivate->nsfunc_d;
  nfunc = nfunc_a+nfunc_c+nfunc_d;
  bs1 = nfunc_c/hole_k;
  bs2 = nfunc_d+nfunc_a;
  coeff = pkv_GetScratchMemf ( nfunc );
  grad  = pkv_GetScratchMemf ( nfunc );
  asize = pkn_Block1ArraySize ( hole_k, bs1, bs2 );
  hessian = pkv_GetScratchMemf ( asize );
  if ( !coeff || !grad || !hessian ) {
    PKV_SIGNALERROR ( LIB_EGHOLE, ERRCODE_2, ERRMSG_2 );
    goto failure;
  }
  memset ( coeff, 0, nfunc*sizeof(float) );
  if ( !_g2h_ComputeSplNLFuncGradHessianf ( domain, nlsprivate, coeff,
                                            &func, grad, hessian ) )
    goto failure;
  f = fopen ( "g2hessianf.txt", "w+" );
  amat = sprivate->SAMat;
  for ( i = 0; i < nfunc; i++ )
    for ( j = 0; j <= i; j++ ) {
      pos1 = pkn_Block1FindElemPos ( hole_k, bs1, bs2, i, j );
      pos2 = pkn_Block2FindElemPos ( hole_k, nfunc_c/hole_k, nfunc_d/hole_k,
                                     nfunc_a, i, j );
      if ( pos1 >= 0 && pos2 >= 0 )
        if ( amat[pos2] || hessian[pos1] )
          fprintf ( f, "%4d,%4d: %12.5f %12.5f %14.7f\n", i, j,
                    0.5*amat[pos2], hessian[pos1], 0.5*amat[pos2]-hessian[pos1] );
    }
  fclose ( f );
  printf ( "%s\n", "g2hessianf.txt" );

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*_TestHessian*/
#endif
/* ////////////////////////////////////////////////////////////////////////// */
boolean g2h_NLSplFillHolef ( GHoleDomainf *domain, const point3f *hole_cp,  
                     float *acoeff, void *usrptr,
                     void (*outpatch) ( int n, int lknu, const float *knu,
                                        int m, int lknv, const float *knv,
                                        const point3f *cp, void *usrptr ) )
{
  typedef void outscf ( int n, int lknu, const float *knu,
                        int m, int lknv, const float *knv,
                        const float *cp, void *usrptr );

  void  *sp;
  G2HNLSPrivatef     *nlsprivate;
  G2HNLPrivatef      *nlprivate;
  G2HolePrivateRecf  *privateG2;
  G2HoleSPrivateRecf *sprivate;
  int   hole_k, nfunc_a, nfunc_c, nfunc_d;
  float *fc00;

  sp = pkv_GetScratchMemTop ();
  privateG2 = domain->privateG2;
  sprivate  = domain->SprivateG2;
  if ( !sprivate )
    goto failure;
  if ( !sprivate->SLMat )
    if ( !g2h_DecomposeSplMatrixf ( domain ) )
      goto failure;
  nlsprivate = _g2h_AllocSplNLPrf ( domain );
  if ( !nlsprivate )
    goto failure;
  _g2h_nlprivf = nlprivate = &nlsprivate->nlpr;
  hole_k  = domain->hole_k;
  nfunc_a = privateG2->nfunc_a;
  nfunc_c = sprivate->nsfunc_c;
  nfunc_d = sprivate->nsfunc_d;
  if ( !_g2h_InitSplNLprf ( domain, nlsprivate, 0 ) )
    goto failure;

  if ( !_g2h_ComputeNLNormalf ( domain, nlprivate, hole_cp ) )
    goto failure;
  g2h_ReflectVectorsf ( 12*hole_k+1, hole_cp, nlprivate->rhole_cp );
  nlprivate->auxc = 0;
  if ( !g2h_SplFillHolef ( domain, 3, (float*)nlprivate->rhole_cp,
                           (float*)nlprivate->acoeff, NULL, g2h_splnloutpatchf ) )
    goto failure;

  if ( !_g2h_TabSplNLBasisFunctionsf ( domain, nlsprivate ) )
    goto failure;

#ifdef DEBUG_HESSIAN
_TestHessian ( domain, nlsprivate );
#endif

  if ( !g2h_SplNLNewtonf ( domain, nlsprivate ) )
    goto failure;

  g2h_ReflectVectorsf ( nfunc_a+nfunc_c+nfunc_d, nlprivate->acoeff,
                        nlprivate->acoeff );
  if ( acoeff )
    memcpy ( acoeff, nlprivate->acoeff,
             (nfunc_a+nfunc_c+nfunc_d)*sizeof(vector3f) );

  fc00 = pkv_GetScratchMemf ( (G2_CROSSDEGSUM+6)*2*hole_k*3 );
  if ( !fc00 ) {
    domain->error_code = G2H_ERROR_NO_SCRATCH_MEMORY;
    goto failure;
  }
  if ( !_g2h_SetSplRightSidef ( domain, 3, (float*)hole_cp, fc00, NULL ) )
    goto failure;

  if ( !_g2h_OutputSplPatchesf ( domain, 3, (float*)nlprivate->acoeff, fc00,
                                 usrptr, (outscf*)outpatch ) )
    goto failure;

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*g2h_NLSplFillHolef*/

/* ///////////////////////////////////////////////////////////////////////// */
static boolean g2h_SplNLConstrNewtonf ( GHoleDomainf *domain,
                                        G2HNLSPrivatef *nlsprivate, int nconstr,
                                        float *SCmat )
{
#define EPSF 2.0e-4
  void    *sp;
  G2HolePrivateRecf  *privateG2;
  G2HoleSPrivateRecf *sprivate;
  G2HNLPrivatef      *nlpr;
  int     itn, jtn, ktn, hole_k, nfunc_a, nfunc_c, nfunc_d, nfacd, nfunc;
  int     bs1, bs2, bs3, asize, esize, diagblsize, sideblsize, esideblsize;
  float   func, *grad, *coeff, *hessian, *E22, *cE22,
          *cT, *aa, *D1, *y, *y1, *M, *f, *hkk, *hki, *E22kk, *E22ki;
  float   func0, func1, aux, gn, gn0 = 0.0, gn1, dco, dyn, dyn1;
  int     i, j;
  boolean positive;

  sp = pkv_GetScratchMemTop ();
  hole_k = domain->hole_k;
  privateG2 = domain->privateG2;
  sprivate  = domain->SprivateG2;
  nlpr      = &nlsprivate->nlpr;
  nfunc_a = privateG2->nfunc_a;
  nfunc_c = sprivate->nsfunc_c;
  nfunc_d = sprivate->nsfunc_d;
  nfacd = nfunc_a+nfunc_c+nfunc_d;
  nfunc = nfacd-nconstr;

  bs1 = nfunc_c/hole_k;
  bs2 = nfunc_d+nfunc_a;
  bs3 = bs2-nconstr;
  diagblsize  = (bs1*(bs1+1))/2;
  sideblsize  = bs1*bs2;
  esideblsize = bs1*bs3;
  coeff   = pkv_GetScratchMemf ( nfacd );
  grad    = pkv_GetScratchMemf ( nfacd );
  asize   = pkn_Block1ArraySize ( hole_k, bs1, bs2 );
  esize   = pkn_Block1ArraySize ( hole_k, bs1, bs3 );
  hessian = pkv_GetScratchMemf ( asize );
  E22     = pkv_GetScratchMemf ( esize );
  cE22    = pkv_GetScratchMemf ( esize );
  cT      = pkv_GetScratchMemf ( bs2*nconstr );
  aa      = pkv_GetScratchMemf ( 2*nconstr );
  D1      = pkv_GetScratchMemf ( (nconstr*(nconstr+1))/2 );
  y       = pkv_GetScratchMemf ( nfacd );
  y1      = pkv_GetScratchMemf ( nfacd );
  M       = pkv_GetScratchMemf ( (bs2*(bs2+1))/2 );
  f       = pkv_GetScratchMemf ( nfacd );
  if ( !coeff || !grad || !hessian || !E22 || !cE22 || !cT || !aa ||
       !D1 || !y || !y1 || !M || !f )
    goto failure;
  hkk = &hessian[pkn_Block1FindBlockPos ( hole_k, bs1, bs2, hole_k, hole_k )];
  hki = &hessian[pkn_Block1FindBlockPos ( hole_k, bs1, bs2, hole_k, 0 )];
  E22kk = &E22[pkn_Block1FindBlockPos ( hole_k, bs1, bs3, hole_k, hole_k)];
  E22ki = &E22[pkn_Block1FindBlockPos ( hole_k, bs1, bs3, hole_k, 0 )];

        /* setup the initial point */
  pkv_Selectf ( nfacd, 1, 3, 1, &nlpr->acoeff[0].z, coeff );

        /* step 1: decompose the constraint equations matrix */
  pkv_TransposeMatrixf ( nconstr, bs2, nfacd, &SCmat[nfunc_c], nconstr, cT );
  pkn_QRDecomposeMatrixf ( bs2, nconstr, cT, aa );
  for ( i = 0; i < nconstr; i++ )
    for ( j = i; j < nconstr; j++ )
      D1[pkn_SymMatIndex(i,j)] = cT[i*nconstr+j];

        /* Newton iterations */
  func0 = 1.0e+38;
  for ( itn = ktn = 0; ; itn++ ) {
          /* step 2 */
    memcpy ( y, coeff, nfacd*sizeof(float) );
    pkn_multiReflectVectorf ( bs2, nconstr, cT, aa, 1, 1, &y[nfunc_c] );
    pkn_LowerTrMatrixSolvef ( nconstr, D1, 1, 3, &nlpr->rhole_cp[12*hole_k+1].z,
                              1, &y[nfunc_c] );
    for ( i = nfunc_c; i < nfunc_c+nconstr; i++ )
      y[i] = -y[i];

          /* step 3 */
    if ( !_g2h_ComputeSplNLFuncGradHessianf ( domain, nlsprivate, coeff,
                              &func, grad, hessian ) )
      goto failure;

    if ( !pkn_ComputeQTSQf ( bs2, hkk, nconstr, cT, aa, M ) )
      goto failure;
    for ( i = 0; i < bs3; i++ )
      for ( j = i; j < bs3; j++ )
        E22kk[pkn_SymMatIndex(i,j)] = M[pkn_SymMatIndex(nconstr+i,nconstr+j)];
    for ( i = 0; i < hole_k; i++ )
      pkn_multiReflectVectorf ( bs2, nconstr, cT, aa, bs1, bs1,
                                &hki[i*sideblsize] );
    pkv_Selectf ( hole_k, esideblsize, sideblsize, esideblsize,
                  &hki[nconstr*bs1], E22ki );
    memcpy ( E22, hessian, hole_k*diagblsize*sizeof(float) );
    memcpy ( cE22, E22, esize*sizeof(float) );

          /* step 4 */
    memcpy ( f, grad, nfacd*sizeof(float) );
    pkn_multiReflectVectorf ( bs2, nconstr, cT, aa, 1, 1, &f[nfunc_c] );
    memmove ( &f[nfunc_c], &f[nfunc_c+nconstr], bs3*sizeof(float) );
    gn = (float)sqrt ( pkn_ScalarProductf ( nfunc, f, f ) );
    if ( itn == 0 )
      gn0 = gn;

          /* step 5 */
    if ( (positive = pkn_Block1CholeskyDecompMf ( hole_k, bs1, bs3, E22 )) ) {
      pkn_Block1LowerTrMSolvef ( hole_k, bs1, bs3, E22, 1, 1, f );
      pkn_Block1UpperTrMSolvef ( hole_k, bs1, bs3, E22, 1, 1, f );
    }
    else {

printf ( "!" );

      if ( !pkn_Block1SymMatrixMultf ( hole_k, bs1, bs3, cE22, 1, 1, f, 1, y1 ) )
        goto failure;
      aux = (float)pkn_ScalarProductf ( nfunc, f, y1 );
      if ( aux <= 0.0 || aux < EPSF*gn )
        goto failure;
      pkn_MultMatrixNumf ( 1, nfunc, 1, f, gn/aux, 1, f );
    }

          /* step 6 */
    dyn = (float)sqrt ( pkn_ScalarProductf ( nfunc, f, f ) );
    memmove ( &f[nfunc_c+nconstr], &f[nfunc_c], bs3*sizeof(float) );
    dco = (float)sqrt ( pkn_ScalarProductf ( nfacd, coeff, coeff ) );
    for ( aux = 1.0; aux >= EPSF; aux *= 0.5 ) {
      for ( i = 0; i < nfunc_c; i++ )
        y1[i] = y[i]+aux*f[i];
      memcpy ( &y1[nfunc_c], &y[nfunc_c], nconstr*sizeof(float) );
      for ( i = nfunc_c+nconstr; i < nfacd; i++ )
        y1[i] = y[i]+aux*f[i];
      pkn_multiInvReflectVectorf ( bs2, nconstr, cT, aa, 1, 1, &y1[nfunc_c] );
      func1 = _g2h_ComputeSplNLFuncf ( domain, nlsprivate, y1 );
      if ( func1 < func )
        break;
    }
    memcpy ( coeff, y1, nfacd*sizeof(float) );
    func = func1;
    ktn ++;

    if ( positive && aux > 0.1 ) {
        /* Now the Hessian matrix is positive-definite; */
        /* as it is expensive to compute, we try to make some */
        /* extra iterations with the same Hessian. */

      for ( jtn = 0; jtn < 10; jtn++ ) {
            /* step 2' */
        memcpy ( y, coeff, nfacd*sizeof(float) );
        pkn_multiReflectVectorf ( bs2, nconstr, cT, aa, 1, 1, &y[nfunc_c] );
        pkn_LowerTrMatrixSolvef ( nconstr, D1, 1, 3, &nlpr->rhole_cp[12*hole_k+1].z,
                                  1, &y[nfunc_c] );
        for ( i = nfunc_c; i < nfunc_c+nconstr; i++ )
          y[i] = -y[i];
            /* step 3' */
        _g2h_ComputeSplNLFuncGradf ( domain, nlsprivate, coeff, &func0, grad );
            /* step 4' */
        memcpy ( f, grad, nfacd*sizeof(float) );
        pkn_multiReflectVectorf ( bs2, nconstr, cT, aa, 1, 1, &f[nfunc_c] );
        memmove ( &f[nfunc_c], &f[nfunc_c+nconstr], bs3*sizeof(float) );
        gn1 = (float)sqrt ( pkn_ScalarProductf ( nfunc, f, f ) );
            /* step 5' */
        pkn_Block1LowerTrMSolvef ( hole_k, bs1, bs3, E22, 1, 1, f );
        pkn_Block1UpperTrMSolvef ( hole_k, bs1, bs3, E22, 1, 1, f );
            /* step 6' */
        dyn1 = (float)sqrt ( pkn_ScalarProductf ( nfunc, f, f ) );
        memmove ( &f[nfunc_c+nconstr], &f[nfunc_c], bs3*sizeof(float) );
        for ( i = 0; i < nfunc_c; i++ )
          y1[i] = y[i]+f[i];
        memcpy ( &y1[nfunc_c], &y[nfunc_c], nconstr*sizeof(float) );
        for ( i = nfunc_c+nconstr; i < nfacd; i++ )
          y1[i] = y[i]+f[i];
        pkn_multiInvReflectVectorf ( bs2, nconstr, cT, aa, 1, 1, &y1[nfunc_c] );

        func1 = _g2h_ComputeSplNLFuncf ( domain, nlsprivate, y1 );
        if ( func1 >= func0 || gn1 >= 0.5*gn )
          break;
        memcpy ( coeff, y1, nfacd*sizeof(float) );
        func = func1;
        ktn ++;
        gn = gn1;
        aux = 1.0;
        dco = (float)sqrt ( pkn_ScalarProductf ( nfacd, coeff, coeff ) );
        dyn = dyn1;
      }
    }
    if ( _g2h_StopItf ( itn, gn0, gn, dco, dyn, aux ) )
      break;
  }

printf ( "itn = %d, ktn = %d\n", itn+1, ktn );
/*
printf ( "func = %f\n", func1 );
*/
  pkv_Selectf ( nfacd, 1, 1, 3, coeff, &nlpr->acoeff[0].z );

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
#undef EPSF
} /*g2h_SplNLConstrNewtonf*/

boolean g2h_NLSplFillHoleConstrf ( GHoleDomainf *domain, const point3f *hole_cp,
                     int nconstr, const vector3f *constr,
                     float *acoeff, void *usrptr,
                     void (*outpatch) ( int n, int lknu, const float *knu,
                                        int m, int lknv, const float *knv,
                                        const point3f *cp, void *usrptr ) )
{
  typedef void outscf ( int n, int lknu, const float *knu,
                        int m, int lknv, const float *knv,
                        const float *cp, void *usrptr );

  void  *sp;
  G2HNLSPrivatef     *nlsprivate;
  G2HNLPrivatef      *nlprivate;
  G2HolePrivateRecf  *privateG2;
  G2HoleSPrivateRecf *sprivate;
  int   hole_k, nfunc_a, nfunc_c, nfunc_d, nfacd;
  float *fc00;

  sp = pkv_GetScratchMemTop ();
  privateG2 = domain->privateG2;
  sprivate = domain->SprivateG2;
  if ( !sprivate )
    goto failure;
  if ( !sprivate->SLMat )
    if ( !g2h_DecomposeSplMatrixf ( domain ) )
      goto failure;
  nlsprivate = _g2h_AllocSplNLPrf ( domain );
  if ( !nlsprivate )
    goto failure;
  _g2h_nlprivf = nlprivate = &nlsprivate->nlpr;
  if ( !_g2h_InitSplNLprf ( domain, nlsprivate, nconstr ) )
    goto failure;
  hole_k = domain->hole_k;
  nfunc_a = privateG2->nfunc_a;
  nfunc_c = sprivate->nsfunc_c;
  nfunc_d = sprivate->nsfunc_d;
  nfacd = nfunc_a+nfunc_c+nfunc_d;

  if ( !_g2h_ComputeNLNormalf ( domain, nlprivate, hole_cp ) )
    goto failure;
  g2h_ReflectVectorsf ( 12*hole_k+1, hole_cp, nlprivate->rhole_cp );
  g2h_ReflectVectorsf ( nconstr, constr, &nlprivate->rhole_cp[12*hole_k+1] );

  nlprivate->auxc = 0;
  if ( !g2h_SplFillHoleConstrf ( domain, 3, (float*)nlprivate->rhole_cp,
                    nconstr, (float*)&nlprivate->rhole_cp[12*hole_k+1],
                    (float*)nlprivate->acoeff, NULL, g2h_splnloutpatchf ) )
    goto failure;

  if ( !_g2h_TabSplNLBasisFunctionsf ( domain, nlsprivate ) )
    goto failure;

  if ( !g2h_SplNLConstrNewtonf ( domain, nlsprivate, nconstr, sprivate->SCmat ) )
    goto failure;

  g2h_ReflectVectorsf ( nfacd, nlprivate->acoeff, nlprivate->acoeff );
  if ( acoeff )
    memcpy ( acoeff, nlprivate->acoeff, nfacd*sizeof(vector3f) );

  fc00 = pkv_GetScratchMemf ( (G2_CROSSDEGSUM+6)*2*hole_k*3 );
  if ( !fc00 )
    goto failure;
  if ( !_g2h_SetSplRightSidef ( domain, 3, (float*)hole_cp, fc00, NULL ) )
    goto failure;

  if ( !_g2h_OutputSplPatchesf ( domain, 3, (float*)nlprivate->acoeff, fc00,
                                 usrptr, (outscf*)outpatch ) )
    goto failure;

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*g2h_NLSplFillHoleConstrf*/

/* ///////////////////////////////////////////////////////////////////////// */
static boolean _g2h_ReflectSplAltConstrMatrixf ( GHoleDomainf *domain,
                                          vector3f *reflv, float *RACmat )
{
  void   *sp;
  G2HolePrivateRecf  *privateG2;
  G2HoleSPrivateRecf *sprivate;
  int    i, j, nconstr, nfunc_a, nfunc_c, nfunc_d, nbf;
  float *buf;

  sp = pkv_GetScratchMemTop ();
  privateG2 = domain->privateG2;
  sprivate = domain->SprivateG2;
  nfunc_a = privateG2->nfunc_a;
  nfunc_c = sprivate->nsfunc_c;
  nfunc_d = sprivate->nsfunc_d;
  nbf = nfunc_a+nfunc_c+nfunc_d;
  nconstr = sprivate->splnaconstr;
  memcpy ( RACmat, sprivate->ASCmat, nconstr*nbf*3*sizeof(float) );
  buf = pkv_GetScratchMemf ( nbf );
  if ( !buf ) {
    domain->error_code = G2H_ERROR_NO_SCRATCH_MEMORY;
    goto failure;
  }
  for ( i = 0; i < nconstr; i++ ) {
    pkn_MultMatrixNumf ( 1, nbf, 0, RACmat, reflv->x, 0, buf );
    pkn_AddMatrixMf ( 1, nbf, 0, buf, 0, &RACmat[nbf], reflv->y, 0, buf );
    pkn_AddMatrixMf ( 1, nbf, 0, buf, 0, &RACmat[2*nbf], reflv->z, 0, buf );
    for ( j = 0; j < nbf; j++ ) {
      RACmat[j]            -= (float)(2.0*buf[j]*reflv->x);
      RACmat[nbf+j]   -= (float)(2.0*buf[j]*reflv->y);
      RACmat[2*nbf+j] -= (float)(2.0*buf[j]*reflv->z);
    }
    RACmat = &RACmat[3*nbf];
  }
  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*_g2h_ReflectSplAltConstrMatrixf*/

boolean g2h_NLSplFillHoleAltConstrf ( GHoleDomainf *domain, const point3f *hole_cp,
                     int naconstr, const float *constr,
                     float *acoeff, void *usrptr,
                     void (*outpatch) ( int n, int lknu, const float *knu,
                                        int m, int lknv, const float *knv,
                                        const point3f *cp, void *usrptr ) )
{
  typedef void outscf ( int n, int lknu, const float *knu,
                        int m, int lknv, const float *knv,
                        const float *cp, void *usrptr );

  void    *sp;
  G2HNLPrivatef      *nlprivate;
  G2HNLSPrivatef     *nlsprivate;
  G2HolePrivateRecf  *privateG2;
  G2HoleSPrivateRecf *sprivate;
  float  *fc00, *cmat, *saveconstr, *newconstr, *rsconstr;
  int     hole_k, nfunc_a, nfunc_c, nfunc_d, nbf, i, j;
  boolean restore;

  sp = pkv_GetScratchMemTop ();
  restore = false;
  privateG2 = domain->privateG2;
  sprivate = domain->SprivateG2;
  if ( !sprivate )
    goto failure;
  if ( !sprivate->SLMat )
    if ( !g2h_DecomposeSplMatrixf ( domain ) )
      goto failure;
  nlsprivate = _g2h_AllocSplNLPrf ( domain );
  if ( !nlsprivate ) {
    domain->error_code = G2H_ERROR_NO_SCRATCH_MEMORY;
    goto failure;
  }
  _g2h_nlprivf = nlprivate = &nlsprivate->nlpr;
  if ( !_g2h_InitSplNLprf ( domain, nlsprivate, naconstr ) )
    goto failure;
  if ( naconstr != sprivate->splnaconstr || sprivate->splacdim != 3 ) {
    domain->error_code = G2H_ERROR_INCONSISTENT_CONSTR;
    goto failure;
  }
  hole_k = domain->hole_k;
  nfunc_a = privateG2->nfunc_a;
  nfunc_c = sprivate->nsfunc_c;
  nfunc_d = sprivate->nsfunc_d;
  nbf = nfunc_a+nfunc_c+nfunc_d;
  saveconstr = pkv_GetScratchMemf ( 7*nbf*naconstr );
  if ( !saveconstr ) {
    domain->error_code = G2H_ERROR_NO_SCRATCH_MEMORY;
    goto failure;
  }
  newconstr = &saveconstr[3*nbf*naconstr];
  cmat = &newconstr[3*nbf*naconstr];
  memcpy ( saveconstr, sprivate->ASCmat, nbf*naconstr*3*sizeof(float) );
  memcpy ( newconstr, sprivate->ASCmat, nbf*naconstr*3*sizeof(float) );
  if ( !_g2h_ComputeNLNormalf ( domain, nlprivate, hole_cp ) )
    goto failure;
  g2h_ReflectVectorsf ( 12*hole_k+1, hole_cp, nlprivate->rhole_cp );
  if ( !_g2h_ReflectSplAltConstrMatrixf ( domain, &nlprivate->reflv, newconstr ) )
    goto failure;

  restore = true;
  if ( !g2h_SetSplAltConstraintMatrixf ( domain, 3, naconstr, newconstr ) )
    goto failure;
  nlprivate->auxc = 0;
  if ( !g2h_SplFillHoleAltConstrf ( domain, 3, (float*)nlprivate->rhole_cp,
                    naconstr, constr,
                    (float*)nlprivate->acoeff, NULL, g2h_splnloutpatchf ) )
    goto failure;

  pkv_Selectf ( naconstr, nbf, 3*nbf, nbf, &newconstr[2*nbf], cmat );
  rsconstr = &nlprivate->rhole_cp[12*hole_k+1].z;
  pkv_Selectf ( naconstr, 1, 1, 3, constr, rsconstr );
  for ( j = 0; j < naconstr; j++ )
    for ( i = 0; i < nbf; i++ )
      rsconstr[3*j] += newconstr[3*j*nbf+i]*nlprivate->acoeff[i].x +
                       newconstr[(3*j+1)*nbf+i]*nlprivate->acoeff[i].y;

  if ( !_g2h_TabSplNLBasisFunctionsf ( domain, nlsprivate ) )
    goto failure;

  if ( !g2h_SplNLConstrNewtonf ( domain, nlsprivate, naconstr, cmat ) )
    goto failure;

  g2h_ReflectVectorsf ( nbf, nlprivate->acoeff, nlprivate->acoeff );
  if ( acoeff )
    memcpy ( acoeff, nlprivate->acoeff, nbf*sizeof(vector3f) );

  fc00 = pkv_GetScratchMemf ( (G2_CROSSDEGSUM+6)*2*hole_k*3 );
  if ( !fc00 )
    goto failure;
  if ( !_g2h_SetSplRightSidef ( domain, 3, (float*)hole_cp, fc00, NULL ) )
    goto failure;

  if ( !_g2h_OutputSplPatchesf ( domain, 3, (float*)nlprivate->acoeff, fc00,
                                 usrptr, (outscf*)outpatch ) )
    goto failure;

  g2h_SetSplAltConstraintMatrixf ( domain, 3, naconstr, saveconstr );
  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  if ( restore )
    g2h_SetSplAltConstraintMatrixf ( domain, 3, naconstr, saveconstr );
  pkv_SetScratchMemTop ( sp );
  return false;
} /*g2h_NLSplFillHoleAltConstrf*/

