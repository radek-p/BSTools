
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2008, 2015                            */
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

/* ///////////////////////////////////////////////////////////////////////// */
static G1HNLSPrivated *_g1hq2_AllocSplNLPrd ( GHoleDomaind *domain )
{
  G1HolePrivateRecd  *privateG1;
  G1HoleSPrivateRecd *sprivateG1;
  G1HNLSPrivated     *nlspr;
  int     hole_k, nk, m1, m2, rr, nkn, psize;
  int     nfunc_a, nfunc_c, nfunc_d;
  int     *fkn, *lkn, *cfuncvi, *dfuncvi;
  double  *tkn;
  int     i, j, k, fn, nr, nc, njcurves, jtabsize;
  int     i0, i1, j0, j1, nzc;
  boolean jumpC, jumpD;

  if ( (nlspr = pkv_GetScratchMem ( sizeof(G1HNLSPrivated) )) ) {
    privateG1 = domain->privateG1;
    sprivateG1 = domain->SprivateG1;
    hole_k = domain->hole_k;
    nk = sprivateG1->nk;
    m1 = sprivateG1->m1;
    m2 = sprivateG1->m2;
    nfunc_a = privateG1->nfunc_a;
    nfunc_c = sprivateG1->nsfunc_c;
    nfunc_d = sprivateG1->nsfunc_d;
    rr = G1H_FINALDEG-3+nk*m2;
    nkn = nlspr->nkn = (nk+1)*G1_QUAD_FACTOR;
    tkn = nlspr->tkn = pkv_GetScratchMemd ( nkn );
    fkn = nlspr->fkn = pkv_GetScratchMem ( 2*rr*sizeof(int) );
    nlspr->cb = pkv_GetScratchMemd ( 4*rr*nkn );
    cfuncvi = nlspr->cfuncvi = pkv_GetScratchMem ( nfunc_c*sizeof(int) );
    dfuncvi = nlspr->dfuncvi = pkv_GetScratchMem ( nfunc_d*sizeof(int) );
    if ( !tkn || !fkn || !nlspr->cb || !cfuncvi || !dfuncvi )
      return NULL;
    lkn = nlspr->lkn = &fkn[rr];
    nlspr->cbt = &nlspr->cb[rr*nkn];
    nlspr->cbtt = &nlspr->cbt[rr*nkn];
    nlspr->cbttt = &nlspr->cbtt[rr*nkn];
    _gh_PrepareTabKnotsd ( nkn, privateG1->opt_quad, tkn );
    psize = sprivateG1->lastfpknot-G1H_FINALDEG;
    nlspr->psize = psize*psize;

        /* compute the numbers of quadrature knots for each C block function */
    _g1h_TabBSFuncDer3d ( G1H_FINALDEG,
            sprivateG1->lastcknot, sprivateG1->cknots,
            2, G1H_FINALDEG+nk*m2-2, nkn, tkn, nlspr->fkn, nlspr->lkn,
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
    rr = 2*nk*m1;  /* == nfunc_d/hole_k */
    for ( fn = 0;  fn < rr;  fn++ ) {
      _g1h_FuncDSuppd ( hole_k, nk, m1, fn, 0, &nzc, &i0, &i1, &j0, &j1 );
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

        /* allocate the arrays for the jump samples */
    nlspr->jumpC = jumpC = (boolean)(G1H_FINALDEG-m2 < 2);
    nlspr->jumpD = jumpD = (boolean)(G1_CROSS00DEG-G1_BF01DEG-m1 < 1);
    if ( jumpC || jumpD )
      nlspr->njcurves = (short)(njcurves = 3+2*nk);
    else
      nlspr->njcurves = (short)(njcurves = 3);
    nlspr->nlpr.ctang = pkv_GetScratchMem ( njcurves*nkn*hole_k*sizeof(vector2d) );
    k = njcurves*(nfunc_a+1)*hole_k;
    nlspr->jcfs = nkn*k;
    k += (njcurves+1)*nfunc_c;
    nlspr->jdfs = nkn*k;
    k += (2*njcurves+1)*nfunc_d;
    k *= nkn;
    nlspr->jtabsize = jtabsize = 5*k;
    nlspr->nlpr.cpsiu = pkv_GetScratchMemd ( jtabsize );
    if ( !nlspr->nlpr.ctang || !nlspr->nlpr.cpsiu )
      return NULL;
    nlspr->nlpr.cpsiv = &nlspr->nlpr.cpsiu[k];
    nlspr->nlpr.cpsiuu = &nlspr->nlpr.cpsiv[k];
    nlspr->nlpr.cpsiuv = &nlspr->nlpr.cpsiuu[k];
    nlspr->nlpr.cpsivv = &nlspr->nlpr.cpsiuv[k];
    return nlspr;
  }
  else
    return NULL;
} /*_g1hq2_AllocSplNLPrd*/

static boolean _g1hq2_InitSplNLprd ( GHoleDomaind *domain,
                                     G1HNLSPrivated *nlspr,
                                     int nconstr )
{
  G1HolePrivateRecd  *privateG1;
  G1HoleSPrivateRecd *sprivateG1;
  G1HNLPrivated *nlpr;
  int hole_k, nfunc_a, nfunc_c, nfunc_d, nk;
  int nfunc, nkn, nkn2, ftabsize;

  privateG1 = domain->privateG1;
  sprivateG1 = domain->SprivateG1;
  hole_k = domain->hole_k;
  nlpr = &nlspr->nlpr;
  nfunc_a = privateG1->nfunc_a;
  nfunc_c = sprivateG1->nsfunc_c;
  nfunc_d = sprivateG1->nsfunc_d;
  nfunc = nfunc_a+nfunc_c+nfunc_d;
  ftabsize = nlspr->ftabsize;
  nk = sprivateG1->nk;
  nkn = (nk+1)*G1_QUAD_FACTOR;
  nkn2 = nkn*nkn;
  nlpr->auxc = 0;
  nlpr->nldi = pkv_GetScratchMem ( nlspr->psize*hole_k*sizeof(point3d) );
  nlpr->acoeff = pkv_GetScratchMem ( nfunc*sizeof(vector3d) );
  nlpr->rhole_cp = pkv_GetScratchMem ( (12*hole_k+1+nconstr)*sizeof(vector3d) );
  nlpr->diu = pkv_GetScratchMem ( 9*nkn2*hole_k*sizeof(vector2d) );
  nlpr->jac = pkv_GetScratchMemd ( nkn2*hole_k );
  nlpr->psiu = pkv_GetScratchMemd ( 9*ftabsize );

  if ( !nlpr->nldi || !nlpr->acoeff || !nlpr->rhole_cp ||
       !nlpr->diu || !nlpr->jac || !nlpr->psiu ) {
    domain->error_code = G1H_ERROR_NO_SCRATCH_MEMORY;
    return false;
  }

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

/* arrays for second derivative jumps are allocated by _g1hq2_AllocSplNLPrd */

  return true;
} /*_g1hq2_InitSplNLprd*/

static boolean _g1hq2_TabSplNLBasisFunctionsd ( GHoleDomaind *domain,
                                                G1HNLSPrivated *nlspr )
{
  void *sp;
  G1HNLPrivated      *nlpr;
  GHolePrivateRecd   *privateG;
  G1HolePrivateRecd  *privateG1;
  G1HoleSPrivateRecd *sprivate;
  int      hole_k, nk, m1, m2;
  int      nfunc_a, nfunc_b, nfunc_c, nfabc;
  int      i, j, k, l, f, rr, fn, nkn, nkn2, kN, kNQ2,
           i0, i1, j0, j1, ik, jk, nzc, sn, dn, ii, jj, fN;
  int      *fkn, *lkn;
  point2d  ddi;
  vector2d *di, *diu, *div, *diuu, *diuv, *divv, *diuuu, *diuuv, *diuvv, *divvv;
  double   *fc00, *fc01, *fc10, *fc11, *fd00, *fd01, *fd10, *fd11,
           *bc00, *bc01, *bc10, *bc11, *bd00, *bd01, *bd10, *bd11;
  double   *psiu, *psiv, *psiuu, *psiuv, *psivv,
           *psiuuu, *psiuuv, *psiuvv, *psivvv;
  double   *cb, *cbt, *cbtt, *cbttt;
  int      lastomcknot, lastpvknot, lastfpknot, disize;
  double   *omcknots, *pvknots, *fcomc, *fpv, *fpu, *tkn;
  double   *hfunc, *dhfunc, *ddhfunc, *dddhfunc;
  double   p, pu, puu, puuu, q, qv, qvv, qvvv, pq[9], bvz;
  double   p0, p0u, p0uu, p0uuu, p1, p1u, p1uu, p1uuu;
  unsigned char *bfcpn;
  double   wsp[4];

  sp = pkv_GetScratchMemTop ();
  hole_k = domain->hole_k;
  nlpr = &nlspr->nlpr;
  if ( !_g1hq2_FindNLDomainDiameterd ( domain, nlpr ) )
    goto failure;
  privateG  = domain->privateG;
  privateG1 = domain->privateG1;
  sprivate  = domain->SprivateG1;
  nk = sprivate->nk;
  m1 = sprivate->m1;
  m2 = sprivate->m2;
  nkn = nlspr->nkn;
  nkn2 = nkn*nkn;
  tkn = nlspr->tkn;
  nfunc_a = privateG1->nfunc_a;
  nfunc_b = privateG1->nfunc_b;
  nfunc_c = sprivate->nsfunc_c;
  nfabc = nfunc_a+nfunc_b+nfunc_c;
  hfunc = pkv_GetScratchMemd ( 16*nkn );
  lastfpknot = sprivate->lastfpknot;
  disize = lastfpknot-G1H_FINALDEG;
  disize *= disize;
  di = pkv_GetScratchMem ( disize*sizeof(point2d) );
  if ( !hfunc || !di ) {
    domain->error_code = G1H_ERROR_NO_SCRATCH_MEMORY;
    goto failure;
  }
  dhfunc = &hfunc[4*nkn];
  ddhfunc = &dhfunc[4*nkn];
  dddhfunc = &ddhfunc[4*nkn];
  if ( !mbs_TabCubicHFuncDer3d ( 0.0, 1.0, nkn, tkn, hfunc, dhfunc, ddhfunc, dddhfunc ) )
    goto failure;

        /* compute the Jacobian */
  for ( k = kN = kNQ2 = 0;
        k < hole_k;
        k++, kN += disize, kNQ2 += nkn2 ) {
    pkv_Selectd ( disize, 2, 3, 2, &nlpr->nldi[kN].x, &di[0].x );
    for ( i = l = 0;  i < nkn;  i++ )
          /* this is not an optimal algorithm, to be improved */
      for ( j = 0;  j < nkn;  j++, l++ ) {
        mbs_deBoorDer3Pd ( G1H_FINALDEG, lastfpknot, sprivate->fpknots,
                   G1H_FINALDEG, lastfpknot, sprivate->fpknots, 2,
                   2*(lastfpknot-G1H_FINALDEG), &di[0].x,
                   tkn[i], tkn[j], (double*)&ddi,
                   (double*)&nlpr->diu[kNQ2+l], (double*)&nlpr->div[kNQ2+l],
                   (double*)&nlpr->diuu[kNQ2+l], (double*)&nlpr->diuv[kNQ2+l],
                   (double*)&nlpr->divv[kNQ2+l], (double*)&nlpr->diuuu[kNQ2+l],
                   (double*)&nlpr->diuuv[kNQ2+l], (double*)&nlpr->diuvv[kNQ2+l],
                   (double*)&nlpr->divvv[kNQ2+l] );
        nlpr->jac[kNQ2+l] = (double)det2d ( &nlpr->diu[kNQ2+l],
                                           &nlpr->div[kNQ2+l] );
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
    for ( k = kNQ2 = 0;  k < hole_k;  k++, fN += nkn2, kNQ2 += nkn2 ) {
      _g1h_GetBFAPatchCurvesd ( domain, f, k, &fc00, &fc01, &fd00, &fd01 );
      _g1hq2_TabNLDer0d ( nkn, tkn, hfunc, dhfunc, ddhfunc, dddhfunc,
                          &nlpr->diu[kNQ2], &nlpr->div[kNQ2],
                          &nlpr->diuu[kNQ2], &nlpr->diuv[kNQ2],
                          &nlpr->divv[kNQ2], &nlpr->diuuu[kNQ2],
                          &nlpr->diuuv[kNQ2], &nlpr->diuvv[kNQ2],
                          &nlpr->divvv[kNQ2],
                          fc00, fc01, fd00, fd01,
                          &nlpr->psiu[fN], &nlpr->psiv[fN],
                          &nlpr->psiuu[fN], &nlpr->psiuv[fN], &nlpr->psivv[fN],
                          &nlpr->psiuuu[fN], &nlpr->psiuuv[fN],
                          &nlpr->psiuvv[fN], &nlpr->psivvv[fN]);
    }
  }
        /* the fixed combination of the block B functions */
  bc00 = pkv_GetScratchMemd ( 2*(G1_CROSSDEGSUM+4)*hole_k );
  if ( !bc00 ) {
    domain->error_code = G1H_ERROR_NO_SCRATCH_MEMORY;
    goto failure;
  }
  bc01 = &bc00[(G1_CROSS00DEG+1)*hole_k];  bc10 = &bc01[(G1_CROSS01DEG+1)*hole_k];
  bc11 = &bc10[(G1_CROSS10DEG+1)*hole_k];  bd00 = &bc11[(G1_CROSS11DEG+1)*hole_k];
  bd01 = &bd00[(G1_CROSS00DEG+1)*hole_k];  bd10 = &bd01[(G1_CROSS01DEG+1)*hole_k];
  bd11 = &bd10[(G1_CROSS10DEG+1)*hole_k];

  memset ( bc00, 0, 2*(G1_CROSSDEGSUM+4)*hole_k*sizeof(double) );
  bfcpn = privateG->bfcpn;
  for ( f = 0; f < nfunc_b; f++ ) {
        /* find the Coons representation of the constant part of the solution */
    bvz = nlpr->rhole_cp[bfcpn[f]].z;
    for ( k = 0; k < hole_k; k++ ) {
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
    _g1hq2_TabNLDerd ( nkn, tkn, hfunc, dhfunc, ddhfunc, dddhfunc,
             &nlpr->diu[kNQ2], &nlpr->div[kNQ2],
             &nlpr->diuu[kNQ2], &nlpr->diuv[kNQ2], &nlpr->divv[kNQ2],
             &nlpr->diuuu[kNQ2], &nlpr->diuuv[kNQ2],
             &nlpr->diuvv[kNQ2], &nlpr->divvv[kNQ2],
             &bc00[(G1_CROSS00DEG+1)*k], &bc01[(G1_CROSS01DEG+1)*k],
             &bc10[(G1_CROSS10DEG+1)*k], &bc11[(G1_CROSS11DEG+1)*k],
             &bd00[(G1_CROSS00DEG+1)*k], &bd01[(G1_CROSS01DEG+1)*k],
             &bd10[(G1_CROSS10DEG+1)*k], &bd11[(G1_CROSS11DEG+1)*k],
             &nlpr->psiu[fN], &nlpr->psiv[fN],
             &nlpr->psiuu[fN], &nlpr->psiuv[fN], &nlpr->psivv[fN],
             &nlpr->psiuuu[fN], &nlpr->psiuuv[fN],
             &nlpr->psiuvv[fN], &nlpr->psivvv[fN] );
  }
        /* now the block C functions */
  rr    = G1H_FINALDEG-3+nk*m2;
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
            _g2h_TensDer3d ( p, pu, puu, puuu, q, qv, qvv, qvvv, pq );
            if ( !_pkn_Comp2iDerivatives3d ( diu[dn].x, diu[dn].y, div[dn].x, div[dn].y,
                    diuu[dn].x, diuu[dn].y, diuv[dn].x, diuv[dn].y,
                    divv[dn].x, divv[dn].y, diuuu[dn].x, diuuu[dn].y,
                    diuuv[dn].x, diuuv[dn].y, diuvv[dn].x, diuvv[dn].y,
                    divvv[dn].x, divvv[dn].y,
                    1, &pq[0], &pq[1], &pq[2], &pq[3], &pq[4], &pq[5], &pq[6],
                    &pq[7], &pq[8],
                    &psiu[sn], &psiv[sn], &psiuu[sn], &psiuv[sn], &psivv[sn],
                    &psiuuu[sn], &psiuuv[sn], &psiuvv[sn], &psivvv[sn], wsp ) )
              goto failure;
          }
        }
      }
    }
  }
        /* and the block D functions */
  lastomcknot = sprivate->lastomcknot;
  omcknots = sprivate->omcknots;
  lastpvknot = sprivate->lastpvknot;
  pvknots = sprivate->pvknots;
  fcomc = pkv_GetScratchMemd ( lastomcknot+2*lastpvknot-
                               (G1_CROSS00DEG+2*G1_CROSS01DEG) );
  if ( !fcomc ) {
    domain->error_code = G1H_ERROR_NO_SCRATCH_MEMORY;
    goto failure;
  }
  fpv = &fcomc[lastomcknot-G1_CROSS00DEG];
  fpu = &fpv[lastpvknot-G1_CROSS01DEG];
  rr = 2*nk*m1;  /* == nfunc_d/hole_k */
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
      _g1h_FuncDSuppd ( hole_k, nk, m1, fn, i, &nzc, &i0, &i1, &j0, &j1 );
      _g1h_GetSplDBasisCrossDerd ( domain, nfabc+fn, i, fcomc, fpv, fpu );
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
          mbs_deBoorDer3C1d ( G1_CROSS00DEG, lastomcknot, omcknots, fcomc,
                              tkn[ik], &p0, &p0u, &p0uu, &p0uuu );
        mbs_deBoorDer3C1d ( G1_CROSS01DEG, lastpvknot, pvknots, fpv,
                            tkn[ik], &p1, &p1u, &p1uu, &p1uuu );
        for ( jk = 0, jj = 4*j0;  jk <= j1;  jk++, jj += 4, sn++, dn++ ) {
          if ( nzc == 0 ) {
            pq[0] = p0u*hfunc[jj]   + p1u*hfunc[jj+2];
            pq[1] = p0*dhfunc[jj]   + p1*dhfunc[jj+2];
            pq[2] = p0uu*hfunc[jj]  + p1uu*hfunc[jj+2];
            pq[3] = p0u*dhfunc[jj]  + p1u*dhfunc[jj+2];
            pq[4] = p0*ddhfunc[jj]  + p1*ddhfunc[jj+2];
            pq[5] = p0uuu*hfunc[jj] + p1uuu*hfunc[jj+2];
            pq[6] = p0uu*dhfunc[jj] + p1uu*dhfunc[jj+2];
            pq[7] = p0u*ddhfunc[jj] + p1u*ddhfunc[jj+2];
            pq[8] = p0*dddhfunc[jj] + p1*dddhfunc[jj+2];
          }
          else {
            pq[0] = p1u*hfunc[jj+2];
            pq[1] = p1*dhfunc[jj+2];
            pq[2] = p1uu*hfunc[jj+2];
            pq[3] = p1u*dhfunc[jj+2];
            pq[4] = p1*ddhfunc[jj+2];
            pq[5] = p1uuu*hfunc[jj+2];
            pq[6] = p1uu*dhfunc[jj+2];
            pq[7] = p1u*ddhfunc[jj+2];
            pq[8] = p1*dddhfunc[jj+2];
          }
          if ( !_pkn_Comp2iDerivatives3d ( diu[dn].x, diu[dn].y, div[dn].x, div[dn].y,
                     diuu[dn].x, diuu[dn].y, diuv[dn].x, diuv[dn].y,
                     divv[dn].x, divv[dn].y, diuuu[dn].x, diuuu[dn].y,
                     diuuv[dn].x, diuuv[dn].y, diuvv[dn].x, diuvv[dn].y,
                     divvv[dn].x, divvv[dn].y,
                     1, &pq[0], &pq[1], &pq[2], &pq[3], &pq[4], &pq[5], &pq[6],
                     &pq[7], &pq[8],
                     &psiu[sn], &psiv[sn], &psiuu[sn], &psiuv[sn], &psivv[sn],
                     &psiuuu[sn], &psiuuv[sn], &psiuvv[sn], &psivvv[sn], wsp ) )
            goto failure;
        }
      }
    }
  }
          /* and then in Omega_{i-1}*/
  for ( ii = hole_k-1, i = fn = 0;  i < hole_k;  ii = i++ ) {
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
      _g1h_FuncDSuppd ( hole_k, nk, m1, fn, ii, &nzc, &i0, &i1, &j0, &j1 );
      _g1h_GetSplDBasisCrossDerd ( domain, nfabc+fn, i, fcomc, fpv, fpu );
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
          mbs_deBoorDer3C1d ( G1_CROSS00DEG, lastomcknot, omcknots, fcomc,
                              tkn[jk], &p0, &p0u, &p0uu, &p0uuu );
        mbs_deBoorDer3C1d ( G1_CROSS01DEG, lastpvknot, pvknots, fpu,
                            tkn[jk], &p1, &p1u, &p1uu, &p1uuu );
        for ( ik = i0, jj = 4*i0, sn = jk-j0, dn = i0*nkn+jk;
              ik <= i1;
              ik++, jj += 4, sn += (j1-j0+1), dn += nkn ) {
          if ( nzc == 0 ) {
            pq[0] = p0*dhfunc[jj]   + p1*dhfunc[jj+2];
            pq[1] = p0u*hfunc[jj]   + p1u*hfunc[jj+2];
            pq[2] = p0*ddhfunc[jj]  + p1*ddhfunc[jj+2];
            pq[3] = p0u*dhfunc[jj]  + p1u*dhfunc[jj+2];
            pq[4] = p0uu*hfunc[jj]  + p1uu*hfunc[jj+2];
            pq[5] = p0*dddhfunc[jj] + p1*dddhfunc[jj+2];
            pq[6] = p0u*ddhfunc[jj] + p1u*ddhfunc[jj+2];
            pq[7] = p0uu*dhfunc[jj] + p1uu*dhfunc[jj+2];
            pq[8] = p0uuu*hfunc[jj] + p1uuu*hfunc[jj+2];
          }
          else {
            pq[0] = p1*dhfunc[jj+2];
            pq[1] = p1u*hfunc[jj+2];
            pq[2] = p1*ddhfunc[jj+2];
            pq[3] = p1u*dhfunc[jj+2];
            pq[4] = p1uu*hfunc[jj+2];
            pq[5] = p1*dddhfunc[jj+2];
            pq[6] = p1u*ddhfunc[jj+2];
            pq[7] = p1uu*dhfunc[jj+2];
            pq[8] = p1uuu*hfunc[jj+2];
          }
          if ( !_pkn_Comp2iDerivatives3d ( diu[dn].x, diu[dn].y, div[dn].x, div[dn].y,
                     diuu[dn].x, diuu[dn].y, diuv[dn].x, diuv[dn].y,
                     divv[dn].x, divv[dn].y, diuuu[dn].x, diuuu[dn].y,
                     diuuv[dn].x, diuuv[dn].y, diuvv[dn].x, diuvv[dn].y,
                     divvv[dn].x, divvv[dn].y,
                     1, &pq[0], &pq[1], &pq[2], &pq[3], &pq[4], &pq[5], &pq[6],
                     &pq[7], &pq[8],
                     &psiu[sn], &psiv[sn], &psiuu[sn], &psiuv[sn], &psivv[sn],
                     &psiuuu[sn], &psiuuv[sn], &psiuvv[sn], &psivvv[sn], wsp ) )
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
} /*_g1hq2_TabSplNLBasisFunctionsd*/

static void _g1hq2_TabSplCJumpd ( int njcurves, int nk, int nkn,
             int fkni, int lkni, const double *tbsi, const double *tbsti,
             const double *tbstti,
             const double *atbsi, const double *atbsti,
             const double *atbstt0i, const double *atbstt1i,
             int fknj, int lknj, const double *tbsj, const double *tbstj,
             const double *tbsttj,
             const double *atbsj, const double *atbstj,
             const double *atbstt0j, const double *atbstt1j,
             double *ctrd, double *ctrdu0,
             double *cpsiu, double *cpsiv,
             double *cpsiuu, double *cpsiuv, double *cpsivv )
{
  int    i, j;
  double der[5];
  double *A11, *A21, *A22;

        /* curve v == 0 */
  if ( atbstt1j[0] )
    for ( i = fkni; i <= lkni; i++ ) {
      _g1h_TensDer2d ( tbsi[i], tbsti[i], tbstti[i],
                       atbsj[0], atbstj[0], atbstt1j[0], der );
      A11 = &ctrd[38*i];  A21 = &A11[4];  A22 = &A21[6];
      cpsiu[i] = A11[0]*der[0] + A11[1]*der[1];
      cpsiv[i] = A11[2]*der[0] + A11[3]*der[1];
      cpsiuu[i] = A21[0]*der[0] + A21[1]*der[1] +
                  A22[0]*der[2] + A22[1]*der[3] + A22[2]*der[4];
      cpsiuv[i] = A21[2]*der[0] + A21[3]*der[1] +
                  A22[3]*der[2] + A22[4]*der[3] + A22[5]*der[4];
      cpsivv[i] = A21[4]*der[0] + A21[5]*der[1] +
                  A22[6]*der[2] + A22[7]*der[3] + A22[8]*der[4];
    }
        /* curve v == 1 */
  if ( atbstt0j[nk+1] )
    for ( i = fkni; i <= lkni; i++ ) {
      _g1h_TensDer2d ( tbsi[i], tbsti[i], tbstti[i],
                       atbsj[nk+1], atbstj[nk+1], atbstt0j[nk+1], der );
      A11 = &ctrd[38*(nkn+i)];  A21 = &A11[4];  A22 = &A21[6];
      cpsiu[nkn+i] = -(A11[0]*der[0] + A11[1]*der[1]);
      cpsiv[nkn+i] = -(A11[2]*der[0] + A11[3]*der[1]);
      cpsiuu[nkn+i] = -(A21[0]*der[0] + A21[1]*der[1] +
                        A22[0]*der[2] + A22[1]*der[3] + A22[2]*der[4]);
      cpsiuv[nkn+i] = -(A21[2]*der[0] + A21[3]*der[1] +
                        A22[3]*der[2] + A22[4]*der[3] + A22[5]*der[4]);
      cpsivv[nkn+i] = -(A21[4]*der[0] + A21[5]*der[1] +
                        A22[6]*der[2] + A22[7]*der[3] + A22[8]*der[4]);
    }
        /* curve u == 0 */
  if ( atbstt1i[0] )
    for ( i = fknj; i <= lknj; i++ ) {
      _g1h_TensDer2d ( atbsi[0], atbsti[0], atbstt1i[0],
                       tbsj[i], tbstj[i], tbsttj[i], der );
      A11 = &ctrdu0[38*i+19];  A21 = &A11[4];  A22 = &A21[6];
      cpsiu[3*nkn+i] = -(A11[0]*der[0] + A11[1]*der[1]);
      cpsiv[3*nkn+i] = -(A11[2]*der[0] + A11[3]*der[1]);
      cpsiuu[3*nkn+i] = -(A21[0]*der[0] + A21[1]*der[1] +
                          A22[0]*der[2] + A22[1]*der[3] + A22[2]*der[4]);
      cpsiuv[3*nkn+i] = -(A21[2]*der[0] + A21[3]*der[1] +
                          A22[3]*der[2] + A22[4]*der[3] + A22[5]*der[4]);
      cpsivv[3*nkn+i] = -(A21[4]*der[0] + A21[5]*der[1] +
                          A22[6]*der[2] + A22[7]*der[3] + A22[8]*der[4]);
    }
        /* curve u == 1 */
  if ( atbstt0i[nk+1] )
    for ( i = fknj; i <= lknj; i++ ) {
      _g1h_TensDer2d ( atbsi[nk+1], atbsti[nk+1], atbstt0i[nk+1],
                       tbsj[i], tbstj[i], tbsttj[i], der );
      A11 = &ctrd[38*(2*nkn+i)];  A21 = &A11[4];  A22 = &A21[6];
      cpsiu[2*nkn+i] = -(A11[0]*der[0] + A11[1]*der[1]);
      cpsiv[2*nkn+i] = -(A11[2]*der[0] + A11[3]*der[1]);
      cpsiuu[2*nkn+i] = -(A21[0]*der[0] + A21[1]*der[1] +
                          A22[0]*der[2] + A22[1]*der[3] + A22[2]*der[4]);
      cpsiuv[2*nkn+i] = -(A21[2]*der[0] + A21[3]*der[1] +
                          A22[3]*der[2] + A22[4]*der[3] + A22[5]*der[4]);
      cpsivv[2*nkn+i] = -(A21[4]*der[0] + A21[5]*der[1] +
                          A22[6]*der[2] + A22[7]*der[3] + A22[8]*der[4]);
    }

  if ( njcurves > 3 ) {
        /* curves inside Omega_i */
          /* curves v == const */
    for ( j = 1; j <= nk; j++ )
      if ( atbsj[j] || atbstj[j] || atbstt0j[j] || atbstt1j[j] )
        for ( i = fkni; i <= lkni; i++ ) {
          _g1h_TensDer2d ( tbsi[i], tbsti[i], tbstti[i],
                           atbsj[j], atbstj[j], atbstt1j[j], der );
          A11 = &ctrd[38*((2+j)*nkn+i)];  A21 = &A11[4];  A22 = &A21[6];
          cpsiu[(3+j)*nkn+i] = A11[0]*der[0] + A11[1]*der[1];
          cpsiv[(3+j)*nkn+i] = A11[2]*der[0] + A11[3]*der[1];
          cpsiuu[(3+j)*nkn+i] = A21[0]*der[0] + A21[1]*der[1] +
                                A22[0]*der[2] + A22[1]*der[3] + A22[2]*der[4];
          cpsiuv[(3+j)*nkn+i] = A21[2]*der[0] + A21[3]*der[1] +
                                A22[3]*der[2] + A22[4]*der[3] + A22[5]*der[4];
          cpsivv[(3+j)*nkn+i] = A21[4]*der[0] + A21[5]*der[1] +
                                A22[6]*der[2] + A22[7]*der[3] + A22[8]*der[4];
          _g1h_TensDer2d ( tbsi[i], tbsti[i], tbstti[i],
                           atbsj[j], atbstj[j], atbstt0j[j], der );
          A11 = &ctrd[38*((2+j)*nkn+i)+19];  A21 = &A11[4];  A22 = &A21[6];
          cpsiuu[(3+j)*nkn+i] -= A21[0]*der[0] + A21[1]*der[1] +
                                 A22[0]*der[2] + A22[1]*der[3] + A22[2]*der[4];
          cpsiuv[(3+j)*nkn+i] -= A21[2]*der[0] + A21[3]*der[1] +
                                 A22[3]*der[2] + A22[4]*der[3] + A22[5]*der[4];
          cpsivv[(3+j)*nkn+i] -= A21[4]*der[0] + A21[5]*der[1] +
                                 A22[6]*der[2] + A22[7]*der[3] + A22[8]*der[4];
        }
          /* curves u == const */
    for ( j = 1; j <= nk; j++ )
      if ( atbsi[j] || atbsti[j] || atbstt0i[j] || atbstt1i[j] )
        for ( i = fknj; i <= lknj; i++ ) {
          _g1h_TensDer2d ( atbsi[j], atbsti[j], atbstt1i[j],
                           tbsj[i], tbstj[i], tbsttj[i], der );
          A11 = &ctrd[38*((2+nk+j)*nkn+i)];  A21 = &A11[4];  A22 = &A21[6];
          cpsiu[(3+nk+j)*nkn+i] = A11[0]*der[0] + A11[1]*der[1];
          cpsiv[(3+nk+j)*nkn+i] = A11[2]*der[0] + A11[3]*der[1];
          cpsiuu[(3+nk+j)*nkn+i] = A21[0]*der[0] + A21[1]*der[1] +
                                   A22[0]*der[2] + A22[1]*der[3] + A22[2]*der[4];
          cpsiuv[(3+nk+j)*nkn+i] = A21[2]*der[0] + A21[3]*der[1] +
                                   A22[3]*der[2] + A22[4]*der[3] + A22[5]*der[4];
          cpsivv[(3+nk+j)*nkn+i] = A21[4]*der[0] + A21[5]*der[1] +
                                   A22[6]*der[2] + A22[7]*der[3] + A22[8]*der[4];
          _g1h_TensDer2d ( atbsi[j], atbsti[j], atbstt0i[j],
                           tbsj[i], tbstj[i], tbsttj[i], der );
          A11 = &ctrd[38*((2+nk+j)*nkn+i)+19];  A21 = &A11[4];  A22 = &A21[6];
          cpsiuu[(3+nk+j)*nkn+i] -= A21[0]*der[0] + A21[1]*der[1] +
                                    A22[0]*der[2] + A22[1]*der[3] + A22[2]*der[4];
          cpsiuv[(3+nk+j)*nkn+i] -= A21[2]*der[0] + A21[3]*der[1] +
                                    A22[3]*der[2] + A22[4]*der[3] + A22[5]*der[4];
          cpsivv[(3+nk+j)*nkn+i] -= A21[4]*der[0] + A21[5]*der[1] +
                                    A22[6]*der[2] + A22[7]*der[3] + A22[8]*der[4];
        }
  }
} /*_g1hq2_TabSplCJumpd*/

static boolean _g1hq2_TabSplDJumpd ( GHoleDomaind *domain,
                  int njcurves, int nk, int m1,
                  int nkn, const double *tkn,
                  const double *hfunc, const double *dhfunc, const double *ddhfunc,
                  const double *ahfunc, const double *adhfunc, const double *addhfunc,
                  int nzc, int fkni, int lkni,
                  int lastomcknot, const double *omcknots, const double *fcomc,
                  int lastpvknot, const double *pvknots,
                  const double *pv, const double *pu,
                  double *trd, double *trdi, double *trdj,
                  double *cpsiu, double *cpsiv,
                  double *cpsiuu, double *cpsiuv, double *cpsivv )
{
  void   *sp;
  int    i, j, ii, jj;
  double *fp, *A11, *A21, *A22, der[6];
  double f00, f01, f02, f02a, f10, f11, f12, f12a;

  sp = pkv_GetScratchMemTop ();
  fp = pkv_GetScratchMemd ( 9*nkn );
  if ( !fp ) {
    domain->error_code = G1H_ERROR_NO_SCRATCH_MEMORY;
    goto failure;
  }

  memset ( fp, 0, 9*nkn*sizeof(double) );
  ii = 9*fkni;
  if ( nzc == 0 )
    mbs_TabBSCurveDer2d ( 1, G1_CROSS00DEG, lastomcknot, omcknots, fcomc,
            lkni-fkni+1, &tkn[fkni], 9, &fp[ii], &fp[ii+1], &fp[ii+2] );
  mbs_TabBSCurveDer2d ( 1, G1_CROSS01DEG, lastpvknot, pvknots, pv,
          lkni-fkni+1, &tkn[fkni], 9, &fp[ii+3], &fp[ii+4], &fp[ii+5] );
  mbs_TabBSCurveDer2d ( 1, G1_CROSS01DEG, lastpvknot, pvknots, pu,
          lkni-fkni+1, &tkn[fkni], 9, &fp[ii+6], &fp[ii+7], &fp[ii+8] );

        /* curves in Omega_k */
          /* 0: v == 0 */
  for ( i = fkni;  i <= lkni;  i++, ii += 9 ) {
    switch ( nzc ) {
  default:
      der[0] = fp[ii+1]*ahfunc[0]  + fp[ii+4]*ahfunc[2];  
      der[1] = fp[ii]*adhfunc[0]   + fp[ii+3]*adhfunc[2]; 
      der[2] = fp[ii+2]*ahfunc[0]  + fp[ii+5]*ahfunc[2];  
      der[3] = fp[ii+1]*adhfunc[0] + fp[ii+4]*adhfunc[2]; 
      der[4] = fp[ii]*addhfunc[0]  + fp[ii+3]*addhfunc[2];
      break;
  case 1:
      der[0] = fp[ii+4]*ahfunc[2];  
      der[1] = fp[ii+3]*adhfunc[2]; 
      der[2] = fp[ii+5]*ahfunc[2];  
      der[3] = fp[ii+4]*adhfunc[2]; 
      der[4] = fp[ii+3]*addhfunc[2];
      break;
    }
    A11 = &trd[38*i];  A21 = &A11[4];  A22 = &A21[6];
    cpsiu[i] = A11[0]*der[0] + A11[1]*der[1];
    cpsiv[i] = A11[2]*der[0] + A11[3]*der[1];
    cpsiuu[i] = A21[0]*der[0] + A21[1]*der[1] +
                A22[0]*der[2] + A22[1]*der[3] + A22[2]*der[4];
    cpsiuv[i] = A21[2]*der[0] + A21[3]*der[1] +
                A22[3]*der[2] + A22[4]*der[3] + A22[5]*der[4];
    cpsivv[i] = A21[4]*der[0] + A21[5]*der[1] +
                A22[6]*der[2] + A22[7]*der[3] + A22[8]*der[4];
    switch ( nzc ) {
  default:
      der[0] = fp[ii]*adhfunc[0]   + fp[ii+6]*adhfunc[2];
      der[1] = fp[ii+1]*ahfunc[0]  + fp[ii+7]*ahfunc[2]; 
      der[2] = fp[ii]*addhfunc[0]  + fp[ii+6]*addhfunc[2];
      der[3] = fp[ii+1]*adhfunc[0] + fp[ii+7]*adhfunc[2]; 
      der[4] = fp[ii+2]*ahfunc[0]  + fp[ii+8]*ahfunc[2];  
      break;
  case 1:
      der[0] = fp[ii+6]*adhfunc[2];
      der[1] = fp[ii+7]*ahfunc[2]; 
      der[2] = fp[ii+6]*addhfunc[2];
      der[3] = fp[ii+7]*adhfunc[2]; 
      der[4] = fp[ii+8]*ahfunc[2];  
      break;
    }
    A11 = &trd[38*i+19];  A21 = &A11[4];  A22 = &A21[6];
    cpsiuu[i] -= A21[0]*der[0] + A21[1]*der[1] +
                 A22[0]*der[2] + A22[1]*der[3] + A22[2]*der[4];
    cpsiuv[i] -= A21[2]*der[0] + A21[3]*der[1] +
                 A22[3]*der[2] + A22[4]*der[3] + A22[5]*der[4];
    cpsivv[i] -= A21[4]*der[0] + A21[5]*der[1] +
                 A22[6]*der[2] + A22[7]*der[3] + A22[8]*der[4];
  }
          /* 1: v == 1 */
  if ( njcurves > 3 )
    j = 4*(nk+1);
  else
    j = 4;
  for ( i = fkni, ii = 9*fkni;  i <= lkni;  i++, ii += 9 ) {
    switch ( nzc ) {
  default:
      der[0] = fp[ii+1]*ahfunc[j]  + fp[ii+4]*ahfunc[j+2];
      der[1] = fp[ii]*adhfunc[j]   + fp[ii+3]*adhfunc[j+2];
      der[2] = fp[ii+2]*ahfunc[j]  + fp[ii+5]*ahfunc[j+2]; 
      der[3] = fp[ii+1]*adhfunc[j] + fp[ii+4]*adhfunc[j+2];
      der[4] = fp[ii]*addhfunc[j]  + fp[ii+3]*addhfunc[j+2];
      break;
  case 1:
      der[0] = fp[ii+4]*ahfunc[j+2];
      der[1] = fp[ii+3]*adhfunc[j+2];
      der[2] = fp[ii+5]*ahfunc[j+2]; 
      der[3] = fp[ii+4]*adhfunc[j+2];
      der[4] = fp[ii+3]*addhfunc[j+2];
      break;
    }
    A11 = &trd[38*(nkn+i)];  A21 = &A11[4];  A22 = &A21[6];
    cpsiu[nkn+i] = A11[0]*der[0] + A11[1]*der[1];
    cpsiv[nkn+i] = A11[2]*der[0] + A11[3]*der[1];
    cpsiuu[nkn+i] = -(A21[0]*der[0] + A21[1]*der[1] +
                      A22[0]*der[2] + A22[1]*der[3] + A22[2]*der[4]);
    cpsiuv[nkn+i] = -(A21[2]*der[0] + A21[3]*der[1] +
                      A22[3]*der[2] + A22[4]*der[3] + A22[5]*der[4]);
    cpsivv[nkn+i] = -(A21[4]*der[0] + A21[5]*der[1] +
                      A22[6]*der[2] + A22[7]*der[3] + A22[8]*der[4]);
  }
          /* 2: u == 1 */
  if ( fcomc[lastomcknot-G1_CROSS00DEG-3] || pv[lastpvknot-G1_CROSS01DEG-3] ) {
    if ( nzc == 0 )
      mbs_deBoorDer2C1d ( G1_CROSS00DEG, lastomcknot, omcknots, fcomc, 1.0,
                          &f00, &f01, &f02 );
    mbs_deBoorDer2C1d ( G1_CROSS01DEG, lastpvknot, pvknots, pv, 1.0,
                        &f10, &f11, &f12 );
    for ( i = ii = 0;  i < nkn;  i++, ii += 4 ) {
      switch ( nzc ) {
    default:
        der[0] = f01*hfunc[ii]   + f11*hfunc[ii+2];  
        der[1] = f00*dhfunc[ii]  + f10*dhfunc[ii+2]; 
        der[2] = f02*hfunc[ii]   + f12*hfunc[ii+2];  
        der[3] = f01*dhfunc[ii]  + f11*dhfunc[ii+2]; 
        der[4] = f00*ddhfunc[ii] + f10*ddhfunc[ii+2];
        break;
    case 1:
        der[0] = f11*hfunc[ii+2];  
        der[1] = f10*dhfunc[ii+2]; 
        der[2] = f12*hfunc[ii+2];  
        der[3] = f11*dhfunc[ii+2]; 
        der[4] = f10*ddhfunc[ii+2];
        break;
      }
      A11 = &trd[38*(2*nkn+i)];  A21 = &A11[4];  A22 = &A21[6];
      cpsiu[2*nkn+i] = A11[0]*der[0] + A11[1]*der[1];
      cpsiv[2*nkn+i] = A11[2]*der[0] + A11[3]*der[1];
      cpsiuu[2*nkn+i] = -(A21[0]*der[0] + A21[1]*der[1] +
                          A22[0]*der[2] + A22[1]*der[3] + A22[2]*der[4]);
      cpsiuv[2*nkn+i] = -(A21[2]*der[0] + A21[3]*der[1] +
                          A22[3]*der[2] + A22[4]*der[3] + A22[5]*der[4]);
      cpsivv[2*nkn+i] = -(A21[4]*der[0] + A21[5]*der[1] +
                          A22[6]*der[2] + A22[7]*der[3] + A22[8]*der[4]);
    }
  }
          /* 3: u == 0 */
  if ( fcomc[2] || pv[2] ) {
    if ( nzc == 0 )
      mbs_deBoorDer2C1d ( G1_CROSS00DEG, lastomcknot, omcknots, fcomc, 0.0,
                          &f00, &f01, &f02 );
    mbs_deBoorDer2C1d ( G1_CROSS01DEG, lastpvknot, pvknots, pv, 0.0,
                        &f10, &f11, &f12 );
    for ( i = ii = 0;  i < nkn;  i++, ii += 4 ) {
      switch ( nzc ) {
    default:
        der[0] = f01*hfunc[ii]   + f11*hfunc[ii+2];  
        der[1] = f00*dhfunc[ii]  + f10*dhfunc[ii+2]; 
        der[2] = f02*hfunc[ii]   + f12*hfunc[ii+2];  
        der[3] = f01*dhfunc[ii]  + f11*dhfunc[ii+2]; 
        der[4] = f00*ddhfunc[ii] + f10*ddhfunc[ii+2];
        break;
    case 1:
        der[0] = f11*hfunc[ii+2];  
        der[1] = f10*dhfunc[ii+2]; 
        der[2] = f12*hfunc[ii+2];  
        der[3] = f11*dhfunc[ii+2]; 
        der[4] = f10*ddhfunc[ii+2];
        break;
      }
      A11 = &trdi[38*i+19];  A21 = &A11[4];  A22 = &A21[6];
      cpsiu[3*nkn+i] = A11[0]*der[0] + A11[1]*der[1];
      cpsiv[3*nkn+i] = A11[2]*der[0] + A11[3]*der[1];
      cpsiuu[3*nkn+i] = -(A21[0]*der[0] + A21[1]*der[1] +
                          A22[0]*der[2] + A22[1]*der[3] + A22[2]*der[4]);
      cpsiuv[3*nkn+i] = -(A21[2]*der[0] + A21[3]*der[1] +
                          A22[3]*der[2] + A22[4]*der[3] + A22[5]*der[4]);
      cpsivv[3*nkn+i] = -(A21[4]*der[0] + A21[5]*der[1] +
                          A22[6]*der[2] + A22[7]*der[3] + A22[8]*der[4]);
    }
  }
  if ( njcurves > 3 ) {
          /* 6+j: v == const */
    for ( j = 1, jj = 4;  j <= nk;  j++, jj += 4 ) {
      for ( i = fkni, ii = 9*i;  i <= lkni;  i++, ii += 9 ) {
        switch ( nzc ) {
      default:
          der[0] = fp[ii+1]*ahfunc[jj]  + fp[ii+4]*ahfunc[jj+2];
          der[1] = fp[ii]*adhfunc[jj]   + fp[ii+3]*adhfunc[jj+2];
          der[2] = fp[ii+2]*ahfunc[jj]  + fp[ii+5]*ahfunc[jj+2]; 
          der[3] = fp[ii+1]*adhfunc[jj] + fp[ii+4]*adhfunc[jj+2];
          der[4] = fp[ii]*addhfunc[jj]  + fp[ii+3]*addhfunc[jj+2];
          break;
      case 1:   
          der[0] = fp[ii+4]*ahfunc[jj+2];
          der[1] = fp[ii+3]*adhfunc[jj+2];
          der[2] = fp[ii+5]*ahfunc[jj+2]; 
          der[3] = fp[ii+4]*adhfunc[jj+2];
          der[4] = fp[ii+3]*addhfunc[jj+2];
          break;
        }
        A11 = &trd[38*((2+j)*nkn+i)];  A21 = &A11[4];  A22 = &A21[6];
        cpsiu[(6+j)*nkn+i] = A11[0]*der[0] + A11[1]*der[1];
        cpsiv[(6+j)*nkn+i] = A11[2]*der[0] + A11[3]*der[1];
        cpsiuu[(6+j)*nkn+i] = (A21[0]-A21[0+19])*der[0] + (A21[1]-A21[1+19])*der[1] +
                    (A22[0]-A22[0+19])*der[2] + (A22[1]-A22[1+19])*der[3] +
                    (A22[2]-A22[2+19])*der[4];
        cpsiuv[(6+j)*nkn+i] = (A21[2]-A21[2+19])*der[0] + (A21[3]-A21[3+19])*der[1] +
                    (A22[3]-A22[3+19])*der[2] + (A22[4]-A22[4+19])*der[3] +
                    (A22[5]-A22[5+19])*der[4];
        cpsivv[(6+j)*nkn+i] = (A21[4]-A21[4+19])*der[0] + (A21[5]-A21[5+19])*der[1] +
                    (A22[6]-A22[6+19])*der[2] + (A22[7]-A22[7+19])*der[3] +
                    (A22[8]-A22[8+19])*der[4];
      }
    }
          /* 6+nk+j: u == const */
    for ( j = 1; j <= nk; j++ ) {
      jj = G1_CROSS00DEG-1+j*m1;
      if ( nzc == 0 ) {
        mbs_deBoorDer2C1d ( G1_CROSS00DEG, lastomcknot, omcknots, fcomc,
                            omcknots[jj], &f00, &f01, &f02 );
        mbs_deBoorDer2C1d ( G1_CROSS00DEG, jj+G1_CROSS00DEG, omcknots, fcomc,
                            omcknots[jj], &f00, &f01, &f02a );
      }
      jj = G1_CROSS01DEG-3+j*(m1+G1_BF01DEG);
      mbs_deBoorDer2C1d ( G1_CROSS01DEG, lastpvknot, pvknots, pv,
                          pvknots[jj], &f10, &f11, &f12 );
      mbs_deBoorDer2C1d ( G1_CROSS01DEG, jj+G1_CROSS01DEG, pvknots, pv,
                          pvknots[jj], &f10, &f11, &f12a );
      for ( i = ii = 0;  i < nkn;  i++, ii += 4 ) {
        switch ( nzc ) {
      default:
          der[0] = f01*hfunc[ii]   + f11*hfunc[ii+2];  
          der[1] = f00*dhfunc[ii]  + f10*dhfunc[ii+2]; 
          der[2] = f02*hfunc[ii]   + f12*hfunc[ii+2];  
          der[5] = f02a*hfunc[ii]  + f12a*hfunc[ii+2];  
          der[3] = f01*dhfunc[ii]  + f11*dhfunc[ii+2]; 
          der[4] = f00*ddhfunc[ii] + f10*ddhfunc[ii+2];
          break;
      case 1:
          der[0] = f11*hfunc[ii+2];  
          der[1] = f10*dhfunc[ii+2]; 
          der[2] = f12*hfunc[ii+2];  
          der[5] = f12a*hfunc[ii+2];
          der[3] = f11*dhfunc[ii+2]; 
          der[4] = f10*ddhfunc[ii+2];
          break;
        }
        A11 = &trd[38*((2+nk+j)*nkn+i)];  A21 = &A11[4];  A22 = &A21[6];
        cpsiu[(6+nk+j)*nkn+i] = A11[0]*der[0] + A11[1]*der[1];
        cpsiv[(6+nk+j)*nkn+i] = A11[2]*der[0] + A11[3]*der[1];
        cpsiuu[(6+nk+j)*nkn+i] = A21[0]*der[0] + A21[1]*der[1] +
                              A22[0]*der[2] + A22[1]*der[3] + A22[2]*der[4];
        cpsiuv[(6+nk+j)*nkn+i] = A21[2]*der[0] + A21[3]*der[1] +
                              A22[3]*der[2] + A22[4]*der[3] + A22[5]*der[4];
        cpsivv[(6+nk+j)*nkn+i] = A21[4]*der[0] + A21[5]*der[1] +
                              A22[6]*der[2] + A22[7]*der[3] + A22[8]*der[4];
        A11 = &trd[38*((2+nk+j)*nkn+i)+19];  A21 = &A11[4];  A22 = &A21[6];
        cpsiuu[(6+nk+j)*nkn+i] -= A21[0]*der[0] + A21[1]*der[1] +
                               A22[0]*der[5] + A22[1]*der[3] + A22[2]*der[4];
        cpsiuv[(6+nk+j)*nkn+i] -= A21[2]*der[0] + A21[3]*der[1] +
                               A22[3]*der[5] + A22[4]*der[3] + A22[5]*der[4];
        cpsivv[(6+nk+j)*nkn+i] -= A21[4]*der[0] + A21[5]*der[1] +
                               A22[6]*der[5] + A22[7]*der[3] + A22[8]*der[4];
      }
    }
  }
        /* curves in Omega_{k-1} */
          /* 4: v == 0 */
  if ( nzc == 0 )
    mbs_deBoorDer2C1d ( G1_CROSS00DEG, lastomcknot, omcknots, fcomc,
                        0.0, &f00, &f01, &f02 );
  mbs_deBoorDer2C1d ( G1_CROSS01DEG, lastpvknot, pvknots, pu,
                      0.0, &f10, &f11, &f12 );
  for ( i = ii = 0;  i < nkn;  i++, ii += 4 ) {
    switch ( nzc ) {
  default:
      der[0] = f00*dhfunc[ii]  + f10*dhfunc[ii+2]; 
      der[1] = f01*hfunc[ii]   + f11*hfunc[ii+2];  
      der[2] = f00*ddhfunc[ii] + f10*ddhfunc[ii+2];
      der[3] = f01*dhfunc[ii]  + f11*dhfunc[ii+2];
      der[4] = f02*hfunc[ii]   + f12*hfunc[ii+2];
      break;
  case 1:
      der[0] = f10*dhfunc[ii+2]; 
      der[1] = f11*hfunc[ii+2];  
      der[2] = f10*ddhfunc[ii+2];
      der[3] = f11*dhfunc[ii+2];
      der[4] = f12*hfunc[ii+2];
      break;
    }
    A11 = &trdj[38*i];  A21 = &A11[4];  A22 = &A21[6];
    cpsiu[4*nkn+i] = A11[0]*der[0] + A11[1]*der[1];
    cpsiv[4*nkn+i] = A11[2]*der[0] + A11[3]*der[1];
    cpsiuu[4*nkn+i] = A21[0]*der[0] + A21[1]*der[1] +
                      A22[0]*der[2] + A22[1]*der[3] + A22[2]*der[4];
    cpsiuv[4*nkn+i] = A21[2]*der[0] + A21[3]*der[1] +
                      A22[3]*der[2] + A22[4]*der[3] + A22[5]*der[4];
    cpsivv[4*nkn+i] = A21[4]*der[0] + A21[5]*der[1] +
                      A22[6]*der[2] + A22[7]*der[3] + A22[8]*der[4];
  }
          /* 5: v == 1 */
  if ( nzc == 0 )
    mbs_deBoorDer2C1d ( G1_CROSS00DEG, lastomcknot, omcknots, fcomc,
                        1.0, &f00, &f01, &f02 );
  mbs_deBoorDer2C1d ( G1_CROSS01DEG, lastpvknot, pvknots, pu,
                      1.0, &f10, &f11, &f12 );
  for ( i = ii = 0;  i < nkn;  i++, ii += 4 ) {
    switch ( nzc ) {
  default:
      der[0] = f00*dhfunc[ii]  + f10*dhfunc[ii+2]; 
      der[1] = f01*hfunc[ii]   + f11*hfunc[ii+2];  
      der[2] = f00*ddhfunc[ii] + f10*ddhfunc[ii+2];
      der[3] = f01*dhfunc[ii]  + f11*dhfunc[ii+2];
      der[4] = f02*hfunc[ii]   + f12*hfunc[ii+2];
      break;
  case 1:
      der[0] = f10*dhfunc[ii+2]; 
      der[1] = f11*hfunc[ii+2];  
      der[2] = f10*ddhfunc[ii+2];
      der[3] = f11*dhfunc[ii+2];
      der[4] = f12*hfunc[ii+2];
      break;
    }
    A11 = &trdj[38*(nkn+i)];  A21 = &A11[4];  A22 = &A21[6];
    cpsiu[5*nkn+i] = A11[0]*der[0] + A11[1]*der[1];
    cpsiv[5*nkn+i] = A11[2]*der[0] + A11[3]*der[1];
    cpsiuu[5*nkn+i] = -(A21[0]*der[0] + A21[1]*der[1] +
                        A22[0]*der[2] + A22[1]*der[3] + A22[2]*der[4]);
    cpsiuv[5*nkn+i] = -(A21[2]*der[0] + A21[3]*der[1] +
                        A22[3]*der[2] + A22[4]*der[3] + A22[5]*der[4]);
    cpsivv[5*nkn+i] = -(A21[4]*der[0] + A21[5]*der[1] +
                        A22[6]*der[2] + A22[7]*der[3] + A22[8]*der[4]);
  }
          /* 6: u == 1 */
  if ( njcurves > 3 )
    jj = 4*(nk+1);
  else
    jj = 4;
  for ( i = fkni, ii = 9*fkni;  i <= lkni;  i++, ii += 9 ) {
    switch ( nzc ) {
  default:
      der[0] = fp[ii]*adhfunc[jj]   + fp[ii+6]*adhfunc[jj+2];
      der[1] = fp[ii+1]*ahfunc[jj]  + fp[ii+7]*ahfunc[jj+2]; 
      der[2] = fp[ii]*addhfunc[jj]  + fp[ii+6]*addhfunc[jj+2];
      der[3] = fp[ii+1]*adhfunc[jj] + fp[ii+7]*adhfunc[jj+2]; 
      der[4] = fp[ii+2]*ahfunc[jj]  + fp[ii+8]*ahfunc[jj+2];  
      break;
  case 1:
      der[0] = fp[ii+6]*adhfunc[jj+2];
      der[1] = fp[ii+7]*ahfunc[jj+2];
      der[2] = fp[ii+6]*addhfunc[jj+2];
      der[3] = fp[ii+7]*adhfunc[jj+2];
      der[4] = fp[ii+8]*ahfunc[jj+2];
      break;
    }
    A11 = &trdj[38*(2*nkn+i)];  A21 = &A11[4];  A22 = &A21[6];
    cpsiu[6*nkn+i] = A11[0]*der[0] + A11[1]*der[1];
    cpsiv[6*nkn+i] = A11[2]*der[0] + A11[3]*der[1];
    cpsiuu[6*nkn+i] = -(A21[0]*der[0] + A21[1]*der[1] +
                        A22[0]*der[2] + A22[1]*der[3] + A22[2]*der[4]);
    cpsiuv[6*nkn+i] = -(A21[2]*der[0] + A21[3]*der[1] +
                        A22[3]*der[2] + A22[4]*der[3] + A22[5]*der[4]);
    cpsivv[6*nkn+i] = -(A21[4]*der[0] + A21[5]*der[1] +
                        A22[6]*der[2] + A22[7]*der[3] + A22[8]*der[4]);
  }
          /* 6+2*nk+j: v == const */
  if ( njcurves > 3 ) {
    for ( j = 1; j <= nk;  j++ ) {
      jj = G1_CROSS00DEG-1+j*m1;
      if ( nzc == 0 ) {
        mbs_deBoorDer2C1d ( G1_CROSS00DEG, lastomcknot, omcknots, fcomc,
                            omcknots[jj], &f00, &f01, &f02 );
        mbs_deBoorDer2C1d ( G1_CROSS00DEG, jj+G1_CROSS00DEG, omcknots, fcomc,
                            omcknots[jj], &f00, &f01, &f02a );
      }
      jj = G1_CROSS01DEG-3+j*(m1+G1_BF01DEG);
      mbs_deBoorDer2C1d ( G1_CROSS01DEG, lastpvknot, pvknots, pu,
                          pvknots[jj], &f10, &f11, &f12 );
      mbs_deBoorDer2C1d ( G1_CROSS01DEG, jj+G1_CROSS01DEG, pvknots, pu,
                          pvknots[jj], &f10, &f11, &f12a );
      for ( i = ii = 0;  i < nkn;  i++, ii += 4 ) {
        switch ( nzc ) {
      default:
          der[0] = f00*dhfunc[ii]  + f10*dhfunc[ii+2]; 
          der[1] = f01*hfunc[ii]   + f11*hfunc[ii+2];  
          der[2] = f00*ddhfunc[ii] + f10*ddhfunc[ii+2];
          der[3] = f01*dhfunc[ii]  + f11*dhfunc[ii+2];
          der[4] = f02*hfunc[ii]   + f12*hfunc[ii+2];
          der[5] = f02a*hfunc[ii]  + f12a*hfunc[ii+2];
          break;
      case 1:
          der[0] = f10*dhfunc[ii+2]; 
          der[1] = f11*hfunc[ii+2];  
          der[2] = f10*ddhfunc[ii+2];
          der[3] = f11*dhfunc[ii+2];
          der[4] = f12*hfunc[ii+2];
          der[5] = f12a*hfunc[ii+2];
          break;
        }
        A11 = &trdj[38*((2+j)*nkn+i)];  A21 = &A11[4];  A22 = &A21[6];
        cpsiu[(6+2*nk+j)*nkn+i] = A11[0]*der[0] + A11[1]*der[1];
        cpsiv[(6+2*nk+j)*nkn+i] = A11[2]*der[0] + A11[3]*der[1];
        cpsiuu[(6+2*nk+j)*nkn+i] = A21[0]*der[0] + A21[1]*der[1] +
                                   A22[0]*der[2] + A22[1]*der[3] + A22[2]*der[4];
        cpsiuv[(6+2*nk+j)*nkn+i] = A21[2]*der[0] + A21[3]*der[1] +
                                   A22[3]*der[2] + A22[4]*der[3] + A22[5]*der[4];
        cpsivv[(6+2*nk+j)*nkn+i] = A21[4]*der[0] + A21[5]*der[1] +
                                   A22[6]*der[2] + A22[7]*der[3] + A22[8]*der[4];
        A11 = &trdj[38*((2+j)*nkn+i)+19];  A21 = &A11[4];  A22 = &A21[6];
        cpsiuu[(6+2*nk+j)*nkn+i] -= A21[0]*der[0] + A21[1]*der[1] +
                                    A22[0]*der[2] + A22[1]*der[3] + A22[2]*der[5];
        cpsiuv[(6+2*nk+j)*nkn+i] -= A21[2]*der[0] + A21[3]*der[1] +
                                    A22[3]*der[2] + A22[4]*der[3] + A22[5]*der[5];
        cpsivv[(6+2*nk+j)*nkn+i] -= A21[4]*der[0] + A21[5]*der[1] +
                                    A22[6]*der[2] + A22[7]*der[3] + A22[8]*der[5];
      }
    }
          /* 6+3*nk+j: u == const */
    for ( j = 1, jj = 4;  j <= nk;  j++, jj += 4 ) {
      for ( i = fkni, ii = 9*i;  i <= lkni;  i++, ii += 9 ) {
        switch ( nzc ) {
      default:
          der[0] = fp[ii]*adhfunc[jj]   + fp[ii+6]*adhfunc[jj+2];
          der[1] = fp[ii+1]*ahfunc[jj]  + fp[ii+7]*ahfunc[jj+2]; 
          der[2] = fp[ii]*addhfunc[jj]  + fp[ii+6]*addhfunc[jj+2];
          der[3] = fp[ii+1]*adhfunc[jj] + fp[ii+7]*adhfunc[jj+2]; 
          der[4] = fp[ii+2]*ahfunc[jj]  + fp[ii+8]*ahfunc[jj+2];  
          break;
      case 1:
          der[0] = fp[ii+6]*adhfunc[jj+2];
          der[1] = fp[ii+7]*ahfunc[jj+2];
          der[2] = fp[ii+6]*addhfunc[jj+2];
          der[3] = fp[ii+7]*adhfunc[jj+2];
          der[4] = fp[ii+8]*ahfunc[jj+2];
          break;
        }
        A11 = &trdj[38*((2+nk+j)*nkn+i)];  A21 = &A11[4];  A22 = &A21[6];
        cpsiu[(6+3*nk+j)*nkn+i] = A11[0]*der[0] + A11[1]*der[1];
        cpsiv[(6+3*nk+j)*nkn+i] = A11[2]*der[0] + A11[3]*der[1];
        cpsiuu[(6+3*nk+j)*nkn+i] = (A21[0]-A21[0+19])*der[0] + (A21[1]-A21[1+19])*der[1] +
                 (A22[0]-A22[0+19])*der[2] + (A22[1]-A22[1+19])*der[3] +
                 (A22[2]-A22[2+19])*der[4];
        cpsiuv[(6+3*nk+j)*nkn+i] = (A21[2]-A21[2+19])*der[0] + (A21[3]-A21[3+19])*der[1] +
                 (A22[3]-A22[3+19])*der[2] + (A22[4]-A22[4+19])*der[3] +
                 (A22[5]-A22[5+19])*der[4];
        cpsivv[(6+3*nk+j)*nkn+i] = (A21[4]-A21[4+19])*der[0] + (A21[5]-A21[5+19])*der[1] +
                 (A22[6]-A22[6+19])*der[2] + (A22[7]-A22[7+19])*der[3] +
                 (A22[8]-A22[8+19])*der[4];
      }
    }
  }

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*_g1hq2_TabSplDJumpd*/

static boolean _g1hq2_TabSplNLBasisFJumpsd ( GHoleDomaind *domain,
                                             G1HNLSPrivated *nlspr )
{
  void *sp;
  G1HNLPrivated      *nlpr;
  GHolePrivateRecd   *privateG;
  G1HolePrivateRecd  *privateG1;
  G1HoleSPrivateRecd *sprivate;
  int      hole_k, nk, m1, m2;
  int      nfunc_a, nfunc_b, nfunc_c;
  int      nkn, nakn, dipitch, disize, lastfpknot, lastcknot;
  double   *fpknots, *cknots, *tkn, *hfunc, *dhfunc, *ddhfunc;
  double   *aknots, *ahfunc, *adhfunc, *addhfunc;
  point2d  *di, *sicp;
  double   *tabeu, *tabev, *tabeuu, *tabeuv, *tabevv,
           *tabfu, *tabfv, *tabfuu, *tabfuv, *tabfvv;
  double   *bc00, *bc01, *bc10, *bc11, *bd00, *bd01, *bd10, *bd11,
           *ec00, *ec01, *ed00, *ed01,
           *fc00, *fc01, *fc10, *fc11, *fd00, *fd01, *fd10, *fd11;
  double   *ctrd, *A11, *A21, *A22;
  double   uv, bvz;
  int      njcurves, jtabsize, rr, m, bstsize;
  int      k, ik, kk, ki, kN, kNQ2, i, j, l, fn, fN, jcfs, jdfs;
  vector2d cdi, cdiu, cdiv, cdiuu, cdiuv, cdivv;
  unsigned char *bfcpn;
  boolean  jumpC, jumpD;
  int      *fkn, *lkn;
  double   *tbs, *tbst, *tbstt, *atbs, *atbst, *atbstt0, *atbstt1;
  int      i00, i01, j00, j01, nzc, lastomcknot, lastpvknot;
  double   *fcomc, *pv, *pu, *omcknots, *pvknots;

  sp = pkv_GetScratchMemTop ();
  hole_k = domain->hole_k;
  nlpr = &nlspr->nlpr;
  privateG  = domain->privateG;
  privateG1 = domain->privateG1;
  sprivate  = domain->SprivateG1;
  jtabsize = nlspr->jtabsize/5;
  nk = sprivate->nk;
  m1 = sprivate->m1;
  m2 = sprivate->m2;
  nkn = nlspr->nkn; 
  tkn = nlspr->tkn; 
  nfunc_a = privateG1->nfunc_a;
  nfunc_b = privateG1->nfunc_b;
  nfunc_c = sprivate->nsfunc_c;
  hfunc = pkv_GetScratchMemd ( 16*nkn );
  lastcknot = sprivate->lastcknot;
  cknots = sprivate->cknots;
  lastfpknot = sprivate->lastfpknot;
  fpknots = sprivate->fpknots;
  dipitch = lastfpknot-G1H_FINALDEG; 
  disize = dipitch*dipitch;
  dipitch *= 2;
  njcurves = nlspr->njcurves;
  di = pkv_GetScratchMem ( disize*sizeof(point2d) );
  sicp = pkv_GetScratchMem ( 16*sizeof(point2d) );
  jumpC = nlspr->jumpC;
  jumpD = nlspr->jumpD;
  if ( jumpC || jumpD )
    nakn = 2+nk;
  else
    nakn = 2;
  tabeu = pkv_GetScratchMemd ( 10*nakn*nkn );
  aknots = pkv_GetScratchMemd ( 13*nakn );
  bc00 = pkv_GetScratchMemd ( 2*(G1_CROSSDEGSUM+4)*hole_k );
  ctrd = pkv_GetScratchMemd ( njcurves*nkn*38*hole_k );
  if ( !hfunc || !di || !sicp || !tabeu || !aknots || !bc00 || !ctrd ) {
    domain->error_code = G1H_ERROR_NO_SCRATCH_MEMORY;
    goto failure;
  }

  tabev = &tabeu[nakn*nkn];    tabeuu = &tabev[nakn*nkn];
  tabeuv = &tabeuu[nakn*nkn];  tabevv = &tabeuv[nakn*nkn];
  tabfu = &tabevv[nakn*nkn];   tabfv = &tabfu[nakn*nkn];
  tabfuu = &tabfv[nakn*nkn];   tabfuv = &tabfuu[nakn*nkn];
  tabfvv = &tabfuv[nakn*nkn];

  dhfunc = &hfunc[4*nkn];
  ddhfunc = &dhfunc[4*nkn];
  ahfunc = &aknots[nakn];
  adhfunc = &ahfunc[4*nakn];
  addhfunc = &adhfunc[4*nakn];

  aknots[0] = 0.0;
  if ( nakn == 2 )
    aknots[1] = 1.0;
  else {
        /* to avoid rounding errors, the knots are copied */
    for ( i = j = 1; i <= lastfpknot; i++ )
      if ( fpknots[i] > fpknots[i-1] )
        aknots[j++] = fpknots[i];
  }

  memset ( ctrd, 0, njcurves*nkn*38*hole_k*sizeof(double) );
  for ( k = kN = kNQ2 = 0, kk = 1;
        k < hole_k;
        k++, kk = (k+1) % hole_k, kN += disize, kNQ2 += njcurves*nkn ) {
            /* get the spline domain patch */
    pkv_Selectd ( disize, 2, 3, 2, &nlpr->nldi[kN], di );
    for ( i = 0; i < nkn; i++ ) {
      mbs_deBoorDer2Pd ( G1H_FINALDEG, lastfpknot, fpknots,
                         G1H_FINALDEG, lastfpknot, fpknots,
                         2, dipitch, &di->x, tkn[i], 0.0,
                         &cdi.x, &cdiu.x, &cdiv.x, &cdiuu.x, &cdiuv.x, &cdivv.x );
      nlpr->ctang[kNQ2+i] = cdiu;
      if ( !_g1hq2_SetupCTrdd ( &cdiu, &cdiv, &cdiuu, &cdiuv, &cdivv,
                                &ctrd[38*(kNQ2+i)] ) )
        goto failure;
      mbs_deBoorDer2Pd ( G1H_FINALDEG, lastfpknot, fpknots,
                         G1H_FINALDEG, lastfpknot, fpknots,
                         2, dipitch, &di->x, tkn[i], 1.0,
                         &cdi.x, &cdiu.x, &cdiv.x, &cdiuu.x, &cdiuv.x, &cdivv.x );
      nlpr->ctang[kNQ2+nkn+i] = cdiu;
      if ( !_g1hq2_SetupCTrdd ( &cdiu, &cdiv, &cdiuu, &cdiuv, &cdivv,
                                &ctrd[38*(kNQ2+nkn+i)] ) )
        goto failure;
      mbs_deBoorDer2Pd ( G1H_FINALDEG, lastfpknot, fpknots,
                         G1H_FINALDEG, lastfpknot, fpknots,
                         2, dipitch, &di->x, 1.0, tkn[i],
                         &cdi.x, &cdiu.x, &cdiv.x, &cdiuu.x, &cdiuv.x, &cdivv.x );
      nlpr->ctang[kNQ2+2*nkn+i] = cdiv;
      if ( !_g1hq2_SetupCTrdd ( &cdiu, &cdiv, &cdiuu, &cdiuv, &cdivv,
                                &ctrd[38*(kNQ2+2*nkn+i)] ) )
        goto failure;
      mbs_deBoorDer2Pd ( G1H_FINALDEG, lastfpknot, fpknots,
                         G1H_FINALDEG, lastfpknot, fpknots,
                         2, dipitch, &di->x, 0.0, tkn[i],
                         &cdi.x, &cdiu.x, &cdiv.x, &cdiuu.x, &cdiuv.x, &cdivv.x );
      if ( !_g1hq2_SetupCTrdd ( &cdiu, &cdiv, &cdiuu, &cdiuv, &cdivv,
                                &ctrd[38*(kk*njcurves*nkn+i)+19] ) )
        goto failure;
    }
    if ( njcurves > 3 ) {
        /* coefficients for the curves inside Omega_i */
          /* v == const */
      for ( j = 1; j <= nk; j++ ) {
        uv = aknots[j];
        for ( i = 0; i < nkn; i++ ) {
          mbs_deBoorDer2Pd ( G1H_FINALDEG, lastfpknot, fpknots,
                    G1H_FINALDEG, 2*G1H_FINALDEG+1, &fpknots[j*(G1H_FINALDEG-1)],
                    2, dipitch, &di[j*(G1H_FINALDEG-1)].x, tkn[i], uv,
                    &cdi.x, &cdiu.x, &cdiv.x, &cdiuu.x, &cdiuv.x, &cdivv.x );
          nlpr->ctang[kNQ2+(2+j)*nkn+i] = cdiu;
          if ( !_g1hq2_SetupCTrdd ( &cdiu, &cdiv, &cdiuu, &cdiuv, &cdivv,
                                    &ctrd[38*(kNQ2+(2+j)*nkn+i)] ) )
            goto failure;
          mbs_deBoorDer2Pd ( G1H_FINALDEG, lastfpknot, fpknots,
                    G1H_FINALDEG, 2*G1H_FINALDEG+1, &fpknots[(j-1)*(G1H_FINALDEG-1)],
                    2, dipitch, &di[(j-1)*(G1H_FINALDEG-1)].x, tkn[i], uv,
                    &cdi.x, &cdiu.x, &cdiv.x, &cdiuu.x, &cdiuv.x, &cdivv.x );
          if ( !_g1hq2_SetupCTrdd ( &cdiu, &cdiv, &cdiuu, &cdiuv, &cdivv,
                                    &ctrd[38*(kNQ2+(2+j)*nkn+i)+19] ) )
            goto failure;
        }
      }
          /* u == const */
      for ( j = 1; j <= nk; j++ ) {
        uv = aknots[j];
        for ( i = 0; i < nkn; i++ ) {
          mbs_deBoorDer2Pd ( G1H_FINALDEG, 2*G1H_FINALDEG+1,
                    &fpknots[j*(G1H_FINALDEG-1)],
                    G1H_FINALDEG, lastfpknot, fpknots, 2, dipitch,
                    &di[j*(G1H_FINALDEG-1)*(lastfpknot-G1H_FINALDEG)].x,
                    uv, tkn[i],
                    &cdi.x, &cdiu.x, &cdiv.x, &cdiuu.x, &cdiuv.x, &cdivv.x );
          nlpr->ctang[kNQ2+(2+nk+j)*nkn+i] = cdiv;
          if ( !_g1hq2_SetupCTrdd ( &cdiu, &cdiv, &cdiuu, &cdiuv, &cdivv,
                                    &ctrd[38*(kNQ2+(2+nk+j)*nkn+i)] ) )
            goto failure;
          mbs_deBoorDer2Pd ( G1H_FINALDEG, 2*G1H_FINALDEG+1,
                    &fpknots[(j-1)*(G1H_FINALDEG-1)],
                    G1H_FINALDEG, lastfpknot, fpknots, 2, dipitch,
                    &di[(j-1)*(G1H_FINALDEG-1)*(lastfpknot-G1H_FINALDEG)].x,
                    uv, tkn[i],
                    &cdi.x, &cdiu.x, &cdiv.x, &cdiuu.x, &cdiuv.x, &cdivv.x );
          nlpr->ctang[kNQ2+(2+nk+j)*nkn+i] = cdiv;
          if ( !_g1hq2_SetupCTrdd ( &cdiu, &cdiv, &cdiuu, &cdiuv, &cdivv,
                                    &ctrd[38*(kNQ2+(2+nk+j)*nkn+i)+19] ) )
            goto failure;
        }
      }
    }
        /* coefficients for the domain surrounding patches */
    if ( !_g1hq2_FindDomSurrndPatchd ( domain, nlpr, kk, 1, sicp ) )
      goto failure;
    for ( i = 0; i < nkn; i++ ) {
      if ( !mbs_BCHornerDer2Pd ( 3, 3, 2, (double*)sicp,
                                 0.0, tkn[i], &cdi.x, &cdiu.x, &cdiv.x,
                                 &cdiuu.x, &cdiuv.x, &cdivv.x ) )
        goto failure;
      if ( !_g1hq2_SetupCTrdd ( &cdiu, &cdiv, &cdiuu, &cdiuv, &cdivv,
                                &ctrd[38*(kNQ2+nkn+i)+19] ) )
        goto failure;
    }
    if ( !_g1hq2_FindDomSurrndPatchd ( domain, nlpr, k, 2, sicp ) )
      goto failure;
    for ( i = 0; i < nkn; i++ ) {
      if ( !mbs_BCHornerDer2Pd ( 3, 3, 2, (double*)sicp,
                                 0.0, tkn[i], &cdi.x, &cdiu.x, &cdiv.x,
                                 &cdiuu.x, &cdiuv.x, &cdivv.x ) )
        goto failure;
      if ( !_g1hq2_SetupCTrdd ( &cdiu, &cdiv, &cdiuu, &cdiuv, &cdivv,
                                &ctrd[38*(kNQ2+2*nkn+i)+19] ) )
        goto failure;
    }
  }

        /* now the computation of the second derivative jumps of the */
        /* basis functions */
  memset ( nlpr->cpsiu, 0, jtabsize*sizeof(double) );
  memset ( nlpr->cpsiv, 0, jtabsize*sizeof(double) );
  memset ( nlpr->cpsiuu, 0, jtabsize*sizeof(double) );
  memset ( nlpr->cpsiuv, 0, jtabsize*sizeof(double) );
  memset ( nlpr->cpsivv, 0, jtabsize*sizeof(double) );
  if ( !mbs_TabCubicHFuncDer2d ( 0.0, 1.0, nakn, aknots, ahfunc, adhfunc, addhfunc ) )
    goto failure;
          /* block A functions */
  for ( fn = 0; fn < nfunc_a; fn++ ) {
    _g1h_GetBFAPatchCurvesd ( domain, fn, hole_k-1, &ec00, &ec01, &ed00, &ed01 );
    for ( k = 0;  k < hole_k;  k++ ) {
      _g1h_GetBFAPatchCurvesd ( domain, fn, k, &fc00, &fc01, &fd00, &fd01 );
      if ( !mbs_TabBezC1Coons0Der2d ( 1, nkn, tkn, hfunc, dhfunc, ddhfunc,
                   nakn, aknots, ahfunc, adhfunc, addhfunc,
                   G1_CROSS00DEG, fc00, G1_CROSS01DEG, fc01,
                   G1_CROSS00DEG, fd00, G1_CROSS01DEG, fd01,
                   NULL, tabfu, tabfv, tabfuu, tabfuv, tabfvv ) )
        goto failure;
      if ( !mbs_TabBezC1Coons0Der2d ( 1, 1, aknots, ahfunc, adhfunc, addhfunc,
                   nkn, tkn, hfunc, dhfunc, ddhfunc,
                   G1_CROSS00DEG, ec00, G1_CROSS01DEG, ec01,
                   G1_CROSS00DEG, ed00, G1_CROSS01DEG, ed01,
                   NULL, tabeu, tabev, tabeuu, tabeuv, tabevv ) )
        goto failure;
      for ( i = j = 0, fN = njcurves*(fn*hole_k+k)*nkn;
            i < nkn;
            i++, j += nakn, fN++ ) {
        A11 = &ctrd[38*(njcurves*k*nkn+i)];  A21 = &A11[4];  A22 = &A21[6];
        nlpr->cpsiu[fN] = A11[0]*tabfu[j] + A11[1]*tabfv[j];
        nlpr->cpsiv[fN] = A11[2]*tabfu[j] + A11[3]*tabfv[j];
        nlpr->cpsiuu[fN] = A21[0]*tabfu[j] + A21[1]*tabfv[j] +
                           A22[0]*tabfuu[j] + A22[1]*tabfuv[j] + A22[2]*tabfvv[j];
        nlpr->cpsiuv[fN] = A21[2]*tabfu[j] + A21[3]*tabfv[j] +
                           A22[3]*tabfuu[j] + A22[4]*tabfuv[j] + A22[5]*tabfvv[j];
        nlpr->cpsivv[fN] = A21[4]*tabfu[j] + A21[5]*tabfv[j] +
                           A22[6]*tabfuu[j] + A22[7]*tabfuv[j] + A22[8]*tabfvv[j];
        A11 = &ctrd[38*(njcurves*k*nkn+i)+19];  A21 = &A11[4];  A22 = &A21[6];
        nlpr->cpsiuu[fN] -= A21[0]*tabeu[i] + A21[1]*tabev[i] +
                            A22[0]*tabeuu[i] + A22[1]*tabeuv[i] + A22[2]*tabevv[i];
        nlpr->cpsiuv[fN] -= A21[2]*tabeu[i] + A21[3]*tabev[i] +
                            A22[3]*tabeuu[i] + A22[4]*tabeuv[i] + A22[5]*tabevv[i];
        nlpr->cpsivv[fN] -= A21[4]*tabeu[i] + A21[5]*tabev[i] +
                            A22[6]*tabeuu[i] + A22[7]*tabeuv[i] + A22[8]*tabevv[i];
      }
      for ( i = 0, j = nakn-1, fN = (njcurves*(fn*hole_k+k)+1)*nkn;
            i < nkn;
            i++, j += nakn, fN++ ) {
        A11 = &ctrd[38*((njcurves*k+1)*nkn+i)];  A21 = &A11[4];  A22 = &A21[6];
        nlpr->cpsiu[fN] = A11[0]*tabfu[j] + A11[1]*tabfv[j];
        nlpr->cpsiv[fN] = A11[2]*tabfu[j] + A11[3]*tabfv[j];
        nlpr->cpsiuu[fN] = -(A21[0]*tabfu[j] + A21[1]*tabfv[j] +
                           A22[0]*tabfuu[j] + A22[1]*tabfuv[j] + A22[2]*tabfvv[j]);
        nlpr->cpsiuv[fN] = -(A21[2]*tabfu[j] + A21[3]*tabfv[j] +
                           A22[3]*tabfuu[j] + A22[4]*tabfuv[j] + A22[5]*tabfvv[j]);
        nlpr->cpsivv[fN] = -(A21[4]*tabfu[j] + A21[5]*tabfv[j] +
                           A22[6]*tabfuu[j] + A22[7]*tabfuv[j] + A22[8]*tabfvv[j]);
      }
      if ( njcurves > 3 ) {
          /* compute the derivatives and jumps inside Omega_i */
            /* v == const */
        for ( j = 1; j <= nk; j++ )
          for ( i = 0, fN = (njcurves*(fn*hole_k+k)+2+j)*nkn, l = j;
                i < nkn;
                i++, fN++, l += nakn ) {
            A11 = &ctrd[38*((njcurves*k+2+j)*nkn+i)];
            A21 = &A11[4];  A22 = &A21[6];
            nlpr->cpsiu[fN] = A11[0]*tabfu[l] + A11[1]*tabfv[l];
            nlpr->cpsiv[fN] = A11[2]*tabfu[l] + A11[3]*tabfv[l];
            nlpr->cpsiuu[fN] = (A21[0]-A21[0+19])*tabfu[l] + (A21[1]-A21[1+19])*tabfv[l] +
                       (A22[0]-A22[0+19])*tabfuu[l] + (A22[1]-A22[1+19])*tabfuv[l] +
                       (A22[2]-A22[2+19])*tabfvv[l];
            nlpr->cpsiuv[fN] = (A21[2]-A21[2+19])*tabfu[l] + (A21[3]-A21[3+19])*tabfv[l] +
                       (A22[3]-A22[3+19])*tabfuu[l] + (A22[4]-A22[4+19])*tabfuv[l] +
                       (A22[5]-A22[5+19])*tabfvv[l];
            nlpr->cpsivv[fN] = (A21[4]-A21[4+19])*tabfu[l] + (A21[5]-A21[5+19])*tabfv[l] +
                       (A22[6]-A22[6+19])*tabfuu[l] + (A22[7]-A22[7+19])*tabfuv[l] +
                       (A22[8]-A22[8+19])*tabfvv[l];
          }
      }
      if ( !mbs_TabBezC1Coons0Der2d ( 1, nakn-1, &aknots[1],
                  &ahfunc[4], &adhfunc[4], &addhfunc[4],
                  nkn, tkn, hfunc, dhfunc, ddhfunc,
                  G1_CROSS00DEG, fc00, G1_CROSS01DEG, fc01,
                  G1_CROSS00DEG, fd00, G1_CROSS01DEG, fd01,
                  NULL, tabfu, tabfv, tabfuu, tabfuv, tabfvv ) )
        goto failure;
      for ( i = 0, fN = (njcurves*(fn*hole_k+k)+2)*nkn, l = (nakn-2)*nkn;
            i < nkn;
            i++, fN++, l++ ) {
        A11 = &ctrd[38*((njcurves*k+2)*nkn+i)];  A21 = &A11[4];  A22 = &A21[6];
        nlpr->cpsiu[fN] = A11[0]*tabfu[l] + A11[1]*tabfv[l];
        nlpr->cpsiv[fN] = A11[2]*tabfu[l] + A11[3]*tabfv[l];
        nlpr->cpsiuu[fN] = -(A21[0]*tabfu[l] + A21[1]*tabfv[l] +
                           A22[0]*tabfuu[l] + A22[1]*tabfuv[l] + A22[2]*tabfvv[l]);
        nlpr->cpsiuv[fN] = -(A21[2]*tabfu[l] + A21[3]*tabfv[l] +
                           A22[3]*tabfuu[l] + A22[4]*tabfuv[l] + A22[5]*tabfvv[l]);
        nlpr->cpsivv[fN] = -(A21[4]*tabfu[l] + A21[5]*tabfv[l] +
                           A22[6]*tabfuu[l] + A22[7]*tabfuv[l] + A22[8]*tabfvv[l]);
      }
      if ( njcurves > 3 ) {
          /* compute the derivatives and jumps inside Omega_i */
            /* u == const */
        for ( j = 1; j <= nk; j++ )
          for ( i = 0, fN = (njcurves*(fn*hole_k+k)+2+nk+j)*nkn, l = (j-1)*nkn;
                i < nkn;
                i++, fN++, l++ ) {
            A11 = &ctrd[38*((njcurves*k+2+nk+j)*nkn+i)];
            A21 = &A11[4];  A22 = &A21[6];
            nlpr->cpsiu[fN] = A11[0]*tabfu[l] + A11[1]*tabfv[l];
            nlpr->cpsiv[fN] = A11[2]*tabfu[l] + A11[3]*tabfv[l];
            nlpr->cpsiuu[fN] = (A21[0]-A21[0+19])*tabfu[l] + (A21[1]-A21[1+19])*tabfv[l] +
                       (A22[0]-A22[0+19])*tabfuu[l] + (A22[1]-A22[1+19])*tabfuv[l] +
                       (A22[2]-A22[2+19])*tabfvv[l];
            nlpr->cpsiuv[fN] = (A21[2]-A21[2+19])*tabfu[l] + (A21[3]-A21[3+19])*tabfv[l] +
                       (A22[3]-A22[3+19])*tabfuu[l] + (A22[4]-A22[4+19])*tabfuv[l] +
                       (A22[5]-A22[5+19])*tabfvv[l];
            nlpr->cpsivv[fN] = (A21[4]-A21[4+19])*tabfu[l] + (A21[5]-A21[5+19])*tabfv[l] +
                       (A22[6]-A22[6+19])*tabfuu[l] + (A22[7]-A22[7+19])*tabfuv[l] +
                       (A22[8]-A22[8+19])*tabfvv[l];
          }
      }
      ec00 = fc00;  ec01 = fc01;  ed00 = fd00;  ed01 = fd01;
    }
  }
          /* the fixed linear combination of the block B functions */
  bc01 = &bc00[(G1_CROSS00DEG+1)*hole_k];  bc10 = &bc01[(G1_CROSS01DEG+1)*hole_k];
  bc11 = &bc10[(G1_CROSS10DEG+1)*hole_k];  bd00 = &bc11[(G1_CROSS11DEG+1)*hole_k];
  bd01 = &bd00[(G1_CROSS00DEG+1)*hole_k];  bd10 = &bd01[(G1_CROSS01DEG+1)*hole_k];
  bd11 = &bd10[(G1_CROSS10DEG+1)*hole_k];
  memset ( bc00, 0, 2*(G1_CROSSDEGSUM+4)*hole_k*sizeof(double) );
  bfcpn = privateG->bfcpn;
  for ( fn = 0; fn < nfunc_b; fn++ ) {
    bvz = nlpr->rhole_cp[bfcpn[fn]].z;
    for ( k = 0; k < hole_k; k++ ) {
      _g1h_GetBFBPatchCurvesd ( domain, fn, k, &fc00, &fc01, &fc10, &fc11,
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
  for ( k = 0, ik = hole_k-1;
        k < hole_k;
        ik = k++ ) {
    if ( !mbs_TabBezC1CoonsDer2d ( 1, nkn, tkn, hfunc, dhfunc, ddhfunc,
                nakn, aknots, ahfunc, adhfunc, addhfunc,
                G1_CROSS00DEG, &bc00[(G1_CROSS00DEG+1)*k],
                G1_CROSS01DEG, &bc01[(G1_CROSS01DEG+1)*k],
                G1_CROSS10DEG, &bc10[(G1_CROSS10DEG+1)*k],
                G1_CROSS11DEG, &bc11[(G1_CROSS11DEG+1)*k],
                G1_CROSS00DEG, &bd00[(G1_CROSS00DEG+1)*k],
                G1_CROSS01DEG, &bd01[(G1_CROSS01DEG+1)*k],
                G1_CROSS10DEG, &bd10[(G1_CROSS10DEG+1)*k],
                G1_CROSS11DEG, &bd11[(G1_CROSS11DEG+1)*k],
                NULL, tabfu, tabfv, tabfuu, tabfuv, tabfvv ) )
      goto failure;
    if ( !mbs_TabBezC1CoonsDer2d ( 1, 1, aknots, ahfunc, adhfunc, addhfunc,
                nkn, tkn, hfunc, dhfunc, ddhfunc,
                G1_CROSS00DEG, &bc00[(G1_CROSS00DEG+1)*ik],
                G1_CROSS01DEG, &bc01[(G1_CROSS01DEG+1)*ik],
                G1_CROSS10DEG, &bc10[(G1_CROSS10DEG+1)*ik],
                G1_CROSS11DEG, &bc11[(G1_CROSS11DEG+1)*ik],
                G1_CROSS00DEG, &bd00[(G1_CROSS00DEG+1)*ik],
                G1_CROSS01DEG, &bd01[(G1_CROSS01DEG+1)*ik],
                G1_CROSS10DEG, &bd10[(G1_CROSS10DEG+1)*ik],
                G1_CROSS11DEG, &bd11[(G1_CROSS11DEG+1)*ik],
                NULL, tabeu, tabev, tabeuu, tabeuv, tabevv ) )
      goto failure;
    for ( i = j = 0, fN = njcurves*(nfunc_a*hole_k+k)*nkn;
          i < nkn;
          i++, j += nakn, fN++ ) {
      A11 = &ctrd[38*(njcurves*k*nkn+i)];  A21 = &A11[4];  A22 = &A21[6];
      nlpr->cpsiu[fN] = A11[0]*tabfu[j] + A11[1]*tabfv[j];
      nlpr->cpsiv[fN] = A11[2]*tabfu[j] + A11[3]*tabfv[j];
      nlpr->cpsiuu[fN] = A21[0]*tabfu[j] + A21[1]*tabfv[j] +
                         A22[0]*tabfuu[j] + A22[1]*tabfuv[j] + A22[2]*tabfvv[j];
      nlpr->cpsiuv[fN] = A21[2]*tabfu[j] + A21[3]*tabfv[j] +
                         A22[3]*tabfuu[j] + A22[4]*tabfuv[j] + A22[5]*tabfvv[j];
      nlpr->cpsivv[fN] = A21[4]*tabfu[j] + A21[5]*tabfv[j] +
                         A22[6]*tabfuu[j] + A22[7]*tabfuv[j] + A22[8]*tabfvv[j];
      A11 = &ctrd[38*(njcurves*k*nkn+i)+19];  A21 = &A11[4];  A22 = &A21[6];
      nlpr->cpsiuu[fN] -= A21[0]*tabeu[i] + A21[1]*tabev[i] +
                          A22[0]*tabeuu[i] + A22[1]*tabeuv[i] + A22[2]*tabevv[i];
      nlpr->cpsiuv[fN] -= A21[2]*tabeu[i] + A21[3]*tabev[i] +
                          A22[3]*tabeuu[i] + A22[4]*tabeuv[i] + A22[5]*tabevv[i];
      nlpr->cpsivv[fN] -= A21[4]*tabeu[i] + A21[5]*tabev[i] +
                          A22[6]*tabeuu[i] + A22[7]*tabeuv[i] + A22[8]*tabevv[i];
    }
    for ( i = 0, j = nakn-1, fN = (njcurves*(nfunc_a*hole_k+k)+1)*nkn;
          i < nkn;
          i++, j += nakn, fN++ ) {
      A11 = &ctrd[38*((njcurves*k+1)*nkn+i)];  A21 = &A11[4];  A22 = &A21[6];
      nlpr->cpsiu[fN] = A11[0]*tabfu[j] + A11[1]*tabfv[j];
      nlpr->cpsiv[fN] = A11[2]*tabfu[j] + A11[3]*tabfv[j];
      nlpr->cpsiuu[fN] = -(A21[0]*tabfu[j] + A21[1]*tabfv[j] +
                          A22[0]*tabfuu[j] + A22[1]*tabfuv[j] + A22[2]*tabfvv[j]);
      nlpr->cpsiuv[fN] = -(A21[2]*tabfu[j] + A21[3]*tabfv[j] +
                          A22[3]*tabfuu[j] + A22[4]*tabfuv[j] + A22[5]*tabfvv[j]);
      nlpr->cpsivv[fN] = -(A21[4]*tabfu[j] + A21[5]*tabfv[j] +
                          A22[6]*tabfuu[j] + A22[7]*tabfuv[j] + A22[8]*tabfvv[j]);
    }
    if ( njcurves > 3 ) {
        /* compute the derivatives and jumps inside Omega_i */
          /* v == const */
      for ( j = 1; j <= nk; j++ )
        for ( i = 0, fN = (njcurves*(nfunc_a*hole_k+k)+2+j)*nkn, l = j;
              i < nkn;
              i++, fN++, l += nakn ) {
          A11 = &ctrd[38*((njcurves*k+2+j)*nkn+i)];
          A21 = &A11[4];  A22 = &A21[6];
          nlpr->cpsiu[fN] = A11[0]*tabfu[l] + A11[1]*tabfv[l];
          nlpr->cpsiv[fN] = A11[2]*tabfu[l] + A11[3]*tabfv[l];
          nlpr->cpsiuu[fN] = (A21[0]-A21[0+19])*tabfu[l] + (A21[1]-A21[1+19])*tabfv[l] +
                     (A22[0]-A22[0+19])*tabfuu[l] + (A22[1]-A22[1+19])*tabfuv[l] +
                     (A22[2]-A22[2+19])*tabfvv[l];
          nlpr->cpsiuv[fN] = (A21[2]-A21[2+19])*tabfu[l] + (A21[3]-A21[3+19])*tabfv[l] +
                     (A22[3]-A22[3+19])*tabfuu[l] + (A22[4]-A22[4+19])*tabfuv[l] +
                     (A22[5]-A22[5+19])*tabfvv[l];
          nlpr->cpsivv[fN] = (A21[4]-A21[4+19])*tabfu[l] + (A21[5]-A21[5+19])*tabfv[l] +
                     (A22[6]-A22[6+19])*tabfuu[l] + (A22[7]-A22[7+19])*tabfuv[l] +
                     (A22[8]-A22[8+19])*tabfvv[l];
      }
    }
    if ( !mbs_TabBezC1CoonsDer2d ( 1, nakn-1, &aknots[1],
                &ahfunc[4], &adhfunc[4], &addhfunc[4],
                nkn, tkn, hfunc, dhfunc, ddhfunc,
                G1_CROSS00DEG, &bc00[(G1_CROSS00DEG+1)*k],
                G1_CROSS01DEG, &bc01[(G1_CROSS01DEG+1)*k],
                G1_CROSS10DEG, &bc10[(G1_CROSS10DEG+1)*k],
                G1_CROSS11DEG, &bc11[(G1_CROSS11DEG+1)*k],
                G1_CROSS00DEG, &bd00[(G1_CROSS00DEG+1)*k],
                G1_CROSS01DEG, &bd01[(G1_CROSS01DEG+1)*k],
                G1_CROSS10DEG, &bd10[(G1_CROSS10DEG+1)*k],
                G1_CROSS11DEG, &bd11[(G1_CROSS11DEG+1)*k],
                NULL, tabfu, tabfv, tabfuu, tabfuv, tabfvv ) )
      goto failure;
    for ( i = 0, fN = (njcurves*(nfunc_a*hole_k+k)+2)*nkn, l = (nakn-2)*nkn;
          i < nkn;
          i++, fN++, l++ ) {
      A11 = &ctrd[38*((njcurves*k+2)*nkn+i)];  A21 = &A11[4];  A22 = &A21[6];
      nlpr->cpsiu[fN] = A11[0]*tabfu[l] + A11[1]*tabfv[l];
      nlpr->cpsiv[fN] = A11[2]*tabfu[l] + A11[3]*tabfv[l];
      nlpr->cpsiuu[fN] -= A21[0]*tabfu[l] + A21[1]*tabfv[l] +
                          A22[0]*tabfuu[l] + A22[1]*tabfuv[l] + A22[2]*tabfvv[l];
      nlpr->cpsiuv[fN] -= A21[2]*tabfu[l] + A21[3]*tabfv[l] +
                          A22[3]*tabfuu[l] + A22[4]*tabfuv[l] + A22[5]*tabfvv[l];
      nlpr->cpsivv[fN] -= A21[4]*tabfu[l] + A21[5]*tabfv[l] +
                          A22[6]*tabfuu[l] + A22[7]*tabfuv[l] + A22[8]*tabfvv[l];
    }
    if ( njcurves > 3 ) {
        /* compute the derivatives and jumps inside Omega_i */
          /* u == const */
      for ( j = 1; j <= nk; j++ )
        for ( i = 0, fN = (njcurves*(nfunc_a*hole_k+k)+2+nk+j)*nkn, l = (j-1)*nkn;
              i < nkn;
              i++, fN++, l++ ) {
          A11 = &ctrd[38*((njcurves*k+2+nk+j)*nkn+i)];
          A21 = &A11[4];  A22 = &A21[6];
          nlpr->cpsiu[fN] = A11[0]*tabfu[l] + A11[1]*tabfv[l];
          nlpr->cpsiv[fN] = A11[2]*tabfu[l] + A11[3]*tabfv[l];
          nlpr->cpsiuu[fN] = (A21[0]-A21[0+19])*tabfu[l] + (A21[1]-A21[1+19])*tabfv[l] +
                     (A22[0]-A22[0+19])*tabfuu[l] + (A22[1]-A22[1+19])*tabfuv[l] +
                     (A22[2]-A22[2+19])*tabfvv[l];
          nlpr->cpsiuv[fN] = (A21[2]-A21[2+19])*tabfu[l] + (A21[3]-A21[3+19])*tabfv[l] +
                     (A22[3]-A22[3+19])*tabfuu[l] + (A22[4]-A22[4+19])*tabfuv[l] +
                     (A22[5]-A22[5+19])*tabfvv[l];
          nlpr->cpsivv[fN] = (A21[4]-A21[4+19])*tabfu[l] + (A21[5]-A21[5+19])*tabfv[l] +
                     (A22[6]-A22[6+19])*tabfuu[l] + (A22[7]-A22[7+19])*tabfuv[l] +
                     (A22[8]-A22[8+19])*tabfvv[l];
        }
    }
  }
          /* now deal with the derivatives outside Omega, */
          /* the arrays bc00 .. bd11 may be reused */
  memset ( bc00, 0, 2*(G1_CROSSDEGSUM+4)*hole_k*sizeof(double) );
  for ( fn = 0; fn < nfunc_b; fn++ ) {
    bvz = nlpr->rhole_cp[bfcpn[fn]].z;
    for ( k = 0, kk = 1;  k < hole_k;  k++, kk = (k+1) % hole_k ) {
      gh_GetDomSurrndBFuncd ( domain, fn, kk, 1, (double*)di );
      pkn_AddMatrixMd ( 1, 16, 0, &bc00[k*2*16], 0, (double*)di, bvz,
                        0, &bc00[k*2*16] );
      gh_GetDomSurrndBFuncd ( domain, fn, k, 2, (double*)di );
      pkn_AddMatrixMd ( 1, 16, 0, &bc00[(k*2+1)*16], 0, (double*)di, bvz,
                        0, &bc00[(k*2+1)*16] );
    }
  }  
  for ( k = 0; k < hole_k; k++ ) {
    for ( i = 0, fN = (njcurves*(nfunc_a*hole_k+k)+1)*nkn;  i < nkn;  i++, fN++ ) {
      if ( !mbs_BCHornerDer2Pd ( 3, 3, 1, &bc00[k*2*16], 0.0, tkn[i],
                                 &tabeu[1], tabeu, tabev, tabeuu, tabeuv, tabevv ) )
        goto failure;
      A11 = &ctrd[38*((njcurves*k+1)*nkn+i)+19];  A21 = &A11[4];  A22 = &A21[6];
      nlpr->cpsiuu[fN] += A21[0]*tabeu[0] + A21[1]*tabev[0] +
                          A22[0]*tabeuu[0] + A22[1]*tabeuv[0] + A22[2]*tabevv[0];
      nlpr->cpsiuv[fN] += A21[2]*tabeu[0] + A21[3]*tabev[0] +
                          A22[3]*tabeuu[0] + A22[4]*tabeuv[0] + A22[5]*tabevv[0];
      nlpr->cpsivv[fN] += A21[4]*tabeu[0] + A21[5]*tabev[0] +
                          A22[6]*tabeuu[0] + A22[7]*tabeuv[0] + A22[8]*tabevv[0];
    }
    for ( i = 0, fN = (njcurves*(nfunc_a*hole_k+k)+2)*nkn;  i < nkn;  i++, fN++ ) {
      if ( !mbs_BCHornerDer2Pd ( 3, 3, 1, &bc00[(k*2+1)*16], 0.0, tkn[i],
                                 &tabeu[1], tabeu, tabev, tabeuu, tabeuv, tabevv ) )
        goto failure;
      A11 = &ctrd[38*((njcurves*k+2)*nkn+i)+19];  A21 = &A11[4];  A22 = &A21[6];
      nlpr->cpsiuu[fN] += A21[0]*tabeu[0] + A21[1]*tabev[0] +
                          A22[0]*tabeuu[0] + A22[1]*tabeuv[0] + A22[2]*tabevv[0];
      nlpr->cpsiuv[fN] += A21[2]*tabeu[0] + A21[3]*tabev[0] +
                          A22[3]*tabeuu[0] + A22[4]*tabeuv[0] + A22[5]*tabevv[0];
      nlpr->cpsivv[fN] += A21[4]*tabeu[0] + A21[5]*tabev[0] +
                          A22[6]*tabeuu[0] + A22[7]*tabeuv[0] + A22[8]*tabevv[0];
    }
  }  
          /* block C functions */
  rr = G1H_FINALDEG-3+nk*m2;
  m = rr*nkn;
  tbs = pkv_GetScratchMemd ( 3*m );
  bstsize = rr*(nk+2);
  atbs = pkv_GetScratchMemd ( 4*bstsize );
  fkn = pkv_GetScratchMemi ( 2*rr );
  if ( !tbs || !atbs || !fkn ) {
    domain->error_code = G1H_ERROR_NO_SCRATCH_MEMORY;
    goto failure;
  }
  tbst = &tbs[m];  tbstt = &tbst[m];
  atbst = &atbs[bstsize];  atbstt0 = &atbst[bstsize];
  atbstt1 = &atbstt0[bstsize];
  lkn = &fkn[rr];
  if ( !_g1h_TabBSFuncDer2d ( G1H_FINALDEG, lastcknot, cknots, 2,
                 G1H_FINALDEG+nk*m2-2, nkn, tkn, fkn, lkn,
                 tbs, tbst, tbstt ) )
    goto failure;
  if ( !_g1h_TabBSFuncDer2Jd ( rr, G1H_FINALDEG, nk, m2, lastcknot, cknots,
                               atbs, atbst, atbstt0, atbstt1 ) )
    goto failure;
  for ( k = fn = 0, kk = 1, jcfs = nlspr->jcfs;
        k < hole_k;
        k++, kk = (k+1) % hole_k ) {
    for ( i = 0; i < rr; i++ )
      for ( j = 0;  j < rr;  j++, fn++, jcfs += (njcurves+1)*nkn )
        _g1hq2_TabSplCJumpd ( njcurves, nk, nkn,
             fkn[i], lkn[i], &tbs[i*nkn], &tbst[i*nkn], &tbstt[i*nkn],
             &atbs[i*(nk+2)], &atbst[i*(nk+2)],
             &atbstt0[i*(nk+2)], &atbstt1[i*(nk+2)],
             fkn[j], lkn[j], &tbs[j*nkn], &tbst[j*nkn], &tbstt[j*nkn],
             &atbs[j*(nk+2)], &atbst[j*(nk+2)],
             &atbstt0[j*(nk+2)], &atbstt1[j*(nk+2)],
             &ctrd[38*k*njcurves*nkn],
             &ctrd[38*kk*njcurves*nkn],
             &nlpr->cpsiu[jcfs], &nlpr->cpsiv[jcfs],
             &nlpr->cpsiuu[jcfs], &nlpr->cpsiuv[jcfs], &nlpr->cpsivv[jcfs] );
  }
          /* block D functions */
  lastomcknot = sprivate->lastomcknot;
  omcknots = sprivate->omcknots;
  lastpvknot = sprivate->lastpvknot;
  pvknots = sprivate->pvknots;
  fcomc = bc00;
  pv = &fcomc[lastomcknot-G1_CROSS00DEG];
  pu = &pv[lastpvknot-G1_CROSS01DEG];
  rr = 2*nk*m1;
  for ( k = fn = 0, kk = 1, ki = hole_k-1, jdfs = nlspr->jdfs;
        k < hole_k;
        ki = k++, kk = (k+1) % hole_k )
    for ( i = 0;  i < rr;  i++, fn++, jdfs += (2*njcurves+1)*nkn ) {
      _g1h_GetSplDBasisCrossDerd ( domain, nfunc_a+nfunc_b+nfunc_c+fn, k,
                                   fcomc, pv, pu );
      _g1h_FuncDSuppd ( hole_k, nk, m1, fn, k,
                        &nzc, &i00, &i01, &j00, &j01 );
      if ( !_g1hq2_TabSplDJumpd ( domain, njcurves, nk, m1,
                nkn, tkn, hfunc, dhfunc, ddhfunc,          
                ahfunc, adhfunc, addhfunc,
                nzc, i00, i01,
                lastomcknot, omcknots, fcomc, lastpvknot, pvknots, pv, pu,
                &ctrd[38*k*njcurves*nkn], &ctrd[38*kk*njcurves*nkn],
                &ctrd[38*ki*njcurves*nkn], 
                &nlpr->cpsiu[jdfs], &nlpr->cpsiv[jdfs],
                &nlpr->cpsiuu[jdfs], &nlpr->cpsiuv[jdfs], &nlpr->cpsivv[jdfs] ) )
        goto failure;
    }

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*_g1hq2_TabSplNLBasisFJumpsd*/

/* ///////////////////////////////////////////////////////////////////////// */
static double _g1hq2_ComputeSplNLFuncd ( GHoleDomaind *domain,
                                        G1HNLSPrivated *nlspr,
                                        const double *coeff,
                                        double C )
{
  G1HolePrivateRecd  *privateG1;
  G1HoleSPrivateRecd *sprivate;
  G1HNLPrivated      *nlpr;
  G2HNLFuncd         f;
  G1Q2HNLFuncd       cf;
  int   hole_k, nfunc_a, nfunc_c, nfunc_d, nfcd, njcurves;
  double c, funct, cfunct;
  int   nk, m1, m2, nkn, nkn2, fi, rrc, rrd, bs1, jcfs, jdfs;
  int   i, j, k, ki, kl, kk, kn, knot, fn,
        fic, fjc, i0, i1, j0, j1, nzc, nc, knc, ii;
  int   *fkn, *lkn, *cfuncvi, *dfuncvi;

  privateG1 = domain->privateG1;
  sprivate = domain->SprivateG1;
  nlpr = &nlspr->nlpr;
  hole_k = domain->hole_k;
  nfunc_a = privateG1->nfunc_a;
  nfunc_c = sprivate->nsfunc_c;
  nfunc_d = sprivate->nsfunc_d;
  nfcd = nfunc_c+nfunc_d;
  nk   = sprivate->nk;
  m1   = sprivate->m1;
  m2   = sprivate->m2;
  nkn  = nlspr->nkn;
  nkn2 = nkn*nkn;
  fkn = nlspr->fkn;
  lkn = nlspr->lkn;
  cfuncvi = nlspr->cfuncvi;
  dfuncvi = nlspr->dfuncvi;

  rrc = G1H_FINALDEG-3+nk*m2;
  bs1 = rrc*rrc;
  rrd = 2*nk*m1;
  funct = cfunct = 0.0;
        /* integration in the area Omega */
  for ( k = knot = 0, kk = 1;
        k < hole_k;
        k++, kk = (k+1) % hole_k ) {

        /* evaluate the function and its derivatives */
    for ( i = 0; i < nkn; i++ )
      for ( j = 0;  j < nkn; j++, knot++ ) {
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
          _g1h_FuncDSuppd ( hole_k, nk, m1, fn, k, 
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
          _g1h_FuncDSuppd ( hole_k, nk, m1, fn, k, 
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
        _g2h_IntFunc1ad ( &f, &funct );
      }
  }
  funct /= (double)nkn2;

        /* integration along the curves */
  njcurves = nlspr->njcurves;
  jcfs = nlspr->jcfs;
  jdfs = nlspr->jdfs;
  for ( k = 0, ki = hole_k-1, kl = 1;
        k < hole_k;
        ki = k++, kl = (k+1) % hole_k ) {
          /* first, integration along the curves which always must be */
          /* taken into account, i.e. boundaries of Omega_i */
    for ( nc = knc = 0, knot = k*njcurves*nkn;  nc < njcurves;  nc++ )
      for ( kn = 0;  kn < nkn;  kn++, knc++, knot++ ) {
        /* evaluate the function and its derivatives */
        fi = (nfunc_a*hole_k+k)*njcurves*nkn + knc;
        cf.pu = nlpr->cpsiu[fi];     cf.pv = nlpr->cpsiv[fi];
        cf.jpuu = nlpr->cpsiuu[fi];  cf.jpuv = nlpr->cpsiuv[fi];
        cf.jpvv = nlpr->cpsivv[fi];
        /* add the linear combination of the block A functions */
        for ( i = 0; i < nfunc_a; i++ ) {
          fi = (i*hole_k+k)*njcurves*nkn + knc;
          c = coeff[nfunc_c+nfunc_d+i];
          cf.pu -= c*nlpr->cpsiu[fi];     cf.pv -= c*nlpr->cpsiv[fi];
          cf.jpuu -= c*nlpr->cpsiuu[fi];  cf.jpuv -= c*nlpr->cpsiuv[fi];
          cf.jpvv -= c*nlpr->cpsivv[fi];
        }
        /* add the linear combination of the block C functions */
          /* the Omega_k block */
        for ( i = 0, ii = k*bs1;  i < bs1;  i++, ii++ ) {
          if ( nc < 3 )
            fi = jcfs + (ii*(njcurves+1)+nc)*nkn + kn;
          else
            fi = jcfs + (ii*(njcurves+1)+nc+1)*nkn + kn;
          c = coeff[ii];
          cf.pu -= c*nlpr->cpsiu[fi];     cf.pv -= c*nlpr->cpsiv[fi];
          cf.jpuu -= c*nlpr->cpsiuu[fi];  cf.jpuv -= c*nlpr->cpsiuv[fi];
          cf.jpvv -= c*nlpr->cpsivv[fi];
        }
          /* the Omega_{k-1} block */
        if ( nc == 0 )
          for ( i = 0, ii = ki*bs1;  i < rrc;  i++, ii++ ) {
            fi = jcfs + (ii*(njcurves+1)+3)*nkn + kn;
            c = coeff[ii];
            cf.pu -= c*nlpr->cpsiu[fi];     cf.pv -= c*nlpr->cpsiv[fi];
            cf.jpuu -= c*nlpr->cpsiuu[fi];  cf.jpuv -= c*nlpr->cpsiuv[fi];
            cf.jpvv -= c*nlpr->cpsivv[fi];
          }
        /* add the linear combination of the block D functions */
          /* the Omega_k block */
        for ( i = 0, ii = k*rrd+i;  i < rrd;  i++, ii++ ) {
          if ( nc < 3 )
            fi = jdfs + (ii*(2*njcurves+1)+nc)*nkn + kn;
          else
            fi = jdfs + (ii*(2*njcurves+1)+nc+4)*nkn + kn;
          c = coeff[nfunc_c+ii];
          cf.pu -= c*nlpr->cpsiu[fi];     cf.pv -= c*nlpr->cpsiv[fi];
          cf.jpuu -= c*nlpr->cpsiuu[fi];  cf.jpuv -= c*nlpr->cpsiuv[fi];
          cf.jpvv -= c*nlpr->cpsivv[fi];
        }
          /* the Omega_{k-1} block */
        if ( nc == 0 )
          for ( i = 0, ii = ki*rrd+i;  i < rrd;  i++, ii++ ) {
            fi = jdfs + (ii*(2*njcurves+1)+3)*nkn + kn;
            c = coeff[nfunc_c+ii];
            cf.pu -= c*nlpr->cpsiu[fi];     cf.pv -= c*nlpr->cpsiv[fi];
            cf.jpuu -= c*nlpr->cpsiuu[fi];  cf.jpuv -= c*nlpr->cpsiuv[fi];
            cf.jpvv -= c*nlpr->cpsivv[fi];
          }
          /* the Omega_{k+1} block */
        for ( i = 0, ii = kl*rrd+i;  i < rrd;  i++, ii++ ) {
          if ( nc < 3 )
            fi = jdfs + (ii*(2*njcurves+1)+nc+4)*nkn + kn;
          else
            fi = jdfs + (ii*(2*njcurves+1)+2*nk+nc+4)*nkn + kn;
          c = coeff[nfunc_c+ii];
          cf.pu -= c*nlpr->cpsiu[fi];     cf.pv -= c*nlpr->cpsiv[fi];
          cf.jpuu -= c*nlpr->cpsiuu[fi];  cf.jpuv -= c*nlpr->cpsiuv[fi];
          cf.jpvv -= c*nlpr->cpsivv[fi];
        }

        cf.tang = nlpr->ctang[knot];
        /* integrate the functional value */
        _g1hq2_IntFunc1ad ( &cf, &cfunct );
      }
    }
  cfunct /= (double)nkn;

  return funct + C*cfunct;
} /*_g1hq2_ComputeSplNLFuncd*/

static boolean _g1hq2_ComputeSplNLFuncGradd ( GHoleDomaind *domain,
                                              G1HNLSPrivated *nlspr,
                                              const double *coeff,
                                              double C,
                                              double *func, double *grad )
{
  void   *sp;
  G1HolePrivateRecd  *privateG1;
  G1HoleSPrivateRecd *sprivate; 
  G1HNLPrivated      *nlpr;
  G2HNLFuncd         f;
  G1Q2HNLFuncd       cf;
  double c, funct, cfunct, *cgrad;
  int    hole_k, nfunc_a, nfunc_c, nfunc_d, nfcd, nfacd, njcurves;
  int    nk, m1, m2, nkn, nkn2, fi, rrc, rrd, bs1, jcfs, jdfs,
         fic, fjc, fn, i0, i1, j0, j1, nzc;
  int    *fkn, *lkn, *cfuncvi, *dfuncvi;
  int    i, j, k, kk, ki, kl, kn, knot, nc, knc, ii;

  sp = pkv_GetScratchMemTop ();
  privateG1 = domain->privateG1;
  sprivate = domain->SprivateG1;
  nlpr = &nlspr->nlpr;
  hole_k = domain->hole_k;
  nfunc_a = privateG1->nfunc_a;
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

  rrc = G1H_FINALDEG-3+nk*m2;
  bs1 = rrc*rrc;
  rrd = 2*nk*m1;

  cgrad = pkv_GetScratchMemd ( nfacd );
  if ( !cgrad ) {
    domain->error_code = G1H_ERROR_NO_SCRATCH_MEMORY;
    goto failure;
  }
  funct = cfunct = 0.0;
  memset ( grad, 0, nfacd*sizeof(double) );
  memset ( cgrad, 0, nfacd*sizeof(double) );

        /* integration in the area Omega */
  for ( k = knot = 0, kk = 1;
        k < hole_k;
        k++, kk = (k+1) % hole_k ) {

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
          _g1h_FuncDSuppd ( hole_k, nk, m1, fn, k, 
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
          _g1h_FuncDSuppd ( hole_k, nk, m1, fn, k, 
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
        _g2h_IntFunc1bd ( &f, &funct );

/* 2. Integration of the functional gradient */
          /* derivatives with respect to the block A coefficients */
        for ( fn = 0; fn < nfunc_a; fn++ ) {
          fi = fn*nkn2*hole_k + knot;
          f.psiu = nlpr->psiu[fi];      f.psiv = nlpr->psiv[fi];
          f.psiuu = nlpr->psiuu[fi];    f.psiuv = nlpr->psiuv[fi];
          f.psivv = nlpr->psivv[fi];    f.psiuuu = nlpr->psiuuu[fi];
          f.psiuuv = nlpr->psiuuv[fi];  f.psiuvv = nlpr->psiuvv[fi];
          f.psivvv = nlpr->psivvv[fi];
          _g2h_IntFunc2bd ( &f, &grad[nfcd+fn] );
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
                _g2h_IntFunc2bd ( &f, &grad[fn] );
              }
          /* derivatives with respect to the block D coefficients */
        for ( fic = 0; fic < rrd; fic++ ) {
          fn = k*rrd + fic;   /* function number */
          _g1h_FuncDSuppd ( hole_k, nk, m1, fn, k, 
                            &nzc, &i0, &i1, &j0, &j1 );
          if ( i >= i0 && i <= i1 && j >= j0 && j <= j1 ) {
            fi = dfuncvi[fn] + (i-i0)*(j1-j0+1) + (j-j0);  
            f.psiu = nlpr->psiu[fi];      f.psiv = nlpr->psiv[fi];
            f.psiuu = nlpr->psiuu[fi];    f.psiuv = nlpr->psiuv[fi];
            f.psivv = nlpr->psivv[fi];    f.psiuuu = nlpr->psiuuu[fi];
            f.psiuuv = nlpr->psiuuv[fi];  f.psiuvv = nlpr->psiuvv[fi];
            f.psivvv = nlpr->psivvv[fi];
            _g2h_IntFunc2bd ( &f, &grad[nfunc_c+fn] );
          }
          fn = kk*rrd + fic;  /* function number */
          _g1h_FuncDSuppd ( hole_k, nk, m1, fn, k, 
                            &nzc, &i0, &i1, &j0, &j1 );
          if ( i >= i0 && i <= i1 && j >= j0 && j <= j1 ) {
            fi = dfuncvi[fn] + (i1-i0+1)*(j1-j0+1) + (i-i0)*(j1-j0+1) + (j-j0);
            f.psiu = nlpr->psiu[fi];      f.psiv = nlpr->psiv[fi];
            f.psiuu = nlpr->psiuu[fi];    f.psiuv = nlpr->psiuv[fi];
            f.psivv = nlpr->psivv[fi];    f.psiuuu = nlpr->psiuuu[fi];
            f.psiuuv = nlpr->psiuuv[fi];  f.psiuvv = nlpr->psiuvv[fi];
            f.psivvv = nlpr->psivvv[fi];
            _g2h_IntFunc2bd ( &f, &grad[nfunc_c+fn] );
          }
        }  
      }
  }
  funct /= (double)nkn2;
  pkn_MultMatrixNumd ( 1, nfacd, 0, grad, 1.0/(double)nkn2, 0, grad );

        /* integration along the curves */
  njcurves = nlspr->njcurves;
  jcfs = nlspr->jcfs;
  jdfs = nlspr->jdfs;
  for ( k = 0, ki = hole_k-1, kl = 1;
        k < hole_k;
        ki = k++, kl = (k+1) % hole_k ) {
          /* first, integration along the curves which always must be */
          /* taken into account, i.e. boundaries of Omega_i */
    for ( nc = knc = 0, knot = k*njcurves*nkn;  nc < njcurves;  nc++)
      for ( kn = 0;  kn < nkn;  kn++, knc++, knot++ ) {
        /* evaluate the function and its derivatives */
        fi = (nfunc_a*hole_k+k)*njcurves*nkn + knc;
        cf.pu = nlpr->cpsiu[fi];     cf.pv = nlpr->cpsiv[fi];
        cf.jpuu = nlpr->cpsiuu[fi];  cf.jpuv = nlpr->cpsiuv[fi];
        cf.jpvv = nlpr->cpsivv[fi];
        /* add the linear combination of the block A functions */
        for ( i = 0; i < nfunc_a; i++ ) {
          fi = (i*hole_k+k)*njcurves*nkn + knc;
          c = coeff[nfunc_c+nfunc_d+i];
          cf.pu -= c*nlpr->cpsiu[fi];     cf.pv -= c*nlpr->cpsiv[fi];
          cf.jpuu -= c*nlpr->cpsiuu[fi];  cf.jpuv -= c*nlpr->cpsiuv[fi];
          cf.jpvv -= c*nlpr->cpsivv[fi];
        }
        /* add the linear combination of the block C functions */
          /* the Omega_k block */
        for ( i = 0, ii = k*bs1;  i < bs1;  i++, ii++ ) {
          if ( nc < 3 )
            fi = jcfs + (ii*(njcurves+1)+nc)*nkn + kn;
          else
            fi = jcfs + (ii*(njcurves+1)+nc+1)*nkn + kn;
          c = coeff[ii];
          cf.pu -= c*nlpr->cpsiu[fi];     cf.pv -= c*nlpr->cpsiv[fi];
          cf.jpuu -= c*nlpr->cpsiuu[fi];  cf.jpuv -= c*nlpr->cpsiuv[fi];
          cf.jpvv -= c*nlpr->cpsivv[fi];
        }
          /* the Omega_{k-1} block */
        if ( nc == 0 )
          for ( i = 0, ii = ki*bs1;  i < rrc;  i++, ii++ ) {
            fi = jcfs + (ii*(njcurves+1)+3)*nkn + kn;
            c = coeff[ii];
            cf.pu -= c*nlpr->cpsiu[fi];     cf.pv -= c*nlpr->cpsiv[fi];
            cf.jpuu -= c*nlpr->cpsiuu[fi];  cf.jpuv -= c*nlpr->cpsiuv[fi];
            cf.jpvv -= c*nlpr->cpsivv[fi];
          }
        /* add the linear combination of the block D functions */
          /* the Omega_k block */
        for ( i = 0, ii = k*rrd+i;  i < rrd;  i++, ii++ ) {
          if ( nc < 3 )
            fi = jdfs + (ii*(2*njcurves+1)+nc)*nkn + kn;
          else
            fi = jdfs + (ii*(2*njcurves+1)+nc+4)*nkn + kn;
          c = coeff[nfunc_c+ii];
          cf.pu -= c*nlpr->cpsiu[fi];     cf.pv -= c*nlpr->cpsiv[fi];
          cf.jpuu -= c*nlpr->cpsiuu[fi];  cf.jpuv -= c*nlpr->cpsiuv[fi];
          cf.jpvv -= c*nlpr->cpsivv[fi];
        }
          /* the Omega_{k-1} block */
        if ( nc == 0 )
          for ( i = 0, ii = ki*rrd+i;  i < rrd;  i++, ii++ ) {
            fi = jdfs + (ii*(2*njcurves+1)+3)*nkn + kn;
            c = coeff[nfunc_c+ii];
            cf.pu -= c*nlpr->cpsiu[fi];     cf.pv -= c*nlpr->cpsiv[fi];
            cf.jpuu -= c*nlpr->cpsiuu[fi];  cf.jpuv -= c*nlpr->cpsiuv[fi];
            cf.jpvv -= c*nlpr->cpsivv[fi];
          }
          /* the Omega_{k+1} block */
        for ( i = 0, ii = kl*rrd+i;  i < rrd;  i++, ii++ ) {
          if ( nc < 3 )
            fi = jdfs + (ii*(2*njcurves+1)+nc+4)*nkn + kn;
          else
            fi = jdfs + (ii*(2*njcurves+1)+2*nk+nc+4)*nkn + kn;
          c = coeff[nfunc_c+ii];
          cf.pu -= c*nlpr->cpsiu[fi];     cf.pv -= c*nlpr->cpsiv[fi];
          cf.jpuu -= c*nlpr->cpsiuu[fi];  cf.jpuv -= c*nlpr->cpsiuv[fi];
          cf.jpvv -= c*nlpr->cpsivv[fi];
        }

        cf.tang = nlpr->ctang[knot];
        /* integrate the functional value */
        _g1hq2_IntFunc1bd ( &cf, &cfunct );
        /* integrate the gradient */
          /* block A functions */
        for ( i = 0; i < nfunc_a; i++ ) {
          fi = (i*hole_k+k)*njcurves*nkn + knc;
          ii = nfunc_c+nfunc_d+i;
          cf.psiu = nlpr->cpsiu[fi];     cf.psiv = nlpr->cpsiv[fi];
          cf.jpsiuu = nlpr->cpsiuu[fi];  cf.jpsiuv = nlpr->cpsiuv[fi];
          cf.jpsivv = nlpr->cpsivv[fi];
          _g1hq2_IntFunc2bd ( &cf, &cgrad[ii] );
        }
          /* block C functions */
            /* the Omega_k block */
        for ( i = 0, ii = k*bs1;  i < bs1;  i++, ii++ ) {
          if ( nc < 3 )
            fi = jcfs + (ii*(njcurves+1)+nc)*nkn + kn;
          else
            fi = jcfs + (ii*(njcurves+1)+nc+1)*nkn + kn;
          cf.psiu = nlpr->cpsiu[fi];     cf.psiv = nlpr->cpsiv[fi];
          cf.jpsiuu = nlpr->cpsiuu[fi];  cf.jpsiuv = nlpr->cpsiuv[fi];
          cf.jpsivv = nlpr->cpsivv[fi];
          _g1hq2_IntFunc2bd ( &cf, &cgrad[ii] );
        }
            /* the Omega_{k+1} block */
        if ( nc == 0 )
          for ( i = 0, ii = ki*bs1;  i < rrc;  i++, ii++ ) {
            fi = jcfs + (ii*(njcurves+1)+3)*nkn + kn;
            cf.psiu = nlpr->cpsiu[fi];     cf.psiv = nlpr->cpsiv[fi];
            cf.jpsiuu = nlpr->cpsiuu[fi];  cf.jpsiuv = nlpr->cpsiuv[fi];
            cf.jpsivv = nlpr->cpsivv[fi];
            _g1hq2_IntFunc2bd ( &cf, &cgrad[ii] );
          }
          /* block D functions */
            /* the Omega_k block */
        for ( i = 0, ii = k*rrd+i;  i < rrd;  i++, ii++ ) {
          if ( nc < 3 )
            fi = jdfs + (ii*(2*njcurves+1)+nc)*nkn + kn;
          else
            fi = jdfs + (ii*(2*njcurves+1)+nc+4)*nkn + kn;
          cf.psiu = nlpr->cpsiu[fi];     cf.psiv = nlpr->cpsiv[fi];
          cf.jpsiuu = nlpr->cpsiuu[fi];  cf.jpsiuv = nlpr->cpsiuv[fi];
          cf.jpsivv = nlpr->cpsivv[fi];
          _g1hq2_IntFunc2bd ( &cf, &cgrad[nfunc_c+ii] );
        }
            /* the Omega_{k-1} block */
        if ( nc == 0 )
          for ( i = 0, ii = ki*rrd+i;  i < rrd;  i++, ii++ ) {
            fi = jdfs + (ii*(2*njcurves+1)+3)*nkn + kn;
            cf.psiu = nlpr->cpsiu[fi];     cf.psiv = nlpr->cpsiv[fi];
            cf.jpsiuu = nlpr->cpsiuu[fi];  cf.jpsiuv = nlpr->cpsiuv[fi];
            cf.jpsivv = nlpr->cpsivv[fi];
            _g1hq2_IntFunc2bd ( &cf, &cgrad[nfunc_c+ii] );
          }
            /* the Omega_{k+1} block */
        for ( i = 0, ii = kl*rrd+i;  i < rrd;  i++, ii++ ) {
          if ( nc < 3 )
            fi = jdfs + (ii*(2*njcurves+1)+nc+4)*nkn + kn;
          else
            fi = jdfs + (ii*(2*njcurves+1)+2*nk+nc+4)*nkn + kn;
          cf.psiu = nlpr->cpsiu[fi];     cf.psiv = nlpr->cpsiv[fi];
          cf.jpsiuu = nlpr->cpsiuu[fi];  cf.jpsiuv = nlpr->cpsiuv[fi];
          cf.jpsivv = nlpr->cpsivv[fi];
          _g1hq2_IntFunc2bd ( &cf, &cgrad[nfunc_c+ii] );
        }
      }
    }
  cfunct /= (double)nkn;

  *func = funct + C*cfunct;
  pkn_AddMatrixMd ( 1, nfacd, 0, grad, 0, cgrad, C/(double)nkn, 0, grad );

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*_g1hq2_ComputeSplNLFuncGradd*/

static boolean _g1hq2_ComputeSplNLFuncGradHessiand ( GHoleDomaind *domain,
                                G1HNLSPrivated *nlspr,
                                const double *coeff,
                                double C,
                                double *func, double *grad, double *hessian,
                                double *chessian )
{
  void *sp;
  G1HolePrivateRecd  *privateG1;
  G1HoleSPrivateRecd *sprivate; 
  G1HNLPrivated      *nlpr;
  G2HNLFuncd         f;
  G1Q2HNLFuncd       cf;
  int   hole_k, nfunc_a, nfunc_c, nfunc_d, nfcd, nfacd, njcurves;
  int   hsize;
  double *Li, *Bi, *BiLT, *Di, *hii, *hki, *hkk, *hkj, *hij, *hjj, *hh;
  double c, funct, cfunct, *cgrad;
  int   *fkn, *lkn, *cfuncvi, *dfuncvi;
  int   nk, m1, m2, nkn, nkn2, fi, fj, rrc, rrd, bs1, bs2,
        fic, fjc, fn, fim, fjm, fm, pos,
        i00, i01, j00, j01, i10, i11, j10, j11, nzc;
  int   k, kk, ki, kl, kn, knot, i, j, ii, jj, nc, knc, jcfs, jdfs;

  sp = pkv_GetScratchMemTop ();
  privateG1 = domain->privateG1;
  sprivate = domain->SprivateG1;
  nlpr = &nlspr->nlpr;
  hole_k = domain->hole_k;
  nfunc_a = privateG1->nfunc_a;
  nfunc_c = sprivate->nsfunc_c;
  nfunc_d = sprivate->nsfunc_d;
  nfcd  = nfunc_c+nfunc_d;
  nfacd = nfcd+nfunc_a;   
  Li = pkv_GetScratchMemd ( 8*nfacd );
  cgrad = pkv_GetScratchMemd ( nfacd );
  if ( !Li || !cgrad ) {
    domain->error_code = G1H_ERROR_NO_SCRATCH_MEMORY;
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

  rrc = G1H_FINALDEG-3+nk*m2;
  bs1 = rrc*rrc;
  rrd = 2*nk*m1;
  bs2 = nfunc_a+nfunc_d+bs1;
  hsize = pkn_Block3ArraySize ( hole_k-1, bs1, bs2 );
  hkk = &hessian[pkn_Block3FindBlockPos(hole_k-1,bs1,bs2,hole_k-1,hole_k-1)];

  funct = cfunct = 0.0;
  memset ( grad, 0, nfacd*sizeof(double) );
  memset ( cgrad, 0, nfacd*sizeof(double) );
  memset ( hessian, 0, hsize*sizeof(double) );
  memset ( chessian, 0, hsize*sizeof(double) );

        /* integration in the area Omega */
  for ( k = knot = 0, kk = 1;
        k < hole_k;
        k++, kk = (k+1) % hole_k ) {
    hii = &hessian[pkn_Block3FindBlockPos(hole_k-1,bs1,bs2,k,k)];
    hki = &hessian[pkn_Block3FindBlockPos(hole_k-1,bs1,bs2,hole_k-1,k)];

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
          _g1h_FuncDSuppd ( hole_k, nk, m1, fn, k, 
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
          _g1h_FuncDSuppd ( hole_k, nk, m1, fn, k, 
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
        _g2h_IntFunc1cd ( &f, &funct );

/* 2. Integration of the functional gradient */
          /* derivatives with respect to the block A coefficients */
        for ( fn = 0; fn < nfunc_a; fn++ ) {
          fi = fn*nkn2*hole_k + knot;
          f.psiu = nlpr->psiu[fi];      f.psiv = nlpr->psiv[fi];
          f.psiuu = nlpr->psiuu[fi];    f.psiuv = nlpr->psiuv[fi];
          f.psivv = nlpr->psivv[fi];    f.psiuuu = nlpr->psiuuu[fi];
          f.psiuuv = nlpr->psiuuv[fi];  f.psiuvv = nlpr->psiuvv[fi];
          f.psivvv = nlpr->psivvv[fi];
          _g2h_IntFunc2cd ( &f, (vector2d*)&Li[2*(nfcd+fn)],
                            &Bi[3*(nfcd+fn)], (vector2d*)&BiLT[2*(nfcd+fn)],
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
                _g2h_IntFunc2cd ( &f, (vector2d*)&Li[2*fn],
                                  &Bi[3*fn], (vector2d*)&BiLT[2*fn],
                                  &Di[fn], &grad[fn] );
              }
          /* derivatives with respect to the block D coefficients */
        for ( fic = 0; fic < rrd; fic++ ) {
          fn = k*rrd + fic;   /* function number */
          _g1h_FuncDSuppd ( hole_k, nk, m1, fn, k,
                            &nzc, &i00, &i01, &j00, &j01 );
          if ( i >= i00 && i <= i01 && j >= j00 && j <= j01 ) {
            fi = dfuncvi[fn] + (i-i00)*(j01-j00+1) + (j-j00);
            f.psiu = nlpr->psiu[fi];      f.psiv = nlpr->psiv[fi];
            f.psiuu = nlpr->psiuu[fi];    f.psiuv = nlpr->psiuv[fi];
            f.psivv = nlpr->psivv[fi];    f.psiuuu = nlpr->psiuuu[fi];
            f.psiuuv = nlpr->psiuuv[fi];  f.psiuvv = nlpr->psiuvv[fi];
            f.psivvv = nlpr->psivvv[fi];
            _g2h_IntFunc2cd ( &f, (vector2d*)&Li[2*(nfunc_c+fn)],
                              &Bi[3*(nfunc_c+fn)], (vector2d*)&BiLT[2*(nfunc_c+fn)],
                              &Di[nfunc_c+fn], &grad[nfunc_c+fn] );
          }
          fn = kk*rrd + fic;  /* function number */
          _g1h_FuncDSuppd ( hole_k, nk, m1, fn, k,
                            &nzc, &i00, &i01, &j00, &j01 );
          if ( i >= i00 && i <= i01 && j >= j00 && j <= j01 ) {
            fi = dfuncvi[fn] + (i01-i00+1)*(j01-j00+1) +
                               (i-i00)*(j01-j00+1) + (j-j00);
            f.psiu = nlpr->psiu[fi];      f.psiv = nlpr->psiv[fi];
            f.psiuu = nlpr->psiuu[fi];    f.psiuv = nlpr->psiuv[fi];
            f.psivv = nlpr->psivv[fi];    f.psiuuu = nlpr->psiuuu[fi];
            f.psiuuv = nlpr->psiuuv[fi];  f.psiuvv = nlpr->psiuvv[fi];
            f.psivvv = nlpr->psivvv[fi];
            _g2h_IntFunc2cd ( &f, (vector2d*)&Li[2*(nfunc_c+fn)],
                              &Bi[3*(nfunc_c+fn)], (vector2d*)&BiLT[2*(nfunc_c+fn)],
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
            pos = pkn_SymMatIndex ( bs1+nfunc_d+fn, bs1+nfunc_d+fm );
            _g2h_IntFunc3cd ( &f, (vector2d*)&Li[2*(nfcd+fn)],
                              (vector2d*)&Li[2*(nfcd+fm)],
                              (vector2d*)&BiLT[2*(nfcd+fn)],
                              (vector2d*)&BiLT[2*(nfcd+fm)],
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
                  if ( k < hole_k-1 )
                    hh = &hki[(bs1+nfunc_d+fn)*bs1 + fic*rrc + fjc];
                  else
                    hh = &hkk[pkn_SymMatIndex( bs1+nfunc_d+fn, fic*rrc+fjc )];
                  _g2h_IntFunc3cd ( &f, (vector2d*)&Li[2*(nfcd+fn)],
                                    (vector2d*)&Li[2*fm],
                                    (vector2d*)&BiLT[2*(nfcd+fn)],
                                    (vector2d*)&BiLT[2*fm],
                                    Di[nfcd+fn], Di[fm], hh );
                }
          /* A x D */
            /* block D functions related with Gamma_i */
          for ( fic = 0, fm = k*rrd;  fic < rrd;  fic++, fm++ ) {
            _g1h_FuncDSuppd ( hole_k, nk, m1, fm, k,
                              &nzc, &i10, &i11, &j10, &j11 );
            if ( i >= i10 && i <= i11 && j >= j10 && j <= j11 ) {
              fj = dfuncvi[fm] + (i-i10)*(j11-j10+1) + j-j10;
              f.psju = nlpr->psiu[fj];      f.psjv = nlpr->psiv[fj];
              f.psjuu = nlpr->psiuu[fj];    f.psjuv = nlpr->psiuv[fj];
              f.psjvv = nlpr->psivv[fj];    f.psjuuu = nlpr->psiuuu[fj];
              f.psjuuv = nlpr->psiuuv[fj];  f.psjuvv = nlpr->psiuvv[fj];
              f.psjvvv = nlpr->psivvv[fj];
              pos = pkn_SymMatIndex ( bs1+nfunc_d+fn, bs1+fm );
              _g2h_IntFunc3cd ( &f, (vector2d*)&Li[2*(nfcd+fn)],
                              (vector2d*)&Li[2*(nfunc_c+fm)],
                              (vector2d*)&BiLT[2*(nfcd+fn)],
                              (vector2d*)&BiLT[2*(nfunc_c+fm)],
                              Di[nfcd+fn], Di[nfunc_c+fm], &hkk[pos] );
            }
          }
            /* block D functions related with Gamma_{i+1} */
          for ( fic = 0, fm = kk*rrd;  fic < rrd;  fic++, fm++ ) {
            _g1h_FuncDSuppd ( hole_k, nk, m1, fm, k,
                              &nzc, &i10, &i11, &j10, &j11 );
            if ( i >= i10 && i <= i11 && j >= j10 && j <= j11 ) {
              fj = dfuncvi[fm] + (i11-i10+1)*(j11-j10+1) + (i-i10)*(j11-j10+1) + j-j10;
              f.psju = nlpr->psiu[fj];      f.psjv = nlpr->psiv[fj];
              f.psjuu = nlpr->psiuu[fj];    f.psjuv = nlpr->psiuv[fj];
              f.psjvv = nlpr->psivv[fj];    f.psjuuu = nlpr->psiuuu[fj];
              f.psjuuv = nlpr->psiuuv[fj];  f.psjuvv = nlpr->psiuvv[fj];
              f.psjvvv = nlpr->psivvv[fj];
              pos = pkn_SymMatIndex ( bs1+nfunc_d+fn, bs1+fm );
              _g2h_IntFunc3cd ( &f, (vector2d*)&Li[2*(nfcd+fn)],
                              (vector2d*)&Li[2*(nfunc_c+fm)],
                              (vector2d*)&BiLT[2*(nfcd+fn)],
                              (vector2d*)&BiLT[2*(nfunc_c+fm)],
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
                          _g2h_IntFunc3cd ( &f, (vector2d*)&Li[2*fn],
                              (vector2d*)&Li[2*fm],
                              (vector2d*)&BiLT[2*fn],
                              (vector2d*)&BiLT[2*fm],
                              Di[fn], Di[fm], &hii[pos] );
                        }
                        else goto cont_c;
                      }
cont_c: ;
          /* C x D */
            /* block D functions related with Gamma_i */
                for ( fim = 0, fm = k*rrd;  fim < rrd;  fim++, fm++ ) {
                  _g1h_FuncDSuppd ( hole_k, nk, m1, fm, k,
                                    &nzc, &i10, &i11, &j10, &j11 );
                  if ( i >= i10 && i <= i11 && j >= j10 && j <= j11 ) {
                    fj = dfuncvi[fm] + (i-i10)*(j11-j10+1) + j-j10;
                    f.psju = nlpr->psiu[fj];      f.psjv = nlpr->psiv[fj];
                    f.psjuu = nlpr->psiuu[fj];    f.psjuv = nlpr->psiuv[fj];
                    f.psjvv = nlpr->psivv[fj];    f.psjuuu = nlpr->psiuuu[fj];
                    f.psjuuv = nlpr->psiuuv[fj];  f.psjuvv = nlpr->psiuvv[fj];
                    f.psjvvv = nlpr->psivvv[fj];
                    if ( k < hole_k-1 )
                      hh = &hki[bs1*(bs1+fm) + fic*rrc + fjc];
                    else
                      hh = &hkk[pkn_SymMatIndex( bs1+fm, fic*rrc+fjc )];
                    _g2h_IntFunc3cd ( &f, (vector2d*)&Li[2*fn],
                              (vector2d*)&Li[2*(nfunc_c+fm)],  
                              (vector2d*)&BiLT[2*fn],
                              (vector2d*)&BiLT[2*(nfunc_c+fm)],
                              Di[fn], Di[nfunc_c+fm], hh );
                  }
                }
            /* block D functions related with Gamma_{i+1} */
                for ( fim = 0, fm = kk*rrd;  fim < rrd;  fim++, fm++ ) {
                  _g1h_FuncDSuppd ( hole_k, nk, m1, fm, k,
                                    &nzc, &i10, &i11, &j10, &j11 );
                  if ( i >= i10 && i <= i11 && j >= j10 && j <= j11 ) {
                    fj = dfuncvi[fm] + (i11-i10+1)*(j11-j10+1) + (i-i10)*(j11-j10+1) + j-j10;
                    f.psju = nlpr->psiu[fj];      f.psjv = nlpr->psiv[fj];
                    f.psjuu = nlpr->psiuu[fj];    f.psjuv = nlpr->psiuv[fj];
                    f.psjvv = nlpr->psivv[fj];    f.psjuuu = nlpr->psiuuu[fj];
                    f.psjuuv = nlpr->psiuuv[fj];  f.psjuvv = nlpr->psiuvv[fj];
                    f.psjvvv = nlpr->psivvv[fj];
                    if ( k < hole_k-1 )
                      hh = &hki[bs1*(bs1+fm) + fic*rrc + fjc];
                    else
                      hh = &hkk[pkn_SymMatIndex( bs1+fm, fic*rrc+fjc )];
                    _g2h_IntFunc3cd ( &f, (vector2d*)&Li[2*fn],
                              (vector2d*)&Li[2*(nfunc_c+fm)],
                              (vector2d*)&BiLT[2*fn],
                              (vector2d*)&BiLT[2*(nfunc_c+fm)],
                              Di[fn], Di[nfunc_c+fm], hh );
                  }
                }
              }
          /* D x D */
        for ( fic = 0, fn = k*rrd;  fic < rrd;  fic++, fn++ ) {
          _g1h_FuncDSuppd ( hole_k, nk, m1, fn, k,
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
                _g1h_FuncDSuppd ( hole_k, nk, m1, fm, k,
                                  &nzc, &i10, &i11, &j10, &j11 );
                if ( i >= i10 && i <= i11 && j >= j10 && j <= j11 ) {
                  fj = dfuncvi[fm] + (i-i10)*(j11-j10+1) + j-j10;
                  f.psju = nlpr->psiu[fj];      f.psjv = nlpr->psiv[fj];
                  f.psjuu = nlpr->psiuu[fj];    f.psjuv = nlpr->psiuv[fj];
                  f.psjvv = nlpr->psivv[fj];    f.psjuuu = nlpr->psiuuu[fj];
                  f.psjuuv = nlpr->psiuuv[fj];  f.psjuvv = nlpr->psiuvv[fj];
                  f.psjvvv = nlpr->psivvv[fj];
                  pos = pkn_SymMatIndex ( bs1+fn, bs1+fm );
                  _g2h_IntFunc3cd ( &f, (vector2d*)&Li[2*(nfunc_c+fn)],
                              (vector2d*)&Li[2*(nfunc_c+fm)],
                              (vector2d*)&BiLT[2*(nfunc_c+fn)],
                              (vector2d*)&BiLT[2*(nfunc_c+fm)],
                              Di[nfunc_c+fn], Di[nfunc_c+fm], &hkk[pos] );
                }
              }
            for ( fim = 0, fm = kk*rrd;  fim < rrd;  fim++, fm++ )
              if ( fm <= fn ) {
                _g1h_FuncDSuppd ( hole_k, nk, m1, fm, k,
                                  &nzc, &i10, &i11, &j10, &j11 );
                if ( i >= i10 && i <= i11 && j >= j10 && j <= j11 ) {
                  fj = dfuncvi[fm] + (i11-i10+1)*(j11-j10+1) + (i-i10)*(j11-j10+1) + j-j10;
                  f.psju = nlpr->psiu[fj];      f.psjv = nlpr->psiv[fj];
                  f.psjuu = nlpr->psiuu[fj];    f.psjuv = nlpr->psiuv[fj];
                  f.psjvv = nlpr->psivv[fj];    f.psjuuu = nlpr->psiuuu[fj];
                  f.psjuuv = nlpr->psiuuv[fj];  f.psjuvv = nlpr->psiuvv[fj];
                  f.psjvvv = nlpr->psivvv[fj];
                  pos = pkn_SymMatIndex ( bs1+fn, bs1+fm );
                  _g2h_IntFunc3cd ( &f, (vector2d*)&Li[2*(nfunc_c+fn)],
                              (vector2d*)&Li[2*(nfunc_c+fm)],
                              (vector2d*)&BiLT[2*(nfunc_c+fn)],
                              (vector2d*)&BiLT[2*(nfunc_c+fm)],
                              Di[nfunc_c+fn], Di[nfunc_c+fm], &hkk[pos] );
                }
              }
          }
        }
        for ( fic = 0, fn = kk*rrd;  fic < rrd;  fic++, fn++ ) {
          _g1h_FuncDSuppd ( hole_k, nk, m1, fn, k,
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
                _g1h_FuncDSuppd ( hole_k, nk, m1, fm, k,
                                  &nzc, &i10, &i11, &j10, &j11 );
                if ( i >= i10 && i <= i11 && j >= j10 && j <= j11 ) {
                  fj = dfuncvi[fm] + (i-i10)*(j11-j10+1) + j-j10;
                  f.psju = nlpr->psiu[fj];      f.psjv = nlpr->psiv[fj];
                  f.psjuu = nlpr->psiuu[fj];    f.psjuv = nlpr->psiuv[fj];
                  f.psjvv = nlpr->psivv[fj];    f.psjuuu = nlpr->psiuuu[fj];
                  f.psjuuv = nlpr->psiuuv[fj];  f.psjuvv = nlpr->psiuvv[fj];
                  f.psjvvv = nlpr->psivvv[fj];
                  pos = pkn_SymMatIndex ( bs1+fn, bs1+fm );
                  _g2h_IntFunc3cd ( &f, (vector2d*)&Li[2*(nfunc_c+fn)],
                              (vector2d*)&Li[2*(nfunc_c+fm)],
                              (vector2d*)&BiLT[2*(nfunc_c+fn)],
                              (vector2d*)&BiLT[2*(nfunc_c+fm)],
                              Di[nfunc_c+fn], Di[nfunc_c+fm], &hkk[pos] );
                }
              }
            for ( fim = 0, fm = kk*rrd;  fim < rrd;  fim++, fm++ )
              if ( fm <= fn ) {
                _g1h_FuncDSuppd ( hole_k, nk, m1, fm, k,
                                  &nzc, &i10, &i11, &j10, &j11 );
                if ( i >= i10 && i <= i11 && j >= j10 && j <= j11 ) {
                  fj = dfuncvi[fm] + (i11-i10+1)*(j11-j10+1) + (i-i10)*(j11-j10+1) + j-j10;
                  f.psju = nlpr->psiu[fj];      f.psjv = nlpr->psiv[fj];
                  f.psjuu = nlpr->psiuu[fj];    f.psjuv = nlpr->psiuv[fj];
                  f.psjvv = nlpr->psivv[fj];    f.psjuuu = nlpr->psiuuu[fj];
                  f.psjuuv = nlpr->psiuuv[fj];  f.psjuvv = nlpr->psiuvv[fj];
                  f.psjvvv = nlpr->psivvv[fj];
                  pos = pkn_SymMatIndex ( bs1+fn, bs1+fm );
                  _g2h_IntFunc3cd ( &f, (vector2d*)&Li[2*(nfunc_c+fn)],
                              (vector2d*)&Li[2*(nfunc_c+fm)],
                              (vector2d*)&BiLT[2*(nfunc_c+fn)],
                              (vector2d*)&BiLT[2*(nfunc_c+fm)],
                              Di[nfunc_c+fn], Di[nfunc_c+fm], &hkk[pos] );
                }
              }
          }
        }
     }
  }
  funct /= (double)nkn2;
  pkn_MultMatrixNumd ( 1, nfacd, 0, grad, 1.0/(double)nkn2, 0, grad );
  pkn_MultMatrixNumd ( 1, hsize, 0, hessian, 1.0/(double)nkn2, 0, hessian );

        /* integration along the curves */
  njcurves = nlspr->njcurves;
  hkk = &chessian[pkn_Block3FindBlockPos(hole_k-1,bs1,bs2,hole_k-1,hole_k-1)];
  jcfs = nlspr->jcfs;
  jdfs = nlspr->jdfs;
  for ( k = 0, ki = hole_k-1, kl = 1;
        k < hole_k;
        ki = k++, kl = (k+1) % hole_k ) {
    hii = &chessian[pkn_Block3FindBlockPos(hole_k-1,bs1,bs2,k,k)];
    hjj = &chessian[pkn_Block3FindBlockPos(hole_k-1,bs1,bs2,ki,ki)];
    hki = &chessian[pkn_Block3FindBlockPos(hole_k-1,bs1,bs2,hole_k-1,k)];
    if ( k > 0 ) {
      hij = &chessian[pkn_Block3FindBlockPos(hole_k-1,bs1,bs2,k,ki)];
      hkj = &chessian[pkn_Block3FindBlockPos(hole_k-1,bs1,bs2,hole_k-1,ki)];
    }
    else
      hij = hkj = &chessian[pkn_Block3FindBlockPos(hole_k-1,bs1,bs2,hole_k-1,0)];
    for ( nc = knc = 0, knot = k*njcurves*nkn;  nc < njcurves;  nc++ )
      for ( kn = 0;  kn < nkn;  kn++, knc++, knot++ ) {
        /* evaluate the function and its derivatives */
        fi = (nfunc_a*hole_k+k)*njcurves*nkn + knc;
        cf.pu = nlpr->cpsiu[fi];     cf.pv = nlpr->cpsiv[fi];
        cf.jpuu = nlpr->cpsiuu[fi];  cf.jpuv = nlpr->cpsiuv[fi];
        cf.jpvv = nlpr->cpsivv[fi];
        /* add the linear combination of the block A functions */
        for ( i = 0; i < nfunc_a; i++ ) {
          fi = (i*hole_k+k)*njcurves*nkn + knc;
          c = coeff[nfunc_c+nfunc_d+i];
          cf.pu -= c*nlpr->cpsiu[fi];     cf.pv -= c*nlpr->cpsiv[fi];
          cf.jpuu -= c*nlpr->cpsiuu[fi];  cf.jpuv -= c*nlpr->cpsiuv[fi];
          cf.jpvv -= c*nlpr->cpsivv[fi];
        }
        /* add the linear combination of the block C functions */
          /* the Omega_k block */
        for ( i = 0, ii = k*bs1;  i < bs1;  i++, ii++ ) {
          if ( nc < 3 )
            fi = jcfs + (ii*(njcurves+1)+nc)*nkn + kn;
          else
            fi = jcfs + (ii*(njcurves+1)+nc+1)*nkn + kn;
          c = coeff[ii];
          cf.pu -= c*nlpr->cpsiu[fi];     cf.pv -= c*nlpr->cpsiv[fi];
          cf.jpuu -= c*nlpr->cpsiuu[fi];  cf.jpuv -= c*nlpr->cpsiuv[fi];
          cf.jpvv -= c*nlpr->cpsivv[fi];
        }
          /* the Omega_{k-1} block */
        if ( nc == 0 )
          for ( i = 0, ii = ki*bs1;  i < rrc;  i++, ii++ ) {
            fi = jcfs + (ii*(njcurves+1)+3)*nkn + kn;
            c = coeff[ii];
            cf.pu -= c*nlpr->cpsiu[fi];     cf.pv -= c*nlpr->cpsiv[fi];
            cf.jpuu -= c*nlpr->cpsiuu[fi];  cf.jpuv -= c*nlpr->cpsiuv[fi];
            cf.jpvv -= c*nlpr->cpsivv[fi];
          }
        /* add the linear combination of the block D functions */
          /* the Omega_k block */
        for ( i = 0, ii = k*rrd+i;  i < rrd;  i++, ii++ ) {
          if ( nc < 3 )
            fi = jdfs + (ii*(2*njcurves+1)+nc)*nkn + kn;
          else
            fi = jdfs + (ii*(2*njcurves+1)+nc+4)*nkn + kn;
          c = coeff[nfunc_c+ii];
          cf.pu -= c*nlpr->cpsiu[fi];     cf.pv -= c*nlpr->cpsiv[fi];
          cf.jpuu -= c*nlpr->cpsiuu[fi];  cf.jpuv -= c*nlpr->cpsiuv[fi];
          cf.jpvv -= c*nlpr->cpsivv[fi];
        }
          /* the Omega_{k-1} block */
        if ( nc == 0 )
          for ( i = 0, ii = ki*rrd+i;  i < rrd;  i++, ii++ ) {
            fi = jdfs + (ii*(2*njcurves+1)+3)*nkn + kn;
            c = coeff[nfunc_c+ii];
            cf.pu -= c*nlpr->cpsiu[fi];     cf.pv -= c*nlpr->cpsiv[fi];
            cf.jpuu -= c*nlpr->cpsiuu[fi];  cf.jpuv -= c*nlpr->cpsiuv[fi];
            cf.jpvv -= c*nlpr->cpsivv[fi];
          }
          /* the Omega_{k+1} block */
        for ( i = 0, ii = kl*rrd+i;  i < rrd;  i++, ii++ ) {
          if ( nc < 3 )
            fi = jdfs + (ii*(2*njcurves+1)+nc+4)*nkn + kn;
          else
            fi = jdfs + (ii*(2*njcurves+1)+2*nk+nc+4)*nkn + kn;
          c = coeff[nfunc_c+ii];
          cf.pu -= c*nlpr->cpsiu[fi];     cf.pv -= c*nlpr->cpsiv[fi];
          cf.jpuu -= c*nlpr->cpsiuu[fi];  cf.jpuv -= c*nlpr->cpsiuv[fi];
          cf.jpvv -= c*nlpr->cpsivv[fi];
        }

        cf.tang = nlpr->ctang[knot];
        /* integrate the functional value */
        _g1hq2_IntFunc1cd ( &cf, &cfunct );

        /* integrate the gradient */
          /* block A functions */
        for ( i = 0; i < nfunc_a; i++ ) {
          fi = (i*hole_k+k)*njcurves*nkn + knc;
          ii = nfunc_c+nfunc_d+i;
          cf.psiu = nlpr->cpsiu[fi];     cf.psiv = nlpr->cpsiv[fi];
          cf.jpsiuu = nlpr->cpsiuu[fi];  cf.jpsiuv = nlpr->cpsiuv[fi];
          cf.jpsivv = nlpr->cpsivv[fi];
          _g1hq2_IntFunc2cd ( &cf, &Di[ii], &Bi[ii], &cgrad[ii] );
        }
          /* block C functions */
            /* the Omega_k block */
        for ( i = 0, ii = k*bs1;  i < bs1;  i++, ii++ ) {
          if ( nc < 3 )
            fi = jcfs + (ii*(njcurves+1)+nc)*nkn + kn;
          else
            fi = jcfs + (ii*(njcurves+1)+nc+1)*nkn + kn;
          cf.psiu = nlpr->cpsiu[fi];     cf.psiv = nlpr->cpsiv[fi];
          cf.jpsiuu = nlpr->cpsiuu[fi];  cf.jpsiuv = nlpr->cpsiuv[fi];
          cf.jpsivv = nlpr->cpsivv[fi];
          _g1hq2_IntFunc2cd ( &cf, &Di[ii], &Bi[ii], &cgrad[ii] );
        }
            /* the Omega_{k+1} block */
        if ( nc == 0 )
          for ( i = 0, ii = ki*bs1;  i < rrc;  i++, ii++ ) {
            fi = jcfs + (ii*(njcurves+1)+3)*nkn + kn;
            cf.psiu = nlpr->cpsiu[fi];     cf.psiv = nlpr->cpsiv[fi];
            cf.jpsiuu = nlpr->cpsiuu[fi];  cf.jpsiuv = nlpr->cpsiuv[fi];
            cf.jpsivv = nlpr->cpsivv[fi];
            _g1hq2_IntFunc2cd ( &cf, &Di[ii], &Bi[ii], &cgrad[ii] );
          }
          /* block D functions */
            /* the Omega_k block */
        for ( i = 0, ii = k*rrd+i;  i < rrd;  i++, ii++ ) {
          if ( nc < 3 )
            fi = jdfs + (ii*(2*njcurves+1)+nc)*nkn + kn;
          else
            fi = jdfs + (ii*(2*njcurves+1)+nc+4)*nkn + kn;
          cf.psiu = nlpr->cpsiu[fi];     cf.psiv = nlpr->cpsiv[fi];
          cf.jpsiuu = nlpr->cpsiuu[fi];  cf.jpsiuv = nlpr->cpsiuv[fi];
          cf.jpsivv = nlpr->cpsivv[fi];
          _g1hq2_IntFunc2cd ( &cf, &Di[nfunc_c+ii], &Bi[nfunc_c+ii],
                              &cgrad[nfunc_c+ii] );
        }
            /* the Omega_{k-1} block */
        if ( nc == 0 )
          for ( i = 0, ii = ki*rrd+i;  i < rrd;  i++, ii++ ) {
            fi = jdfs + (ii*(2*njcurves+1)+3)*nkn + kn;
            cf.psiu = nlpr->cpsiu[fi];     cf.psiv = nlpr->cpsiv[fi];
            cf.jpsiuu = nlpr->cpsiuu[fi];  cf.jpsiuv = nlpr->cpsiuv[fi];
            cf.jpsivv = nlpr->cpsivv[fi];
            _g1hq2_IntFunc2cd ( &cf, &Di[nfunc_c+ii], &Bi[nfunc_c+ii],
                                &cgrad[nfunc_c+ii] );
          }
            /* the Omega_{k+1} block */
        for ( i = 0, ii = kl*rrd+i;  i < rrd;  i++, ii++ ) {
          if ( nc < 3 )
            fi = jdfs + (ii*(2*njcurves+1)+nc+4)*nkn + kn;
          else
            fi = jdfs + (ii*(2*njcurves+1)+2*nk+nc+4)*nkn + kn;
          cf.psiu = nlpr->cpsiu[fi];     cf.psiv = nlpr->cpsiv[fi];
          cf.jpsiuu = nlpr->cpsiuu[fi];  cf.jpsiuv = nlpr->cpsiuv[fi];
          cf.jpsivv = nlpr->cpsivv[fi];
          _g1hq2_IntFunc2cd ( &cf, &Di[nfunc_c+ii], &Bi[nfunc_c+ii],
                              &cgrad[nfunc_c+ii] );
        }

        /* integrate the Hessian */
        for ( i = 0; i < nfunc_a; i++ ) {
          fi = i*njcurves*nkn*hole_k + k*njcurves*nkn + knc;
          ii = nfunc_c+nfunc_d+i;
          cf.psiu = nlpr->cpsiu[fi];     cf.psiv = nlpr->cpsiv[fi];
          cf.jpsiuu = nlpr->cpsiuu[fi];  cf.jpsiuv = nlpr->cpsiuv[fi];
          cf.jpsivv = nlpr->cpsivv[fi];
          /* A x A */
          for ( j = 0; j <= i; j++ ) {
            fj = j*njcurves*nkn*hole_k + k*njcurves*nkn + knc;
            jj = nfunc_c+nfunc_d+j;
            cf.psju = nlpr->cpsiu[fj];     cf.psjv = nlpr->cpsiv[fj];
            cf.jpsjuu = nlpr->cpsiuu[fj];  cf.jpsjuv = nlpr->cpsiuv[fj];
            cf.jpsjvv = nlpr->cpsivv[fj];
            _g1hq2_IntFunc3cd ( &cf, Di[ii], Bi[ii], Di[jj], Bi[jj],
                                &hkk[pkn_SymMatIndex(bs1+nfunc_d+i,bs1+nfunc_d+j)] );
          }
          /* A x C */
            /* C - block Omega_k */
          for ( j = 0, jj = k*bs1;  j < bs1;  j++, jj++ ) {
            if ( nc < 3 )
              fj = jcfs + (jj*(njcurves+1)+nc)*nkn + kn;
            else
              fj = jcfs + (jj*(njcurves+1)+nc+1)*nkn + kn;
            cf.psju = nlpr->cpsiu[fj];     cf.psjv = nlpr->cpsiv[fj];
            cf.jpsjuu = nlpr->cpsiuu[fj];  cf.jpsjuv = nlpr->cpsiuv[fj];
            cf.jpsjvv = nlpr->cpsivv[fj];
            if ( k < hole_k-1 )
              hh = &hki[(bs1+nfunc_d+i)*bs1+j];
            else
              hh = &hkk[pkn_SymMatIndex(bs1+nfunc_d+i,j)];
            _g1hq2_IntFunc3cd ( &cf, Di[ii], Bi[ii], Di[jj], Bi[jj], hh );
          }
            /* C - block Omega_{k-1} */
          if ( nc == 0 )
            for ( j = 0, jj = ki*bs1;  j < rrc;  j++, jj++ ) {
              fj = jcfs + (jj*(njcurves+1)+3)*nkn + kn;
              cf.psju = nlpr->cpsiu[fj];     cf.psjv = nlpr->cpsiv[fj];
              cf.jpsjuu = nlpr->cpsiuu[fj];  cf.jpsjuv = nlpr->cpsiuv[fj];
              cf.jpsjvv = nlpr->cpsivv[fj];
              if ( ki < hole_k-1 )
                hh = &hkj[(bs1+nfunc_d+i)*bs1+j];
              else
                hh = &hkk[pkn_SymMatIndex(bs1+nfunc_d+i,j)];
              _g1hq2_IntFunc3cd ( &cf, Di[ii], Bi[ii], Di[jj], Bi[jj], hh );
            }
          /* A x D */
            /* D - block Omega_k */
          for ( j = 0, jj = k*rrd+j;  j < rrd;  j++, jj++ ) {
            if ( nc < 3 )
              fj = jdfs + (jj*(2*njcurves+1)+nc)*nkn + kn;
            else
              fj = jdfs + (jj*(2*njcurves+1)+nc+4)*nkn + kn;
            cf.psju = nlpr->cpsiu[fj];     cf.psjv = nlpr->cpsiv[fj];
            cf.jpsjuu = nlpr->cpsiuu[fj];  cf.jpsjuv = nlpr->cpsiuv[fj];
            cf.jpsjvv = nlpr->cpsivv[fj];
            hh = &hkk[pkn_SymMatIndex(bs1+nfunc_d+i,bs1+jj)];
            _g1hq2_IntFunc3cd ( &cf, Di[ii], Bi[ii],
                                Di[nfunc_c+jj], Bi[nfunc_c+jj], hh );
          }
            /* D - block Omega_{k-1} */
          if ( nc == 0 )
            for ( j = 0, jj = ki*rrd;  j < 2;  j++, jj++ ) {
              fj = jdfs + (jj*(2*njcurves+1)+3)*nkn + kn;
              cf.psju = nlpr->cpsiu[fj];     cf.psjv = nlpr->cpsiv[fj];
              cf.jpsjuu = nlpr->cpsiuu[fj];  cf.jpsjuv = nlpr->cpsiuv[fj];
              cf.jpsjvv = nlpr->cpsivv[fj];
              hh = &hkk[pkn_SymMatIndex(bs1+nfunc_d+i,bs1+jj)];
              _g1hq2_IntFunc3cd ( &cf, Di[ii], Bi[ii],
                                  Di[nfunc_c+jj], Bi[nfunc_c+jj], hh );
            }
            /* D - block Omega_{k+1} */
          for ( j = 0, jj = kl*rrd+j;  j < rrd;  j++, jj++ ) {
            if ( nc < 3 )
              fj = jdfs + (jj*(2*njcurves+1)+nc+4)*nkn + kn;
            else
              fj = jdfs + (jj*(2*njcurves+1)+2*nk+nc+4)*nkn + kn;
            cf.psju = nlpr->cpsiu[fj];     cf.psjv = nlpr->cpsiv[fj];
            cf.jpsjuu = nlpr->cpsiuu[fj];  cf.jpsjuv = nlpr->cpsiuv[fj];
            cf.jpsjvv = nlpr->cpsivv[fj];
            hh = &hkk[pkn_SymMatIndex(bs1+nfunc_d+i,bs1+jj)];
            _g1hq2_IntFunc3cd ( &cf, Di[ii], Bi[ii],
                                Di[nfunc_c+jj], Bi[nfunc_c+jj], hh );
          }
        }
        /* C - block Omega_k */
        for ( i = 0, ii = k*bs1;  i < bs1;  i++, ii++ ) {
          if ( nc < 3 )
            fi = jcfs + (ii*(njcurves+1)+nc)*nkn + kn;
          else
            fi = jcfs + (ii*(njcurves+1)+nc+1)*nkn + kn;
          cf.psiu = nlpr->cpsiu[fi];     cf.psiv = nlpr->cpsiv[fi];
          cf.jpsiuu = nlpr->cpsiuu[fi];  cf.jpsiuv = nlpr->cpsiuv[fi];
          cf.jpsivv = nlpr->cpsivv[fi];
          /* C x C */
          for ( j = 0, jj = k*bs1;  j <= i;  j++, jj++ ) {
            if ( nc < 3 )
              fj = jcfs + (jj*(njcurves+1)+nc)*nkn + kn;
            else
              fj = jcfs + (jj*(njcurves+1)+nc+1)*nkn + kn;
            cf.psju = nlpr->cpsiu[fj];     cf.psjv = nlpr->cpsiv[fj];
            cf.jpsjuu = nlpr->cpsiuu[fj];  cf.jpsjuv = nlpr->cpsiuv[fj];
            cf.jpsjvv = nlpr->cpsivv[fj];
            _g1hq2_IntFunc3cd ( &cf, Di[ii], Bi[ii], Di[jj], Bi[jj],
                                &hii[pkn_SymMatIndex(i,j)] );
          }
          if ( nc == 0 )
            for ( j = 0, jj = ki*bs1;  j < rrc;  j++, jj++ ) {
              fj = jcfs + (jj*(njcurves+1)+3)*nkn + kn;
              cf.psju = nlpr->cpsiu[fj];     cf.psjv = nlpr->cpsiv[fj];
              cf.jpsjuu = nlpr->cpsiuu[fj];  cf.jpsjuv = nlpr->cpsiuv[fj];
              cf.jpsjvv = nlpr->cpsivv[fj];
              if ( k == 0 )
                hh = &hkj[j*bs1+i];
              else
                hh = &hij[i*bs1+j];
              _g1hq2_IntFunc3cd ( &cf, Di[ii], Bi[ii], Di[jj], Bi[jj], hh );
            }
          /* C x D */
            /* D - block Omega_k */
          for ( j = 0, jj = k*rrd+j;  j < rrd;  j++, jj++ ) {
            if ( nc < 3 )
              fj = jdfs + (jj*(2*njcurves+1)+nc)*nkn + kn;
            else
              fj = jdfs + (jj*(2*njcurves+1)+nc+4)*nkn + kn;
            cf.psju = nlpr->cpsiu[fj];     cf.psjv = nlpr->cpsiv[fj];
            cf.jpsjuu = nlpr->cpsiuu[fj];  cf.jpsjuv = nlpr->cpsiuv[fj];
            cf.jpsjvv = nlpr->cpsivv[fj];
            if ( k < hole_k-1 )
              hh = &hki[(bs1+jj)*bs1+i];
            else
              hh = &hkk[pkn_SymMatIndex(bs1+jj,i)];
            _g1hq2_IntFunc3cd ( &cf, Di[ii], Bi[ii],
                                Di[nfunc_c+jj], Bi[nfunc_c+jj], hh );
          }
            /* D - block Omega_{k-1} */
          if ( nc == 0 )
            for ( j = 0, jj = ki*rrd;  j < 2;  j++, jj++ ) {
              fj = jdfs + (jj*(2*njcurves+1)+3)*nkn + kn;
              cf.psju = nlpr->cpsiu[fj];     cf.psjv = nlpr->cpsiv[fj];
              cf.jpsjuu = nlpr->cpsiuu[fj];  cf.jpsjuv = nlpr->cpsiuv[fj];
              cf.jpsjvv = nlpr->cpsivv[fj];
              if ( k < hole_k-1 )
                hh = &hki[(bs1+jj)*bs1+i];
              else
                hh = &hkk[pkn_SymMatIndex(bs1+jj,i)];
              _g1hq2_IntFunc3cd ( &cf, Di[ii], Bi[ii],
                                  Di[nfunc_c+jj], Bi[nfunc_c+jj], hh );
            }
            /* D - block Omega_{k+1} */
          for ( j = 0, jj = kl*rrd+j;  j < rrd;  j++, jj++ ) {
            if ( nc < 3 )
              fj = jdfs + (jj*(2*njcurves+1)+nc+4)*nkn + kn;
            else
              fj = jdfs + (jj*(2*njcurves+1)+2*nk+nc+4)*nkn + kn;
            cf.psju = nlpr->cpsiu[fj];     cf.psjv = nlpr->cpsiv[fj];
            cf.jpsjuu = nlpr->cpsiuu[fj];  cf.jpsjuv = nlpr->cpsiuv[fj];
            cf.jpsjvv = nlpr->cpsivv[fj];
            if ( k < hole_k-1 )
              hh = &hki[(bs1+jj)*bs1+i];
            else
              hh = &hkk[pkn_SymMatIndex(bs1+jj,i)];
            _g1hq2_IntFunc3cd ( &cf, Di[ii], Bi[ii],
                                Di[nfunc_c+jj], Bi[nfunc_c+jj], hh );
          }
        }
        /* C - block Omega_{k-1} */
        if ( nc == 0 )
          for ( i = 0, ii = ki*bs1;  i < rrc;  i++, ii++ ) {
            fi = jcfs+(ii*(njcurves+1)+3)*nkn + kn;
            cf.psiu = nlpr->cpsiu[fi];     cf.psiv = nlpr->cpsiv[fi];
            cf.jpsiuu = nlpr->cpsiuu[fi];  cf.jpsiuv = nlpr->cpsiuv[fi];
            cf.jpsivv = nlpr->cpsivv[fi];
          /* C x C, C - block Omega_{k-1} */  
            for ( j = 0, jj = ki*bs1;  j <= i;  j++, jj++ ) {
              fj = jcfs+(jj*(njcurves+1)+3)*nkn + kn;
              cf.psju = nlpr->cpsiu[fj];     cf.psjv = nlpr->cpsiv[fj];
              cf.jpsjuu = nlpr->cpsiuu[fj];  cf.jpsjuv = nlpr->cpsiuv[fj];
              cf.jpsjvv = nlpr->cpsivv[fj];
              _g1hq2_IntFunc3cd ( &cf, Di[ii], Bi[ii], Di[jj], Bi[jj],
                                  &hjj[pkn_SymMatIndex(i,j)] );
            }
          /* C x D, D - block Omega_k */
            for ( j = 0, jj = k*rrd+j;  j < rrd;  j++, jj++ ) {
              fj = jdfs + jj*(2*njcurves+1)*nkn + kn;
              cf.psju = nlpr->cpsiu[fj];     cf.psjv = nlpr->cpsiv[fj];
              cf.jpsjuu = nlpr->cpsiuu[fj];  cf.jpsjuv = nlpr->cpsiuv[fj];
              cf.jpsjvv = nlpr->cpsivv[fj];
              if ( ki < hole_k-1 )
                hh = &hkj[(bs1+jj)*bs1+i];
              else
                hh = &hkk[pkn_SymMatIndex(bs1+jj,i)];
              _g1hq2_IntFunc3cd ( &cf, Di[ii], Bi[ii],
                                  Di[nfunc_c+jj], Bi[nfunc_c+jj], hh );
            }
          /* C x D, D - block Omega_{k-1} */
            for ( j = 0, jj = ki*rrd+j;  j < rrd;  j++, jj++ ) {
              fj = jdfs + (jj*(2*njcurves+1)+3)*nkn + kn;
              cf.psju = nlpr->cpsiu[fj];     cf.psjv = nlpr->cpsiv[fj];
              cf.jpsjuu = nlpr->cpsiuu[fj];  cf.jpsjuv = nlpr->cpsiuv[fj];
              cf.jpsjvv = nlpr->cpsivv[fj];
              if ( ki < hole_k-1 )
                hh = &hkj[(bs1+jj)*bs1+i];
              else
                hh = &hkk[pkn_SymMatIndex(bs1+jj,i)];
              _g1hq2_IntFunc3cd ( &cf, Di[ii], Bi[ii],
                                  Di[nfunc_c+jj], Bi[nfunc_c+jj], hh );
            }
          /* C x D, D - block Omega_{k+1} */
            for ( j = 0, jj = kl*rrd+j;  j < 2;  j++, jj++ ) {
              fj = jdfs + (jj*(2*njcurves+1)+4)*nkn + kn;
              cf.psju = nlpr->cpsiu[fj];     cf.psjv = nlpr->cpsiv[fj];
              cf.jpsjuu = nlpr->cpsiuu[fj];  cf.jpsjuv = nlpr->cpsiuv[fj];
              cf.jpsjvv = nlpr->cpsivv[fj];
              if ( ki < hole_k-1 )
                hh = &hkj[(bs1+jj)*bs1+i];
              else
                hh = &hkk[pkn_SymMatIndex(bs1+jj,i)];
              _g1hq2_IntFunc3cd ( &cf, Di[ii], Bi[ii],
                                  Di[nfunc_c+jj], Bi[nfunc_c+jj], hh );
            }
          }
          /* D x D */
            /* D - block Omega_k */
        for ( i = 0, ii = k*rrd+i;  i < rrd;  i++, ii++ ) {
          if ( nc < 3 )
            fi = jdfs + (ii*(2*njcurves+1)+nc)*nkn + kn;
          else
            fi = jdfs + (ii*(2*njcurves+1)+nc+4)*nkn + kn;
          cf.psiu = nlpr->cpsiu[fi];     cf.psiv = nlpr->cpsiv[fi];
          cf.jpsiuu = nlpr->cpsiuu[fi];  cf.jpsiuv = nlpr->cpsiuv[fi];
          cf.jpsivv = nlpr->cpsivv[fi];
              /* D - block Omega_k */
          for ( j = 0, jj = k*rrd+j;  j <= i;  j++, jj++ ) {
            if ( nc < 3 )
              fj = jdfs + (jj*(2*njcurves+1)+nc)*nkn + kn;
            else
              fj = jdfs + (jj*(2*njcurves+1)+nc+4)*nkn + kn;
            cf.psju = nlpr->cpsiu[fj];     cf.psjv = nlpr->cpsiv[fj];
            cf.jpsjuu = nlpr->cpsiuu[fj];  cf.jpsjuv = nlpr->cpsiuv[fj];
            cf.jpsjvv = nlpr->cpsivv[fj];
            hh = &hkk[pkn_SymMatIndex(bs1+ii,bs1+jj)];
            _g1hq2_IntFunc3cd ( &cf, Di[nfunc_c+ii], Bi[nfunc_c+ii],
                                Di[nfunc_c+jj], Bi[nfunc_c+jj], hh );
          }
              /* D - block Omega_{k-1} */
          if ( nc == 0 )
            for ( j = 0, jj = ki*rrd;  j < 2;  j++, jj++ ) {
              fj = jdfs + (jj*(2*njcurves+1)+3)*nkn + kn;
              cf.psju = nlpr->cpsiu[fj];     cf.psjv = nlpr->cpsiv[fj];
              cf.jpsjuu = nlpr->cpsiuu[fj];  cf.jpsjuv = nlpr->cpsiuv[fj];
              cf.jpsjvv = nlpr->cpsivv[fj];
              hh = &hkk[pkn_SymMatIndex(bs1+ii,bs1+jj)];
              _g1hq2_IntFunc3cd ( &cf, Di[nfunc_c+ii], Bi[nfunc_c+ii],
                                  Di[nfunc_c+jj], Bi[nfunc_c+jj], hh );
            }
              /* D - block Omega_{k+1} */
          for ( j = 0, jj = kl*rrd+j;  j < rrd;  j++, jj++ ) {
            if ( nc < 3 )
              fj = jdfs + (jj*(2*njcurves+1)+nc+4)*nkn + kn;
            else
              fj = jdfs + (jj*(2*njcurves+1)+2*nk+nc+4)*nkn + kn;
            cf.psju = nlpr->cpsiu[fj];     cf.psjv = nlpr->cpsiv[fj];
            cf.jpsjuu = nlpr->cpsiuu[fj];  cf.jpsjuv = nlpr->cpsiuv[fj];
            cf.jpsjvv = nlpr->cpsivv[fj];
            hh = &hkk[pkn_SymMatIndex(bs1+ii,bs1+jj)];
            _g1hq2_IntFunc3cd ( &cf, Di[nfunc_c+ii], Bi[nfunc_c+ii],
                                Di[nfunc_c+jj], Bi[nfunc_c+jj], hh );
          }
        }
            /* D - block Omega_{k-1} */
        if ( nc == 0 )
          for ( i = 0, ii = ki*rrd;  i < 2;  i++, ii++ ) {
            fi = jdfs + (ii*(2*njcurves+1)+3)*nkn + kn;
            cf.psiu = nlpr->cpsiu[fi];     cf.psiv = nlpr->cpsiv[fi];
            cf.jpsiuu = nlpr->cpsiuu[fi];  cf.jpsiuv = nlpr->cpsiuv[fi];
            cf.jpsivv = nlpr->cpsivv[fi];
              /* D - block Omega_{k-1} */
            for ( j = 0, jj = ki*rrd;  j <= i;  j++, jj++ ) {
              fj = jdfs + (jj*(2*njcurves+1)+3)*nkn + kn;
              cf.psju = nlpr->cpsiu[fj];     cf.psjv = nlpr->cpsiv[fj];
              cf.jpsjuu = nlpr->cpsiuu[fj];  cf.jpsjuv = nlpr->cpsiuv[fj];
              cf.jpsjvv = nlpr->cpsivv[fj];
              hh = &hkk[pkn_SymMatIndex(bs1+ii,bs1+jj)];
              _g1hq2_IntFunc3cd ( &cf, Di[nfunc_c+ii], Bi[nfunc_c+ii],
                                  Di[nfunc_c+jj], Bi[nfunc_c+jj], hh );
            }
              /* D - block Omega_{k+1} */
            for ( j = 0, jj = kl*rrd+j;  j < rrd;  j++, jj++ ) {
              fj = jdfs + (jj*(2*njcurves+1)+nc+4)*nkn + kn;
              cf.psju = nlpr->cpsiu[fj];     cf.psjv = nlpr->cpsiv[fj];
              cf.jpsjuu = nlpr->cpsiuu[fj];  cf.jpsjuv = nlpr->cpsiuv[fj];
              cf.jpsjvv = nlpr->cpsivv[fj];
              hh = &hkk[pkn_SymMatIndex(bs1+ii,bs1+jj)];
              _g1hq2_IntFunc3cd ( &cf, Di[nfunc_c+ii], Bi[nfunc_c+ii],
                                  Di[nfunc_c+jj], Bi[nfunc_c+jj], hh );
            }
          }
            /* D - block Omega_{k+1} */
        for ( i = 0, ii = kl*rrd+i;  i < rrd;  i++, ii++ ) {
          if ( nc < 3 )
            fi = jdfs + (ii*(2*njcurves+1)+nc+4)*nkn + kn;
          else
            fi = jdfs + (ii*(2*njcurves+1)+2*nk+nc+4)*nkn + kn;
          cf.psiu = nlpr->cpsiu[fi];     cf.psiv = nlpr->cpsiv[fi];
          cf.jpsiuu = nlpr->cpsiuu[fi];  cf.jpsiuv = nlpr->cpsiuv[fi];
          cf.jpsivv = nlpr->cpsivv[fi];
              /* D - block Omega_{k+1} */
          for ( j = 0, jj = kl*rrd+j;  j <= i;  j++, jj++ ) {
            if ( nc < 3 )
              fj = jdfs + (jj*(2*njcurves+1)+nc+4)*nkn + kn;
            else
              fj = jdfs + (jj*(2*njcurves+1)+2*nk+nc+4)*nkn + kn;
            cf.psju = nlpr->cpsiu[fj];     cf.psjv = nlpr->cpsiv[fj];
            cf.jpsjuu = nlpr->cpsiuu[fj];  cf.jpsjuv = nlpr->cpsiuv[fj];
            cf.jpsjvv = nlpr->cpsivv[fj];
            hh = &hkk[pkn_SymMatIndex(bs1+ii,bs1+jj)];
            _g1hq2_IntFunc3cd ( &cf, Di[nfunc_c+ii], Bi[nfunc_c+ii],
                                Di[nfunc_c+jj], Bi[nfunc_c+jj], hh );
          }
        }
      }
    }
  cfunct /= (double)nkn;
  *func = funct + C*cfunct;
  pkn_AddMatrixMd ( 1, nfacd, 0, grad, 0, cgrad, C/(double)nkn, 0, grad );
  pkn_AddMatrixMd ( 1, hsize, 0, hessian, 0, chessian, C/(double)nkn,
                    0, hessian );

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*_g1hq2_ComputeSplNLFuncGradHessiand*/

static boolean g1hq2_SplNLNewtond ( GHoleDomaind *domain,
                                    G1HNLSPrivated *nlsprivate )
{
#define EPSF 2.0e-4
  void *sp;
  G1HolePrivateRecd  *privateG1;
  G1HoleSPrivateRecd *sprivate;
  G1HNLPrivated      *nlpr;
  int     itn, jtn, ktn, hole_k, nfunc_a, nfunc_c, nfunc_d, nfunc;
  int     bs1, bs2, asize;
  int     nk, m2;
  double  func, *grad, *coeff, *dcoeff, *hessian, *chess;
  double  func0, func1, aux, gn, gn0 = 0.0, gn1, dco, dyn;
  double  C; 
  boolean positive;

  sp = pkv_GetScratchMemTop ();
  hole_k = domain->hole_k;
  privateG1 = domain->privateG1;
  sprivate = domain->SprivateG1;
  nlpr = &nlsprivate->nlpr;
  nfunc_a = privateG1->nfunc_a;
  nfunc_c = sprivate->nsfunc_c;
  nfunc_d = sprivate->nsfunc_d;
  nfunc = nfunc_a+nfunc_c+nfunc_d;
  nk = sprivate->nk;
  m2 = sprivate->m2;
  C = (double)(5+nk*m2)*sprivate->C1s/nlpr->ddiam;

  coeff = pkv_GetScratchMemd ( 3*nfunc );
  bs1 = nfunc_c/hole_k;
  bs2 = bs1+nfunc_d+nfunc_a;
  asize = pkn_Block3ArraySize ( hole_k-1, bs1, bs2 );
  hessian = pkv_GetScratchMemd ( 2*asize );
  if ( !coeff || !hessian ) {
    domain->error_code = G1H_ERROR_NO_SCRATCH_MEMORY;
    goto failure;
  }
  dcoeff = &coeff[nfunc];
  grad   = &dcoeff[nfunc];
  chess  = &hessian[asize];

        /* setup the initial point */
  pkv_Selectd ( nfunc, 1, 3, 1, &nlpr->acoeff[0].z, coeff );

        /* Newton method iterations */
  for ( itn = ktn = 0; ; itn++ ) {
    if ( !_g1hq2_ComputeSplNLFuncGradHessiand ( domain, nlsprivate, coeff, C,
                                                &func, grad, hessian, chess ) )
      goto failure;
    gn = (double)sqrt ( pkn_ScalarProductd ( nfunc, grad, grad ) );
    if ( itn == 0 ) {
printf ( "func = %f, gn0 = %f\n", func, gn );
      gn0 = gn;
    }
    memcpy ( chess, hessian, asize*sizeof(double) );
    if ( (positive = pkn_Block3CholeskyDecompMd ( hole_k-1, bs1, bs2,
                                                  hessian ) ) ) {
      pkn_Block3LowerTrMSolved ( hole_k-1, bs1, bs2, hessian, 1, 1, grad );
      pkn_Block3UpperTrMSolved ( hole_k-1, bs1, bs2, hessian, 1, 1, grad );
    }
    else {
printf ( "! " );

      pkn_Block3SymMatrixMultd ( hole_k-1, bs1, bs2, chess,
                                 1, 1, grad, 1, dcoeff );
      aux = (double)pkn_ScalarProductd ( nfunc, grad, dcoeff );
      if ( gn < 0.0 || aux < EPSF*gn ) {
        domain->error_code = G1H_ERROR_NL_MINIMIZATION;
        goto failure;
      }
      pkn_MultMatrixNumd ( 1, nfunc, 0, grad, gn/aux, 0, grad );
    }
    dco = (double)sqrt ( pkn_ScalarProductd ( nfunc, coeff, coeff ) );
    dyn = (double)sqrt ( pkn_ScalarProductd ( nfunc, grad, grad ) );

    for ( aux = 1.0; aux > EPSF; aux *= 0.5 ) {
      pkn_AddMatrixMd ( 1, nfunc, 0, coeff, 0, grad, aux, 0, dcoeff );
      func1 = _g1hq2_ComputeSplNLFuncd ( domain, nlsprivate, dcoeff, C );
      if ( func1 < func )
        break;
    }

    memcpy ( coeff, dcoeff, nfunc*sizeof(double) );
    func = func1;
    ktn ++;

    if ( positive && aux > 0.1 ) {
      /* with the positive-definite Hessian we try to mke some */
      /* extra iterations */

      for ( jtn = 0; jtn < 10; jtn ++ ) {

printf ( "+" );

        if ( !_g1hq2_ComputeSplNLFuncGradd ( domain, nlsprivate, coeff,
                                             C, &func0, grad ) )
          goto failure;
        gn1 = (double)sqrt ( pkn_ScalarProductd ( nfunc, grad, grad ) );
        pkn_Block3LowerTrMSolved ( hole_k-1, bs1, bs2, hessian, 1, 1, grad );
        pkn_Block3UpperTrMSolved ( hole_k-1, bs1, bs2, hessian, 1, 1, grad );
        pkn_AddMatrixd ( 1, nfunc, 0, coeff, 0, grad, 0, dcoeff );
        func1 = _g1hq2_ComputeSplNLFuncd ( domain, nlsprivate, dcoeff, C );
        if ( func1 >= func0 || gn1 >= 0.5*gn )
          break;

printf ( "    func = %f, gn = %f\n", func1, gn1 );

        memcpy ( coeff, dcoeff, nfunc*sizeof(double) );
        func = func1;
        ktn ++;
        gn = gn1;
        aux = 1.0;
        dco = (double)sqrt ( pkn_ScalarProductd ( nfunc, coeff, coeff ) );
        dyn = (double)sqrt ( pkn_ScalarProductd ( nfunc, grad, grad ) );
      }
    }

    if ( _g1h_StopItd ( itn, gn0, gn, dco, dyn, aux ) )
      break;
  }

printf ( "itn = %d, ktn = %d, func = %f, gn = %f\n", itn+1, ktn, func1, gn );

  pkv_Selectd ( nfunc, 1, 1, 3, coeff, &nlpr->acoeff[0].z );

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
#undef EPSF
} /*g1hq2_SplNLNewtond*/

/* ///////////////////////////////////////////////////////////////////////// */
boolean g1h_Q2NLSplFillHoled ( GHoleDomaind *domain,
                    const point3d *hole_cp,
                    double *acoeff, void *usrptr,
                    void (*outpatch) ( int n, int lknu, const double *knu,
                                       int m, int lknv, const double *knv,
                                       const point3d *cp, void *usrptr ) )
{
  typedef void outscf ( int n, int lknu, const double *knu,
                        int m, int lknv, const double *knv,
                        const double *cp, void *usrptr );

  void *sp;
  G1HNLSPrivated     *nlsprivate;
  G1HNLPrivated      *nlprivate;
  G1HolePrivateRecd  *privateG1;
  G1HoleSPrivateRecd *sprivate;
  int   hole_k, nfunc_a, nfunc_c, nfunc_d;
  double *fc00;

  sp = pkv_GetScratchMemTop ();
  privateG1 = domain->privateG1;
  sprivate = domain->SprivateG1;
  if ( !sprivate ) {
    domain->error_code = G1H_ERROR_SPLINE_BASIS_NOT_READY;
    goto failure;
  }
  if ( !sprivate->Q2SAMat )
    if ( !g1h_Q2SplComputeFormMatrixd ( domain ) )
      goto failure;
  if ( !sprivate->Q2SLMat )
    if ( !g1h_Q2SplDecomposeMatrixd ( domain ) )
      goto failure;
  nlsprivate = _g1hq2_AllocSplNLPrd ( domain );
  if ( !nlsprivate )
    goto failure;
  _g1h_nlprivd = nlprivate = &nlsprivate->nlpr;
  hole_k = domain->hole_k;
  nfunc_a = privateG1->nfunc_a;
  nfunc_c = sprivate->nsfunc_c;
  nfunc_d = sprivate->nsfunc_d;
  if ( !_g1hq2_InitSplNLprd ( domain, nlsprivate, 0 ) )
    goto failure;

  if ( !_g1h_ComputeNLNormald ( domain, nlprivate, hole_cp ) )
    goto failure;
  g1h_ReflectVectorsd ( 12*hole_k+1, hole_cp, nlprivate->rhole_cp );
  nlprivate->auxc = 0;
  if ( !g1h_Q2SplFillHoled ( domain, 3, (double*)nlprivate->rhole_cp,
                             (double*)nlprivate->acoeff, NULL,
                             g1h_splnloutpatchd ) )
    goto failure;

  if ( !_g1hq2_TabSplNLBasisFunctionsd ( domain, nlsprivate ) )
    goto failure;
  if ( !_g1hq2_TabSplNLBasisFJumpsd ( domain, nlsprivate ) )
    goto failure;

  if ( !g1hq2_SplNLNewtond ( domain, nlsprivate ) )
    goto failure;

  g1h_ReflectVectorsd ( nfunc_a+nfunc_c+nfunc_d, nlprivate->acoeff,
                        nlprivate->acoeff );
  if ( acoeff )
    memcpy ( acoeff, nlprivate->acoeff,
             (nfunc_a+nfunc_c+nfunc_d)*sizeof(vector3d) );

  fc00 = pkv_GetScratchMemd ( (G1_CROSSDEGSUM+4)*2*hole_k*3 );
  if ( !fc00 ) {
    domain->error_code = G1H_ERROR_NO_SCRATCH_MEMORY;
    goto failure;
  }
  if ( !_g1h_SetSplRightSided ( domain, 3, (double*)hole_cp,
                                sprivate->Q2SBMat, fc00, NULL ) )
    goto failure;

  if ( !_g1h_OutputSplPatchesd ( domain, 3, (double*)nlprivate->acoeff, fc00,
                                 usrptr, (outscf*)outpatch ) )
    goto failure;

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*g1h_Q2NLSplFillHoled*/

/* ///////////////////////////////////////////////////////////////////////// */
static boolean g1hq2_SplNLConstrNewtond ( GHoleDomaind *domain,
                                          G1HNLSPrivated *nlsprivate,
                                          int nconstr, double *SCmat )
{
#define EPSF 2.0e-4
  void *sp;
  G1HolePrivateRecd  *privateG1;
  G1HoleSPrivateRecd *sprivate;
  G1HNLPrivated      *nlpr;
  int    itn, jtn, ktn, hole_k, nfunc_a, nfunc_c, nfunc_d, nfacd, nfunc, nk, m2;
  int    bs1, bs2, bs3, asize, esize, diagblsize, subdiagblsize, sideblsize,
         esideblsize;
  double func, *grad, *coeff, *hessian, *E22, *cE22,
         *cT, *aa, *D1, *y, *y1, *M, *f, *hkk, *hij, *hki, *E22kk, *E22ij, *E22ki;
  double C, func0, func1, aux, gn, gn0 = 0.0, gn1, dco, dyn, dyn1;
  int    i, j;
  boolean positive;

  sp = pkv_GetScratchMemTop ();
  hole_k = domain->hole_k;
  privateG1 = domain->privateG1;
  sprivate  = domain->SprivateG1;
  nlpr      = &nlsprivate->nlpr;
  nfunc_a = privateG1->nfunc_a;
  nfunc_c = sprivate->nsfunc_c;
  nfunc_d = sprivate->nsfunc_d;
  nfacd = nfunc_a+nfunc_c+nfunc_d;
  nfunc = nfacd-nconstr;
  nk = sprivate->nk;
  m2 = sprivate->m2;
  C = (double)(5+nk*m2)*sprivate->C1s/nlpr->ddiam;

  bs1 = nfunc_c/hole_k;
  bs2 = bs1+nfunc_d+nfunc_a;
  bs3 = bs2-nconstr;
  diagblsize    = (bs1*(bs1+1))/2;
  subdiagblsize = bs1*bs1;
  sideblsize    = bs1*bs2;
  esideblsize   = bs1*bs3;
  asize = pkn_Block3ArraySize ( hole_k-1, bs1, bs2 );
  esize = pkn_Block3ArraySize ( hole_k-1, bs1, bs3 );
  coeff   = pkv_GetScratchMemd ( nfacd );
  grad    = pkv_GetScratchMemd ( nfacd );
  hessian = pkv_GetScratchMemd ( asize );
  E22     = pkv_GetScratchMemd ( asize );
  cE22    = pkv_GetScratchMemd ( esize );
  cT      = pkv_GetScratchMemd ( bs2*nconstr );
  aa      = pkv_GetScratchMemd ( 2*nconstr );
  D1      = pkv_GetScratchMemd ( (nconstr*(nconstr+1))/2 );
  y       = pkv_GetScratchMemd ( nfacd );
  y1      = pkv_GetScratchMemd ( nfacd );
  M       = pkv_GetScratchMemd ( (bs2*(bs2+1))/2 );
  f       = pkv_GetScratchMemd ( nfacd );
  if ( !coeff || !grad || !hessian || !E22 || !cE22 ||
       !cT || !aa || !D1 || !y || !y1 || !M || !f ) {
    domain->error_code = G1H_ERROR_NO_SCRATCH_MEMORY;
    goto failure;
  }
  hij = &hessian[pkn_Block3FindBlockPos ( hole_k-1, bs1, bs2, 1, 0 )];
  hkk = &hessian[pkn_Block3FindBlockPos ( hole_k-1, bs1, bs2,
                                          hole_k-1, hole_k-1 )];
  hki = &hessian[pkn_Block3FindBlockPos ( hole_k-1, bs1, bs2, hole_k-1, 0 )];
  E22ij = &E22[pkn_Block3FindBlockPos ( hole_k-1, bs1, bs3, 1, 0 )];
  E22kk = &E22[pkn_Block3FindBlockPos ( hole_k-1, bs1, bs3,
                                        hole_k-1, hole_k-1 )];
  E22ki = &E22[pkn_Block3FindBlockPos ( hole_k-1, bs1, bs3, hole_k-1, 0 )];

        /* setup the initial point */
  pkv_Selectd ( nfacd, 1, 3, 1, &nlpr->acoeff[0].z, coeff );

        /* step 1: decompose the constraint equations matrix */
  pkv_TransposeMatrixd ( nconstr, bs2, nfacd, &SCmat[nfacd-bs2],
                         nconstr, cT );
  pkn_QRDecomposeMatrixd ( bs2, nconstr, cT, aa );
  for ( i = 0; i < nconstr; i++ )
    for ( j = i; j < nconstr; j++ )
      D1[pkn_SymMatIndex(i,j)] = cT[i*nconstr+j];

        /* Newton iterations */
  for ( itn = ktn = 0; ; itn++ ) {
          /* step 2 */
    memcpy ( y, coeff, nfacd*sizeof(double) );
    pkn_multiReflectVectord ( bs2, nconstr, cT, aa, 1, 1, &y[nfacd-bs2] );
    pkn_LowerTrMatrixSolved ( nconstr, D1, 1, 3, &nlpr->rhole_cp[12*hole_k+1].z,
                              1, &y[nfacd-bs2] );
    for ( i = nfacd-bs2; i < nfacd-bs2+nconstr; i++ )
      y[i] = -y[i];

          /* step 3 */
    if ( !_g1hq2_ComputeSplNLFuncGradHessiand ( domain, nlsprivate, coeff, C,
                                                &func, grad, hessian, E22 ) )
      goto failure;

    if ( !pkn_ComputeQTSQd ( bs2, hkk, nconstr, cT, aa, M ) )
      goto failure;
    for ( i = 0; i < bs3; i++ )
      for ( j = i; j < bs3; j++ )
        E22kk[pkn_SymMatIndex(i,j)] = M[pkn_SymMatIndex(nconstr+i,nconstr+j)];
    for ( i = 0; i < hole_k-1; i++ )
      pkn_multiReflectVectord ( bs2, nconstr, cT, aa, bs1, bs1,
                                &hki[i*sideblsize] );
    pkv_Selectd ( hole_k-1, esideblsize, sideblsize, esideblsize,
                  &hki[nconstr*bs1], E22ki );
    memcpy ( E22, hessian, (hole_k-1)*diagblsize*sizeof(double) );
    memcpy ( E22ij, hij, (hole_k-2)*subdiagblsize*sizeof(double) );
    memcpy ( cE22, E22, esize*sizeof(double) );
          /* step 4 */
    memcpy ( f, grad, nfacd*sizeof(double) );
    pkn_multiReflectVectord ( bs2, nconstr, cT, aa, 1, 1, &f[nfacd-bs2] );
    memmove ( &f[nfacd-bs2], &f[nfacd-bs3], bs3*sizeof(double) );
    gn = sqrt ( pkn_ScalarProductd ( nfunc, f, f ) );
    if ( itn == 0 )
      gn0 = gn;
          /* step 5 */
    if ( (positive = pkn_Block3CholeskyDecompMd ( hole_k-1, bs1, bs3, E22 ) ) ) {
      pkn_Block3LowerTrMSolved ( hole_k-1, bs1, bs3, E22, 1, 1, f );
      pkn_Block3UpperTrMSolved ( hole_k-1, bs1, bs3, E22, 1, 1, f );
    }
    else {

printf ( "! " );

      pkn_Block3SymMatrixMultd ( hole_k-1, bs1, bs3, cE22, 1, 1, f, 1, y1 );
      aux = pkn_ScalarProductd ( nfunc, f, y1 );
      if ( aux <= 0.0 || aux < EPSF*gn )
        goto failure;
      pkn_MultMatrixNumd ( 1, nfunc, 1, f, gn/aux, 1, f );
    }
          /* step 6 */
    dyn = sqrt ( pkn_ScalarProductd ( nfunc, f, f ) );
    memmove ( &f[nfacd-bs3], &f[nfacd-bs2], bs3*sizeof(double) );
    dco = sqrt ( pkn_ScalarProductd ( nfacd, coeff, coeff ) );
    for ( aux = 1.0; aux >= EPSF; aux *= 0.5 ) {
      for ( i = 0; i < nfacd-bs2; i++ )
        y1[i] = y[i]+aux*f[i];
      memcpy ( &y1[nfacd-bs2], &y[nfacd-bs2], nconstr*sizeof(double) );
      for ( i = nfacd-bs3; i < nfacd; i++ )
        y1[i] = y[i]+aux*f[i];
      pkn_multiInvReflectVectord ( bs2, nconstr, cT, aa, 1, 1, &y1[nfacd-bs2] );
      func1 = _g1hq2_ComputeSplNLFuncd ( domain, nlsprivate, y1, C );
      if ( func1 < func )
        break;
    }
    memcpy ( coeff, y1, nfacd*sizeof(double) );
    func = func1;
    ktn ++;

    if ( positive && aux > 0.1 ) {
        /* Now the Hessian is positive-definite; */
        /* as it is expensive to compute, we try to make some */
        /* extra iterations with the same Hessian */
      for ( jtn = 0; jtn < 10; jtn++ ) {
            /* step 2' */
        memcpy ( y, coeff, nfacd*sizeof(double) );
        pkn_multiReflectVectord ( bs3, nconstr, cT, aa, 1, 1, &y[nfacd-bs2] );
        pkn_LowerTrMatrixSolved ( nconstr, D1, 1, 3, &nlpr->rhole_cp[12*hole_k+1].z,
                                  1, &y[nfacd-bs2] );
        for ( i = nfacd-bs2; i < nfacd-bs3; i++ )
          y[i] = -y[i];
            /* step 3' */
        _g1hq2_ComputeSplNLFuncGradd ( domain, nlsprivate, coeff, C, &func0, grad );
            /* step 4' */
        memcpy ( f, grad, nfacd*sizeof(double) );
        pkn_multiReflectVectord ( bs2, nconstr, cT, aa, 1, 1, &f[nfacd-bs2] );
        memmove ( &f[nfacd-bs2], &f[nfacd-bs3], bs3*sizeof(double) );
        gn1 = sqrt ( pkn_ScalarProductd ( nfunc, f, f ) );
            /* step 5' */
        pkn_Block3LowerTrMSolved ( hole_k-1, bs1, bs3, E22, 1, 1, f );
        pkn_Block3UpperTrMSolved ( hole_k-1, bs1, bs3, E22, 1, 1, f );
            /* step 6' */
        dyn1 = sqrt ( pkn_ScalarProductd ( nfunc, f, f ) );
        memmove ( &f[nfacd-bs3], &f[nfacd-bs2], bs3*sizeof(double) );
        for ( i = 0; i < nfacd-bs2; i++ )
          y1[i] = y[i]+f[i];
        memcpy ( &y1[nfacd-bs2], &y[nfacd-bs2], nconstr*sizeof(double) );
        for ( i = nfacd-bs3; i < nfacd; i++ )
          y1[i] = y[i]+f[i];
        pkn_multiReflectVectord ( bs2, nconstr, cT, aa, 1, 1, &y1[nfacd-bs2] );

        func1 = _g1hq2_ComputeSplNLFuncd ( domain, nlsprivate, y1, C );
        if ( func1 > func0 || gn1 >= 0.5*gn )
          break;

printf ( "    func = %f, gn = %f\n", func1, gn );

        memcpy ( coeff, y1, nfacd*sizeof(double) );
        func = func1;
        ktn ++;
        gn = gn1;
        aux = 1.0;
        dco = sqrt ( pkn_ScalarProductd ( nfacd, coeff, coeff ) );
        dyn = dyn1;
      }
    }
    if ( _g1h_StopItd ( itn, gn0, gn, dco, dyn, aux ) )
      break;
  }

printf ( "itn = %d, ktn = %d\n", itn+1, ktn );
printf ( "func = %f\n", func1 );

  pkv_Selectf ( nfacd, 1, 1, 3, coeff, &nlpr->acoeff[0].z );

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
#undef EPSF
} /*g1hq2_SplNLConstrNewtond*/

boolean g1h_Q2NLSplFillHoleConstrd ( GHoleDomaind *domain, const point3d *hole_cp,
                     int nconstr, const vector3d *constr,
                     double *acoeff, void *usrptr,
                     void (*outpatch) ( int n, int lknu, const double *knu,
                                        int m, int lknv, const double *knv,
                                        const point3d *cp, void *usrptr ) )
{
  typedef void outscf ( int n, int lknu, const double *knu,
                        int m, int lknv, const double *knv,
                        const double *cp, void *usrptr );

  void   *sp;
  G1HNLSPrivated     *nlsprivate;
  G1HNLPrivated      *nlprivate;
  G1HolePrivateRecd  *privateG1;
  G1HoleSPrivateRecd *sprivate;
  int    hole_k, nfunc_a, nfunc_c, nfunc_d;
  double *fc00;

  sp = pkv_GetScratchMemTop ();
  privateG1 = domain->privateG1;
  sprivate = domain->SprivateG1;
  if ( !sprivate ) {
    domain->error_code = G1H_ERROR_SPLINE_BASIS_NOT_READY;
    goto failure;
  }
  if ( !sprivate->Q2SAMat )
    if ( !g1h_Q2SplComputeFormMatrixd ( domain ) )
      goto failure;
  if ( !sprivate->Q2SLMat )
    if ( !g1h_Q2SplDecomposeMatrixd ( domain ) )
      goto failure;
  nlsprivate = _g1hq2_AllocSplNLPrd ( domain );
  if ( !nlsprivate )
    goto failure;
  _g1h_nlprivd = nlprivate = &nlsprivate->nlpr;
  hole_k = domain->hole_k;
  nfunc_a = privateG1->nfunc_a;
  nfunc_c = sprivate->nsfunc_c;
  nfunc_d = sprivate->nsfunc_d;
  if ( !_g1hq2_InitSplNLprd ( domain, nlsprivate, nconstr ) )
    goto failure;

  if ( !_g1h_ComputeNLNormald ( domain, nlprivate, hole_cp ) )
    goto failure;
  g1h_ReflectVectorsd ( 12*hole_k+1, hole_cp, nlprivate->rhole_cp );
  g1h_ReflectVectorsd ( nconstr, constr, &nlprivate->rhole_cp[12*hole_k+1] );

  nlprivate->auxc = 0;
  if ( !g1h_Q2SplFillHoleConstrd ( domain, 3, (double*)nlprivate->rhole_cp,
                    nconstr, (double*)&nlprivate->rhole_cp[12*hole_k+1],
                    (double*)nlprivate->acoeff, NULL, g1h_splnloutpatchd ) )
    goto failure;

  if ( !_g1hq2_TabSplNLBasisFunctionsd ( domain, nlsprivate ) )
    goto failure;
  if ( !_g1hq2_TabSplNLBasisFJumpsd ( domain, nlsprivate ) )
    goto failure;

  if ( !g1hq2_SplNLConstrNewtond ( domain, nlsprivate,
                                   nconstr, sprivate->SCmat ) )
    goto failure;

  g1h_ReflectVectorsd ( nfunc_a+nfunc_c+nfunc_d, nlprivate->acoeff,
                        nlprivate->acoeff );
  if ( acoeff )
    memcpy ( acoeff, nlprivate->acoeff,
             (nfunc_a+nfunc_c+nfunc_d)*sizeof(vector3d) );

  fc00 = pkv_GetScratchMemd ( (G1_CROSSDEGSUM+4)*2*hole_k*3 );
  if ( !fc00 ) {
    domain->error_code = G1H_ERROR_NO_SCRATCH_MEMORY;
    goto failure;
  }
  if ( !_g1h_SetSplRightSided ( domain, 3, (double*)hole_cp,
                                sprivate->Q2SBMat, fc00, NULL ) )
    goto failure;

  if ( !_g1h_OutputSplPatchesd ( domain, 3, (double*)nlprivate->acoeff, fc00,
                                 usrptr, (outscf*)outpatch ) )
    goto failure;

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*g1h_Q2NLSplFillHoleConstrd*/

/* ///////////////////////////////////////////////////////////////////////// */
boolean g1h_Q2NLSplFillHoleAltConstrd ( GHoleDomaind *domain, const point3d *hole_cp,
                     int naconstr, const double *constr,
                     double *acoeff, void *usrptr,
                     void (*outpatch) ( int n, int lknu, const double *knu,
                                        int m, int lknv, const double *knv,
                                        const point3d *cp, void *usrptr ) )
{
  typedef void outscf ( int n, int lknu, const double *knu,
                        int m, int lknv, const double *knv,
                        const double *cp, void *usrptr );

  void    *sp;
  G1HNLPrivated      *nlprivate;
  G1HNLSPrivated     *nlsprivate;
  G1HolePrivateRecd  *privateG1;
  G1HoleSPrivateRecd *sprivate;
  double  *fc00, *cmat, *saveconstr, *newconstr, *rsconstr;
  int     hole_k, nfunc_a, nfunc_c, nfunc_d, nbf, i, j;
  boolean restore;

  sp = pkv_GetScratchMemTop ();
  restore = false;
  privateG1 = domain->privateG1;
  sprivate = domain->SprivateG1;
  if ( !sprivate )
    goto failure;
  if ( !sprivate->Q2SAMat )
    if ( !g1h_Q2SplComputeFormMatrixd ( domain ) )
      goto failure;
  if ( !sprivate->Q2SLMat )
    if ( !g1h_Q2SplDecomposeMatrixd ( domain ) )
      goto failure;
  nlsprivate = _g1hq2_AllocSplNLPrd ( domain );
  if ( !nlsprivate ) {
    domain->error_code = G1H_ERROR_NO_SCRATCH_MEMORY;
    goto failure;
  }
  _g1h_nlprivd = nlprivate = &nlsprivate->nlpr;
  if ( !_g1hq2_InitSplNLprd ( domain, nlsprivate, naconstr ) )
    goto failure;
  if ( naconstr != sprivate->splnaconstr || sprivate->splacdim != 3 ) {
    domain->error_code = G1H_ERROR_INCONSISTENT_CONSTR;
    goto failure;
  }
  hole_k = domain->hole_k;
  nfunc_a = privateG1->nfunc_a;
  nfunc_c = sprivate->nsfunc_c;
  nfunc_d = sprivate->nsfunc_d;
  nbf = nfunc_a+nfunc_c+nfunc_d;
  saveconstr = pkv_GetScratchMemd ( 7*nbf*naconstr );
  if ( !saveconstr ) {
    domain->error_code = G1H_ERROR_NO_SCRATCH_MEMORY;
    goto failure;
  }
  newconstr = &saveconstr[3*nbf*naconstr];
  cmat = &newconstr[3*nbf*naconstr];
  memcpy ( saveconstr, sprivate->ASCmat, nbf*naconstr*3*sizeof(double) );
  memcpy ( newconstr, sprivate->ASCmat, nbf*naconstr*3*sizeof(double) );
  if ( !_g1h_ComputeNLNormald ( domain, nlprivate, hole_cp ) )
    goto failure;
  g1h_ReflectVectorsd ( 12*hole_k+1, hole_cp, nlprivate->rhole_cp );
  if ( !_g1h_ReflectSplAltConstrMatrixd ( domain, &nlprivate->reflv, newconstr ) )
    goto failure;

  restore = true;
  if ( !g1h_SetSplAltConstraintMatrixd ( domain, 3, naconstr, newconstr ) )
    goto failure;
  nlprivate->auxc = 0;
  if ( !g1h_Q2SplFillHoleAltConstrd ( domain, 3, (double*)nlprivate->rhole_cp,
                      naconstr, constr,
                      (double*)nlprivate->acoeff, NULL, g1h_splnloutpatchd ) )
    goto failure;

  pkv_Selectd ( naconstr, nbf, 3*nbf, nbf, &newconstr[2*nbf], cmat );
  rsconstr = &nlprivate->rhole_cp[12*hole_k+1].z;
  pkv_Selectd ( naconstr, 1, 1, 3, constr, rsconstr );
  for ( j = 0; j < naconstr; j++ )
    for ( i = 0; i < nbf; i++ )
      rsconstr[3*j] += newconstr[3*j*nbf+i]*nlprivate->acoeff[i].x +
                       newconstr[(3*j+1)*nbf+i]*nlprivate->acoeff[i].y;

  if ( !_g1hq2_TabSplNLBasisFunctionsd ( domain, nlsprivate ) )
    goto failure;
  if ( !_g1hq2_TabSplNLBasisFJumpsd ( domain, nlsprivate ) )
    goto failure;

  if ( !g1hq2_SplNLConstrNewtond ( domain, nlsprivate, naconstr, cmat ) )
    goto failure;

  g1h_ReflectVectorsd ( nbf, nlprivate->acoeff, nlprivate->acoeff );
  if ( acoeff )
    memcpy ( acoeff, nlprivate->acoeff, nbf*sizeof(vector3d) );

  fc00 = pkv_GetScratchMemd ( (G1_CROSSDEGSUM+4)*2*hole_k*3 );
  if ( !fc00 )
    goto failure;
  if ( !_g1h_SetSplRightSided ( domain, 3, (double*)hole_cp,
                                sprivate->SBMat, fc00, NULL ) )
    goto failure;

  if ( !_g1h_OutputSplPatchesd ( domain, 3, (double*)nlprivate->acoeff, fc00,
                                 usrptr, (outscf*)outpatch ) )
    goto failure;

  g1h_SetSplAltConstraintMatrixd ( domain, 3, naconstr, saveconstr );
  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  if ( restore )
    g1h_SetSplAltConstraintMatrixd ( domain, 3, naconstr, saveconstr );
  pkv_SetScratchMemTop ( sp );
  return false;
} /*g1h_Q2NLSplFillHoleAltConstrd*/

