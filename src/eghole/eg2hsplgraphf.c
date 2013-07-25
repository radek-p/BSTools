
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2007, 2009                            */
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

/* ////////////////////////////////////////////////////////////////////////// */
void g2h_DrawSplBasFuncNumf ( GHoleDomainf *domain,
                        int *nfunc_a, int *nfunc_b, int *nfunc_c, int *nfunc_d )
{
  G2HolePrivateRecf  *privateG2;
  G2HoleSPrivateRecf *sprivate;

  privateG2 = domain->privateG2;
  sprivate = domain->SprivateG2;
  *nfunc_a = privateG2->nfunc_a;
  *nfunc_b = privateG2->nfunc_b;
  *nfunc_c = sprivate->nsfunc_c;
  *nfunc_d = sprivate->nsfunc_d;
} /*g2h_DrawSplBasFuncNumf*/

void g2h_DrawSplBasAuxPatchesf ( GHoleDomainf *domain, int fn,
               void (*drawpatch) ( int n, int lknu, const float *knu,
                                   int m, int lknv, const float *knv,
                                   const point3f *cp ) )
{
  void     *sp;
  G2HolePrivateRecf  *privateG2;
  G2HoleSPrivateRecf *sprivate;
  int      i, hole_k, degu, lknu, lkn, nrows, nzc;
  int      nfunc_a, nfunc_b, nfunc_c, nfunc_d;
  vector2f *auxp, *omc, *omcd, *omcdd;
  float    *bezknots, *omcknots, *auxpknots;
  float    *fcomc, *fcomcd, *fcomcdd, *fcaux;
  point3f  *basp;

  sp = pkv_GetScratchMemTop ();
  hole_k = domain->hole_k;
  privateG2 = domain->privateG2;
  sprivate = domain->SprivateG2;
  if ( !sprivate )
    return;
  nfunc_a = privateG2->nfunc_a;
  nfunc_b = privateG2->nfunc_b;
  nfunc_c = sprivate->nsfunc_c;
  nfunc_d = sprivate->nsfunc_d;
  if ( fn < 0 || fn >= nfunc_a+nfunc_b+nfunc_c+nfunc_d )
    goto finish;
  bezknots = sprivate->bezknots;
  omcknots = sprivate->omcknots;
  lkn = sprivate->lastomcknot;
  omc = privateG2->omc;
  omcd = &omc[(G2H_OMCDEG+1)*hole_k];
  omcdd = &omcd[G2H_OMCDEG*hole_k];

      /* determine the "u" knot sequence for the auxiliary basis function patch */
  degu = 0;
  if ( !mbs_FindBSCommonKnotSequencef ( &degu, &lknu, &auxpknots, 3,
                                        G2H_OMCDEG, lkn, omcknots,
                                        G2H_OMCDEG-1, lkn-2, &omcknots[1],
                                        G2H_OMCDEG-2, lkn-4, &omcknots[2] ) )
    goto finish;
      /* now the knot sequence is in the array on the scratch memory stack */
  nrows = lknu-degu;
  auxp = pkv_GetScratchMem ( 3*nrows*sizeof(point2f) );
  basp = pkv_GetScratchMem ( 3*nrows*sizeof(point3f) );
  fcomc = pkv_GetScratchMemf ( 4*nrows );
  if ( !auxp || !basp || !fcomc )
    goto finish;
  fcomcd = &fcomc[nrows];
  fcomcdd = &fcomcd[nrows];
  fcaux = &fcomcdd[nrows];

  memset ( basp, 0, 3*nrows*sizeof(point3f) );
  for ( i = 0; i < hole_k; i++ ) {
        /* construct the auxiliary domain patch B-spline representation */
    mbs_AdjustBSCRepC2f ( G2H_OMCDEG, 2*G2H_OMCDEG+1,
              &bezknots[G2H_FINALDEG-G2H_OMCDEG], omc,
              degu, lknu, auxpknots, auxp );
    mbs_AdjustBSCRepC2f ( G2H_OMCDEG-1, 2*G2H_OMCDEG-1,
              &bezknots[G2H_FINALDEG-G2H_OMCDEG+1], omcd,
              degu, lknu, auxpknots, &auxp[nrows] );
    mbs_AdjustBSCRepC2f ( G2H_OMCDEG-2, 2*G2H_OMCDEG-3,
              &bezknots[G2H_FINALDEG-G2H_OMCDEG+2], omcdd,
              degu, lknu, auxpknots, &auxp[2*nrows] );
    pkv_Selectf ( nrows, 2, 2, 3, auxp, basp );
    pkn_AddMatrixMf ( nrows, 2, 3, &basp[0].x, 2, &auxp[nrows].x, 0.5,
                      3, &basp[nrows].x );
    pkn_AddMatrixf ( nrows, 2, 3, &basp[0].x, 2, &auxp[nrows].x,
                     3, &basp[2*nrows].x );
    pkn_AddMatrixMf ( nrows, 2, 3, &basp[2*nrows].x, 2, &auxp[2*nrows].x, 0.5,
                      3, &basp[2*nrows].x );

        /* construct the basis function auxiliary patch representation */
        /* - nonzero only for the functions of the group A, B and D */
    if ( fn < nfunc_a ) {               /* to be written */
      /* ************ */
    }
    else if ( fn < nfunc_a+nfunc_b ) {  /* to be written */
      /* ************ */
    }
    else if ( fn < nfunc_a+nfunc_b+nfunc_c )  /* zero */
      ;
    else {
      _g2h_GetSplDBasisAuxpf ( domain, fn, i, &nzc, fcomc, fcomcd, fcomcdd );
      mbs_AdjustBSCRepC1f ( G2H_OMCDEG, lkn, omcknots, fcomc,
                            degu, lknu, auxpknots, fcaux );
      pkv_Selectf ( nrows, 1, 1, 3, fcaux, &basp[0].z );
      mbs_AdjustBSCRepC1f ( G2H_OMCDEG-1, lkn-2, &omcknots[1], fcomcd,
                            degu, lknu, auxpknots, fcaux );
      pkn_AddMatrixMf ( nrows, 1, 3, &basp[0].z, 1, fcaux, 0.5,
                        3, &basp[nrows].z );
      pkn_AddMatrixf ( nrows, 1, 3, &basp[0].z, 1, fcaux, 3, &basp[2*nrows].z );
      mbs_AdjustBSCRepC1f ( G2H_OMCDEG-2, lkn-4, &omcknots[2], fcomcdd,
                            degu, lknu, auxpknots, fcaux );
      pkn_AddMatrixMf ( nrows, 1, 3, &basp[2*nrows].z, 1, fcaux, 0.5,
                        3, &basp[2*nrows].z );
    }
        /* output it */
    drawpatch ( 2, 5, &bezknots[G2H_FINALDEG-2], degu, lknu, auxpknots, basp );

    omc += (G2H_OMCDEG+1);
    omcd += G2H_OMCDEG;
    omcdd += (G2H_OMCDEG-1);
  }

finish:
  pkv_SetScratchMemTop ( sp );
} /*g2h_DrawBasAuxPatchesf*/

void g2h_DrawSplBasFunctionf ( GHoleDomainf *domain, int fn,
             void (*drawpatch) ( int n, int lknu, const float *knu,
                                 int m, int lknv, const float *knv,
                                 const point3f *cp ) )
{
  void     *sp;
  G2HolePrivateRecf  *privateG2;
  G2HoleSPrivateRecf *sprivate;
  int      hole_k, nfunc_a, nfunc_b, nfunc_c, nfunc_d, bfn;
  int      i, j, n, m, k0, k1, k2;
  vector2f *c00, *c01, *c02, *c10, *c11, *c12,
           *d00, *d01, *d02, *d10, *d11, *d12;
  float    *fc00, *fc01, *fc02, *fc10, *fc11, *fc12,
           *fd00, *fd01, *fd02, *fd10, *fd11, *fd12;
  point2f  *di, *aux;
  point3f  *cp;
  int      degu, degv, lastomcknot, lastpvknot, lastpvvknot,
           lastcknot, lastfpknot, lastukn, lastvkn, ncp;
  float    *bezknots, *omcknots, *pvknots, *pvvknots, *cknots;
  float    *fcomc, *pv, *pvv, *pu, *puu, *auxbfp, *bfp, *uknots, *vknots;
  float    zero[2] = {0.0,0.0};

  sp = pkv_GetScratchMemTop ();
  hole_k = domain->hole_k;
  privateG2 = domain->privateG2;
  sprivate = domain->SprivateG2;
  if ( !sprivate )
    return;
  nfunc_a = privateG2->nfunc_a;
  nfunc_b = privateG2->nfunc_b;
  nfunc_c = sprivate->nsfunc_c;
  nfunc_d = sprivate->nsfunc_d;
  if ( fn < 0 || fn >= nfunc_a+nfunc_b+nfunc_c+nfunc_d )
    goto finish;

  bezknots = sprivate->bezknots;
  cknots = sprivate->cknots;
  lastcknot = sprivate->lastcknot;
  omcknots = sprivate->omcknots;
  lastomcknot = sprivate->lastomcknot;
  pvknots = sprivate->pvknots;
  lastpvknot = sprivate->lastpvknot;
  pvvknots = sprivate->pvvknots;
  lastpvvknot = sprivate->lastpvvknot;
  lastfpknot = sprivate->lastfpknot;
  
  n = lastfpknot-G2H_FINALDEG;
  m = G2H_FINALDEG+1;
  di = (point2f*)pkv_GetScratchMem ( m*m*sizeof(point2f) );
  aux = (point2f*)pkv_GetScratchMem ( n*m*sizeof(point2f) );
  cp = (point3f*)pkv_GetScratchMem ( n*n*sizeof(point3f) );
  auxbfp = pkv_GetScratchMemf ( 6*(lastpvvknot-G2_CROSS02DEG) );
  bfp = pkv_GetScratchMemf ( (lastpvvknot-G2_CROSS02DEG)*(lastpvvknot-G2_CROSS02DEG) );
  fcomc = pkv_GetScratchMemf ( lastomcknot+2*(lastpvknot+lastpvvknot)-
                               (G2_CROSS00DEG+2*(G2_CROSS01DEG+G2_CROSS02DEG)) );
  uknots = pkv_GetScratchMemf ( 2*(lastfpknot+1) );
  if ( !di || !aux || !cp || !auxbfp || !bfp || !fcomc || !uknots )
    goto finish;
  pv = &fcomc[lastomcknot-G2_CROSS00DEG];
  pvv = &pv[lastpvknot-G2_CROSS01DEG];
  pu = &pvv[lastpvvknot-G2_CROSS02DEG];
  puu = &pu[lastpvknot-G2_CROSS01DEG];
  vknots = &uknots[lastfpknot+1];

  for ( i = 0; i < hole_k; i++ ) {
        /* get the domain patch */
    _g2h_GetDiPatchCurvesf ( domain, i,
                             &c00, &c01, &c02, &c10, &c11, &c12,
                             &d00, &d01, &d02, &d10, &d11, &d12 );
    mbs_BezC2CoonsToBezf ( 2,
       G2_CROSS00DEG, (float*)c00, G2_CROSS01DEG, (float*)c01, G2_CROSS02DEG, (float*)c02,
       3, (float*)c10, G2_CROSS11DEG, (float*)c11, G2_CROSS12DEG, (float*)c12,
       G2_CROSS00DEG, (float*)d00, G2_CROSS01DEG, (float*)d01, G2_CROSS02DEG, (float*)d02,
       3, (float*)d10, G2_CROSS11DEG, (float*)d11, G2_CROSS12DEG, (float*)d12,
       &degu, &degv, (float*)di );

        /* get the basis function patch */
    if ( fn < nfunc_a ) {                       /* block A */
      _g2h_GetBFAPatchCurvesf ( domain, fn, i,
                                &fc00, &fc01, &fc02, &fd00, &fd01, &fd02 );
      mbs_BezC2CoonsToBezf ( 1,
                    G2_CROSS00DEG, fc00, G2_CROSS01DEG, fc01, G2_CROSS02DEG, fc02,
                    1, zero, 1, zero, 1, zero,
                    G2_CROSS00DEG, fd00, G2_CROSS01DEG, fd01, G2_CROSS02DEG, fd02,
                    1, zero, 1, zero, 1, zero,
                    &degu, &degv, bfp );
      ncp = m*m;
      pkv_Selectf ( ncp, 2, 2, 3, (float*)di, cp );
      pkv_Selectf ( ncp, 1, 1, 3, bfp, &cp[0].z );
      drawpatch ( G2H_FINALDEG, 2*G2H_FINALDEG+1, bezknots,
                  G2H_FINALDEG, 2*G2H_FINALDEG+1, bezknots, cp );
    }
    else if ( fn < nfunc_a+nfunc_b ) {          /* block B */
      _g2h_GetBFBPatchCurvesf ( domain, fn-nfunc_a, i,
                                &fc00, &fc01, &fc02, &fc10, &fc11, &fc12,
                                &fd00, &fd01, &fd02, &fd10, &fd11, &fd12 );
      mbs_BezC2CoonsToBezf ( 1,
                    G2_CROSS00DEG, fc00, G2_CROSS01DEG, fc01, G2_CROSS02DEG, fc02,
                    G2_CROSS10DEG, fc10, G2_CROSS11DEG, fc11, G2_CROSS12DEG, fc12,
                    G2_CROSS00DEG, fd00, G2_CROSS01DEG, fd01, G2_CROSS02DEG, fd02,
                    G2_CROSS10DEG, fd10, G2_CROSS11DEG, fd11, G2_CROSS12DEG, fd12,
                    &degu, &degv, bfp );
      ncp = m*m;
      pkv_Selectf ( ncp, 2, 2, 3, (float*)di, cp );
      pkv_Selectf ( ncp, 1, 1, 3, bfp, &cp[0].z );
      drawpatch ( G2H_FINALDEG, 2*G2H_FINALDEG+1, bezknots,
                  G2H_FINALDEG, 2*G2H_FINALDEG+1, bezknots, cp );
    }
    else if ( fn < nfunc_a+nfunc_b+nfunc_c ) {  /* block C */
      n = lastcknot-G2H_FINALDEG;
      m = G2H_FINALDEG+1;
      mbs_multiAdjustBSCRepf ( 1, 2*(degv+1), degu, 2*degu+1,
                               &bezknots[G2H_FINALDEG-degu], 0, (float*)di,
                               G2H_FINALDEG, lastcknot, cknots, 0, (float*)aux );
      mbs_multiAdjustBSCRepf ( n, 2, degv, 2*degv+1,
                               &bezknots[G2H_FINALDEG-degv], 2*(degv+1), (float*)aux,
                               G2H_FINALDEG, lastcknot, cknots, 2*n, (float*)cp );
      pkv_Rearrangef ( n*n, 2, 2, 3, cp );
      for ( j = 0; j < n*n; j++ )
        cp[j].z = 0.0;

      bfn = fn - (nfunc_a+nfunc_b);
      k0 = bfn / sprivate->csize;  /* which Omega_i? */
      if ( k0 == i ) {
        bfn = bfn % sprivate->csize;  /* which function in the block? */
        k1  = bfn / (lastcknot-G2H_FINALDEG-6);
        k2  = bfn % (lastcknot-G2H_FINALDEG-6);
        cp[(k1+3)*n+k2+3].z = 1.0;
      }

      drawpatch ( G2H_FINALDEG, lastcknot, cknots,
                  G2H_FINALDEG, lastcknot, cknots, cp );
    }
    else {                                      /* block D */
      n = lastpvvknot-G2H_FINALDEG;
      m = G2H_FINALDEG+1;
      mbs_multiAdjustBSCRepf ( 1, 2*(degv+1), degu, 2*degu+1,
                       &bezknots[G2H_FINALDEG-degu], 0, (float*)di,
                       G2H_FINALDEG, lastpvvknot, pvvknots, 0, (float*)aux );
      mbs_multiAdjustBSCRepf ( n, 2, degv, 2*degv+1,
                       &bezknots[G2H_FINALDEG-degv], 2*(degv+1), (float*)aux,
                       G2H_FINALDEG, lastpvvknot, pvvknots, 2*n, (float*)cp );
      pkv_Rearrangef ( n*n, 2, 2, 3, cp );
      for ( j = 0; j < n*n; j++ )
        cp[j].z = 0.0;

      bfn = fn - (nfunc_a+nfunc_b+nfunc_c);
      k0 = bfn / sprivate->dsize;
      if ( k0 == i ) {
        _g2h_GetSplDBasisCrossDerf ( domain, fn, i, fcomc, pv, pvv, pu, puu );
        mbs_BSC2CoonsToBSf ( 1, G2_CROSS00DEG, lastomcknot, omcknots, fcomc,
            G2_CROSS01DEG, lastpvknot, pvknots, pv,
            G2_CROSS02DEG, lastpvvknot, pvvknots, pvv,
            1, 3, &bezknots[G2H_FINALDEG-1], zero,
            1, 3, &bezknots[G2H_FINALDEG-1], zero,
            1, 3, &bezknots[G2H_FINALDEG-1], zero,
            1, 3, &bezknots[G2H_FINALDEG-1], zero,
            1, 3, &bezknots[G2H_FINALDEG-1], zero,
            1, 3, &bezknots[G2H_FINALDEG-1], zero,
            1, 3, &bezknots[G2H_FINALDEG-1], zero,
            1, 3, &bezknots[G2H_FINALDEG-1], zero,
            1, 3, &bezknots[G2H_FINALDEG-1], zero,
            &degu, &lastukn, uknots, &degv, &lastvkn, vknots, auxbfp );
        mbs_multiAdjustBSCRepf ( n, 1,
                      degv, lastvkn, vknots, lastvkn-degv, auxbfp,
                      G2H_FINALDEG, lastpvvknot, pvvknots, n, bfp );
        pkv_Selectf ( n*n, 1, 1, 3, bfp, &cp[0].z );
      }
      else if ( k0 == (i+1) % hole_k ) {
        _g2h_GetSplDBasisCrossDerf ( domain, fn, k0, fcomc, pv, pvv, pu, puu );
        mbs_BSC2CoonsToBSf ( 1,
            1, 3, &bezknots[G2H_FINALDEG-1], zero,
            1, 3, &bezknots[G2H_FINALDEG-1], zero,
            1, 3, &bezknots[G2H_FINALDEG-1], zero,
            1, 3, &bezknots[G2H_FINALDEG-1], zero,
            1, 3, &bezknots[G2H_FINALDEG-1], zero,
            1, 3, &bezknots[G2H_FINALDEG-1], zero,
            G2_CROSS00DEG, lastomcknot, omcknots, fcomc,
            G2_CROSS01DEG, lastpvknot, pvknots, pu,
            G2_CROSS02DEG, lastpvvknot, pvvknots, puu,
            1, 3, &bezknots[G2H_FINALDEG-1], zero,
            1, 3, &bezknots[G2H_FINALDEG-1], zero,
            1, 3, &bezknots[G2H_FINALDEG-1], zero,
            &degu, &lastukn, uknots, &degv, &lastvkn, vknots, auxbfp );
        mbs_multiAdjustBSCRepf ( 1, n,
                                 degu, lastukn, uknots, 0, auxbfp,
                                 G2H_FINALDEG, lastpvvknot, pvvknots, n, bfp );
        pkv_Selectf ( n*n, 1, 1, 3, bfp, &cp[0].z );
      }

      drawpatch ( G2H_FINALDEG, lastpvvknot, pvvknots,
                  G2H_FINALDEG, lastpvvknot, pvvknots, cp );
    }
  }

finish:
  pkv_SetScratchMemTop ( sp );
} /*g2h_DrawSplBasFunctionf*/

/* ////////////////////////////////////////////////////////////////////////// */
void g2h_DrawSplBFAomcf ( GHoleDomainf *domain, int fn,
                          void (*drawpoly)(int degree, int lastknot,
                                           const float *knots,
                                           const float *coeff) )
{
  void  *sp;
  G2HoleSPrivateRecf *privateS;
  int   hole_k, i, lkn;
  float *c00, *c01, *c02, *d00, *d01, *d02, *c, *knb, *kns;

  if ( !(privateS = domain->SprivateG2) )
    return;
  lkn = privateS->lastomcknot;
  knb = &privateS->bezknots[G2H_FINALDEG-G2H_OMCDEG];
  kns = privateS->omcknots;
  sp = pkv_GetScratchMemTop ();
  if ( (c = pkv_GetScratchMemf ( lkn-G2H_OMCDEG )) ) {
    hole_k = domain->hole_k;
    for ( i = 0; i < hole_k; i++ ) {
      _g2h_GetBFAPatchCurvesf ( domain, fn, i, &c00, &c01, &c02,
                                &d00, &d01, &d02 );
      mbs_AdjustBSCRepC1f ( G2H_OMCDEG, 2*G2H_OMCDEG+1, knb, c00,
                            G2H_OMCDEG, lkn, kns, c );
      drawpoly ( G2H_OMCDEG, lkn, kns, c );
    }
  }
  pkv_SetScratchMemTop ( sp );
} /*g2h_DrawSplBFAomcf*/

void g2h_DrawSplBFBomcf ( GHoleDomainf *domain, int fn,
                          void (*drawpoly)(int degree, int lastknot,
                                           const float *knots,
                                           const float *coeff) )
{
  void  *sp;
  G2HoleSPrivateRecf *privateS;
  int   hole_k, i, lkn;
  float *c00, *c01, *c02, *c10, *c11, *c12,
        *d00, *d01, *d02, *d10, *d11, *d12;
  float *c, *knb, *kns;

  if ( !(privateS = domain->SprivateG2) )
    return;
  lkn = privateS->lastomcknot;
  knb = &privateS->bezknots[G2H_FINALDEG-G2H_OMCDEG];
  kns = privateS->omcknots;
  sp = pkv_GetScratchMemTop ();
  if ( (c = pkv_GetScratchMemf ( lkn-G2H_OMCDEG )) ) {
    hole_k = domain->hole_k;
    for ( i = 0; i < hole_k; i++ ) {
      _g2h_GetBFBPatchCurvesf ( domain, fn, i,
                                &c00, &c01, &c02, &c10, &c11, &c12,
                                &d00, &d01, &d02, &d10, &d11, &d12 );
      mbs_AdjustBSCRepC1f ( G2H_OMCDEG, 2*G2H_OMCDEG+1, knb, c00,
                            G2H_OMCDEG, lkn, kns, c );
      drawpoly ( G2H_OMCDEG, lkn, kns, c );
    }
  }
  pkv_SetScratchMemTop ( sp );
} /*g2h_DrawSplBFBomcf*/

void g2h_DrawSplBFDomcf ( GHoleDomainf *domain, int fn,
                          void (*drawpoly)(int degree, int lastknot,
                                           const float *knots,
                                           const float *coeff) )
{
  void  *sp;
  G2HolePrivateRecf *privateG2;
  G2HoleSPrivateRecf *privateS;
  int   hole_k, i, lkn, lpvkn, lpvvkn;
  float *kns, *fcomc, *pv, *pvv, *pu, *puu;

  privateG2 = domain->privateG2;
  if ( !(privateS = domain->SprivateG2) )
    return;

  lkn = privateS->lastomcknot;
  lpvkn = privateS->lastpvknot;
  lpvvkn = privateS->lastpvvknot;
  kns = privateS->omcknots;
  sp = pkv_GetScratchMemTop ();
  if ( (fcomc = pkv_GetScratchMemf (lkn-G2H_OMCDEG+
                    2*(lpvkn-G2_CROSS01DEG+lpvvkn-G2_CROSS02DEG))) ) {
    pv = &fcomc[lkn-G2H_OMCDEG];
    pu = &pv[lpvkn-G2_CROSS01DEG];
    pvv = &pu[lpvkn-G2_CROSS01DEG];
    puu = &pvv[lpvvkn-G2_CROSS02DEG];
    hole_k = domain->hole_k;
    fn += privateG2->nfunc_a+privateG2->nfunc_b+privateS->nsfunc_c;
    for ( i = 0; i < hole_k; i++ ) {
      _g2h_GetSplDBasisCrossDerf ( domain, fn, i, fcomc, pv, pvv, pu, puu );
      drawpoly ( G2H_OMCDEG, lkn, kns, fcomc );
    }
  }
  pkv_SetScratchMemTop ( sp );
} /*g2h_DrawSplBFDomcf*/

void g2h_DrawSplFinalSurfBCf ( GHoleDomainf *domain,
                               int spdimen, const float *hole_cp,
                               const float *acoeff,
                               void (*drawcurve)(int degree, int spdimen,
                                             int lastknot, const float *knots,
                                             const float *cp) )
{
  void *sp;
  int  hole_k, nfunc_a, nfunc_b, nfunc_c, nfunc_d;
  GHolePrivateRecf  *privateG;
  G2HolePrivateRecf *privateG2;
  G2HoleSPrivateRecf *privateS;
  unsigned char *bfcpn;
  float *x, *y, *c;
  float *c00, *c01, *c02, *c10, *c11, *c12, *d00, *d01, *d02, *d10, *d11, *d12;
  int   i, j, lkn, lpvvkn;
  float *knb, *kns;

  hole_k = domain->hole_k;
  privateG = domain->privateG;
  privateG2 = domain->privateG2;
  privateS = domain->SprivateG2;
  if ( !privateS )
    return;
  sp = pkv_GetScratchMemTop ();
  nfunc_a = privateG2->nfunc_a;
  nfunc_b = privateG2->nfunc_b;
  nfunc_c = privateS->nsfunc_c;
  nfunc_d = privateS->nsfunc_d;
  bfcpn = privateG->bfcpn;
  lkn = privateS->lastomcknot;
  lpvvkn = privateS->lastpvvknot;
  knb = &privateS->bezknots[G2H_FINALDEG-G2H_OMCDEG];
  kns = privateS->omcknots;
  if ( (x = pkv_GetScratchMemf ( 3*spdimen*(lpvvkn-G2H_OMCDEG) )) ) {
    y = &x[spdimen*(lpvvkn-G2H_OMCDEG)];
    c = &y[spdimen*(lpvvkn-G2H_OMCDEG)];
    for ( i = 0; i < hole_k; i++ ) {
      memset ( x, 0, spdimen*(G2H_OMCDEG+1)*sizeof(float) );
      for ( j = 0; j < nfunc_a; j++ ) {
        _g2h_GetBFAPatchCurvesf ( domain, j, i,
                                  &c00, &c01, &c02, &d00, &d01, &d02 );
        mbs_AdjustBSCRepC1f ( G2H_OMCDEG, 2*G2H_OMCDEG+1, knb, c00,
                              G2H_OMCDEG, lkn, kns, c );
        pkn_MultMatrixf ( lkn-G2H_OMCDEG, 1, 1, c, spdimen, 0,
                          &acoeff[(nfunc_c+nfunc_d+j)*spdimen], spdimen, y );
        pkn_AddMatrixf ( 1, spdimen*(lkn-G2H_OMCDEG), 0, x, 0, y, 0, x );
      }
      for ( j = 0; j < nfunc_d; j++ ) {
        _g2h_GetSplDBasisCrossDerf ( domain, j+nfunc_a+nfunc_b+nfunc_c, i,
                                     y, c, c, c, c );
        pkn_MultMatrixf ( lkn-G2H_OMCDEG, 1, 1, y, spdimen, 0,
                          &acoeff[(nfunc_c+j)*spdimen], spdimen, c );
        pkn_AddMatrixf ( 1, spdimen*(lkn-G2H_OMCDEG), 0, x, 0, c, 0, x );
      }
      for ( j = 0; j < nfunc_b; j++ ) {
        _g2h_GetBFBPatchCurvesf ( domain, j, i,
                                  &c00, &c01, &c02, &c10, &c11, &c12,
                                  &d00, &d01, &d02, &d10, &d11, &d12 );
        mbs_AdjustBSCRepC1f ( G2H_OMCDEG, 2*G2H_OMCDEG+1, knb, c00,
                              G2H_OMCDEG, lkn, kns, c );
        pkn_MultMatrixf ( lkn-G2H_OMCDEG, 1, 1, c, spdimen, 0,
                          &hole_cp[bfcpn[j]], spdimen, y );
        pkn_AddMatrixf ( 1, spdimen*(lkn-G2H_OMCDEG), 0, x, 0, y, 0, x );
      }
      drawcurve ( G2H_OMCDEG, spdimen, lkn, kns, x );
    }
  }
  pkv_SetScratchMemTop ( sp );
} /*g2h_DrawSplFinalSurfBCf*/

/* ///////////////////////////////////////////////////////////////////////// */
void g2h_DrawSplMatricesf ( GHoleDomainf *domain,
                            void (*drawmatrix)( int k, int r, int s, int t,
                            float *A, float *B ) )
{
  void *sp;
  G2HolePrivateRecf  *privateG2;
  G2HoleSPrivateRecf *sprivate;
  float *A, *B;
  int hole_k, nfunc_a, nfunc_b, nfunc_c, nfunc_d, s1, s2;

  sp = pkv_GetScratchMemTop ();
  privateG2 = domain->privateG2;
  sprivate = domain->SprivateG2;
  if ( !sprivate )
    goto finish;
  if ( !sprivate->SAMat || !sprivate->SBMat )
    goto finish;
  hole_k = domain->hole_k;
  nfunc_a = privateG2->nfunc_a;
  nfunc_b = privateG2->nfunc_b;
  nfunc_c = sprivate->nsfunc_c;
  nfunc_d = sprivate->nsfunc_d;
  s1 = pkn_Block2ArraySize ( hole_k, nfunc_c/hole_k, nfunc_d/hole_k, nfunc_a );
  s2 = nfunc_b*(nfunc_a+nfunc_c+nfunc_d);
  A = pkv_GetScratchMemf ( s1 );
  B = pkv_GetScratchMemf ( s2 );
  if ( !A || !B )
    goto finish;
  memcpy ( A, sprivate->SAMat, s1*sizeof(float) );
  memcpy ( B, sprivate->SBMat, s2*sizeof(float) );
  drawmatrix ( hole_k, nfunc_c/hole_k, nfunc_d/hole_k, nfunc_a, A, B );

finish:
  pkv_SetScratchMemTop ( sp );
} /*g2h_DrawSplMatricesf*/

