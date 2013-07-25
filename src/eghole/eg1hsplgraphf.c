
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2008, 2009                            */
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
#include "eg1holef.h"

#include "eg1hprivatef.h"
#include "eg1herror.h"

/* ////////////////////////////////////////////////////////////////////////// */
void g1h_DrawSplBasFuncNumf ( GHoleDomainf *domain,
                        int *nfunc_a, int *nfunc_b, int *nfunc_c, int *nfunc_d )
{
  G1HolePrivateRecf  *privateG1;
  G1HoleSPrivateRecf *sprivate;

  privateG1 = domain->privateG1;
  sprivate = domain->SprivateG1;
  *nfunc_a = privateG1->nfunc_a;
  *nfunc_b = privateG1->nfunc_b;
  *nfunc_c = sprivate->nsfunc_c;
  *nfunc_d = sprivate->nsfunc_d;
} /*g1h_DrawSplBasFuncNumf*/

void g1h_DrawSplBasAuxPatchesf ( GHoleDomainf *domain, int fn,
               void (*drawpatch) ( int n, int lknu, const float *knu,
                                   int m, int lknv, const float *knv,
                                   const point3f *cp ) )
{
  void     *sp;
  G1HolePrivateRecf  *privateG1;
  G1HoleSPrivateRecf *sprivate;
  int      i, hole_k, degu, lknu, lkn, nrows, nzc;
  int      nfunc_a, nfunc_b, nfunc_c, nfunc_d;
  vector2f *auxp, *omc, *omcd;
  float    *bezknots, *omcknots, *auxpknots;
  float    *fcomc, *fcomcd, *fcaux;
  point3f  *basp;

  sp = pkv_GetScratchMemTop ();
  hole_k = domain->hole_k;
  privateG1 = domain->privateG1;
  sprivate = domain->SprivateG1;
  if ( !sprivate )
    return;
  nfunc_a = privateG1->nfunc_a;
  nfunc_b = privateG1->nfunc_b;
  nfunc_c = sprivate->nsfunc_c;
  nfunc_d = sprivate->nsfunc_d;
  if ( fn < 0 || fn >= nfunc_a+nfunc_b+nfunc_c+nfunc_d )
    goto finish;
  bezknots = sprivate->bezknots;
  omcknots = sprivate->omcknots;
  lkn = sprivate->lastomcknot;
  omc = privateG1->omc;
  omcd = &omc[(G1H_OMCDEG+1)*hole_k];

      /* determine the "u" knot sequence for the auxiliary basis function patch */
  degu = 0;
  if ( !mbs_FindBSCommonKnotSequencef ( &degu, &lknu, &auxpknots, 2,
                                        G1H_OMCDEG, lkn, omcknots,
                                        G1H_OMCDEG-1, lkn-2, &omcknots[1] ) )
    goto finish;
      /* now the knot sequence is in the array on the scratch memory stack */
  nrows = lknu-degu;
  auxp = pkv_GetScratchMem ( 2*nrows*sizeof(point2f) );
  basp = pkv_GetScratchMem ( 2*nrows*sizeof(point3f) );
  fcomc = pkv_GetScratchMemf ( 3*nrows );
  if ( !auxp || !basp || !fcomc )
    goto finish;
  fcomcd = &fcomc[nrows];
  fcaux = &fcomcd[nrows];

  memset ( basp, 0, 2*nrows*sizeof(point3f) );
  for ( i = 0; i < hole_k; i++ ) {
        /* construct the auxiliary domain patch B-spline representation */
    mbs_AdjustBSCRepC2f ( G1H_OMCDEG, 2*G1H_OMCDEG+1,
              &bezknots[G1H_FINALDEG-G1H_OMCDEG], omc,
              degu, lknu, auxpknots, auxp );
    mbs_AdjustBSCRepC2f ( G1H_OMCDEG-1, 2*G1H_OMCDEG-1,
              &bezknots[G1H_FINALDEG-G1H_OMCDEG+1], omcd,
              degu, lknu, auxpknots, &auxp[nrows] );
    pkv_Selectf ( nrows, 2, 2, 3, auxp, basp );
    pkn_AddMatrixf ( nrows, 2, 3, &basp[0].x, 2, &auxp[nrows].x,
                      3, &basp[nrows].x );

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
      _g1h_GetSplDBasisAuxpf ( domain, fn, i, &nzc, fcomc, fcomcd );
      mbs_AdjustBSCRepC1f ( G1H_OMCDEG, lkn, omcknots, fcomc,
                            degu, lknu, auxpknots, fcaux );
      pkv_Selectf ( nrows, 1, 1, 3, fcaux, &basp[0].z );
      mbs_AdjustBSCRepC1f ( G1H_OMCDEG-1, lkn-2, &omcknots[1], fcomcd,
                            degu, lknu, auxpknots, fcaux );
      pkn_AddMatrixf ( nrows, 1, 3, &basp[0].z, 1, fcaux,
                        3, &basp[nrows].z );
    }
        /* output it */
    drawpatch ( 1, 3, &bezknots[G1H_FINALDEG-1], degu, lknu, auxpknots, basp );

    omc += (G1H_OMCDEG+1);
    omcd += G1H_OMCDEG;
  }

finish:
  pkv_SetScratchMemTop ( sp );
} /*g1h_DrawBasAuxPatchesf*/

void g1h_DrawSplBasFunctionf ( GHoleDomainf *domain, int fn,
             void (*drawpatch) ( int n, int lknu, const float *knu,
                                 int m, int lknv, const float *knv,
                                 const point3f *cp ) )
{
  void     *sp;
  G1HolePrivateRecf  *privateG1;
  G1HoleSPrivateRecf *sprivate;
  int      hole_k, nfunc_a, nfunc_b, nfunc_c, nfunc_d, bfn;
  int      i, j, n, m, k0, k1, k2;
  vector2f *c00, *c01, *c10, *c11, *d00, *d01, *d10, *d11;
  float    *fc00, *fc01, *fc10, *fc11,
           *fd00, *fd01, *fd10, *fd11;
  point2f  *di, *aux;
  point3f  *cp;
  int      degu, degv, lastomcknot, lastpvknot,
           lastcknot, lastfpknot, lastukn, lastvkn, ncp;
  float    *bezknots, *omcknots, *pvknots, *cknots;
  float    *fcomc, *pv, *pu, *auxbfp, *bfp, *uknots, *vknots;
  float    zero[2] = {0.0,0.0};

  sp = pkv_GetScratchMemTop ();
  hole_k = domain->hole_k;
  privateG1 = domain->privateG1;
  sprivate = domain->SprivateG1;
  if ( !sprivate )
    return;
  nfunc_a = privateG1->nfunc_a;
  nfunc_b = privateG1->nfunc_b;
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
  lastfpknot = sprivate->lastfpknot;
  
  n = lastfpknot-G1H_FINALDEG;
  m = G1H_FINALDEG+1;
  di = (point2f*)pkv_GetScratchMem ( m*m*sizeof(point2f) );
  aux = (point2f*)pkv_GetScratchMem ( n*m*sizeof(point2f) );
  cp = (point3f*)pkv_GetScratchMem ( n*n*sizeof(point3f) );
  auxbfp = pkv_GetScratchMemf ( 6*(lastpvknot-G1_CROSS01DEG) );
  bfp = pkv_GetScratchMemf ( (lastpvknot-G1_CROSS01DEG)*(lastpvknot-G1_CROSS01DEG) );
  fcomc = pkv_GetScratchMemf ( lastomcknot+2*lastpvknot-
                               (G1_CROSS00DEG+2*G1_CROSS01DEG) );
  uknots = pkv_GetScratchMemf ( 2*(lastfpknot+1) );
  if ( !di || !aux || !cp || !auxbfp || !bfp || !fcomc || !uknots )
    goto finish;
  pv = &fcomc[lastomcknot-G1_CROSS00DEG];
  pu = &pv[lastpvknot-G1_CROSS01DEG];
  vknots = &uknots[lastfpknot+1];

  for ( i = 0; i < hole_k; i++ ) {
        /* get the domain patch */
    _g1h_GetDiPatchCurvesf ( domain, i,
                             &c00, &c01, &c10, &c11, &d00, &d01, &d10, &d11 );
    mbs_BezC1CoonsToBezf ( 2,
       G1_CROSS00DEG, (float*)c00, G1_CROSS01DEG, (float*)c01,
       3, (float*)c10, G1_CROSS11DEG, (float*)c11,
       G1_CROSS00DEG, (float*)d00, G1_CROSS01DEG, (float*)d01,
       3, (float*)d10, G1_CROSS11DEG, (float*)d11,
       &degu, &degv, (float*)di );

        /* get the basis function patch */
    if ( fn < nfunc_a ) {                       /* block A */
      _g1h_GetBFAPatchCurvesf ( domain, fn, i,
                                &fc00, &fc01, &fd00, &fd01 );
      mbs_BezC1CoonsToBezf ( 1,
                    G1_CROSS00DEG, fc00, G1_CROSS01DEG, fc01, 1, zero, 1, zero,
                    G1_CROSS00DEG, fd00, G1_CROSS01DEG, fd01, 1, zero, 1, zero,
                    &degu, &degv, bfp );
      ncp = m*m;
      pkv_Selectf ( ncp, 2, 2, 3, (float*)di, cp );
      pkv_Selectf ( ncp, 1, 1, 3, bfp, &cp[0].z );
      drawpatch ( G1H_FINALDEG, 2*G1H_FINALDEG+1, bezknots,
                  G1H_FINALDEG, 2*G1H_FINALDEG+1, bezknots, cp );
    }
    else if ( fn < nfunc_a+nfunc_b ) {          /* block B */
      _g1h_GetBFBPatchCurvesf ( domain, fn-nfunc_a, i,
                                &fc00, &fc01, &fc10, &fc11,
                                &fd00, &fd01, &fd10, &fd11 );
      mbs_BezC1CoonsToBezf ( 1,
                    G1_CROSS00DEG, fc00, G1_CROSS01DEG, fc01,
                    G1_CROSS10DEG, fc10, G1_CROSS11DEG, fc11,
                    G1_CROSS00DEG, fd00, G1_CROSS01DEG, fd01,
                    G1_CROSS10DEG, fd10, G1_CROSS11DEG, fd11,
                    &degu, &degv, bfp );
      ncp = m*m;
      pkv_Selectf ( ncp, 2, 2, 3, (float*)di, cp );
      pkv_Selectf ( ncp, 1, 1, 3, bfp, &cp[0].z );
      drawpatch ( G1H_FINALDEG, 2*G1H_FINALDEG+1, bezknots,
                  G1H_FINALDEG, 2*G1H_FINALDEG+1, bezknots, cp );
    }
    else if ( fn < nfunc_a+nfunc_b+nfunc_c ) {  /* block C */
      n = lastcknot-G1H_FINALDEG;
      m = G1H_FINALDEG+1;
      mbs_multiAdjustBSCRepf ( 1, 2*(degv+1), degu, 2*degu+1,
                               &bezknots[G1H_FINALDEG-degu], 0, (float*)di,
                               G1H_FINALDEG, lastcknot, cknots, 0, (float*)aux );
      mbs_multiAdjustBSCRepf ( n, 2, degv, 2*degv+1,
                               &bezknots[G1H_FINALDEG-degv], 2*(degv+1), (float*)aux,
                               G1H_FINALDEG, lastcknot, cknots, 2*n, (float*)cp );
      pkv_Rearrangef ( n*n, 2, 2, 3, cp );
      for ( j = 0; j < n*n; j++ )
        cp[j].z = 0.0;

      bfn = fn - (nfunc_a+nfunc_b);
      k0 = bfn / sprivate->csize;  /* which Omega_i? */
      if ( k0 == i ) {
        bfn = bfn % sprivate->csize;  /* which function in the block? */
        k1  = bfn / (lastcknot-G1H_FINALDEG-4);
        k2  = bfn % (lastcknot-G1H_FINALDEG-4);
        cp[(k1+2)*n+k2+2].z = 1.0;
      }

      drawpatch ( G1H_FINALDEG, lastcknot, cknots,
                  G1H_FINALDEG, lastcknot, cknots, cp );
    }
    else {                                      /* block D */
      n = lastpvknot-G1H_FINALDEG;
      m = G1H_FINALDEG+1;
      mbs_multiAdjustBSCRepf ( 1, 2*(degv+1), degu, 2*degu+1,
                       &bezknots[G1H_FINALDEG-degu], 0, (float*)di,
                       G1H_FINALDEG, lastpvknot, pvknots, 0, (float*)aux );
      mbs_multiAdjustBSCRepf ( n, 2, degv, 2*degv+1,
                       &bezknots[G1H_FINALDEG-degv], 2*(degv+1), (float*)aux,
                       G1H_FINALDEG, lastpvknot, pvknots, 2*n, (float*)cp );
      pkv_Rearrangef ( n*n, 2, 2, 3, cp );
      for ( j = 0; j < n*n; j++ )
        cp[j].z = 0.0;

      bfn = fn - (nfunc_a+nfunc_b+nfunc_c);
      k0 = bfn / sprivate->dsize;
      if ( k0 == i ) {
        _g1h_GetSplDBasisCrossDerf ( domain, fn, i, fcomc, pv, pu );
        mbs_BSC1CoonsToBSf ( 1, G1_CROSS00DEG, lastomcknot, omcknots, fcomc,
            G1_CROSS01DEG, lastpvknot, pvknots, pv,
            1, 3, &bezknots[G1H_FINALDEG-1], zero,
            1, 3, &bezknots[G1H_FINALDEG-1], zero,
            1, 3, &bezknots[G1H_FINALDEG-1], zero,
            1, 3, &bezknots[G1H_FINALDEG-1], zero,
            1, 3, &bezknots[G1H_FINALDEG-1], zero,
            1, 3, &bezknots[G1H_FINALDEG-1], zero,
            &degu, &lastukn, uknots, &degv, &lastvkn, vknots, auxbfp );
        mbs_multiAdjustBSCRepf ( n, 1,
                      degv, lastvkn, vknots, lastvkn-degv, auxbfp,
                      G1H_FINALDEG, lastpvknot, pvknots, n, bfp );
        pkv_Selectf ( n*n, 1, 1, 3, bfp, &cp[0].z );
      }
      else if ( k0 == (i+1) % hole_k ) {
        _g1h_GetSplDBasisCrossDerf ( domain, fn, k0, fcomc, pv, pu );
        mbs_BSC1CoonsToBSf ( 1,
            1, 3, &bezknots[G1H_FINALDEG-1], zero,
            1, 3, &bezknots[G1H_FINALDEG-1], zero,
            1, 3, &bezknots[G1H_FINALDEG-1], zero,
            1, 3, &bezknots[G1H_FINALDEG-1], zero,
            G1_CROSS00DEG, lastomcknot, omcknots, fcomc,
            G1_CROSS01DEG, lastpvknot, pvknots, pu,
            1, 3, &bezknots[G1H_FINALDEG-1], zero,
            1, 3, &bezknots[G1H_FINALDEG-1], zero,
            &degu, &lastukn, uknots, &degv, &lastvkn, vknots, auxbfp );
        mbs_multiAdjustBSCRepf ( 1, n,
                                 degu, lastukn, uknots, 0, auxbfp,
                                 G1H_FINALDEG, lastpvknot, pvknots, n, bfp );
        pkv_Selectf ( n*n, 1, 1, 3, bfp, &cp[0].z );
      }

      drawpatch ( G1H_FINALDEG, lastpvknot, pvknots,
                  G1H_FINALDEG, lastpvknot, pvknots, cp );
    }
  }

finish:
  pkv_SetScratchMemTop ( sp );
} /*g1h_DrawSplBasFunctionf*/

/* ////////////////////////////////////////////////////////////////////////// */
void g1h_DrawSplBFAomcf ( GHoleDomainf *domain, int fn,
                          void (*drawpoly)(int degree, int lastknot,
                                           const float *knots,
                                           const float *coeff) )
{
  void  *sp;
  G1HoleSPrivateRecf *privateS;
  int   hole_k, i, lkn;
  float *c00, *c01, *d00, *d01, *c, *knb, *kns;

  if ( !(privateS = domain->SprivateG1) )
    return;
  lkn = privateS->lastomcknot;
  knb = &privateS->bezknots[G1H_FINALDEG-G1H_OMCDEG];
  kns = privateS->omcknots;
  sp = pkv_GetScratchMemTop ();
  if ( (c = pkv_GetScratchMemf ( lkn-G1H_OMCDEG )) ) {
    hole_k = domain->hole_k;
    for ( i = 0; i < hole_k; i++ ) {
      _g1h_GetBFAPatchCurvesf ( domain, fn, i, &c00, &c01, &d00, &d01 );
      mbs_AdjustBSCRepC1f ( G1H_OMCDEG, 2*G1H_OMCDEG+1, knb, c00,
                            G1H_OMCDEG, lkn, kns, c );
      drawpoly ( G1H_OMCDEG, lkn, kns, c );
    }
  }
  pkv_SetScratchMemTop ( sp );
} /*g1h_DrawSplBFAomcf*/

void g1h_DrawSplBFBomcf ( GHoleDomainf *domain, int fn,
                          void (*drawpoly)(int degree, int lastknot,
                                           const float *knots,
                                           const float *coeff) )
{
  void  *sp;
  G1HoleSPrivateRecf *privateS;
  int   hole_k, i, lkn;
  float *c00, *c01, *c10, *c11, *d00, *d01, *d10, *d11, *c, *knb, *kns;

  if ( !(privateS = domain->SprivateG1) )
    return;
  lkn = privateS->lastomcknot;
  knb = &privateS->bezknots[G1H_FINALDEG-G1H_OMCDEG];
  kns = privateS->omcknots;
  sp = pkv_GetScratchMemTop ();
  if ( (c = pkv_GetScratchMemf ( lkn-G1H_OMCDEG )) ) {
    hole_k = domain->hole_k;
    for ( i = 0; i < hole_k; i++ ) {
      _g1h_GetBFBPatchCurvesf ( domain, fn, i, &c00, &c01, &c10, &c11,
                                &d00, &d01, &d10, &d11 );
      mbs_AdjustBSCRepC1f ( G1H_OMCDEG, 2*G1H_OMCDEG+1, knb, c00,
                            G1H_OMCDEG, lkn, kns, c );
      drawpoly ( G1H_OMCDEG, lkn, kns, c );
    }
  }
  pkv_SetScratchMemTop ( sp );
} /*g1h_DrawSplBFBomcf*/

void g1h_DrawSplBFDomcf ( GHoleDomainf *domain, int fn,
                          void (*drawpoly)(int degree, int lastknot,
                                           const float *knots,
                                           const float *coeff) )
{
  void  *sp;
  G1HolePrivateRecf *privateG1;
  G1HoleSPrivateRecf *privateS;
  int   hole_k, i, lkn, lpvkn;
  float *kns, *fcomc, *pv, *pu;

  privateG1 = domain->privateG1;
  if ( !(privateS = domain->SprivateG1) )
    return;

  lkn = privateS->lastomcknot;
  lpvkn = privateS->lastpvknot;
  kns = privateS->omcknots;
  sp = pkv_GetScratchMemTop ();
  if ( (fcomc = pkv_GetScratchMemf (lkn-G1H_OMCDEG+2*(lpvkn-G1_CROSS01DEG))) ) {
    pv = &fcomc[lkn-G1H_OMCDEG];
    pu = &pv[lpvkn-G1_CROSS01DEG];
    hole_k = domain->hole_k;
    fn += privateG1->nfunc_a+privateG1->nfunc_b+privateS->nsfunc_c;
    for ( i = 0; i < hole_k; i++ ) {
      _g1h_GetSplDBasisCrossDerf ( domain, fn, i, fcomc, pv, pu );
      drawpoly ( G1H_OMCDEG, lkn, kns, fcomc );
    }
  }
  pkv_SetScratchMemTop ( sp );
} /*g1h_DrawSplBFDomcf*/

void g1h_DrawSplFinalSurfBCf ( GHoleDomainf *domain,
                               int spdimen, const float *hole_cp,
                               const float *acoeff,
                               void (*drawcurve)(int degree, int spdimen,
                                             int lastknot, const float *knots,
                                             const float *cp) )
{
  void *sp;
  int  hole_k, nfunc_a, nfunc_b, nfunc_c, nfunc_d;
  GHolePrivateRecf  *privateG;
  G1HolePrivateRecf *privateG1;
  G1HoleSPrivateRecf *privateS;
  unsigned char *bfcpn;
  float *x, *y, *c;
  float *c00, *c01, *c10, *c11, *d00, *d01, *d10, *d11;
  int   i, j, lkn, lpvkn;
  float *knb, *kns;

  hole_k = domain->hole_k;
  privateG = domain->privateG;
  privateG1 = domain->privateG1;
  privateS = domain->SprivateG1;
  if ( !privateS )
    return;
  sp = pkv_GetScratchMemTop ();
  nfunc_a = privateG1->nfunc_a;
  nfunc_b = privateG1->nfunc_b;
  nfunc_c = privateS->nsfunc_c;
  nfunc_d = privateS->nsfunc_d;
  bfcpn = privateG->bfcpn;
  lkn = privateS->lastomcknot;
  lpvkn = privateS->lastpvknot;
  knb = &privateS->bezknots[G1H_FINALDEG-G1H_OMCDEG];
  kns = privateS->omcknots;
  if ( (x = pkv_GetScratchMemf ( 3*spdimen*(lpvkn-G1H_OMCDEG) )) ) {
    y = &x[spdimen*(lpvkn-G1H_OMCDEG)];
    c = &y[spdimen*(lpvkn-G1H_OMCDEG)];
    for ( i = 0; i < hole_k; i++ ) {
      memset ( x, 0, spdimen*(G1H_OMCDEG+1)*sizeof(float) );
      for ( j = 0; j < nfunc_a; j++ ) {
        _g1h_GetBFAPatchCurvesf ( domain, j, i, &c00, &c01, &d00, &d01 );
        mbs_AdjustBSCRepC1f ( G1H_OMCDEG, 2*G1H_OMCDEG+1, knb, c00,
                              G1H_OMCDEG, lkn, kns, c );
        pkn_MultMatrixf ( lkn-G1H_OMCDEG, 1, 1, c, spdimen, 0,
                          &acoeff[(nfunc_c+nfunc_d+j)*spdimen], spdimen, y );
        pkn_AddMatrixf ( 1, spdimen*(lkn-G1H_OMCDEG), 0, x, 0, y, 0, x );
      }
      for ( j = 0; j < nfunc_d; j++ ) {
        _g1h_GetSplDBasisCrossDerf ( domain, j+nfunc_a+nfunc_b+nfunc_c, i,
                                     y, c, c );
        pkn_MultMatrixf ( lkn-G1H_OMCDEG, 1, 1, y, spdimen, 0,
                          &acoeff[(nfunc_c+j)*spdimen], spdimen, c );
        pkn_AddMatrixf ( 1, spdimen*(lkn-G1H_OMCDEG), 0, x, 0, c, 0, x );
      }
      for ( j = 0; j < nfunc_b; j++ ) {
        _g1h_GetBFBPatchCurvesf ( domain, j, i, &c00, &c01, &c10, &c11,
                                  &d00, &d01, &d10, &d11 );
        mbs_AdjustBSCRepC1f ( G1H_OMCDEG, 2*G1H_OMCDEG+1, knb, c00,
                              G1H_OMCDEG, lkn, kns, c );
        pkn_MultMatrixf ( lkn-G1H_OMCDEG, 1, 1, c, spdimen, 0,
                          &hole_cp[bfcpn[j]], spdimen, y );
        pkn_AddMatrixf ( 1, spdimen*(lkn-G1H_OMCDEG), 0, x, 0, y, 0, x );
      }
      drawcurve ( G1H_OMCDEG, spdimen, lkn, kns, x );
    }
  }
  pkv_SetScratchMemTop ( sp );
} /*g1h_DrawSplFinalSurfBCf*/

/* ///////////////////////////////////////////////////////////////////////// */
void g1h_DrawSplMatricesf ( GHoleDomainf *domain,
                            void (*drawmatrix)( int k, int r, int s, int t,
                            float *A, float *B ) )
{
  void *sp;
  G1HolePrivateRecf  *privateG1;
  G1HoleSPrivateRecf *sprivate;
  float *A, *B;
  int hole_k, nfunc_a, nfunc_b, nfunc_c, nfunc_d, s1, s2;

  sp = pkv_GetScratchMemTop ();
  privateG1 = domain->privateG1;
  sprivate = domain->SprivateG1;
  if ( !sprivate )
    goto finish;
  if ( !sprivate->SAMat || !sprivate->SBMat )
    goto finish;
  hole_k = domain->hole_k;
  nfunc_a = privateG1->nfunc_a;
  nfunc_b = privateG1->nfunc_b;
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
} /*g1h_DrawSplMatricesf*/

