
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
#include "eg2holed.h"

#include "eg2hprivated.h"
#include "eg2herror.h"

/* ////////////////////////////////////////////////////////////////////////// */
void g2h_DrawSplBasFuncNumd ( GHoleDomaind *domain,
                        int *nfunc_a, int *nfunc_b, int *nfunc_c, int *nfunc_d )
{
  G2HolePrivateRecd  *privateG2;
  G2HoleSPrivateRecd *sprivate;

  privateG2 = domain->privateG2;
  sprivate = domain->SprivateG2;
  *nfunc_a = privateG2->nfunc_a;
  *nfunc_b = privateG2->nfunc_b;
  *nfunc_c = sprivate->nsfunc_c;
  *nfunc_d = sprivate->nsfunc_d;
} /*g2h_DrawSplBasFuncNumd*/

void g2h_DrawSplBasAuxPatchesd ( GHoleDomaind *domain, int fn,
               void (*drawpatch) ( int n, int lknu, const double *knu,
                                   int m, int lknv, const double *knv,
                                   const point3d *cp ) )
{
  void     *sp;
  G2HolePrivateRecd  *privateG2;
  G2HoleSPrivateRecd *sprivate;
  int      i, hole_k, degu, lknu, lkn, nrows, nzc;
  int      nfunc_a, nfunc_b, nfunc_c, nfunc_d;
  vector2d *auxp, *omc, *omcd, *omcdd;
  double   *bezknots, *omcknots, *auxpknots;
  double   *fcomc, *fcomcd, *fcomcdd, *fcaux;
  point3d  *basp;

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
  if ( !mbs_FindBSCommonKnotSequenced ( &degu, &lknu, &auxpknots, 3,
                                        G2H_OMCDEG, lkn, omcknots,
                                        G2H_OMCDEG-1, lkn-2, &omcknots[1],
                                        G2H_OMCDEG-2, lkn-4, &omcknots[2] ) )
    goto finish;
      /* now the knot sequence is in the array on the scratch memory stack */
  nrows = lknu-degu;
  auxp = pkv_GetScratchMem ( 3*nrows*sizeof(point2d) );
  basp = pkv_GetScratchMem ( 3*nrows*sizeof(point3d) );
  fcomc = pkv_GetScratchMemd ( 4*nrows );
  if ( !auxp || !basp || !fcomc )
    goto finish;
  fcomcd = &fcomc[nrows];
  fcomcdd = &fcomcd[nrows];
  fcaux = &fcomcdd[nrows];

  memset ( basp, 0, 3*nrows*sizeof(point3d) );
  for ( i = 0; i < hole_k; i++ ) {
        /* construct the auxiliary domain patch B-spline representation */
    mbs_AdjustBSCRepC2d ( G2H_OMCDEG, 2*G2H_OMCDEG+1,
              &bezknots[G2H_FINALDEG-G2H_OMCDEG], omc,
              degu, lknu, auxpknots, auxp );
    mbs_AdjustBSCRepC2d ( G2H_OMCDEG-1, 2*G2H_OMCDEG-1,
              &bezknots[G2H_FINALDEG-G2H_OMCDEG+1], omcd,
              degu, lknu, auxpknots, &auxp[nrows] );
    mbs_AdjustBSCRepC2d ( G2H_OMCDEG-2, 2*G2H_OMCDEG-3,
              &bezknots[G2H_FINALDEG-G2H_OMCDEG+2], omcdd,
              degu, lknu, auxpknots, &auxp[2*nrows] );
    pkv_Selectd ( nrows, 2, 2, 3, auxp, basp );
    pkn_AddMatrixMd ( nrows, 2, 3, &basp[0].x, 2, &auxp[nrows].x, 0.5,
                      3, &basp[nrows].x );
    pkn_AddMatrixd ( nrows, 2, 3, &basp[0].x, 2, &auxp[nrows].x,
                     3, &basp[2*nrows].x );
    pkn_AddMatrixMd ( nrows, 2, 3, &basp[2*nrows].x, 2, &auxp[2*nrows].x, 0.5,
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
      _g2h_GetSplDBasisAuxpd ( domain, fn, i, &nzc, fcomc, fcomcd, fcomcdd );
      mbs_AdjustBSCRepC1d ( G2H_OMCDEG, lkn, omcknots, fcomc,
                            degu, lknu, auxpknots, fcaux );
      pkv_Selectd ( nrows, 1, 1, 3, fcaux, &basp[0].z );
      mbs_AdjustBSCRepC1d ( G2H_OMCDEG-1, lkn-2, &omcknots[1], fcomcd,
                            degu, lknu, auxpknots, fcaux );
      pkn_AddMatrixMd ( nrows, 1, 3, &basp[0].z, 1, fcaux, 0.5,
                        3, &basp[nrows].z );
      pkn_AddMatrixd ( nrows, 1, 3, &basp[0].z, 1, fcaux, 3, &basp[2*nrows].z );
      mbs_AdjustBSCRepC1d ( G2H_OMCDEG-2, lkn-4, &omcknots[2], fcomcdd,
                            degu, lknu, auxpknots, fcaux );
      pkn_AddMatrixMd ( nrows, 1, 3, &basp[2*nrows].z, 1, fcaux, 0.5,
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
} /*g2h_DrawBasAuxPatchesd*/

void g2h_DrawSplBasFunctiond ( GHoleDomaind *domain, int fn,
             void (*drawpatch) ( int n, int lknu, const double *knu,
                                 int m, int lknv, const double *knv,
                                 const point3d *cp ) )
{
  void     *sp;
  G2HolePrivateRecd  *privateG2;
  G2HoleSPrivateRecd *sprivate;
  int      hole_k, nfunc_a, nfunc_b, nfunc_c, nfunc_d, bfn;
  int      i, j, n, m, k0, k1, k2;
  vector2d *c00, *c01, *c02, *c10, *c11, *c12,
           *d00, *d01, *d02, *d10, *d11, *d12;
  double   *fc00, *fc01, *fc02, *fc10, *fc11, *fc12,
           *fd00, *fd01, *fd02, *fd10, *fd11, *fd12;
  point2d  *di, *aux;
  point3d  *cp;
  int      degu, degv, lastomcknot, lastpvknot, lastpvvknot,
           lastcknot, lastfpknot, lastukn, lastvkn, ncp;
  double   *bezknots, *omcknots, *pvknots, *pvvknots, *cknots;
  double   *fcomc, *pv, *pvv, *pu, *puu, *auxbfp, *bfp, *uknots, *vknots;
  double   zero[2] = {0.0,0.0};

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
  di = (point2d*)pkv_GetScratchMem ( m*m*sizeof(point2d) );
  aux = (point2d*)pkv_GetScratchMem ( n*m*sizeof(point2d) );
  cp = (point3d*)pkv_GetScratchMem ( n*n*sizeof(point3d) );
  auxbfp = pkv_GetScratchMemd ( 6*(lastpvvknot-G2_CROSS02DEG) );
  bfp = pkv_GetScratchMemd ( (lastpvvknot-G2_CROSS02DEG)*(lastpvvknot-G2_CROSS02DEG) );
  fcomc = pkv_GetScratchMemd ( lastomcknot+2*(lastpvknot+lastpvvknot)-
                               (G2_CROSS00DEG+2*(G2_CROSS01DEG+G2_CROSS02DEG)) );
  uknots = pkv_GetScratchMemd ( 2*(lastfpknot+1) );
  if ( !di || !aux || !cp || !auxbfp || !bfp || !fcomc || !uknots )
    goto finish;
  pv = &fcomc[lastomcknot-G2_CROSS00DEG];
  pvv = &pv[lastpvknot-G2_CROSS01DEG];
  pu = &pvv[lastpvvknot-G2_CROSS02DEG];
  puu = &pu[lastpvknot-G2_CROSS01DEG];
  vknots = &uknots[lastfpknot+1];

  for ( i = 0; i < hole_k; i++ ) {
        /* get the domain patch */
    _g2h_GetDiPatchCurvesd ( domain, i,
                             &c00, &c01, &c02, &c10, &c11, &c12,
                             &d00, &d01, &d02, &d10, &d11, &d12 );
    mbs_BezC2CoonsToBezd ( 2,
       G2_CROSS00DEG, (double*)c00, G2_CROSS01DEG, (double*)c01, G2_CROSS02DEG, (double*)c02,
       3, (double*)c10, G2_CROSS11DEG, (double*)c11, G2_CROSS12DEG, (double*)c12,
       G2_CROSS00DEG, (double*)d00, G2_CROSS01DEG, (double*)d01, G2_CROSS02DEG, (double*)d02,
       3, (double*)d10, G2_CROSS11DEG, (double*)d11, G2_CROSS12DEG, (double*)d12,
       &degu, &degv, (double*)di );

        /* get the basis function patch */
    if ( fn < nfunc_a ) {                       /* block A */
      _g2h_GetBFAPatchCurvesd ( domain, fn, i,
                                &fc00, &fc01, &fc02, &fd00, &fd01, &fd02 );
      mbs_BezC2CoonsToBezd ( 1,
                    G2_CROSS00DEG, fc00, G2_CROSS01DEG, fc01, G2_CROSS02DEG, fc02,
                    1, zero, 1, zero, 1, zero,
                    G2_CROSS00DEG, fd00, G2_CROSS01DEG, fd01, G2_CROSS02DEG, fd02,
                    1, zero, 1, zero, 1, zero,
                    &degu, &degv, bfp );
      ncp = m*m;
      pkv_Selectd ( ncp, 2, 2, 3, (double*)di, cp );
      pkv_Selectd ( ncp, 1, 1, 3, bfp, &cp[0].z );
      drawpatch ( G2H_FINALDEG, 2*G2H_FINALDEG+1, bezknots,
                  G2H_FINALDEG, 2*G2H_FINALDEG+1, bezknots, cp );
    }
    else if ( fn < nfunc_a+nfunc_b ) {          /* block B */
      _g2h_GetBFBPatchCurvesd ( domain, fn-nfunc_a, i,
                                &fc00, &fc01, &fc02, &fc10, &fc11, &fc12,
                                &fd00, &fd01, &fd02, &fd10, &fd11, &fd12 );
      mbs_BezC2CoonsToBezd ( 1,
                    G2_CROSS00DEG, fc00, G2_CROSS01DEG, fc01, G2_CROSS02DEG, fc02,
                    G2_CROSS10DEG, fc10, G2_CROSS11DEG, fc11, G2_CROSS12DEG, fc12,
                    G2_CROSS00DEG, fd00, G2_CROSS01DEG, fd01, G2_CROSS02DEG, fd02,
                    G2_CROSS10DEG, fd10, G2_CROSS11DEG, fd11, G2_CROSS12DEG, fd12,
                    &degu, &degv, bfp );
      ncp = m*m;
      pkv_Selectd ( ncp, 2, 2, 3, (double*)di, cp );
      pkv_Selectd ( ncp, 1, 1, 3, bfp, &cp[0].z );
      drawpatch ( G2H_FINALDEG, 2*G2H_FINALDEG+1, bezknots,
                  G2H_FINALDEG, 2*G2H_FINALDEG+1, bezknots, cp );
    }
    else if ( fn < nfunc_a+nfunc_b+nfunc_c ) {  /* block C */
      n = lastcknot-G2H_FINALDEG;
      m = G2H_FINALDEG+1;
      mbs_multiAdjustBSCRepd ( 1, 2*(degv+1), degu, 2*degu+1,
                               &bezknots[G2H_FINALDEG-degu], 0, (double*)di,
                               G2H_FINALDEG, lastcknot, cknots, 0, (double*)aux );
      mbs_multiAdjustBSCRepd ( n, 2, degv, 2*degv+1,
                               &bezknots[G2H_FINALDEG-degv], 2*(degv+1), (double*)aux,
                               G2H_FINALDEG, lastcknot, cknots, 2*n, (double*)cp );
      pkv_Rearranged ( n*n, 2, 2, 3, cp );
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
      mbs_multiAdjustBSCRepd ( 1, 2*(degv+1), degu, 2*degu+1,
                       &bezknots[G2H_FINALDEG-degu], 0, (double*)di,
                       G2H_FINALDEG, lastpvvknot, pvvknots, 0, (double*)aux );
      mbs_multiAdjustBSCRepd ( n, 2, degv, 2*degv+1,
                       &bezknots[G2H_FINALDEG-degv], 2*(degv+1), (double*)aux,
                       G2H_FINALDEG, lastpvvknot, pvvknots, 2*n, (double*)cp );
      pkv_Rearranged ( n*n, 2, 2, 3, cp );
      for ( j = 0; j < n*n; j++ )
        cp[j].z = 0.0;

      bfn = fn - (nfunc_a+nfunc_b+nfunc_c);
      k0 = bfn / sprivate->dsize;
      if ( k0 == i ) {
        _g2h_GetSplDBasisCrossDerd ( domain, fn, i, fcomc, pv, pvv, pu, puu );
        mbs_BSC2CoonsToBSd ( 1, G2_CROSS00DEG, lastomcknot, omcknots, fcomc,
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
        mbs_multiAdjustBSCRepd ( n, 1,
                      degv, lastvkn, vknots, lastvkn-degv, auxbfp,
                      G2H_FINALDEG, lastpvvknot, pvvknots, n, bfp );
        pkv_Selectd ( n*n, 1, 1, 3, bfp, &cp[0].z );
      }
      else if ( k0 == (i+1) % hole_k ) {
        _g2h_GetSplDBasisCrossDerd ( domain, fn, k0, fcomc, pv, pvv, pu, puu );
        mbs_BSC2CoonsToBSd ( 1,
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
        mbs_multiAdjustBSCRepd ( 1, n,
                                 degu, lastukn, uknots, 0, auxbfp,
                                 G2H_FINALDEG, lastpvvknot, pvvknots, n, bfp );
        pkv_Selectd ( n*n, 1, 1, 3, bfp, &cp[0].z );
      }

      drawpatch ( G2H_FINALDEG, lastpvvknot, pvvknots,
                  G2H_FINALDEG, lastpvvknot, pvvknots, cp );
    }
  }

finish:
  pkv_SetScratchMemTop ( sp );
} /*g2h_DrawSplBasFunctiond*/

/* ////////////////////////////////////////////////////////////////////////// */
void g2h_DrawSplBFAomcd ( GHoleDomaind *domain, int fn,
                          void (*drawpoly)(int degree, int lastknot,
                                           const double *knots,
                                           const double *coeff) )
{
  void  *sp;
  G2HoleSPrivateRecd *privateS;
  int   hole_k, i, lkn;
  double *c00, *c01, *c02, *d00, *d01, *d02, *c, *knb, *kns;

  if ( !(privateS = domain->SprivateG2) )
    return;
  lkn = privateS->lastomcknot;
  knb = &privateS->bezknots[G2H_FINALDEG-G2H_OMCDEG];
  kns = privateS->omcknots;
  sp = pkv_GetScratchMemTop ();
  if ( (c = pkv_GetScratchMemd ( lkn-G2H_OMCDEG )) ) {
    hole_k = domain->hole_k;
    for ( i = 0; i < hole_k; i++ ) {
      _g2h_GetBFAPatchCurvesd ( domain, fn, i, &c00, &c01, &c02,
                                &d00, &d01, &d02 );
      mbs_AdjustBSCRepC1d ( G2H_OMCDEG, 2*G2H_OMCDEG+1, knb, c00,
                            G2H_OMCDEG, lkn, kns, c );
      drawpoly ( G2H_OMCDEG, lkn, kns, c );
    }
  }
  pkv_SetScratchMemTop ( sp );
} /*g2h_DrawSplBFAomcd*/

void g2h_DrawSplBFBomcd ( GHoleDomaind *domain, int fn,
                          void (*drawpoly)(int degree, int lastknot,
                                           const double *knots,
                                           const double *coeff) )
{
  void  *sp;
  G2HoleSPrivateRecd *privateS;
  int   hole_k, i, lkn;
  double *c00, *c01, *c02, *c10, *c11, *c12,
        *d00, *d01, *d02, *d10, *d11, *d12;
  double *c, *knb, *kns;

  if ( !(privateS = domain->SprivateG2) )
    return;
  lkn = privateS->lastomcknot;
  knb = &privateS->bezknots[G2H_FINALDEG-G2H_OMCDEG];
  kns = privateS->omcknots;
  sp = pkv_GetScratchMemTop ();
  if ( (c = pkv_GetScratchMemd ( lkn-G2H_OMCDEG )) ) {
    hole_k = domain->hole_k;
    for ( i = 0; i < hole_k; i++ ) {
      _g2h_GetBFBPatchCurvesd ( domain, fn, i,
                                &c00, &c01, &c02, &c10, &c11, &c12,
                                &d00, &d01, &d02, &d10, &d11, &d12 );
      mbs_AdjustBSCRepC1d ( G2H_OMCDEG, 2*G2H_OMCDEG+1, knb, c00,
                            G2H_OMCDEG, lkn, kns, c );
      drawpoly ( G2H_OMCDEG, lkn, kns, c );
    }
  }
  pkv_SetScratchMemTop ( sp );
} /*g2h_DrawSplBFBomcd*/

void g2h_DrawSplBFDomcd ( GHoleDomaind *domain, int fn,
                          void (*drawpoly)(int degree, int lastknot,
                                           const double *knots,
                                           const double *coeff) )
{
  void  *sp;
  G2HolePrivateRecd *privateG2;
  G2HoleSPrivateRecd *privateS;
  int   hole_k, i, lkn, lpvkn, lpvvkn;
  double *kns, *fcomc, *pv, *pvv, *pu, *puu;

  privateG2 = domain->privateG2;
  if ( !(privateS = domain->SprivateG2) )
    return;

  lkn = privateS->lastomcknot;
  lpvkn = privateS->lastpvknot;
  lpvvkn = privateS->lastpvvknot;
  kns = privateS->omcknots;
  sp = pkv_GetScratchMemTop ();
  if ( (fcomc = pkv_GetScratchMemd (lkn-G2H_OMCDEG+
                    2*(lpvkn-G2_CROSS01DEG+lpvvkn-G2_CROSS02DEG))) ) {
    pv = &fcomc[lkn-G2H_OMCDEG];
    pu = &pv[lpvkn-G2_CROSS01DEG];
    pvv = &pu[lpvkn-G2_CROSS01DEG];
    puu = &pvv[lpvvkn-G2_CROSS02DEG];
    hole_k = domain->hole_k;
    fn += privateG2->nfunc_a+privateG2->nfunc_b+privateS->nsfunc_c;
    for ( i = 0; i < hole_k; i++ ) {
      _g2h_GetSplDBasisCrossDerd ( domain, fn, i, fcomc, pv, pvv, pu, puu );
      drawpoly ( G2H_OMCDEG, lkn, kns, fcomc );
    }
  }
  pkv_SetScratchMemTop ( sp );
} /*g2h_DrawSplBFDomcd*/

void g2h_DrawSplFinalSurfBCd ( GHoleDomaind *domain,
                               int spdimen, const double *hole_cp,
                               const double *acoeff,
                               void (*drawcurve)(int degree, int spdimen,
                                             int lastknot, const double *knots,
                                             const double *cp) )
{
  void *sp;
  int  hole_k, nfunc_a, nfunc_b, nfunc_c, nfunc_d;
  GHolePrivateRecd  *privateG;
  G2HolePrivateRecd *privateG2;
  G2HoleSPrivateRecd *privateS;
  unsigned char *bfcpn;
  double *x, *y, *c;
  double *c00, *c01, *c02, *c10, *c11, *c12, *d00, *d01, *d02, *d10, *d11, *d12;
  int   i, j, lkn, lpvvkn;
  double *knb, *kns;

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
  if ( (x = pkv_GetScratchMemd ( 3*spdimen*(lpvvkn-G2H_OMCDEG) )) ) {
    y = &x[spdimen*(lpvvkn-G2H_OMCDEG)];
    c = &y[spdimen*(lpvvkn-G2H_OMCDEG)];
    for ( i = 0; i < hole_k; i++ ) {
      memset ( x, 0, spdimen*(G2H_OMCDEG+1)*sizeof(double) );
      for ( j = 0; j < nfunc_a; j++ ) {
        _g2h_GetBFAPatchCurvesd ( domain, j, i,
                                  &c00, &c01, &c02, &d00, &d01, &d02 );
        mbs_AdjustBSCRepC1d ( G2H_OMCDEG, 2*G2H_OMCDEG+1, knb, c00,
                              G2H_OMCDEG, lkn, kns, c );
        pkn_MultMatrixd ( lkn-G2H_OMCDEG, 1, 1, c, spdimen, 0,
                          &acoeff[(nfunc_c+nfunc_d+j)*spdimen], spdimen, y );
        pkn_AddMatrixd ( 1, spdimen*(lkn-G2H_OMCDEG), 0, x, 0, y, 0, x );
      }
      for ( j = 0; j < nfunc_d; j++ ) {
        _g2h_GetSplDBasisCrossDerd ( domain, j+nfunc_a+nfunc_b+nfunc_c, i,
                                     y, c, c, c, c );
        pkn_MultMatrixd ( lkn-G2H_OMCDEG, 1, 1, y, spdimen, 0,
                          &acoeff[(nfunc_c+j)*spdimen], spdimen, c );
        pkn_AddMatrixd ( 1, spdimen*(lkn-G2H_OMCDEG), 0, x, 0, c, 0, x );
      }
      for ( j = 0; j < nfunc_b; j++ ) {
        _g2h_GetBFBPatchCurvesd ( domain, j, i,
                                  &c00, &c01, &c02, &c10, &c11, &c12,
                                  &d00, &d01, &d02, &d10, &d11, &d12 );
        mbs_AdjustBSCRepC1d ( G2H_OMCDEG, 2*G2H_OMCDEG+1, knb, c00,
                              G2H_OMCDEG, lkn, kns, c );
        pkn_MultMatrixd ( lkn-G2H_OMCDEG, 1, 1, c, spdimen, 0,
                          &hole_cp[bfcpn[j]], spdimen, y );
        pkn_AddMatrixd ( 1, spdimen*(lkn-G2H_OMCDEG), 0, x, 0, y, 0, x );
      }
      drawcurve ( G2H_OMCDEG, spdimen, lkn, kns, x );
    }
  }
  pkv_SetScratchMemTop ( sp );
} /*g2h_DrawSplFinalSurfBCd*/

/* ///////////////////////////////////////////////////////////////////////// */
void g2h_DrawSplMatricesd ( GHoleDomaind *domain,
                            void (*drawmatrix)( int k, int r, int s, int t,
                            double *A, double *B ) )
{
  void  *sp;
  G2HolePrivateRecd  *privateG2;
  G2HoleSPrivateRecd *sprivate;
  double *A, *B;
  int    hole_k, nfunc_a, nfunc_b, nfunc_c, nfunc_d, s1, s2;

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
  A = pkv_GetScratchMemd ( s1 );
  B = pkv_GetScratchMemd ( s2 );
  if ( !A || !B )
    goto finish;
  memcpy ( A, sprivate->SAMat, s1*sizeof(double) );
  memcpy ( B, sprivate->SBMat, s2*sizeof(double) );
  drawmatrix ( hole_k, nfunc_c/hole_k, nfunc_d/hole_k, nfunc_a, A, B );

finish:
  pkv_SetScratchMemTop ( sp );
} /*g2h_DrawSplMatricesd*/

