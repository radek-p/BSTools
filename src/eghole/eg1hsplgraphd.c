
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
#include "eg1holed.h"

#include "eg1hprivated.h"
#include "eg1herror.h"

/* ////////////////////////////////////////////////////////////////////////// */
void g1h_DrawSplBasFuncNumd ( GHoleDomaind *domain,
                        int *nfunc_a, int *nfunc_b, int *nfunc_c, int *nfunc_d )
{
  G1HolePrivateRecd  *privateG1;
  G1HoleSPrivateRecd *sprivate;

  privateG1 = domain->privateG1;
  sprivate = domain->SprivateG1;
  *nfunc_a = privateG1->nfunc_a;
  *nfunc_b = privateG1->nfunc_b;
  *nfunc_c = sprivate->nsfunc_c;
  *nfunc_d = sprivate->nsfunc_d;
} /*g1h_DrawSplBasFuncNumd*/

void g1h_DrawSplBasAuxPatchesd ( GHoleDomaind *domain, int fn,
               void (*drawpatch) ( int n, int lknu, const double *knu,
                                   int m, int lknv, const double *knv,
                                   const point3d *cp ) )
{
  void     *sp;
  G1HolePrivateRecd  *privateG1;
  G1HoleSPrivateRecd *sprivate;
  int      i, hole_k, degu, lknu, lkn, nrows, nzc;
  int      nfunc_a, nfunc_b, nfunc_c, nfunc_d;
  vector2d *auxp, *omc, *omcd;
  double   *bezknots, *omcknots, *auxpknots;
  double   *fcomc, *fcomcd, *fcaux;
  point3d  *basp;

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
  if ( !mbs_FindBSCommonKnotSequenced ( &degu, &lknu, &auxpknots, 2,
                                        G1H_OMCDEG, lkn, omcknots,
                                        G1H_OMCDEG-1, lkn-2, &omcknots[1] ) )
    goto finish;
      /* now the knot sequence is in the array on the scratch memory stack */
  nrows = lknu-degu;
  auxp = pkv_GetScratchMem ( 2*nrows*sizeof(point2d) );
  basp = pkv_GetScratchMem ( 2*nrows*sizeof(point3d) );
  fcomc = pkv_GetScratchMemd ( 3*nrows );
  if ( !auxp || !basp || !fcomc )
    goto finish;
  fcomcd = &fcomc[nrows];
  fcaux = &fcomcd[nrows];

  memset ( basp, 0, 2*nrows*sizeof(point3d) );
  for ( i = 0; i < hole_k; i++ ) {
        /* construct the auxiliary domain patch B-spline representation */
    mbs_AdjustBSCRepC2d ( G1H_OMCDEG, 2*G1H_OMCDEG+1,
              &bezknots[G1H_FINALDEG-G1H_OMCDEG], omc,
              degu, lknu, auxpknots, auxp );
    mbs_AdjustBSCRepC2d ( G1H_OMCDEG-1, 2*G1H_OMCDEG-1,
              &bezknots[G1H_FINALDEG-G1H_OMCDEG+1], omcd,
              degu, lknu, auxpknots, &auxp[nrows] );
    pkv_Selectd ( nrows, 2, 2, 3, auxp, basp );
    pkn_AddMatrixd ( nrows, 2, 3, &basp[0].x, 2, &auxp[nrows].x,
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
      _g1h_GetSplDBasisAuxpd ( domain, fn, i, &nzc, fcomc, fcomcd );
      mbs_AdjustBSCRepC1d ( G1H_OMCDEG, lkn, omcknots, fcomc,
                            degu, lknu, auxpknots, fcaux );
      pkv_Selectd ( nrows, 1, 1, 3, fcaux, &basp[0].z );
      mbs_AdjustBSCRepC1d ( G1H_OMCDEG-1, lkn-2, &omcknots[1], fcomcd,
                            degu, lknu, auxpknots, fcaux );
      pkn_AddMatrixd ( nrows, 1, 3, &basp[0].z, 1, fcaux,
                        3, &basp[nrows].z );
    }
        /* output it */
    drawpatch ( 1, 3, &bezknots[G1H_FINALDEG-1], degu, lknu, auxpknots, basp );

    omc += (G1H_OMCDEG+1);
    omcd += G1H_OMCDEG;
  }

finish:
  pkv_SetScratchMemTop ( sp );
} /*g1h_DrawBasAuxPatchesd*/

void g1h_DrawSplBasFunctiond ( GHoleDomaind *domain, int fn,
             void (*drawpatch) ( int n, int lknu, const double *knu,
                                 int m, int lknv, const double *knv,
                                 const point3d *cp ) )
{
  void     *sp;
  G1HolePrivateRecd  *privateG1;
  G1HoleSPrivateRecd *sprivate;
  int      hole_k, nfunc_a, nfunc_b, nfunc_c, nfunc_d, bfn;
  int      i, j, n, m, k0, k1, k2;
  vector2d *c00, *c01, *c10, *c11, *d00, *d01, *d10, *d11;
  double   *fc00, *fc01, *fc10, *fc11,
           *fd00, *fd01, *fd10, *fd11;
  point2d  *di, *aux;
  point3d  *cp;
  int      degu, degv, lastomcknot, lastpvknot,
           lastcknot, lastfpknot, lastukn, lastvkn, ncp;
  double   *bezknots, *omcknots, *pvknots, *cknots;
  double   *fcomc, *pv, *pu, *auxbfp, *bfp, *uknots, *vknots;
  double   zero[2] = {0.0,0.0};

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
  di = (point2d*)pkv_GetScratchMem ( m*m*sizeof(point2d) );
  aux = (point2d*)pkv_GetScratchMem ( n*m*sizeof(point2d) );
  cp = (point3d*)pkv_GetScratchMem ( n*n*sizeof(point3d) );
  auxbfp = pkv_GetScratchMemd ( 6*(lastpvknot-G1_CROSS01DEG) );
  bfp = pkv_GetScratchMemd ( (lastpvknot-G1_CROSS01DEG)*(lastpvknot-G1_CROSS01DEG) );
  fcomc = pkv_GetScratchMemd ( lastomcknot+2*lastpvknot-
                               (G1_CROSS00DEG+2*G1_CROSS01DEG) );
  uknots = pkv_GetScratchMemd ( 2*(lastfpknot+1) );
  if ( !di || !aux || !cp || !auxbfp || !bfp || !fcomc || !uknots )
    goto finish;
  pv = &fcomc[lastomcknot-G1_CROSS00DEG];
  pu = &pv[lastpvknot-G1_CROSS01DEG];
  vknots = &uknots[lastfpknot+1];

  for ( i = 0; i < hole_k; i++ ) {
        /* get the domain patch */
    _g1h_GetDiPatchCurvesd ( domain, i,
                             &c00, &c01, &c10, &c11, &d00, &d01, &d10, &d11 );
    mbs_BezC1CoonsToBezd ( 2,
       G1_CROSS00DEG, (double*)c00, G1_CROSS01DEG, (double*)c01,
       3, (double*)c10, G1_CROSS11DEG, (double*)c11,
       G1_CROSS00DEG, (double*)d00, G1_CROSS01DEG, (double*)d01,
       3, (double*)d10, G1_CROSS11DEG, (double*)d11,
       &degu, &degv, (double*)di );

        /* get the basis function patch */
    if ( fn < nfunc_a ) {                       /* block A */
      _g1h_GetBFAPatchCurvesd ( domain, fn, i,
                                &fc00, &fc01, &fd00, &fd01 );
      mbs_BezC1CoonsToBezd ( 1,
                    G1_CROSS00DEG, fc00, G1_CROSS01DEG, fc01, 1, zero, 1, zero,
                    G1_CROSS00DEG, fd00, G1_CROSS01DEG, fd01, 1, zero, 1, zero,
                    &degu, &degv, bfp );
      ncp = m*m;
      pkv_Selectd ( ncp, 2, 2, 3, (double*)di, cp );
      pkv_Selectd ( ncp, 1, 1, 3, bfp, &cp[0].z );
      drawpatch ( G1H_FINALDEG, 2*G1H_FINALDEG+1, bezknots,
                  G1H_FINALDEG, 2*G1H_FINALDEG+1, bezknots, cp );
    }
    else if ( fn < nfunc_a+nfunc_b ) {          /* block B */
      _g1h_GetBFBPatchCurvesd ( domain, fn-nfunc_a, i,
                                &fc00, &fc01, &fc10, &fc11,
                                &fd00, &fd01, &fd10, &fd11 );
      mbs_BezC1CoonsToBezd ( 1,
                    G1_CROSS00DEG, fc00, G1_CROSS01DEG, fc01,
                    G1_CROSS10DEG, fc10, G1_CROSS11DEG, fc11,
                    G1_CROSS00DEG, fd00, G1_CROSS01DEG, fd01,
                    G1_CROSS10DEG, fd10, G1_CROSS11DEG, fd11,
                    &degu, &degv, bfp );
      ncp = m*m;
      pkv_Selectd ( ncp, 2, 2, 3, (double*)di, cp );
      pkv_Selectd ( ncp, 1, 1, 3, bfp, &cp[0].z );
      drawpatch ( G1H_FINALDEG, 2*G1H_FINALDEG+1, bezknots,
                  G1H_FINALDEG, 2*G1H_FINALDEG+1, bezknots, cp );
    }
    else if ( fn < nfunc_a+nfunc_b+nfunc_c ) {  /* block C */
      n = lastcknot-G1H_FINALDEG;
      m = G1H_FINALDEG+1;
      mbs_multiAdjustBSCRepd ( 1, 2*(degv+1), degu, 2*degu+1,
                               &bezknots[G1H_FINALDEG-degu], 0, (double*)di,
                               G1H_FINALDEG, lastcknot, cknots, 0, (double*)aux );
      mbs_multiAdjustBSCRepd ( n, 2, degv, 2*degv+1,
                               &bezknots[G1H_FINALDEG-degv], 2*(degv+1), (double*)aux,
                               G1H_FINALDEG, lastcknot, cknots, 2*n, (double*)cp );
      pkv_Rearranged ( n*n, 2, 2, 3, cp );
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
      mbs_multiAdjustBSCRepd ( 1, 2*(degv+1), degu, 2*degu+1,
                       &bezknots[G1H_FINALDEG-degu], 0, (double*)di,
                       G1H_FINALDEG, lastpvknot, pvknots, 0, (double*)aux );
      mbs_multiAdjustBSCRepd ( n, 2, degv, 2*degv+1,
                       &bezknots[G1H_FINALDEG-degv], 2*(degv+1), (double*)aux,
                       G1H_FINALDEG, lastpvknot, pvknots, 2*n, (double*)cp );
      pkv_Rearranged ( n*n, 2, 2, 3, cp );
      for ( j = 0; j < n*n; j++ )
        cp[j].z = 0.0;

      bfn = fn - (nfunc_a+nfunc_b+nfunc_c);
      k0 = bfn / sprivate->dsize;
      if ( k0 == i ) {
        _g1h_GetSplDBasisCrossDerd ( domain, fn, i, fcomc, pv, pu );
        mbs_BSC1CoonsToBSd ( 1, G1_CROSS00DEG, lastomcknot, omcknots, fcomc,
            G1_CROSS01DEG, lastpvknot, pvknots, pv,
            1, 3, &bezknots[G1H_FINALDEG-1], zero,
            1, 3, &bezknots[G1H_FINALDEG-1], zero,
            1, 3, &bezknots[G1H_FINALDEG-1], zero,
            1, 3, &bezknots[G1H_FINALDEG-1], zero,
            1, 3, &bezknots[G1H_FINALDEG-1], zero,
            1, 3, &bezknots[G1H_FINALDEG-1], zero,
            &degu, &lastukn, uknots, &degv, &lastvkn, vknots, auxbfp );
        mbs_multiAdjustBSCRepd ( n, 1,
                      degv, lastvkn, vknots, lastvkn-degv, auxbfp,
                      G1H_FINALDEG, lastpvknot, pvknots, n, bfp );
        pkv_Selectd ( n*n, 1, 1, 3, bfp, &cp[0].z );
      }
      else if ( k0 == (i+1) % hole_k ) {
        _g1h_GetSplDBasisCrossDerd ( domain, fn, k0, fcomc, pv, pu );
        mbs_BSC1CoonsToBSd ( 1,
            1, 3, &bezknots[G1H_FINALDEG-1], zero,
            1, 3, &bezknots[G1H_FINALDEG-1], zero,
            1, 3, &bezknots[G1H_FINALDEG-1], zero,
            1, 3, &bezknots[G1H_FINALDEG-1], zero,
            G1_CROSS00DEG, lastomcknot, omcknots, fcomc,
            G1_CROSS01DEG, lastpvknot, pvknots, pu,
            1, 3, &bezknots[G1H_FINALDEG-1], zero,
            1, 3, &bezknots[G1H_FINALDEG-1], zero,
            &degu, &lastukn, uknots, &degv, &lastvkn, vknots, auxbfp );
        mbs_multiAdjustBSCRepd ( 1, n,
                                 degu, lastukn, uknots, 0, auxbfp,
                                 G1H_FINALDEG, lastpvknot, pvknots, n, bfp );
        pkv_Selectd ( n*n, 1, 1, 3, bfp, &cp[0].z );
      }

      drawpatch ( G1H_FINALDEG, lastpvknot, pvknots,
                  G1H_FINALDEG, lastpvknot, pvknots, cp );
    }
  }

finish:
  pkv_SetScratchMemTop ( sp );
} /*g1h_DrawSplBasFunctiond*/

/* ////////////////////////////////////////////////////////////////////////// */
void g1h_DrawSplBFAomcd ( GHoleDomaind *domain, int fn,
                          void (*drawpoly)(int degree, int lastknot,
                                           const double *knots,
                                           const double *coeff) )
{
  void  *sp;
  G1HoleSPrivateRecd *privateS;
  int   hole_k, i, lkn;
  double *c00, *c01, *d00, *d01, *c, *knb, *kns;

  if ( !(privateS = domain->SprivateG1) )
    return;
  lkn = privateS->lastomcknot;
  knb = &privateS->bezknots[G1H_FINALDEG-G1H_OMCDEG];
  kns = privateS->omcknots;
  sp = pkv_GetScratchMemTop ();
  if ( (c = pkv_GetScratchMemd ( lkn-G1H_OMCDEG )) ) {
    hole_k = domain->hole_k;
    for ( i = 0; i < hole_k; i++ ) {
      _g1h_GetBFAPatchCurvesd ( domain, fn, i, &c00, &c01, &d00, &d01 );
      mbs_AdjustBSCRepC1d ( G1H_OMCDEG, 2*G1H_OMCDEG+1, knb, c00,
                            G1H_OMCDEG, lkn, kns, c );
      drawpoly ( G1H_OMCDEG, lkn, kns, c );
    }
  }
  pkv_SetScratchMemTop ( sp );
} /*g1h_DrawSplBFAomcd*/

void g1h_DrawSplBFBomcd ( GHoleDomaind *domain, int fn,
                          void (*drawpoly)(int degree, int lastknot,
                                           const double *knots,
                                           const double *coeff) )
{
  void  *sp;
  G1HoleSPrivateRecd *privateS;
  int   hole_k, i, lkn;
  double *c00, *c01, *c10, *c11, *d00, *d01, *d10, *d11, *c, *knb, *kns;

  if ( !(privateS = domain->SprivateG1) )
    return;
  lkn = privateS->lastomcknot;
  knb = &privateS->bezknots[G1H_FINALDEG-G1H_OMCDEG];
  kns = privateS->omcknots;
  sp = pkv_GetScratchMemTop ();
  if ( (c = pkv_GetScratchMemd ( lkn-G1H_OMCDEG )) ) {
    hole_k = domain->hole_k;
    for ( i = 0; i < hole_k; i++ ) {
      _g1h_GetBFBPatchCurvesd ( domain, fn, i, &c00, &c01, &c10, &c11,
                                &d00, &d01, &d10, &d11 );
      mbs_AdjustBSCRepC1d ( G1H_OMCDEG, 2*G1H_OMCDEG+1, knb, c00,
                            G1H_OMCDEG, lkn, kns, c );
      drawpoly ( G1H_OMCDEG, lkn, kns, c );
    }
  }
  pkv_SetScratchMemTop ( sp );
} /*g1h_DrawSplBFBomcd*/

void g1h_DrawSplBFDomcd ( GHoleDomaind *domain, int fn,
                          void (*drawpoly)(int degree, int lastknot,
                                           const double *knots,
                                           const double *coeff) )
{
  void  *sp;
  G1HolePrivateRecd *privateG1;
  G1HoleSPrivateRecd *privateS;
  int   hole_k, i, lkn, lpvkn;
  double *kns, *fcomc, *pv, *pu;

  privateG1 = domain->privateG1;
  if ( !(privateS = domain->SprivateG1) )
    return;

  lkn = privateS->lastomcknot;
  lpvkn = privateS->lastpvknot;
  kns = privateS->omcknots;
  sp = pkv_GetScratchMemTop ();
  if ( (fcomc = pkv_GetScratchMemd (lkn-G1H_OMCDEG+2*(lpvkn-G1_CROSS01DEG))) ) {
    pv = &fcomc[lkn-G1H_OMCDEG];
    pu = &pv[lpvkn-G1_CROSS01DEG];
    hole_k = domain->hole_k;
    fn += privateG1->nfunc_a+privateG1->nfunc_b+privateS->nsfunc_c;
    for ( i = 0; i < hole_k; i++ ) {
      _g1h_GetSplDBasisCrossDerd ( domain, fn, i, fcomc, pv, pu );
      drawpoly ( G1H_OMCDEG, lkn, kns, fcomc );
    }
  }
  pkv_SetScratchMemTop ( sp );
} /*g1h_DrawSplBFDomcd*/

void g1h_DrawSplFinalSurfBCd ( GHoleDomaind *domain,
                               int spdimen, const double *hole_cp,
                               const double *acoeff,
                               void (*drawcurve)(int degree, int spdimen,
                                             int lastknot, const double *knots,
                                             const double *cp) )
{
  void *sp;
  int  hole_k, nfunc_a, nfunc_b, nfunc_c, nfunc_d;
  GHolePrivateRecd  *privateG;
  G1HolePrivateRecd *privateG1;
  G1HoleSPrivateRecd *privateS;
  unsigned char *bfcpn;
  double *x, *y, *c;
  double *c00, *c01, *c10, *c11, *d00, *d01, *d10, *d11;
  int   i, j, lkn, lpvkn;
  double *knb, *kns;

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
  if ( (x = pkv_GetScratchMemd ( 3*spdimen*(lpvkn-G1H_OMCDEG) )) ) {
    y = &x[spdimen*(lpvkn-G1H_OMCDEG)];
    c = &y[spdimen*(lpvkn-G1H_OMCDEG)];
    for ( i = 0; i < hole_k; i++ ) {
      memset ( x, 0, spdimen*(G1H_OMCDEG+1)*sizeof(double) );
      for ( j = 0; j < nfunc_a; j++ ) {
        _g1h_GetBFAPatchCurvesd ( domain, j, i, &c00, &c01, &d00, &d01 );
        mbs_AdjustBSCRepC1d ( G1H_OMCDEG, 2*G1H_OMCDEG+1, knb, c00,
                              G1H_OMCDEG, lkn, kns, c );
        pkn_MultMatrixd ( lkn-G1H_OMCDEG, 1, 1, c, spdimen, 0,
                          &acoeff[(nfunc_c+nfunc_d+j)*spdimen], spdimen, y );
        pkn_AddMatrixd ( 1, spdimen*(lkn-G1H_OMCDEG), 0, x, 0, y, 0, x );
      }
      for ( j = 0; j < nfunc_d; j++ ) {
        _g1h_GetSplDBasisCrossDerd ( domain, j+nfunc_a+nfunc_b+nfunc_c, i,
                                     y, c, c );
        pkn_MultMatrixd ( lkn-G1H_OMCDEG, 1, 1, y, spdimen, 0,
                          &acoeff[(nfunc_c+j)*spdimen], spdimen, c );
        pkn_AddMatrixd ( 1, spdimen*(lkn-G1H_OMCDEG), 0, x, 0, c, 0, x );
      }
      for ( j = 0; j < nfunc_b; j++ ) {
        _g1h_GetBFBPatchCurvesd ( domain, j, i, &c00, &c01, &c10, &c11,
                                  &d00, &d01, &d10, &d11 );
        mbs_AdjustBSCRepC1d ( G1H_OMCDEG, 2*G1H_OMCDEG+1, knb, c00,
                              G1H_OMCDEG, lkn, kns, c );
        pkn_MultMatrixd ( lkn-G1H_OMCDEG, 1, 1, c, spdimen, 0,
                          &hole_cp[bfcpn[j]], spdimen, y );
        pkn_AddMatrixd ( 1, spdimen*(lkn-G1H_OMCDEG), 0, x, 0, y, 0, x );
      }
      drawcurve ( G1H_OMCDEG, spdimen, lkn, kns, x );
    }
  }
  pkv_SetScratchMemTop ( sp );
} /*g1h_DrawSplFinalSurfBCd*/

/* ///////////////////////////////////////////////////////////////////////// */
void g1h_DrawSplMatricesd ( GHoleDomaind *domain,
                            void (*drawmatrix)( int k, int r, int s, int t,
                            double *A, double *B ) )
{
  void  *sp;
  G1HolePrivateRecd  *privateG1;
  G1HoleSPrivateRecd *sprivate;
  double *A, *B;
  int    hole_k, nfunc_a, nfunc_b, nfunc_c, nfunc_d, s1, s2;

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
  A = pkv_GetScratchMemd ( s1 );
  B = pkv_GetScratchMemd ( s2 );
  if ( !A || !B )
    goto finish;
  memcpy ( A, sprivate->SAMat, s1*sizeof(double) );
  memcpy ( B, sprivate->SBMat, s2*sizeof(double) );
  drawmatrix ( hole_k, nfunc_c/hole_k, nfunc_d/hole_k, nfunc_a, A, B );

finish:
  pkv_SetScratchMemTop ( sp );
} /*g1h_DrawSplMatricesd*/

