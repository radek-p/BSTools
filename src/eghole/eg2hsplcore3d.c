
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

#undef CONST_
#define CONST_

#include "eg2holed.h"
#include "eg2hprivated.h"
#include "eg2herror.h"

/* ///////////////////////////////////////////////////////////////////////// */
boolean g2h_DecomposeSplMatrixd ( GHoleDomaind *domain )
{
  void   *sp;
  G2HolePrivateRecd  *privateG2;
  G2HoleSPrivateRecd *sprivate;
  int    hole_k, nfunc_a, nfunc_c, nfunc_d;
  int    size;
  double *lmat;

  sp = pkv_GetScratchMemTop ();
  privateG2 = domain->privateG2;
  sprivate = domain->SprivateG2;
  if ( !sprivate )
    goto failure;
  if ( !sprivate->SLMat ) {
    if ( !sprivate->SAMat )
      if ( !g2h_ComputeSplFormMatrixd ( domain ) )
        goto failure;
    hole_k = domain->hole_k;
    nfunc_a = privateG2->nfunc_a;
    nfunc_c = sprivate->nsfunc_c;
    nfunc_d = sprivate->nsfunc_d;
    size = pkn_Block2ArraySize ( hole_k, nfunc_c/hole_k, nfunc_d/hole_k,
                                 nfunc_a );
    lmat = sprivate->SLMat = malloc ( size*sizeof(double) );
    if ( !lmat )
      goto failure;
    memcpy ( lmat, sprivate->SAMat, size*sizeof(double) );
    if ( !pkn_Block2CholeskyDecompMd ( hole_k, nfunc_c/hole_k, nfunc_d/hole_k,
                                       nfunc_a, lmat ) ) {
      domain->error_code = G2H_ERROR_NONPOSITIVE_MATRIX;
      goto failure;
    }
  }

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*g2h_DecomposeSplMatrixd*/

boolean _g2h_SetSplRightSided ( GHoleDomaind *domain,
                                int spdimen, CONST_ double *hole_cp,
                                double *fc00, double *b )
{
  void   *sp;
  GHolePrivateRecd   *privateG;
  G2HolePrivateRecd  *privateG2;
  G2HoleSPrivateRecd *sprivate;
  int    hole_k, nfunc_a, nfunc_b, nfunc_c, nfunc_d, nbf;
  double *bmat, *x, *y, *cp;
  unsigned char *bfcpn;
  double *bc00, *bc01, *bc02, *bc10, *bc11, *bc12,
         *bd00, *bd01, *bd02, *bd10, *bd11, *bd12;
  double        *fc01, *fc02, *fc10, *fc11, *fc12,
         *fd00, *fd01, *fd02, *fd10, *fd11, *fd12;
  int    i, j, k, l;

  sp = pkv_GetScratchMemTop ();
  hole_k = domain->hole_k;
  privateG  = domain->privateG;
  privateG2 = domain->privateG2;
  sprivate = domain->SprivateG2;
  if ( !sprivate )
    goto failure;
  nfunc_a = privateG2->nfunc_a;
  nfunc_b = privateG2->nfunc_b;
  nfunc_c = sprivate->nsfunc_c;
  nfunc_d = sprivate->nsfunc_d;
  nbf = nfunc_a+nfunc_c+nfunc_d;
  bfcpn = privateG->bfcpn;
  bmat = sprivate->SBMat;
  x = pkv_GetScratchMemd ( spdimen*nfunc_b );
  y = pkv_GetScratchMemd ( spdimen*(G2H_FINALDEG+1)*(G2H_FINALDEG+1) );
  if ( !x || !y )
    goto failure;

  memset ( fc00, 0, (G2_CROSSDEGSUM+6)*2*hole_k*spdimen*sizeof(double) );
  G2GetFCAddresses ();

  for ( j = 0; j < nfunc_b; j++ ) {
    l = bfcpn[j];
    cp = &hole_cp[l*spdimen];
    memcpy  ( &x[j*spdimen], cp, spdimen*sizeof(double) );
    for ( i = k = 0;  i < hole_k;  i++, k += spdimen ) {
      _g2h_GetBFBPatchCurvesd ( domain, j, i,
                                &bc00, &bc01, &bc02, &bc10, &bc11, &bc12,
                                &bd00, &bd01, &bd02, &bd10, &bd11, &bd12 );

      pkn_MultMatrixd ( G2_CROSS00DEG+1, 1, 1, bc00, spdimen, 0, cp, spdimen, y );
      pkn_AddMatrixd ( 1, spdimen*(G2_CROSS00DEG+1), 0, &fc00[k*(G2_CROSS00DEG+1)],
                       0, y, 0, &fc00[k*(G2_CROSS00DEG+1)] );
      pkn_MultMatrixd ( G2_CROSS01DEG+1, 1, 1, bc01, spdimen, 0, cp, spdimen, y );
      pkn_AddMatrixd ( 1, spdimen*(G2_CROSS01DEG+1), 0, &fc01[k*(G2_CROSS01DEG+1)],
                       0, y, 0, &fc01[k*(G2_CROSS01DEG+1)] );
      pkn_MultMatrixd ( G2_CROSS02DEG+1, 1, 1, bc02, spdimen, 0, cp, spdimen, y );
      pkn_AddMatrixd ( 1, spdimen*(G2_CROSS02DEG+1), 0, &fc02[k*(G2_CROSS02DEG+1)],
                       0, y, 0, &fc02[k*(G2_CROSS02DEG+1)] );
      pkn_MultMatrixd ( G2_CROSS10DEG+1, 1, 1, bc10, spdimen, 0, cp, spdimen, y );
      pkn_AddMatrixd ( 1, spdimen*(G2_CROSS10DEG+1), 0, &fc10[k*(G2_CROSS10DEG+1)],
                       0, y, 0, &fc10[k*(G2_CROSS10DEG+1)] );
      pkn_MultMatrixd ( G2_CROSS11DEG+1, 1, 1, bc11, spdimen, 0, cp, spdimen, y );
      pkn_AddMatrixd ( 1, spdimen*(G2_CROSS11DEG+1), 0, &fc11[k*(G2_CROSS11DEG+1)],
                       0, y, 0, &fc11[k*(G2_CROSS11DEG+1)] );
      pkn_MultMatrixd ( G2_CROSS12DEG+1, 1, 1, bc12, spdimen, 0, cp, spdimen, y );
      pkn_AddMatrixd ( 1, spdimen*(G2_CROSS12DEG+1), 0, &fc12[k*(G2_CROSS12DEG+1)],
                       0, y, 0, &fc12[k*(G2_CROSS12DEG+1)] );

      pkn_MultMatrixd ( G2_CROSS00DEG+1, 1, 1, bd00, spdimen, 0, cp, spdimen, y );
      pkn_AddMatrixd ( 1, spdimen*(G2_CROSS00DEG+1), 0, &fd00[k*(G2_CROSS00DEG+1)],
                       0, y, 0, &fd00[k*(G2_CROSS00DEG+1)] );
      pkn_MultMatrixd ( G2_CROSS01DEG+1, 1, 1, bd01, spdimen, 0, cp, spdimen, y );
      pkn_AddMatrixd ( 1, spdimen*(G2_CROSS01DEG+1), 0, &fd01[k*(G2_CROSS01DEG+1)],
                       0, y, 0, &fd01[k*(G2_CROSS01DEG+1)] );
      pkn_MultMatrixd ( G2_CROSS02DEG+1, 1, 1, bd02, spdimen, 0, cp, spdimen, y );
      pkn_AddMatrixd ( 1, spdimen*(G2_CROSS02DEG+1), 0, &fd02[k*(G2_CROSS02DEG+1)],
                       0, y, 0, &fd02[k*(G2_CROSS02DEG+1)] );
      pkn_MultMatrixd ( G2_CROSS10DEG+1, 1, 1, bd10, spdimen, 0, cp, spdimen, y );
      pkn_AddMatrixd ( 1, spdimen*(G2_CROSS10DEG+1), 0, &fd10[k*(G2_CROSS10DEG+1)],
                       0, y, 0, &fd10[k*(G2_CROSS10DEG+1)] );
      pkn_MultMatrixd ( G2_CROSS11DEG+1, 1, 1, bd11, spdimen, 0, cp, spdimen, y );
      pkn_AddMatrixd ( 1, spdimen*(G2_CROSS11DEG+1), 0, &fd11[k*(G2_CROSS11DEG+1)],
                       0, y, 0, &fd11[k*(G2_CROSS11DEG+1)] );
      pkn_MultMatrixd ( G2_CROSS12DEG+1, 1, 1, bd12, spdimen, 0, cp, spdimen, y );
      pkn_AddMatrixd ( 1, spdimen*(G2_CROSS12DEG+1), 0, &fd12[k*(G2_CROSS12DEG+1)],
                       0, y, 0, &fd12[k*(G2_CROSS12DEG+1)] );
    }
  }

  if ( b )
    pkn_MultMatrixd ( nbf, nfunc_b, nfunc_b, bmat,
                      spdimen, spdimen, x, spdimen, b );

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*_g2h_SetSplRightSided*/

boolean _g2h_OutputSplPatchesd ( GHoleDomaind *domain,
                int spdimen, CONST_ double *x, double *fc00, void *usrptr,
                void (*outpatch) ( int n, int lknu, const double *knu,
                                   int m, int lknv, const double *knv,
                                   const double *cp, void *usrptr ) )
{
  void   *sp;
  G2HolePrivateRecd  *privateG2;
  G2HoleSPrivateRecd *sprivate;
  double *bc00, *bc01, *bc02, *bd00, *bd01, *bd02;
  double        *fc01, *fc02, *fc10, *fc11, *fc12,
         *fd00, *fd01, *fd02, *fd10, *fd11, *fd12;
  double *sfc00, *sfc01, *sfc02, *sfd00, *sfd01, *sfd02;
  double *fcomc, *pv, *pvv, *pu, *puu;
  double *y, *yv, *yvv, *yu, *yuu, *z, *xx, *cp;
  int    hole_k, nk, m1, m2, csiz, csize, dsize;
  int    nfunc_a, nfunc_b, nfunc_c, nfunc_d, nabc;
  int    degu, degv;
  int    lastomcknot, lastpvknot, lastpvvknot, lastcknot, lastfpknot;
  double *omcknots, *pvknots, *pvvknots, *bezknots, *cknots, *fpknots, *aknots;
  int    i, j, k, l, m;


  sp = pkv_GetScratchMemTop ();

  hole_k = domain->hole_k;
  privateG2 = domain->privateG2;
  sprivate = domain->SprivateG2;
  nk = sprivate->nk;
  m1 = sprivate->m1;
  m2 = sprivate->m2;
  csize = sprivate->csize;
  csiz = G2H_FINALDEG-5+nk*m2;  /* this is sqrt(csize) */
  dsize = sprivate->dsize;
  nfunc_a = privateG2->nfunc_a;
  nfunc_b = privateG2->nfunc_b;
  nfunc_c = sprivate->nsfunc_c;
  nfunc_d = sprivate->nsfunc_d;
  nabc = nfunc_a+nfunc_b+nfunc_c;
  lastomcknot = sprivate->lastomcknot;
  lastpvknot  = sprivate->lastpvknot;
  lastpvvknot = sprivate->lastpvvknot;
  lastcknot   = sprivate->lastcknot;
  lastfpknot  = sprivate->lastfpknot;
  bezknots = sprivate->bezknots;
  omcknots = sprivate->omcknots;
  pvknots  = sprivate->pvknots;
  pvvknots = sprivate->pvvknots;
  cknots   = sprivate->cknots;
  fpknots  = sprivate->fpknots;

  aknots = pkv_GetScratchMemd ( 2*(lastfpknot+1) );
  sfc00 = pkv_GetScratchMemd ( (lastomcknot-G2_CROSS00DEG+lastpvknot-G2_CROSS01DEG+
                                lastpvvknot-G2_CROSS02DEG)*spdimen*2*hole_k );
  y = pkv_GetScratchMemd (
          3*spdimen*(lastfpknot-G2H_FINALDEG)*(lastfpknot-G2H_FINALDEG) );
  if ( !aknots || !sfc00 || !y )
    goto failure;

  G2GetFCAddresses ();
  G2GetSFCAddresses ();

  for ( j = 0; j < nfunc_a; j++ ) {
    cp = &x[(nfunc_c+nfunc_d+j)*spdimen];
    for ( i = k = 0;  i < hole_k;  i++, k += spdimen ) {
      _g2h_GetBFAPatchCurvesd ( domain, j, i,
                                &bc00, &bc01, &bc02, &bd00, &bd01, &bd02 );

      pkn_MultMatrixd ( G2_CROSS00DEG+1, 1, 1, bc00, spdimen, 0, cp, spdimen, y );
      pkn_SubtractMatrixd ( 1, spdimen*(G2_CROSS00DEG+1), 0, &fc00[k*(G2_CROSS00DEG+1)],
                            0, y, 0, &fc00[k*(G2_CROSS00DEG+1)] );
      pkn_MultMatrixd ( G2_CROSS01DEG+1, 1, 1, bc01, spdimen, 0, cp, spdimen, y );
      pkn_SubtractMatrixd ( 1, spdimen*(G2_CROSS01DEG+1), 0, &fc01[k*(G2_CROSS01DEG+1)],
                            0, y, 0, &fc01[k*(G2_CROSS01DEG+1)] );
      pkn_MultMatrixd ( G2_CROSS02DEG+1, 1, 1, bc02, spdimen, 0, cp, spdimen, y );
      pkn_SubtractMatrixd ( 1, spdimen*(G2_CROSS02DEG+1), 0, &fc02[k*(G2_CROSS02DEG+1)],
                            0, y, 0, &fc02[k*(G2_CROSS02DEG+1)] );

      pkn_MultMatrixd ( G2_CROSS00DEG+1, 1, 1, bd00, spdimen, 0, cp, spdimen, y );
      pkn_SubtractMatrixd ( 1, spdimen*(G2_CROSS00DEG+1), 0, &fd00[k*(G2_CROSS00DEG+1)],
                            0, y, 0, &fd00[k*(G2_CROSS00DEG+1)] );
      pkn_MultMatrixd ( G2_CROSS01DEG+1, 1, 1, bd01, spdimen, 0, cp, spdimen, y );
      pkn_SubtractMatrixd ( 1, spdimen*(G2_CROSS01DEG+1), 0, &fd01[k*(G2_CROSS01DEG+1)],
                            0, y, 0, &fd01[k*(G2_CROSS01DEG+1)] );
      pkn_MultMatrixd ( G2_CROSS02DEG+1, 1, 1, bd02, spdimen, 0, cp, spdimen, y );
      pkn_SubtractMatrixd ( 1, spdimen*(G2_CROSS02DEG+1), 0, &fd02[k*(G2_CROSS02DEG+1)],
                            0, y, 0, &fd02[k*(G2_CROSS02DEG+1)] );
    }
  }

        /* At this point the term made of the block A and B basis functions */
        /* is represented in the Coons form, by the Bezier curves. The next */
        /* step is to obtain the B-spline representation of the curves and  */
        /* add the term represented by the block D functions.               */

  yv = &y[(lastfpknot-G2H_FINALDEG)*spdimen];
  yvv = &yv[(lastpvknot-G2_CROSS01DEG)*spdimen];
  yu = &yvv[(lastpvvknot-G2_CROSS02DEG)*spdimen];
  yuu = &yu[(lastpvknot-G2_CROSS01DEG)*spdimen];
  fcomc = &yuu[(lastpvvknot-G2_CROSS02DEG)*spdimen];
  pv = &fcomc[lastomcknot-G2_CROSS00DEG];
  pvv = &pv[lastpvknot-G2_CROSS01DEG];
  pu = &pvv[lastpvvknot-G2_CROSS02DEG];
  puu = &pu[lastpvknot-G2_CROSS01DEG];
  z = &puu[lastpvvknot-G2_CROSS02DEG];
  for ( i = 0;  i < hole_k;  i++ ) {
    l = (i+hole_k-1) % hole_k;
    memset ( y, 0, (lastfpknot-G2H_FINALDEG+
             2*(lastpvknot-G2_CROSS01DEG+lastpvvknot-G2_CROSS02DEG))*
             spdimen*sizeof(double) );
    pkv_Selectd ( nk*m1, spdimen, 3*spdimen, spdimen,
                  &x[(nfunc_c+i*dsize)*spdimen], &y[5*spdimen] );
    mbs_multiSubtractBSCurvesd ( 1, spdimen, G2_CROSS00DEG, 2*G2_CROSS00DEG+1,
            &bezknots[G2H_FINALDEG-G2_CROSS00DEG], 0,
            &fc00[i*(G2_CROSS00DEG+1)*spdimen],
            G2_CROSS00DEG, lastomcknot, omcknots, 0, y, &degu, &m, z, 0,
            &sfc00[i*(lastomcknot-G2_CROSS00DEG)*spdimen] );
    memcpy ( &sfd00[l*(lastomcknot-G2_CROSS00DEG)*spdimen],
             &sfc00[i*(lastomcknot-G2_CROSS00DEG)*spdimen],
             (lastomcknot-G2_CROSS00DEG)*spdimen*sizeof(double) );

    for ( j = 0; j < nk*m1; j++ ) {
      xx = &x[(nfunc_c+i*dsize+3*j)*spdimen];

      _g2h_GetSplDBasisCrossDerd ( domain, nabc+i*dsize+3*j,
                                   i, fcomc, pv, pvv, pu, puu );
      pkn_MultMatrixAddd ( lastpvknot-G2_CROSS01DEG, 1, 1, pv,
                           spdimen, 0, xx, spdimen, yv );
      pkn_MultMatrixAddd ( lastpvknot-G2_CROSS01DEG, 1, 1, pu,
                           spdimen, 0, xx, spdimen, yu );
      pkn_MultMatrixAddd ( lastpvvknot-G2_CROSS02DEG, 1, 1, pvv,
                           spdimen, 0, xx, spdimen, yvv );
      pkn_MultMatrixAddd ( lastpvvknot-G2_CROSS02DEG, 1, 1, puu,
                           spdimen, 0, xx, spdimen, yuu );

      _g2h_GetSplDBasisCrossDerd ( domain, nabc+i*dsize+3*j+1,
                                   i, fcomc, pv, pvv, pu, puu );
      pkn_MultMatrixAddd ( lastpvknot-G2_CROSS01DEG, 1, 1, pv,
                           spdimen, 0, &xx[spdimen], spdimen, yv );
      pkn_MultMatrixAddd ( lastpvknot-G2_CROSS01DEG, 1, 1, pu,
                           spdimen, 0, &xx[spdimen], spdimen, yu );
      pkn_MultMatrixAddd ( lastpvvknot-G2_CROSS02DEG, 1, 1, pvv,
                           spdimen, 0, &xx[spdimen], spdimen, yvv );
      pkn_MultMatrixAddd ( lastpvvknot-G2_CROSS02DEG, 1, 1, puu,
                           spdimen, 0, &xx[spdimen], spdimen, yuu );

      _g2h_GetSplDBasisCrossDerd ( domain, nabc+i*dsize+3*j+2,
                                   i, fcomc, pv, pvv, pu, puu );
      pkn_MultMatrixAddd ( lastpvvknot-G2_CROSS02DEG, 1, 1, pvv,
                           spdimen, 0, &xx[2*spdimen], spdimen, yvv );
      pkn_MultMatrixAddd ( lastpvvknot-G2_CROSS02DEG, 1, 1, puu,
                           spdimen, 0, &xx[2*spdimen], spdimen, yuu );
    }
    mbs_multiSubtractBSCurvesd ( 1, spdimen, G2_CROSS01DEG, 2*G2_CROSS01DEG+1,
            &bezknots[G2H_FINALDEG-G2_CROSS01DEG], 0, &fc01[i*(G2_CROSS01DEG+1)*spdimen],
            G2_CROSS01DEG, lastpvknot, pvknots, 0, yv, &degu, &m, z, 0,
            &sfc01[i*(lastpvknot-G2_CROSS01DEG)*spdimen] );
    mbs_multiSubtractBSCurvesd ( 1, spdimen, G2_CROSS02DEG, 2*G2_CROSS02DEG+1,
            &bezknots[G2H_FINALDEG-G2_CROSS02DEG], 0, &fc02[i*(G2_CROSS02DEG+1)*spdimen],
            G2_CROSS02DEG, lastpvvknot, pvvknots, 0, yvv, &degu, &m, z, 0,
            &sfc02[i*(lastpvvknot-G2_CROSS02DEG)*spdimen] );
    mbs_multiSubtractBSCurvesd ( 1, spdimen, G2_CROSS01DEG, 2*G2_CROSS01DEG+1,
            &bezknots[G2H_FINALDEG-G2_CROSS01DEG], 0, &fd01[l*(G2_CROSS01DEG+1)*spdimen],
            G2_CROSS01DEG, lastpvknot, pvknots, 0, yu, &degv, &m, z, 0,
            &sfd01[l*(lastpvknot-G2_CROSS01DEG)*spdimen] );
    mbs_multiSubtractBSCurvesd ( 1, spdimen, G2_CROSS02DEG, 2*G2_CROSS02DEG+1,
            &bezknots[G2H_FINALDEG-G2_CROSS02DEG], 0, &fd02[l*(G2_CROSS02DEG+1)*spdimen],
            G2_CROSS02DEG, lastpvvknot, pvvknots, 0, yuu, &degv, &m, z, 0,
            &sfd02[l*(lastpvvknot-G2_CROSS02DEG)*spdimen] );
  }

        /* The next step is to find the B-spline representation of the     */
        /* Coons patch and add the terms represented by the block C basis  */
        /* functions */

  yu = &y[spdimen*(lastfpknot-G2H_FINALDEG)*(lastfpknot-G2H_FINALDEG)];
  yv = &y[spdimen*(lastfpknot-G2H_FINALDEG)*(lastfpknot-G2H_FINALDEG)];
  for ( i = 0; i < hole_k; i++ ) {
    mbs_BSC2CoonsToBSd ( spdimen,
            G2_CROSS00DEG, lastomcknot, omcknots,
            &sfc00[i*(lastomcknot-G2_CROSS00DEG)*spdimen],
            G2_CROSS01DEG, lastpvknot,  pvknots,
            &sfc01[i*(lastpvknot-G2_CROSS01DEG)*spdimen],
            G2_CROSS02DEG, lastpvvknot, pvvknots,
            &sfc02[i*(lastpvvknot-G2_CROSS02DEG)*spdimen],
            G2_CROSS10DEG, 2*G2_CROSS10DEG+1, &bezknots[G2H_FINALDEG-G2_CROSS10DEG],
            &fc10[i*(G2_CROSS10DEG+1)*spdimen],
            G2_CROSS11DEG, 2*G2_CROSS11DEG+1, &bezknots[G2H_FINALDEG-G2_CROSS11DEG],
            &fc11[i*(G2_CROSS11DEG+1)*spdimen],
            G2_CROSS12DEG, 2*G2_CROSS12DEG+1, &bezknots[G2H_FINALDEG-G2_CROSS12DEG],
            &fc12[i*(G2_CROSS12DEG+1)*spdimen],
            G2_CROSS00DEG, lastomcknot, omcknots,
            &sfd00[i*(lastomcknot-G2_CROSS00DEG)*spdimen],
            G2_CROSS01DEG, lastpvknot,  pvknots,
            &sfd01[i*(lastpvknot-G2_CROSS01DEG)*spdimen],
            G2_CROSS02DEG, lastpvvknot, pvvknots,
            &sfd02[i*(lastpvvknot-G2_CROSS02DEG)*spdimen],
            G2_CROSS10DEG, 2*G2_CROSS10DEG+1, &bezknots[G2H_FINALDEG-G2_CROSS10DEG],
            &fd10[i*(G2_CROSS10DEG+1)*spdimen],
            G2_CROSS11DEG, 2*G2_CROSS11DEG+1, &bezknots[G2H_FINALDEG-G2_CROSS11DEG],
            &fd11[i*(G2_CROSS11DEG+1)*spdimen],
            G2_CROSS12DEG, 2*G2_CROSS12DEG+1, &bezknots[G2H_FINALDEG-G2_CROSS12DEG],
            &fd12[i*(G2_CROSS12DEG+1)*spdimen],
            &degu, &l, aknots, &degv, &m, &aknots[lastfpknot+1], y );

        /* convert A+B+D to the final form */
    mbs_multiAdjustBSCRepd ( 1, (m-degv)*spdimen, degu, l, aknots, 0, y,
                             G2H_FINALDEG, lastfpknot, fpknots, 0, yu );
    mbs_multiAdjustBSCRepd ( lastfpknot-G2H_FINALDEG, spdimen, degv, m,
                             &aknots[lastfpknot+1], (m-degv)*spdimen, yu,
                             G2H_FINALDEG, lastfpknot, fpknots,
                             (lastfpknot-G2H_FINALDEG)*spdimen, y );

        /* now process the terms with the block C basis functions */
    memset ( yu, 0, spdimen*(lastcknot-G2H_FINALDEG)*(lastcknot-G2H_FINALDEG)*
                    sizeof(double) );
    pkv_Selectd ( csiz, csiz*spdimen, csiz*spdimen,
                  (lastcknot-G2H_FINALDEG)*spdimen, &x[i*csize*spdimen],
                  &yu[(3*(lastcknot-G2H_FINALDEG)+3)*spdimen] );

    mbs_multiAdjustBSCRepd ( 1, (lastcknot-G2H_FINALDEG)*spdimen, G2H_FINALDEG,
                             lastcknot, cknots, 0, yu,
                             G2H_FINALDEG, lastfpknot, fpknots, 0, yv );
    mbs_multiAdjustBSCRepd ( lastfpknot-G2H_FINALDEG, spdimen, G2H_FINALDEG,
                             lastcknot, cknots, (lastcknot-G2H_FINALDEG)*spdimen,
                             yv, G2H_FINALDEG, lastfpknot, fpknots,
                             (lastfpknot-G2H_FINALDEG)*spdimen, yu );
    pkn_SubtractMatrixd ( lastfpknot-G2H_FINALDEG-6,
                          (lastfpknot-G2H_FINALDEG-6)*spdimen,
                          (lastfpknot-G2H_FINALDEG)*spdimen,
                          &y[(3*(lastfpknot-G2H_FINALDEG)+3)*spdimen],
                          (lastfpknot-G2H_FINALDEG)*spdimen,
                          &yu[(3*(lastfpknot-G2H_FINALDEG)+3)*spdimen],
                          (lastfpknot-G2H_FINALDEG)*spdimen,
                          &y[(3*(lastfpknot-G2H_FINALDEG)+3)*spdimen] );

    outpatch ( G2H_FINALDEG, lastfpknot, fpknots,
               G2H_FINALDEG, lastfpknot, fpknots, y, usrptr );
  }

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*_g2h_OutputSplPatchesd*/

/* ///////////////////////////////////////////////////////////////////////// */
boolean g2h_SplFillHoled ( GHoleDomaind *domain,
               int spdimen, CONST_ double *hole_cp,
               double *acoeff, void *usrptr,
               void (*outpatch) ( int n, int lknu, const double *knu,
                                  int m, int lknv, const double *knv,
                                  const double *cp, void *usrptr ) )
{ 
  void   *sp;
  G2HolePrivateRecd  *privateG2;
  G2HoleSPrivateRecd *sprivate;
  int    hole_k, nfunc_a, nfunc_c, nfunc_d, nbf;
  double *lmat, *x, *fc00;

  sp = pkv_GetScratchMemTop ();
  hole_k = domain->hole_k;
  privateG2 = domain->privateG2;
  sprivate = domain->SprivateG2;
  if ( !sprivate )
    goto failure;
  nfunc_a = privateG2->nfunc_a;
  nfunc_c = sprivate->nsfunc_c;
  nfunc_d = sprivate->nsfunc_d;
  nbf = nfunc_a+nfunc_c+nfunc_d;
  if ( !sprivate->SLMat )
    if ( !g2h_DecomposeSplMatrixd ( domain ) )
      goto failure;
  lmat = sprivate->SLMat;

  x = pkv_GetScratchMemd ( spdimen*nbf );
  fc00 = pkv_GetScratchMemd ( (G2_CROSSDEGSUM+6)*2*hole_k*spdimen );
  if ( !x || !fc00 )
    goto failure;

  if ( !_g2h_SetSplRightSided ( domain, spdimen, hole_cp, fc00, x ) )
    goto failure;

  pkn_Block2LowerTrMSolved ( hole_k, nfunc_c/hole_k, nfunc_d/hole_k, nfunc_a,
                             lmat, spdimen, spdimen, x );
  pkn_Block2UpperTrMSolved ( hole_k, nfunc_c/hole_k, nfunc_d/hole_k, nfunc_a,
                             lmat, spdimen, spdimen, x );
  if ( acoeff )
    memcpy ( acoeff, x, spdimen*nbf*sizeof(double) );

  if ( !_g2h_OutputSplPatchesd ( domain, spdimen, x, fc00, usrptr, outpatch ) )
    goto failure;

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*g2h_SplFillHoled*/

/* ///////////////////////////////////////////////////////////////////////// */
boolean g2h_SplFillHoleConstrd ( GHoleDomaind *domain,
               int spdimen, CONST_ double *hole_cp,
               int nconstr, CONST_ double *constr,
               double *acoeff, void *usrptr,
               void (*outpatch) ( int n, int lknu, const double *knu,
                                  int m, int lknv, const double *knv,
                                  const double *cp, void *usrptr ) )
{
  void   *sp;
  G2HolePrivateRecd  *privateG2;
  G2HoleSPrivateRecd *sprivate;
  double *fc00, *b, *x, *y, *cmat, *rcmat, *lmat;
  int    hole_k, nfunc_a, nfunc_c, nfunc_d, nbf;

  sp = pkv_GetScratchMemTop ();

  hole_k = domain->hole_k;
  privateG2 = domain->privateG2;
  sprivate = domain->SprivateG2;
  if ( !sprivate->SCmat || !sprivate->SRCmat ||
       sprivate->splnconstr != nconstr ) {
    domain->error_code = G2H_ERROR_UNDEFINED_CONSTR;
    goto failure;
  }
  nfunc_a = privateG2->nfunc_a;
  nfunc_c = sprivate->nsfunc_c;
  nfunc_d = sprivate->nsfunc_d;
  nbf = nfunc_a+nfunc_c+nfunc_d;
  cmat = sprivate->SCmat;
  rcmat = sprivate->SRCmat;
  lmat = sprivate->SLMat;

  b = pkv_GetScratchMemd ( spdimen*nbf );
  x = pkv_GetScratchMemd ( spdimen*nbf );
  fc00 = pkv_GetScratchMemd ( (G2_CROSSDEGSUM+6)*2*hole_k*spdimen );
  if ( !b || !x || !fc00 )
    goto failure;

  if ( !_g2h_SetSplRightSided ( domain, spdimen, hole_cp, fc00, x ) )
    goto failure;

  pkn_Block2LowerTrMSolved ( hole_k, nfunc_c/hole_k, nfunc_d/hole_k, nfunc_a,
                             lmat, spdimen, spdimen, x );
  pkn_Block2UpperTrMSolved ( hole_k, nfunc_c/hole_k, nfunc_d/hole_k, nfunc_a,
                             lmat, spdimen, spdimen, x );
              /* now x is the solution without the constraints and it is time */   
              /* to employ the constraints that have been previously set */
  pkn_MultMatrixd ( nconstr, nbf, nbf, cmat,
                    spdimen, spdimen, x, spdimen, b );
  pkn_AddMatrixd ( 1, spdimen*nconstr, 0, b, 0, constr, 0, b );
              /* solve the system with the Schur matrix */
  pkn_LowerTrMatrixSolved ( nconstr, rcmat, spdimen, spdimen, b, spdimen, b );
  pkn_UpperTrMatrixSolved ( nconstr, rcmat, spdimen, spdimen, b, spdimen, b );
              /* compute the constraints correction */
  y = pkv_GetScratchMemd ( spdimen*nbf );
  if ( !y ) {
    domain->error_code = G2H_ERROR_NO_SCRATCH_MEMORY;
    goto failure;
  }
  pkn_MultTMatrixd ( nconstr, nbf, nbf, cmat,
                     spdimen, spdimen, b, spdimen, y );
  pkn_Block2LowerTrMSolved ( hole_k, nfunc_c/hole_k, nfunc_d/hole_k, nfunc_a,
                             lmat, spdimen, spdimen, y );
  pkn_Block2UpperTrMSolved ( hole_k, nfunc_c/hole_k, nfunc_d/hole_k, nfunc_a,
                             lmat, spdimen, spdimen, y );
  pkn_SubtractMatrixd ( 1, spdimen*nbf, 0, x, 0, y, 0, x );
  pkv_FreeScratchMemd ( spdimen*nbf );
  if ( acoeff )
    memcpy ( acoeff, x, spdimen*nbf*sizeof(double) );

  if ( !_g2h_OutputSplPatchesd ( domain, spdimen, x, fc00, usrptr, outpatch ) )
    goto failure;

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*g2h_SplFillHoleConstrd*/

/* ///////////////////////////////////////////////////////////////////////// */
boolean g2h_SplFillHoleAltConstrd ( GHoleDomaind *domain,
               int spdimen, CONST_ double *hole_cp,
               int naconstr, CONST_ double *constr,
               double *acoeff, void *usrptr,
               void (*outpatch) ( int n, int lknu, const double *knu,
                                  int m, int lknv, const double *knv,
                                  const double *cp, void *usrptr ) )
{
  void   *sp;
  G2HolePrivateRecd  *privateG2;
  G2HoleSPrivateRecd *sprivate;
  double *fc00, *b, *x, *y, *acmat, *arcmat, *lmat;
  int    hole_k, nfunc_a, nfunc_b, nfunc_c, nfunc_d, nbf;

  sp = pkv_GetScratchMemTop ();

  hole_k = domain->hole_k;
  privateG2 = domain->privateG2;
  sprivate  = domain->SprivateG2;
  if ( !sprivate->ASCmat || !sprivate->ASRCmat ||
       sprivate->splnaconstr != naconstr ) {
    domain->error_code = G2H_ERROR_UNDEFINED_CONSTR;
    goto failure;
  }
  nfunc_a = privateG2->nfunc_a;
  nfunc_b = privateG2->nfunc_b;
  nfunc_c = sprivate->nsfunc_c;
  nfunc_d = sprivate->nsfunc_d;
  nbf = nfunc_a+nfunc_c+nfunc_d;
  acmat  = sprivate->ASCmat;
  arcmat = sprivate->ASRCmat;
  lmat = sprivate->SLMat;

  b = pkv_GetScratchMemd ( spdimen*nbf );
  x = pkv_GetScratchMemd ( spdimen*max(nbf,nfunc_b) );
  y = pkv_GetScratchMemd ( spdimen*nbf );
  fc00 = pkv_GetScratchMemd ( (G2_CROSSDEGSUM+6)*2*hole_k*spdimen );
  if ( !b || !x || !y || !fc00 ) {
    domain->error_code = G2H_ERROR_NO_SCRATCH_MEMORY;
    goto failure;
  }

  if ( !_g2h_SetSplRightSided ( domain, spdimen, hole_cp, fc00, x ) )
    goto failure;

  pkn_Block2LowerTrMSolved ( hole_k, nfunc_c/hole_k, nfunc_d/hole_k, nfunc_a,
                             lmat, spdimen, spdimen, x );
  pkn_Block2UpperTrMSolved ( hole_k, nfunc_c/hole_k, nfunc_d/hole_k, nfunc_a,
                             lmat, spdimen, spdimen, x );
              /* now x is the solution without the constraints and it is time */   
              /* to employ the constraints that have been previously set */
  pkv_TransposeMatrixd ( nbf, spdimen, spdimen, x, nbf, y );
  pkn_MultMatrixd ( naconstr, spdimen*nbf, spdimen*nbf, acmat,
                    1, 1, y, 1, b );
  pkn_AddMatrixd ( 1, naconstr, 0, b, 0, constr, 0, b );
              /* solve the system with the Schur matrix */
  pkn_LowerTrMatrixSolved ( naconstr, arcmat, 1, 1, b, 1, b );
  pkn_UpperTrMatrixSolved ( naconstr, arcmat, 1, 1, b, 1, b );
              /* compute the constraints correction */
  pkn_MultTMatrixd ( naconstr, nbf*spdimen, nbf*spdimen, acmat,
                     1, 1, b, 1, y );
  pkv_TransposeMatrixd ( spdimen, nbf, nbf, y, spdimen, b );
  pkn_Block2LowerTrMSolved ( hole_k, nfunc_c/hole_k, nfunc_d/hole_k, nfunc_a,
                             lmat, spdimen, spdimen, b );
  pkn_Block2UpperTrMSolved ( hole_k, nfunc_c/hole_k, nfunc_d/hole_k, nfunc_a,
                             lmat, spdimen, spdimen, b );
  pkn_SubtractMatrixd ( 1, spdimen*nbf, 0, x, 0, b, 0, x );
  if ( acoeff )
    memcpy ( acoeff, x, spdimen*nbf*sizeof(double) );

  if ( !_g2h_OutputSplPatchesd ( domain, spdimen, x, fc00, usrptr, outpatch ) )
    goto failure;

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*g2h_SplFillHoleAltConstrd*/

