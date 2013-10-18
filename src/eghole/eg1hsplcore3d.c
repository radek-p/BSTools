
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

#undef CONST_
#define CONST_

#include "eg1holed.h"
#include "eg1hprivated.h"
#include "eg1herror.h"

/* ///////////////////////////////////////////////////////////////////////// */
boolean g1h_DecomposeSplMatrixd ( GHoleDomaind *domain )
{
  void   *sp;
  G1HolePrivateRecd  *privateG1;
  G1HoleSPrivateRecd *sprivate;
  int    hole_k, nfunc_a, nfunc_c, nfunc_d;
  int    size;
  double *lmat;

  sp = pkv_GetScratchMemTop ();
  privateG1 = domain->privateG1;
  sprivate = domain->SprivateG1;
  if ( !sprivate )
    goto failure;
  if ( !sprivate->SLMat ) {
    if ( !sprivate->SAMat )
      if ( !g1h_ComputeSplFormMatrixd ( domain ) )
        goto failure;
    hole_k = domain->hole_k;
    nfunc_a = privateG1->nfunc_a;
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
      domain->error_code = G1H_ERROR_NONPOSITIVE_MATRIX;
      goto failure;
    }
  }

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*g1h_DecomposeSplMatrixd*/

boolean _g1h_SetSplRightSided ( GHoleDomaind *domain,
                                int spdimen, CONST_ double *hole_cp,
                                const double *bmat,
                                double *fc00, double *b )
{
  void   *sp;
  GHolePrivateRecd   *privateG;
  G1HolePrivateRecd  *privateG1;
  G1HoleSPrivateRecd *sprivate;
  int    hole_k, nfunc_a, nfunc_b, nfunc_c, nfunc_d, nbf;
  double *x, *y, *cp;
  unsigned char *bfcpn;
  double *bc00, *bc01, *bc10, *bc11, *bd00, *bd01, *bd10, *bd11;
  double        *fc01, *fc10, *fc11, *fd00, *fd01, *fd10, *fd11;
  int    i, j, k, l;

  sp = pkv_GetScratchMemTop ();
  hole_k = domain->hole_k;
  privateG  = domain->privateG;
  privateG1 = domain->privateG1;
  sprivate = domain->SprivateG1;
  if ( !sprivate )
    goto failure;
  nfunc_a = privateG1->nfunc_a;
  nfunc_b = privateG1->nfunc_b;
  nfunc_c = sprivate->nsfunc_c;
  nfunc_d = sprivate->nsfunc_d;
  nbf = nfunc_a+nfunc_c+nfunc_d;
  bfcpn = privateG->bfcpn;
  x = pkv_GetScratchMemd ( spdimen*nfunc_b );
  y = pkv_GetScratchMemd ( spdimen*(G1H_FINALDEG+1)*(G1H_FINALDEG+1) );
  if ( !x || !y )
    goto failure;

  memset ( fc00, 0, (G1_CROSSDEGSUM+4)*2*hole_k*spdimen*sizeof(double) );
  G1GetFCAddresses ();

  for ( j = 0; j < nfunc_b; j++ ) {
    l = bfcpn[j];
    cp = &hole_cp[l*spdimen];
    memcpy  ( &x[j*spdimen], cp, spdimen*sizeof(double) );
    for ( i = k = 0;  i < hole_k;  i++, k += spdimen ) {
      _g1h_GetBFBPatchCurvesd ( domain, j, i,
                                &bc00, &bc01, &bc10, &bc11,
                                &bd00, &bd01, &bd10, &bd11 );

      pkn_MultMatrixd ( G1_CROSS00DEG+1, 1, 1, bc00, spdimen, 0, cp, spdimen, y );
      pkn_AddMatrixd ( 1, spdimen*(G1_CROSS00DEG+1), 0, &fc00[k*(G1_CROSS00DEG+1)],
                       0, y, 0, &fc00[k*(G1_CROSS00DEG+1)] );
      pkn_MultMatrixd ( G1_CROSS01DEG+1, 1, 1, bc01, spdimen, 0, cp, spdimen, y );
      pkn_AddMatrixd ( 1, spdimen*(G1_CROSS01DEG+1), 0, &fc01[k*(G1_CROSS01DEG+1)],
                       0, y, 0, &fc01[k*(G1_CROSS01DEG+1)] );
      pkn_MultMatrixd ( G1_CROSS10DEG+1, 1, 1, bc10, spdimen, 0, cp, spdimen, y );
      pkn_AddMatrixd ( 1, spdimen*(G1_CROSS10DEG+1), 0, &fc10[k*(G1_CROSS10DEG+1)],
                       0, y, 0, &fc10[k*(G1_CROSS10DEG+1)] );
      pkn_MultMatrixd ( G1_CROSS11DEG+1, 1, 1, bc11, spdimen, 0, cp, spdimen, y );
      pkn_AddMatrixd ( 1, spdimen*(G1_CROSS11DEG+1), 0, &fc11[k*(G1_CROSS11DEG+1)],
                       0, y, 0, &fc11[k*(G1_CROSS11DEG+1)] );

      pkn_MultMatrixd ( G1_CROSS00DEG+1, 1, 1, bd00, spdimen, 0, cp, spdimen, y );
      pkn_AddMatrixd ( 1, spdimen*(G1_CROSS00DEG+1), 0, &fd00[k*(G1_CROSS00DEG+1)],
                       0, y, 0, &fd00[k*(G1_CROSS00DEG+1)] );
      pkn_MultMatrixd ( G1_CROSS01DEG+1, 1, 1, bd01, spdimen, 0, cp, spdimen, y );
      pkn_AddMatrixd ( 1, spdimen*(G1_CROSS01DEG+1), 0, &fd01[k*(G1_CROSS01DEG+1)],
                       0, y, 0, &fd01[k*(G1_CROSS01DEG+1)] );
      pkn_MultMatrixd ( G1_CROSS10DEG+1, 1, 1, bd10, spdimen, 0, cp, spdimen, y );
      pkn_AddMatrixd ( 1, spdimen*(G1_CROSS10DEG+1), 0, &fd10[k*(G1_CROSS10DEG+1)],
                       0, y, 0, &fd10[k*(G1_CROSS10DEG+1)] );
      pkn_MultMatrixd ( G1_CROSS11DEG+1, 1, 1, bd11, spdimen, 0, cp, spdimen, y );
      pkn_AddMatrixd ( 1, spdimen*(G1_CROSS11DEG+1), 0, &fd11[k*(G1_CROSS11DEG+1)],
                       0, y, 0, &fd11[k*(G1_CROSS11DEG+1)] );
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
} /*_g1h_SetSplRightSided*/

boolean _g1h_OutputSplPatchesd ( GHoleDomaind *domain,
                int spdimen, CONST_ double *x, double *fc00, void *usrptr,
                void (*outpatch) ( int n, int lknu, const double *knu,
                                   int m, int lknv, const double *knv,
                                   const double *cp, void *usrptr ) )
{
  void   *sp;
  G1HolePrivateRecd  *privateG1;
  G1HoleSPrivateRecd *sprivate;
  double *bc00, *bc01, *bd00, *bd01;
  double        *fc01, *fc10, *fc11,
         *fd00, *fd01, *fd10, *fd11;
  double *sfc00, *sfc01, *sfd00, *sfd01;
  double *fcomc, *pv, *pu;
  double *y, *yv, *yu, *z, *xx, *cp;
  int    hole_k, nk, m1, m2, csiz, csize, dsize;
  int    nfunc_a, nfunc_b, nfunc_c, nfunc_d, nabc;
  int    degu, degv;
  int    lastomcknot, lastpvknot, lastcknot, lastfpknot;
  double *omcknots, *pvknots, *bezknots, *cknots, *fpknots, *aknots;
  int    i, j, k, l, m;


  sp = pkv_GetScratchMemTop ();

  hole_k = domain->hole_k;
  privateG1 = domain->privateG1;
  sprivate = domain->SprivateG1;
  nk = sprivate->nk;
  m1 = sprivate->m1;
  m2 = sprivate->m2;
  csize = sprivate->csize;
  csiz = G1H_FINALDEG-3+nk*m2;  /* this is sqrt(csize) */
  dsize = sprivate->dsize;
  nfunc_a = privateG1->nfunc_a;
  nfunc_b = privateG1->nfunc_b;
  nfunc_c = sprivate->nsfunc_c;
  nfunc_d = sprivate->nsfunc_d;
  nabc = nfunc_a+nfunc_b+nfunc_c;
  lastomcknot = sprivate->lastomcknot;
  lastpvknot  = sprivate->lastpvknot;
  lastcknot   = sprivate->lastcknot;
  lastfpknot  = sprivate->lastfpknot;
  bezknots = sprivate->bezknots;
  omcknots = sprivate->omcknots;
  pvknots  = sprivate->pvknots;
  cknots   = sprivate->cknots;
  fpknots  = sprivate->fpknots;

  aknots = pkv_GetScratchMemd ( 2*(lastfpknot+1) );
  sfc00 = pkv_GetScratchMemd ( (lastomcknot-G1_CROSS00DEG+lastpvknot-G1_CROSS01DEG)*
                                spdimen*2*hole_k );
  y = pkv_GetScratchMemd (
          3*spdimen*(lastfpknot-G1H_FINALDEG)*(lastfpknot-G1H_FINALDEG) );
  if ( !aknots || !sfc00 || !y )
    goto failure;

  G1GetFCAddresses ();
  G1GetSFCAddresses ();

  for ( j = 0; j < nfunc_a; j++ ) {
    cp = &x[(nfunc_c+nfunc_d+j)*spdimen];
    for ( i = k = 0;  i < hole_k;  i++, k += spdimen ) {
      _g1h_GetBFAPatchCurvesd ( domain, j, i, &bc00, &bc01, &bd00, &bd01 );

      pkn_MultMatrixd ( G1_CROSS00DEG+1, 1, 1, bc00, spdimen, 0, cp, spdimen, y );
      pkn_SubtractMatrixd ( 1, spdimen*(G1_CROSS00DEG+1), 0, &fc00[k*(G1_CROSS00DEG+1)],
                            0, y, 0, &fc00[k*(G1_CROSS00DEG+1)] );
      pkn_MultMatrixd ( G1_CROSS01DEG+1, 1, 1, bc01, spdimen, 0, cp, spdimen, y );
      pkn_SubtractMatrixd ( 1, spdimen*(G1_CROSS01DEG+1), 0, &fc01[k*(G1_CROSS01DEG+1)],
                            0, y, 0, &fc01[k*(G1_CROSS01DEG+1)] );

      pkn_MultMatrixd ( G1_CROSS00DEG+1, 1, 1, bd00, spdimen, 0, cp, spdimen, y );
      pkn_SubtractMatrixd ( 1, spdimen*(G1_CROSS00DEG+1), 0, &fd00[k*(G1_CROSS00DEG+1)],
                            0, y, 0, &fd00[k*(G1_CROSS00DEG+1)] );
      pkn_MultMatrixd ( G1_CROSS01DEG+1, 1, 1, bd01, spdimen, 0, cp, spdimen, y );
      pkn_SubtractMatrixd ( 1, spdimen*(G1_CROSS01DEG+1), 0, &fd01[k*(G1_CROSS01DEG+1)],
                            0, y, 0, &fd01[k*(G1_CROSS01DEG+1)] );
    }
  }

        /* At this point the term made of the block A and B basis functions */
        /* is represented in the Coons form, by the Bezier curves. The next */
        /* step is to obtain the B-spline representation of the curves and  */
        /* add the term represented by the block D functions.               */

  yv = &y[(lastfpknot-G1H_FINALDEG)*spdimen];
  yu = &yv[(lastpvknot-G1_CROSS01DEG)*spdimen];
  fcomc = &yu[(lastpvknot-G1_CROSS01DEG)*spdimen];
  pv = &fcomc[lastomcknot-G1_CROSS00DEG];
  pu = &pv[lastpvknot-G1_CROSS01DEG];
  z = &pu[lastpvknot-G1_CROSS01DEG];
  for ( i = 0;  i < hole_k;  i++ ) {
    l = (i+hole_k-1) % hole_k;
    memset ( y, 0, (lastfpknot-G1H_FINALDEG+
             2*(lastpvknot-G1_CROSS01DEG))*spdimen*sizeof(double) );
    pkv_Selectd ( nk*m1, spdimen, 2*spdimen, spdimen,
                  &x[(nfunc_c+i*dsize)*spdimen], &y[3*spdimen] );
    if ( !mbs_multiSubtractBSCurvesd ( 1, spdimen, G1_CROSS00DEG, 2*G1_CROSS00DEG+1,
            &bezknots[G1H_FINALDEG-G1_CROSS00DEG], 0,
            &fc00[i*(G1_CROSS00DEG+1)*spdimen],
            G1_CROSS00DEG, lastomcknot, omcknots, 0, y, &degu, &m, z, 0,
            &sfc00[i*(lastomcknot-G1_CROSS00DEG)*spdimen] ) )
      goto failure;
    memcpy ( &sfd00[l*(lastomcknot-G1_CROSS00DEG)*spdimen],
             &sfc00[i*(lastomcknot-G1_CROSS00DEG)*spdimen],
             (lastomcknot-G1_CROSS00DEG)*spdimen*sizeof(double) );

    for ( j = 0; j < nk*m1; j++ ) {
      xx = &x[(nfunc_c+i*dsize+2*j)*spdimen];

      _g1h_GetSplDBasisCrossDerd ( domain, nabc+i*dsize+2*j, i, fcomc, pv, pu );
      pkn_MultMatrixAddd ( lastpvknot-G1_CROSS01DEG, 1, 1, pv,
                           spdimen, 0, xx, spdimen, yv );
      pkn_MultMatrixAddd ( lastpvknot-G1_CROSS01DEG, 1, 1, pu,
                           spdimen, 0, xx, spdimen, yu );

      _g1h_GetSplDBasisCrossDerd ( domain, nabc+i*dsize+2*j+1, i, fcomc, pv, pu );
      pkn_MultMatrixAddd ( lastpvknot-G1_CROSS01DEG, 1, 1, pv,
                           spdimen, 0, &xx[spdimen], spdimen, yv );
      pkn_MultMatrixAddd ( lastpvknot-G1_CROSS01DEG, 1, 1, pu,
                           spdimen, 0, &xx[spdimen], spdimen, yu );
    }
    if ( !mbs_multiSubtractBSCurvesd ( 1, spdimen, G1_CROSS01DEG, 2*G1_CROSS01DEG+1,
            &bezknots[G1H_FINALDEG-G1_CROSS01DEG], 0, &fc01[i*(G1_CROSS01DEG+1)*spdimen],
            G1_CROSS01DEG, lastpvknot, pvknots, 0, yv, &degu, &m, z, 0,
            &sfc01[i*(lastpvknot-G1_CROSS01DEG)*spdimen] ) )
      goto failure;
    if ( !mbs_multiSubtractBSCurvesd ( 1, spdimen, G1_CROSS01DEG, 2*G1_CROSS01DEG+1,
            &bezknots[G1H_FINALDEG-G1_CROSS01DEG], 0, &fd01[l*(G1_CROSS01DEG+1)*spdimen],
            G1_CROSS01DEG, lastpvknot, pvknots, 0, yu, &degv, &m, z, 0,
            &sfd01[l*(lastpvknot-G1_CROSS01DEG)*spdimen] ) )
      goto failure;
  }

        /* The next step is to find the B-spline representation of the     */
        /* Coons patch and add the terms represented by the block C basis  */
        /* functions */

  yu = &y[spdimen*(lastfpknot-G1H_FINALDEG)*(lastfpknot-G1H_FINALDEG)];
  yv = &y[spdimen*(lastfpknot-G1H_FINALDEG)*(lastfpknot-G1H_FINALDEG)];
  for ( i = 0; i < hole_k; i++ ) {
    mbs_BSC1CoonsToBSd ( spdimen,
            G1_CROSS00DEG, lastomcknot, omcknots,
            &sfc00[i*(lastomcknot-G1_CROSS00DEG)*spdimen],
            G1_CROSS01DEG, lastpvknot,  pvknots,
            &sfc01[i*(lastpvknot-G1_CROSS01DEG)*spdimen],
            G1_CROSS10DEG, 2*G1_CROSS10DEG+1, &bezknots[G1H_FINALDEG-G1_CROSS10DEG],
            &fc10[i*(G1_CROSS10DEG+1)*spdimen],
            G1_CROSS11DEG, 2*G1_CROSS11DEG+1, &bezknots[G1H_FINALDEG-G1_CROSS11DEG],
            &fc11[i*(G1_CROSS11DEG+1)*spdimen],
            G1_CROSS00DEG, lastomcknot, omcknots,
            &sfd00[i*(lastomcknot-G1_CROSS00DEG)*spdimen],
            G1_CROSS01DEG, lastpvknot,  pvknots,
            &sfd01[i*(lastpvknot-G1_CROSS01DEG)*spdimen],
            G1_CROSS10DEG, 2*G1_CROSS10DEG+1, &bezknots[G1H_FINALDEG-G1_CROSS10DEG],
            &fd10[i*(G1_CROSS10DEG+1)*spdimen],
            G1_CROSS11DEG, 2*G1_CROSS11DEG+1, &bezknots[G1H_FINALDEG-G1_CROSS11DEG],
            &fd11[i*(G1_CROSS11DEG+1)*spdimen],
            &degu, &l, aknots, &degv, &m, &aknots[lastfpknot+1], y );

        /* convert A+B+D to the final form */
    mbs_multiAdjustBSCRepd ( 1, (m-degv)*spdimen, degu, l, aknots, 0, y,
                             G1H_FINALDEG, lastfpknot, fpknots, 0, yu );
    mbs_multiAdjustBSCRepd ( lastfpknot-G1H_FINALDEG, spdimen, degv, m,
                             &aknots[lastfpknot+1], (m-degv)*spdimen, yu,
                             G1H_FINALDEG, lastfpknot, fpknots,
                             (lastfpknot-G1H_FINALDEG)*spdimen, y );

        /* now process the terms with the block C basis functions */
    memset ( yu, 0, spdimen*(lastcknot-G1H_FINALDEG)*(lastcknot-G1H_FINALDEG)*
                    sizeof(double) );
    pkv_Selectd ( csiz, csiz*spdimen, csiz*spdimen,
                  (lastcknot-G1H_FINALDEG)*spdimen, &x[i*csize*spdimen],
                  &yu[(2*(lastcknot-G1H_FINALDEG)+2)*spdimen] );

    mbs_multiAdjustBSCRepd ( 1, (lastcknot-G1H_FINALDEG)*spdimen, G1H_FINALDEG,
                             lastcknot, cknots, 0, yu,
                             G1H_FINALDEG, lastfpknot, fpknots, 0, yv );
    mbs_multiAdjustBSCRepd ( lastfpknot-G1H_FINALDEG, spdimen, G1H_FINALDEG,
                             lastcknot, cknots, (lastcknot-G1H_FINALDEG)*spdimen,
                             yv, G1H_FINALDEG, lastfpknot, fpknots,
                             (lastfpknot-G1H_FINALDEG)*spdimen, yu );
    pkn_SubtractMatrixd ( lastfpknot-G1H_FINALDEG-4,
                          (lastfpknot-G1H_FINALDEG-4)*spdimen,
                          (lastfpknot-G1H_FINALDEG)*spdimen,
                          &y[(2*(lastfpknot-G1H_FINALDEG)+2)*spdimen],
                          (lastfpknot-G1H_FINALDEG)*spdimen,
                          &yu[(2*(lastfpknot-G1H_FINALDEG)+2)*spdimen],
                          (lastfpknot-G1H_FINALDEG)*spdimen,
                          &y[(2*(lastfpknot-G1H_FINALDEG)+2)*spdimen] );

    outpatch ( G1H_FINALDEG, lastfpknot, fpknots,
               G1H_FINALDEG, lastfpknot, fpknots, y, usrptr );
  }

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*_g1h_OutputSplPatchesd*/

/* ///////////////////////////////////////////////////////////////////////// */
boolean g1h_SplFillHoled ( GHoleDomaind *domain,
               int spdimen, CONST_ double *hole_cp,
               double *acoeff, void *usrptr,
               void (*outpatch) ( int n, int lknu, const double *knu,
                                  int m, int lknv, const double *knv,
                                  const double *cp, void *usrptr ) )
{ 
  void   *sp;
  G1HolePrivateRecd  *privateG1;
  G1HoleSPrivateRecd *sprivate;
  int    hole_k, nfunc_a, nfunc_c, nfunc_d, nbf;
  double *bmat, *lmat, *x, *fc00;

  sp = pkv_GetScratchMemTop ();
  hole_k = domain->hole_k;
  privateG1 = domain->privateG1;
  sprivate = domain->SprivateG1;
  if ( !sprivate )
    goto failure;
  nfunc_a = privateG1->nfunc_a;
  nfunc_c = sprivate->nsfunc_c;
  nfunc_d = sprivate->nsfunc_d;
  nbf = nfunc_a+nfunc_c+nfunc_d;
  if ( !sprivate->SLMat )
    if ( !g1h_DecomposeSplMatrixd ( domain ) )
      goto failure;
  lmat = sprivate->SLMat;
  bmat = sprivate->SBMat;

  x = pkv_GetScratchMemd ( spdimen*nbf );
  fc00 = pkv_GetScratchMemd ( (G1_CROSSDEGSUM+4)*2*hole_k*spdimen );
  if ( !x || !fc00 )
    goto failure;

  if ( !_g1h_SetSplRightSided ( domain, spdimen, hole_cp, bmat, fc00, x ) )
    goto failure;

  pkn_Block2LowerTrMSolved ( hole_k, nfunc_c/hole_k, nfunc_d/hole_k, nfunc_a,
                             lmat, spdimen, spdimen, x );
  pkn_Block2UpperTrMSolved ( hole_k, nfunc_c/hole_k, nfunc_d/hole_k, nfunc_a,
                             lmat, spdimen, spdimen, x );
  if ( acoeff )
    memcpy ( acoeff, x, spdimen*nbf*sizeof(double) );

  if ( !_g1h_OutputSplPatchesd ( domain, spdimen, x, fc00, usrptr, outpatch ) )
    goto failure;

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*g1h_SplFillHoled*/

/* ///////////////////////////////////////////////////////////////////////// */
boolean g1h_SplFillHoleConstrd ( GHoleDomaind *domain,
               int spdimen, CONST_ double *hole_cp,
               int nconstr, CONST_ double *constr,
               double *acoeff, void *usrptr,
               void (*outpatch) ( int n, int lknu, const double *knu,
                                  int m, int lknv, const double *knv,
                                  const double *cp, void *usrptr ) )
{
  void   *sp;
  G1HolePrivateRecd  *privateG1;
  G1HoleSPrivateRecd *sprivate;
  double *fc00, *b, *x, *y, *cmat, *rcmat, *lmat, *bmat;
  int    hole_k, nfunc_a, nfunc_c, nfunc_d, nbf;

  sp = pkv_GetScratchMemTop ();

  hole_k = domain->hole_k;
  privateG1 = domain->privateG1;
  sprivate = domain->SprivateG1;
  if ( !sprivate->SCmat || !sprivate->SRCmat ||
       sprivate->splnconstr != nconstr ) {
    domain->error_code = G1H_ERROR_UNDEFINED_CONSTR;
    goto failure;
  }
  nfunc_a = privateG1->nfunc_a;
  nfunc_c = sprivate->nsfunc_c;
  nfunc_d = sprivate->nsfunc_d;
  nbf = nfunc_a+nfunc_c+nfunc_d;
  cmat = sprivate->SCmat;
  rcmat = sprivate->SRCmat;
  lmat = sprivate->SLMat;
  bmat = sprivate->SBMat;

  b = pkv_GetScratchMemd ( spdimen*nbf );
  x = pkv_GetScratchMemd ( spdimen*nbf );
  fc00 = pkv_GetScratchMemd ( (G1_CROSSDEGSUM+4)*2*hole_k*spdimen );
  if ( !b || !x || !fc00 )
    goto failure;

  if ( !_g1h_SetSplRightSided ( domain, spdimen, hole_cp, bmat, fc00, x ) )
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
    domain->error_code = G1H_ERROR_NO_SCRATCH_MEMORY;
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

  if ( !_g1h_OutputSplPatchesd ( domain, spdimen, x, fc00, usrptr, outpatch ) )
    goto failure;

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*g1h_SplFillHoleConstrd*/

/* ///////////////////////////////////////////////////////////////////////// */
boolean g1h_SplFillHoleAltConstrd ( GHoleDomaind *domain,
               int spdimen, CONST_ double *hole_cp,
               int naconstr, CONST_ double *constr,
               double *acoeff, void *usrptr,
               void (*outpatch) ( int n, int lknu, const double *knu,
                                  int m, int lknv, const double *knv,
                                  const double *cp, void *usrptr ) )
{
  void   *sp;
  G1HolePrivateRecd  *privateG1;
  G1HoleSPrivateRecd *sprivate;
  double *fc00, *b, *x, *y, *acmat, *arcmat, *lmat, *bmat;
  int    hole_k, nfunc_a, nfunc_b, nfunc_c, nfunc_d, nbf;

  sp = pkv_GetScratchMemTop ();

  hole_k = domain->hole_k;
  privateG1 = domain->privateG1;
  sprivate  = domain->SprivateG1;
  if ( !sprivate->ASCmat || !sprivate->ASRCmat ||
       sprivate->splnaconstr != naconstr ) {
    domain->error_code = G1H_ERROR_UNDEFINED_CONSTR;
    goto failure;
  }
  nfunc_a = privateG1->nfunc_a;
  nfunc_b = privateG1->nfunc_b;
  nfunc_c = sprivate->nsfunc_c;
  nfunc_d = sprivate->nsfunc_d;
  nbf = nfunc_a+nfunc_c+nfunc_d;
  acmat  = sprivate->ASCmat;
  arcmat = sprivate->ASRCmat;
  lmat = sprivate->SLMat;
  bmat = sprivate->SBMat;

  b = pkv_GetScratchMemd ( spdimen*nbf );
  x = pkv_GetScratchMemd ( spdimen*max(nbf,nfunc_b) );
  y = pkv_GetScratchMemd ( spdimen*nbf );
  fc00 = pkv_GetScratchMemd ( (G1_CROSSDEGSUM+4)*2*hole_k*spdimen );
  if ( !b || !x || !y || !fc00 ) {
    domain->error_code = G1H_ERROR_NO_SCRATCH_MEMORY;
    goto failure;
  }

  if ( !_g1h_SetSplRightSided ( domain, spdimen, hole_cp, bmat, fc00, x ) )
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

  if ( !_g1h_OutputSplPatchesd ( domain, spdimen, x, fc00, usrptr, outpatch ) )
    goto failure;

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*g1h_SplFillHoleAltConstrd*/

