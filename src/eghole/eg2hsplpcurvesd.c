
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

#undef CONST_
#define CONST_

#include "eg2holed.h"
#include "eg2hprivated.h"
#include "eg2herror.h"


boolean g2h_GetSplFinalPatchCurvesd ( GHoleDomaind *domain, int spdimen,
                                      CONST_ double *hole_cp, double *acoeff,
                                      void (*outcurve) ( int n, int lkn,
                                          const double *kn, const double *cp ) )
{
#define NCTLP (lastomcknot-G2_CROSS00DEG)
  void          *sp;
  int           hole_k, lastomcknot;
  GHolePrivateRecd   *privateG;
  G2HolePrivateRecd  *privateG2;
  G2HoleSPrivateRecd *sprivate;
  int           nfunc_a, nfunc_b, nfunc_c, nfunc_d, nabc, nk, m1;
  double        *bcp, *cp, *y;
  double        *bc00, *bc01, *bc02, *bc10, *bc11, *bc12,
                *bd00, *bd01, *bd02, *bd10, *bd11, *bd12;
  double        *omcknots, *bezknots;
  unsigned char *bfcpn;
  int           i, j, k, l;

  sp = pkv_GetScratchMemTop ();
  hole_k = domain->hole_k;
  privateG  = domain->privateG;
  privateG2 = domain->privateG2;
  sprivate  = domain->SprivateG2;
  if ( !sprivate )
    goto failure;
  nfunc_a = privateG2->nfunc_a;
  nfunc_b = 6*hole_k+1;
  nfunc_c = sprivate->nsfunc_c;
  nfunc_d = sprivate->nsfunc_d;
  nabc = nfunc_a+nfunc_b+nfunc_c;
  nk = sprivate->nk;
  m1 = sprivate->m1;
  lastomcknot = sprivate->lastomcknot;
  omcknots = sprivate->omcknots;
  bezknots = &sprivate->bezknots[G2H_FINALDEG-G2_CROSS00DEG];
  bfcpn = privateG->bfcpn;
  bcp = pkv_GetScratchMemd ( hole_k*spdimen*NCTLP );
  y = pkv_GetScratchMemd ( hole_k*spdimen*NCTLP );
  if ( !bcp || !y )
    goto failure;

  memset ( bcp, 0, hole_k*spdimen*NCTLP*sizeof(double) );
  for ( j = k = 0; j < nfunc_b; j++ ) {
    l = bfcpn[j];
    cp = &hole_cp[l*spdimen];
    for ( i = k = 0;  i < hole_k;  i++, k += spdimen ) {
      _g2h_GetBFBPatchCurvesd ( domain, j, i, &bc00, &bc01, &bc02, &bc10, &bc11, &bc12,
                                &bd00, &bd01, &bd02, &bd10, &bd11, &bd12 );
      pkn_MultMatrixd ( G2_CROSS00DEG+1, 1, 1, bc00, spdimen, 0, cp, spdimen, y );
      pkn_AddMatrixd ( 1, spdimen*(G2_CROSS00DEG+1), 0, &bcp[k*NCTLP],
                       0, y, 0, &bcp[k*NCTLP] );
    }
  }
  for ( j = 0; j < nfunc_a; j++ ) {
    cp = &acoeff[(nfunc_c+nfunc_d+j)*spdimen];
    for ( i = k = 0;  i < hole_k;  i++, k += spdimen ) {
      _g2h_GetBFAPatchCurvesd ( domain, j, i, &bc00, &bc01, &bc02, &bd00, &bd01, &bd02 );
      pkn_MultMatrixd ( G2_CROSS00DEG+1, 1, 1, bc00, spdimen, 0, cp, spdimen, y );
      pkn_SubtractMatrixd ( 1, spdimen*(G2_CROSS00DEG+1), 0, &bcp[k*NCTLP],
                            0, y, 0, &bcp[k*NCTLP] );
    }
  }
        /* convert the Bezier representation to the final B-spline */
  mbs_multiOsloInsertKnotsd ( hole_k, spdimen, G2_CROSS00DEG,
                              2*G2_CROSS00DEG+1, bezknots, spdimen*NCTLP, bcp,
                              lastomcknot, omcknots, spdimen*NCTLP, y );
  memcpy ( bcp, y, hole_k*spdimen*NCTLP*sizeof(double) );
        /* add the linear combination of the block D functions */
  for ( i = 0; i < hole_k; i++ ) {
    for ( j = 0;  j < nk*m1;  j++ ) {
      cp = &acoeff[(nfunc_c+(3*i*nk*m1+j))*spdimen];
      _g2h_GetSplDBasisCrossDerd ( domain, nabc+3*(i*nk*m1+j), i,
              y, &y[NCTLP], &y[NCTLP], &y[NCTLP], &y[NCTLP] );
      pkn_MultMatrixd ( NCTLP, 1, 1, y, spdimen, 0, cp, spdimen, &y[NCTLP] );
      pkn_SubtractMatrixd ( 1, spdimen*NCTLP, 0, &bcp[k*NCTLP],
                            0, y, 0, &bcp[k*NCTLP] );
    }
  }

  for ( i = k = 0;  i < hole_k;  i++, k+= spdimen )
    outcurve ( G2_CROSS00DEG, lastomcknot, omcknots, &bcp[k*NCTLP] );

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
#undef NCTLP
} /*g2h_GetSplFinalPatchCurvesd*/

boolean g2h_GetSplABasisFPatchCurved ( GHoleDomaind *domain, int fn, int i,
                                       int *lkn, double *kn, double *bfpc )
{
  void *sp;
  G2HolePrivateRecd  *privateG2;
  G2HoleSPrivateRecd *sprivate;
  int   hole_k, nfunc_a, lastomcknot;
  double *bc00, *bc01, *bc02, *bd00, *bd01, *bd02, *bezknots;

  sp = pkv_GetScratchMemTop ();
  hole_k = domain->hole_k;
  if ( i < 0 || i >= hole_k )
    goto failure;
  privateG2 = domain->privateG2;
  sprivate = domain->SprivateG2;
  if ( !sprivate )
    goto failure;
  nfunc_a = privateG2->nfunc_a;
  if ( fn < 0 || fn >= nfunc_a )
    goto failure;
  _g2h_GetBFAPatchCurvesd ( domain, fn, i, &bc00, &bc01, &bc02, &bd00, &bd01, &bd02 );
  *lkn = lastomcknot = sprivate->lastomcknot;
  memcpy ( kn, sprivate->omcknots, (lastomcknot+1)*sizeof(double) );
  bezknots = &sprivate->bezknots[G2H_FINALDEG-G2_CROSS00DEG];
  mbs_multiOsloInsertKnotsd ( 1, 1, G2_CROSS00DEG,
                              2*G2_CROSS00DEG+1, bezknots, 0, bc00,
                              lastomcknot, kn, 0, bfpc );
  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*g2h_GetSplABasisFPatchCurved*/

boolean g2h_GetSplBBasisFPatchCurved ( GHoleDomaind *domain, int fn, int i,
                                       int *lkn, double *kn, double *bfpc )
{
  void  *sp;
  G2HoleSPrivateRecd *sprivate;
  int   hole_k, nfunc_b, lastomcknot;
  double *bc00, *bc01, *bc02, *bc10, *bc11, *bc12,
        *bd00, *bd01, *bd02, *bd10, *bd11, *bd12, *bezknots;

  sp = pkv_GetScratchMemTop ();
  hole_k = domain->hole_k;
  if ( i < 0 || i >= hole_k )
    goto failure;
  nfunc_b = 6*hole_k+1;
  if ( fn < 0 || fn >= nfunc_b )
    goto failure;
  sprivate = domain->SprivateG2;
  if ( !sprivate )
    goto failure;
  _g2h_GetBFBPatchCurvesd ( domain, fn, i, &bc00, &bc01, &bc02, &bc10, &bc11, &bc12,
                            &bd00, &bd01, &bd02, &bd10, &bd11, &bd12 );
  *lkn = lastomcknot = sprivate->lastomcknot;
  memcpy ( kn, sprivate->omcknots, (lastomcknot+1)*sizeof(double) );
  bezknots = &sprivate->bezknots[G2H_FINALDEG-G2_CROSS00DEG];
  mbs_multiOsloInsertKnotsd ( 1, 1, G2_CROSS00DEG,
                              2*G2_CROSS00DEG+1, bezknots, 0, bc00,
                              lastomcknot, kn, 0, bfpc );
  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*g2h_GetSplBBasisFPatchCurved*/

boolean g2h_GetSplDBasisFPatchCurved ( GHoleDomaind *domain, int fn, int i,
                                       int *lkn, double *kn, double *bfpc )
{
  void  *sp;
  G2HolePrivateRecd  *privateG2;
  G2HoleSPrivateRecd *sprivate;
  int   hole_k, nfunc_a, nfunc_b, nfunc_c, nfunc_d, nabc;
  int   lastomcknot, lastpvknot, lastpvvknot;
  double *pu, *pv, *puu, *pvv;

  sp = pkv_GetScratchMemTop ();
  hole_k = domain->hole_k;
  if ( i < 0 || i >= hole_k )
    goto failure;
  privateG2 = domain->privateG2;
  sprivate = domain->SprivateG2;
  if ( !sprivate )
    goto failure;
  nfunc_a = privateG2->nfunc_a;
  nfunc_b = 6*hole_k+1;
  nfunc_c = sprivate->nsfunc_c;
  nfunc_d = sprivate->nsfunc_d;
  if ( fn < 0 || fn >= nfunc_d )
    goto failure;
  lastpvknot = sprivate->lastpvknot;
  lastpvvknot = sprivate->lastpvvknot;
  pu = pkv_GetScratchMemd ( 2*(lastpvknot-G2_CROSS01DEG+lastpvvknot-G2_CROSS02DEG) );
  if ( !pu )
    goto failure;
  pv = &pu[lastpvknot-G2_CROSS01DEG];  puu = &pv[lastpvknot-G2_CROSS01DEG];
  pvv = &puu[lastpvvknot-G2_CROSS02DEG];
  nabc = nfunc_a+nfunc_b+nfunc_c;
  *lkn = lastomcknot = sprivate->lastomcknot;
  memcpy ( kn, sprivate->omcknots, (lastomcknot+1)*sizeof(double) );
  _g2h_GetSplDBasisCrossDerd ( domain, nabc+fn, i, bfpc, pv, pvv, pu, puu );

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*g2h_GetSplDBasisFPatchCurved*/

