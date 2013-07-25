
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

#include "eg2holef.h"
#include "eg2hprivatef.h"
#include "eg2herror.h"


boolean g2h_GetSplFinalPatchCurvesf ( GHoleDomainf *domain, int spdimen,
                                      CONST_ float *hole_cp, float *acoeff,
                                      void (*outcurve) ( int n, int lkn,
                                          const float *kn, const float *cp ) )
{
#define NCTLP (lastomcknot-G2_CROSS00DEG)
  void          *sp;
  int           hole_k, lastomcknot;
  GHolePrivateRecf   *privateG;
  G2HolePrivateRecf  *privateG2;
  G2HoleSPrivateRecf *sprivate;
  int           nfunc_a, nfunc_b, nfunc_c, nfunc_d, nabc, nk, m1;
  float         *bcp, *cp, *y;
  float         *bc00, *bc01, *bc02, *bc10, *bc11, *bc12,
                *bd00, *bd01, *bd02, *bd10, *bd11, *bd12;
  float         *omcknots, *bezknots;
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
  bcp = pkv_GetScratchMemf ( hole_k*spdimen*NCTLP );
  y = pkv_GetScratchMemf ( hole_k*spdimen*NCTLP );
  if ( !bcp || !y )
    goto failure;

  memset ( bcp, 0, hole_k*spdimen*NCTLP*sizeof(float) );
  for ( j = k = 0; j < nfunc_b; j++ ) {
    l = bfcpn[j];
    cp = &hole_cp[l*spdimen];
    for ( i = k = 0;  i < hole_k;  i++, k += spdimen ) {
      _g2h_GetBFBPatchCurvesf ( domain, j, i, &bc00, &bc01, &bc02, &bc10, &bc11, &bc12,
                                &bd00, &bd01, &bd02, &bd10, &bd11, &bd12 );
      pkn_MultMatrixf ( G2_CROSS00DEG+1, 1, 1, bc00, spdimen, 0, cp, spdimen, y );
      pkn_AddMatrixf ( 1, spdimen*(G2_CROSS00DEG+1), 0, &bcp[k*NCTLP],
                       0, y, 0, &bcp[k*NCTLP] );
    }
  }
  for ( j = 0; j < nfunc_a; j++ ) {
    cp = &acoeff[(nfunc_c+nfunc_d+j)*spdimen];
    for ( i = k = 0;  i < hole_k;  i++, k += spdimen ) {
      _g2h_GetBFAPatchCurvesf ( domain, j, i, &bc00, &bc01, &bc02, &bd00, &bd01, &bd02 );
      pkn_MultMatrixf ( G2_CROSS00DEG+1, 1, 1, bc00, spdimen, 0, cp, spdimen, y );
      pkn_SubtractMatrixf ( 1, spdimen*(G2_CROSS00DEG+1), 0, &bcp[k*NCTLP],
                            0, y, 0, &bcp[k*NCTLP] );
    }
  }
        /* convert the Bezier representation to the final B-spline */
  mbs_multiOsloInsertKnotsf ( hole_k, spdimen, G2_CROSS00DEG,
                              2*G2_CROSS00DEG+1, bezknots, spdimen*NCTLP, bcp,
                              lastomcknot, omcknots, spdimen*NCTLP, y );
  memcpy ( bcp, y, hole_k*spdimen*NCTLP*sizeof(float) );
        /* add the linear combination of the block D functions */
  for ( i = 0; i < hole_k; i++ ) {
    for ( j = 0;  j < nk*m1;  j++ ) {
      cp = &acoeff[(nfunc_c+(3*i*nk*m1+j))*spdimen];
      _g2h_GetSplDBasisCrossDerf ( domain, nabc+3*(i*nk*m1+j), i,
              y, &y[NCTLP], &y[NCTLP], &y[NCTLP], &y[NCTLP] );
      pkn_MultMatrixf ( NCTLP, 1, 1, y, spdimen, 0, cp, spdimen, &y[NCTLP] );
      pkn_SubtractMatrixf ( 1, spdimen*NCTLP, 0, &bcp[k*NCTLP],
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
} /*g2h_GetSplFinalPatchCurvesf*/

boolean g2h_GetSplABasisFPatchCurvef ( GHoleDomainf *domain, int fn, int i,
                                       int *lkn, float *kn, float *bfpc )
{
  void *sp;
  G2HolePrivateRecf  *privateG2;
  G2HoleSPrivateRecf *sprivate;
  int   hole_k, nfunc_a, lastomcknot;
  float *bc00, *bc01, *bc02, *bd00, *bd01, *bd02, *bezknots;

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
  _g2h_GetBFAPatchCurvesf ( domain, fn, i, &bc00, &bc01, &bc02, &bd00, &bd01, &bd02 );
  *lkn = lastomcknot = sprivate->lastomcknot;
  memcpy ( kn, sprivate->omcknots, (lastomcknot+1)*sizeof(float) );
  bezknots = &sprivate->bezknots[G2H_FINALDEG-G2_CROSS00DEG];
  mbs_multiOsloInsertKnotsf ( 1, 1, G2_CROSS00DEG,
                              2*G2_CROSS00DEG+1, bezknots, 0, bc00,
                              lastomcknot, kn, 0, bfpc );
  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*g2h_GetSplABasisFPatchCurvef*/

boolean g2h_GetSplBBasisFPatchCurvef ( GHoleDomainf *domain, int fn, int i,
                                       int *lkn, float *kn, float *bfpc )
{
  void  *sp;
  G2HoleSPrivateRecf *sprivate;
  int   hole_k, nfunc_b, lastomcknot;
  float *bc00, *bc01, *bc02, *bc10, *bc11, *bc12,
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
  _g2h_GetBFBPatchCurvesf ( domain, fn, i, &bc00, &bc01, &bc02, &bc10, &bc11, &bc12,
                            &bd00, &bd01, &bd02, &bd10, &bd11, &bd12 );
  *lkn = lastomcknot = sprivate->lastomcknot;
  memcpy ( kn, sprivate->omcknots, (lastomcknot+1)*sizeof(float) );
  bezknots = &sprivate->bezknots[G2H_FINALDEG-G2_CROSS00DEG];
  mbs_multiOsloInsertKnotsf ( 1, 1, G2_CROSS00DEG,
                              2*G2_CROSS00DEG+1, bezknots, 0, bc00,
                              lastomcknot, kn, 0, bfpc );
  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*g2h_GetSplBBasisFPatchCurvef*/

boolean g2h_GetSplDBasisFPatchCurvef ( GHoleDomainf *domain, int fn, int i,
                                       int *lkn, float *kn, float *bfpc )
{
  void  *sp;
  G2HolePrivateRecf  *privateG2;
  G2HoleSPrivateRecf *sprivate;
  int   hole_k, nfunc_a, nfunc_b, nfunc_c, nfunc_d, nabc;
  int   lastomcknot, lastpvknot, lastpvvknot;
  float *pu, *pv, *puu, *pvv;

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
  pu = pkv_GetScratchMemf ( 2*(lastpvknot-G2_CROSS01DEG+lastpvvknot-G2_CROSS02DEG) );
  if ( !pu )
    goto failure;
  pv = &pu[lastpvknot-G2_CROSS01DEG];  puu = &pv[lastpvknot-G2_CROSS01DEG];
  pvv = &puu[lastpvvknot-G2_CROSS02DEG];
  nabc = nfunc_a+nfunc_b+nfunc_c;
  *lkn = lastomcknot = sprivate->lastomcknot;
  memcpy ( kn, sprivate->omcknots, (lastomcknot+1)*sizeof(float) );
  _g2h_GetSplDBasisCrossDerf ( domain, nabc+fn, i, bfpc, pv, pvv, pu, puu );

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*g2h_GetSplDBasisFPatchCurvef*/

