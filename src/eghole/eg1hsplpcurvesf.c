
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

#include "eg1holef.h"
#include "eg1hprivatef.h"
#include "eg1herror.h"


boolean g1h_GetSplFinalPatchCurvesf ( GHoleDomainf *domain, int spdimen,
                                      CONST_ float *hole_cp, float *acoeff,
                                      void (*outcurve) ( int n, int lkn,
                                          const float *kn, const float *cp ) )
{
#define NCTLP (lastomcknot-G1_CROSS00DEG)
  void          *sp;
  int           hole_k, lastomcknot;
  GHolePrivateRecf   *privateG;
  G1HolePrivateRecf  *privateG1;
  G1HoleSPrivateRecf *sprivate;
  int           nfunc_a, nfunc_b, nfunc_c, nfunc_d, nabc, nk, m1;
  float         *bcp, *cp, *y;
  float         *bc00, *bc01, *bc10, *bc11, *bd00, *bd01, *bd10, *bd11;
  float         *omcknots, *bezknots;
  unsigned char *bfcpn;
  int           i, j, k, l;

  sp = pkv_GetScratchMemTop ();
  hole_k = domain->hole_k;
  privateG  = domain->privateG;
  privateG1 = domain->privateG1;
  sprivate  = domain->SprivateG1;
  if ( !sprivate )
    goto failure;
  nfunc_a = privateG1->nfunc_a;
  nfunc_b = 6*hole_k+1;
  nfunc_c = sprivate->nsfunc_c;
  nfunc_d = sprivate->nsfunc_d;
  nabc = nfunc_a+nfunc_b+nfunc_c;
  nk = sprivate->nk;
  m1 = sprivate->m1;
  lastomcknot = sprivate->lastomcknot;
  omcknots = sprivate->omcknots;
  bezknots = &sprivate->bezknots[G1H_FINALDEG-G1_CROSS00DEG];
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
      _g1h_GetBFBPatchCurvesf ( domain, j, i, &bc00, &bc01, &bc10, &bc11,
                                &bd00, &bd01, &bd10, &bd11 );
      pkn_MultMatrixf ( G1_CROSS00DEG+1, 1, 1, bc00, spdimen, 0, cp, spdimen, y );
      pkn_AddMatrixf ( 1, spdimen*(G1_CROSS00DEG+1), 0, &bcp[k*NCTLP],
                       0, y, 0, &bcp[k*NCTLP] );
    }
  }
  for ( j = 0; j < nfunc_a; j++ ) {
    cp = &acoeff[(nfunc_c+nfunc_d+j)*spdimen];
    for ( i = k = 0;  i < hole_k;  i++, k += spdimen ) {
      _g1h_GetBFAPatchCurvesf ( domain, j, i, &bc00, &bc01, &bd00, &bd01 );
      pkn_MultMatrixf ( G1_CROSS00DEG+1, 1, 1, bc00, spdimen, 0, cp, spdimen, y );
      pkn_SubtractMatrixf ( 1, spdimen*(G1_CROSS00DEG+1), 0, &bcp[k*NCTLP],
                            0, y, 0, &bcp[k*NCTLP] );
    }
  }
        /* convert the Bezier representation to the final B-spline */
  mbs_multiOsloInsertKnotsf ( hole_k, spdimen, G1_CROSS00DEG,
                              2*G1_CROSS00DEG+1, bezknots, spdimen*NCTLP, bcp,
                              lastomcknot, omcknots, spdimen*NCTLP, y );
  memcpy ( bcp, y, hole_k*spdimen*NCTLP*sizeof(float) );
        /* add the linear combination of the block D functions */
  for ( i = 0; i < hole_k; i++ ) {
    for ( j = 0;  j < nk*m1;  j++ ) {
      cp = &acoeff[(nfunc_c+(2*i*nk*m1+j))*spdimen];
      _g1h_GetSplDBasisCrossDerf ( domain, nabc+2*(i*nk*m1+j), i,
              y, &y[NCTLP], &y[NCTLP] );
      pkn_MultMatrixf ( NCTLP, 1, 1, y, spdimen, 0, cp, spdimen, &y[NCTLP] );
      pkn_SubtractMatrixf ( 1, spdimen*NCTLP, 0, &bcp[k*NCTLP],
                            0, y, 0, &bcp[k*NCTLP] );
    }
  }

  for ( i = k = 0;  i < hole_k;  i++, k+= spdimen )
    outcurve ( G1_CROSS00DEG, lastomcknot, omcknots, &bcp[k*NCTLP] );

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
#undef NCTLP
} /*g1h_GetSplFinalPatchCurvesf*/

boolean g1h_GetSplABasisFPatchCurvef ( GHoleDomainf *domain, int fn, int i,
                                       int *lkn, float *kn, float *bfpc )
{
  void *sp;
  G1HolePrivateRecf  *privateG1;
  G1HoleSPrivateRecf *sprivate;
  int   hole_k, nfunc_a, lastomcknot;
  float *bc00, *bc01, *bd00, *bd01, *bezknots;

  sp = pkv_GetScratchMemTop ();
  hole_k = domain->hole_k;
  if ( i < 0 || i >= hole_k )
    goto failure;
  privateG1 = domain->privateG1;
  sprivate = domain->SprivateG1;
  if ( !sprivate )
    goto failure;
  nfunc_a = privateG1->nfunc_a;
  if ( fn < 0 || fn >= nfunc_a )
    goto failure;
  _g1h_GetBFAPatchCurvesf ( domain, fn, i, &bc00, &bc01, &bd00, &bd01 );
  *lkn = lastomcknot = sprivate->lastomcknot;
  memcpy ( kn, sprivate->omcknots, (lastomcknot+1)*sizeof(float) );
  bezknots = &sprivate->bezknots[G1H_FINALDEG-G1_CROSS00DEG];
  mbs_multiOsloInsertKnotsf ( 1, 1, G1_CROSS00DEG,
                              2*G1_CROSS00DEG+1, bezknots, 0, bc00,
                              lastomcknot, kn, 0, bfpc );
  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*g1h_GetSplABasisFPatchCurvef*/

boolean g1h_GetSplBBasisFPatchCurvef ( GHoleDomainf *domain, int fn, int i,
                                       int *lkn, float *kn, float *bfpc )
{
  void  *sp;
  G1HoleSPrivateRecf *sprivate;
  int   hole_k, nfunc_b, lastomcknot;
  float *bc00, *bc01, *bc10, *bc11, *bd00, *bd01, *bd10, *bd11, *bezknots;

  sp = pkv_GetScratchMemTop ();
  hole_k = domain->hole_k;
  if ( i < 0 || i >= hole_k )
    goto failure;
  nfunc_b = 6*hole_k+1;
  if ( fn < 0 || fn >= nfunc_b )
    goto failure;
  sprivate = domain->SprivateG1;
  if ( !sprivate )
    goto failure;
  _g1h_GetBFBPatchCurvesf ( domain, fn, i, &bc00, &bc01, &bc10, &bc11,
                            &bd00, &bd01, &bd10, &bd11 );
  *lkn = lastomcknot = sprivate->lastomcknot;
  memcpy ( kn, sprivate->omcknots, (lastomcknot+1)*sizeof(float) );
  bezknots = &sprivate->bezknots[G1H_FINALDEG-G1_CROSS00DEG];
  mbs_multiOsloInsertKnotsf ( 1, 1, G1_CROSS00DEG,
                              2*G1_CROSS00DEG+1, bezknots, 0, bc00,
                              lastomcknot, kn, 0, bfpc );
  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*g1h_GetSplBBasisFPatchCurvef*/

boolean g1h_GetSplDBasisFPatchCurvef ( GHoleDomainf *domain, int fn, int i,
                                       int *lkn, float *kn, float *bfpc )
{
  void  *sp;
  G1HolePrivateRecf  *privateG1;
  G1HoleSPrivateRecf *sprivate;
  int   hole_k, nfunc_a, nfunc_b, nfunc_c, nfunc_d, nabc;
  int   lastomcknot, lastpvknot;
  float *pu, *pv;

  sp = pkv_GetScratchMemTop ();
  hole_k = domain->hole_k;
  if ( i < 0 || i >= hole_k )
    goto failure;
  privateG1 = domain->privateG1;
  sprivate = domain->SprivateG1;
  if ( !sprivate )
    goto failure;
  nfunc_a = privateG1->nfunc_a;
  nfunc_b = 6*hole_k+1;
  nfunc_c = sprivate->nsfunc_c;
  nfunc_d = sprivate->nsfunc_d;
  if ( fn < 0 || fn >= nfunc_d )
    goto failure;
  lastpvknot = sprivate->lastpvknot;
  pu = pkv_GetScratchMemf ( 2*(lastpvknot-G1_CROSS01DEG) );
  if ( !pu )
    goto failure;
  pv = &pu[lastpvknot-G1_CROSS01DEG];
  nabc = nfunc_a+nfunc_b+nfunc_c;
  *lkn = lastomcknot = sprivate->lastomcknot;
  memcpy ( kn, sprivate->omcknots, (lastomcknot+1)*sizeof(float) );
  _g1h_GetSplDBasisCrossDerf ( domain, nabc+fn, i, bfpc, pv, pu );

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*g1h_GetSplDBasisFPatchCurvef*/

