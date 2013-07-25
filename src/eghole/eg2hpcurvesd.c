
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


boolean g2h_GetFinalPatchCurvesd ( GHoleDomaind *domain, int spdimen,
                                   CONST_ double *hole_cp, double *acoeff,
                                   void (*outcurve) ( int n, const double *cp ) )
{
  void          *sp;
  int           hole_k;
  GHolePrivateRecd  *privateG;
  G2HolePrivateRecd *privateG2;
  int           nfunc_a, nfunc_b;
  double        *bcp, *cp, *y;
  double        *bc00, *bc01, *bc02, *bc10, *bc11, *bc12,
                *bd00, *bd01, *bd02, *bd10, *bd11, *bd12;
  unsigned char *bfcpn;
  int           i, j, k, l;

  sp = pkv_GetScratchMemTop ();
  hole_k = domain->hole_k;
  privateG  = domain->privateG;
  privateG2 = domain->privateG2;
  nfunc_a = privateG2->nfunc_a;
  nfunc_b = 6*hole_k+1;
  bfcpn = privateG->bfcpn;
  bcp = pkv_GetScratchMemd ( hole_k*spdimen*(G2_CROSS00DEG+1) );
  y = pkv_GetScratchMemd ( spdimen*(G2_CROSS00DEG+1) );
  if ( !bcp || !y )
    goto failure;

  memset ( bcp, 0, hole_k*spdimen*(G2_CROSS00DEG+1)*sizeof(double) );
  for ( j = 0; j < nfunc_b; j++ ) {
    l = bfcpn[j];
    cp = &hole_cp[l*spdimen];
    for ( i = k = 0;  i < hole_k;  i++, k += spdimen ) {
      _g2h_GetBFBPatchCurvesd ( domain, j, i, &bc00, &bc01, &bc02, &bc10, &bc11, &bc12,
                                &bd00, &bd01, &bd02, &bd10, &bd11, &bd12 );
      pkn_MultMatrixd ( G2_CROSS00DEG+1, 1, 1, bc00, spdimen, 0, cp, spdimen, y );
      pkn_AddMatrixd ( 1, spdimen*(G2_CROSS00DEG+1), 0, &bcp[k*(G2_CROSS00DEG+1)],
                       0, y, 0, &bcp[k*(G2_CROSS00DEG+1)] );
    }
  }
  for ( j = 0; j < nfunc_a; j++ ) {
    cp = &acoeff[j*spdimen];
    for ( i = k = 0;  i < hole_k;  i++, k += spdimen ) {
      _g2h_GetBFAPatchCurvesd ( domain, j, i, &bc00, &bc01, &bc02, &bd00, &bd01, &bd02 );
      pkn_MultMatrixd ( G2_CROSS00DEG+1, 1, 1, bc00, spdimen, 0, cp, spdimen, y );
      pkn_SubtractMatrixd ( 1, spdimen*(G2_CROSS00DEG+1), 0, &bcp[k*(G2_CROSS00DEG+1)],
                            0, y, 0, &bcp[k*(G2_CROSS00DEG+1)] );
    }
  }

  for ( i = k = 0;  i < hole_k;  i++, k+= spdimen )
    outcurve ( G2_CROSS00DEG, &bcp[k*(G2_CROSS00DEG+1)] );

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*g2h_GetFinalPatchCurvesd*/

boolean g2h_GetExtFinalPatchCurvesd ( GHoleDomaind *domain, int spdimen,
                                      CONST_ double *hole_cp, double *acoeff,
                                      void (*outcurve) ( int n, const double *cp ) )
{
  int hole_k, nfunc_c;

  hole_k = domain->hole_k;
  nfunc_c = hole_k*G2_DBDIM;
  return g2h_GetFinalPatchCurvesd ( domain, spdimen, hole_cp,
                                    &acoeff[spdimen*nfunc_c], outcurve );
} /*g2h_GetExtFinalPatchCurvesd*/

boolean g2h_GetABasisFPatchCurved ( GHoleDomaind *domain, int fn, int i,
                                    double *bfpc )
{
  G2HolePrivateRecd *privateG2;
  int    hole_k, nfunc_a;
  double *bc00, *bc01, *bc02, *bd00, *bd01, *bd02;

  hole_k = domain->hole_k;
  if ( i < 0 || i >= hole_k )
    return false;
  privateG2 = domain->privateG2;
  nfunc_a = privateG2->nfunc_a;
  if ( fn < 0 || fn >= nfunc_a )
    return false;
  _g2h_GetBFAPatchCurvesd ( domain, fn, i, &bc00, &bc01, &bc02,
                            &bd00, &bd01, &bd02 );
  memcpy ( bfpc, bc00, (G2_CROSS00DEG+1)*sizeof(double) );
  return true;
} /*g2h_GetABasisFPatchCurved*/

boolean g2h_GetBBasisFPatchCurved ( GHoleDomaind *domain, int fn, int i,
                                    double *bfpc )
{
  int    hole_k, nfunc_b;
  double *bc00, *bc01, *bc02, *bc10, *bc11, *bc12,
         *bd00, *bd01, *bd02, *bd10, *bd11, *bd12;

  hole_k = domain->hole_k;
  if ( i < 0 || i >= hole_k )
    return false;
  nfunc_b = 6*hole_k+1;
  if ( fn < 0 || fn >= nfunc_b )
    return false;
  _g2h_GetBFBPatchCurvesd ( domain, fn, i, &bc00, &bc01, &bc02, &bc10, &bc11, &bc12,
                            &bd00, &bd01, &bd02, &bd10, &bd11, &bd12 );
  memcpy ( bfpc, bc00, (G2_CROSS00DEG+1)*sizeof(double) );
  return true;
} /*g2h_GetBBasisFPatchCurved*/

