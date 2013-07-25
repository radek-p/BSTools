
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2005, 2009                            */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdio.h>

#include "pkvaria.h"
#include "pknum.h"
#include "pkgeom.h"
#include "multibs.h"

#include "msgpool.h"

/* ////////////////////////////////////////// */
/* multiple knots insertion/removal with use  */
/* of the Oslo algorithm                      */

void mbs_multiOsloInsertKnotsf ( int ncurves, int spdimen, int degree,
                                 int inlastknot, const float *inknots,
                                 int inpitch, float *inctlpoints,
                                 int outlastknot, const float *outknots,
                                 int outpitch, float *outctlpoints )
{
  void          *stp;
  int           k, nc, nr, as;
  bandm_profile *aprof;
  float         *aa;

  if ( inlastknot == outlastknot ) { /* nothing really to do, just copy data */
     pkv_Selectf ( ncurves, spdimen*(inlastknot-degree), inpitch, outpitch,
                   inctlpoints, outctlpoints );
    return;
  }

  stp = pkv_GetScratchMemTop ();

/* Construct the Oslo matrix */
  nr = outlastknot-degree;
  nc = inlastknot-degree;
  aprof = (bandm_profile*)pkv_GetScratchMem ( (nc+1)*sizeof(bandm_profile) );
  if ( !aprof ) {
    pkv_SignalError ( LIB_MULTIBS, 16, ERRMSG_0 );
    exit ( 1 );
  }
  as = mbs_BuildOsloMatrixProfilef ( degree, inlastknot, inknots,
                                     outlastknot, outknots, aprof );
  aa = pkv_GetScratchMemf ( as );
  if ( !aa ) {
    pkv_SignalError ( LIB_MULTIBS, 17, ERRMSG_0 );
    exit ( 1 );
  }
  mbs_BuildOsloMatrixf ( degree, inlastknot, inknots, outknots, aprof, aa );

/* Multiply the vectors of control points by the Oslo matrix. */
  for ( k = 0; k < ncurves; k++ )
    pkn_multiBandmMultVectorf ( nr, nc, aprof, aa, spdimen,
                                &inctlpoints[k*inpitch],
                                &outctlpoints[k*outpitch] );

  pkv_SetScratchMemTop ( stp );
} /*mbs_multiOsloInsertKnotsf*/

void mbs_multiOsloRemoveKnotsLSQf ( int ncurves, int spdimen, int degree,
                                    int inlastknot, const float *inknots,
                                    int inpitch, float *inctlpoints,
                                    int outlastknot, const float *outknots,
                                    int outpitch, float *outctlpoints )
{
  void          *stp;
  int           k, nc, nr, as, qs, rs;
  bandm_profile *aprof, *qprof, *rprof;
  float         *aa, *qa, *ra;

  if ( inlastknot == outlastknot ) { /* nothing really to do, just copy data */
     pkv_Selectf ( ncurves, spdimen*(inlastknot-degree), inpitch, outpitch,
                   inctlpoints, outctlpoints );
    return;
  }

   stp = pkv_GetScratchMemTop ();

/* Construct the Oslo matrix */
  nr = inlastknot-degree;
  nc = outlastknot-degree;
  aprof = (bandm_profile*)pkv_GetScratchMem ( (nc+1)*sizeof(bandm_profile) );
  if ( !aprof ) {
    pkv_SignalError ( LIB_MULTIBS, 18, ERRMSG_0 );
    exit ( 1 );
  }
  as = mbs_BuildOsloMatrixProfilef ( degree, outlastknot, outknots,
                                     inlastknot, inknots, aprof );
  as = max ( as, spdimen*(inlastknot-degree) );
                     /* greater size, as later it is used for something else */
  aa = pkv_GetScratchMemf ( as );
  if ( !aa ) {
    pkv_SignalError ( LIB_MULTIBS, 19, ERRMSG_0 );
    exit ( 1 );
  }
  mbs_BuildOsloMatrixf ( degree, outlastknot, outknots, inknots, aprof, aa );

/* Decompose the Oslo matrix into the orthogonal and upper triangular     */
/* factors. The orthogonal matrix is represented as a composition of      */
/* Householder reflections. All the matrices are stored as band matrices. */
  pkn_BandmFindQRMSizes ( nc, aprof, &qs, &rs );
  qprof = (bandm_profile*)pkv_GetScratchMem ( (nc+1)*sizeof(bandm_profile) );
  rprof = (bandm_profile*)pkv_GetScratchMem ( (nc+1)*sizeof(bandm_profile) );
  qa = pkv_GetScratchMemf ( qs );
  ra = pkv_GetScratchMemf ( rs );
  if ( !qprof || !rprof || !qa || !ra ) {
    pkv_SignalError ( LIB_MULTIBS, 20, ERRMSG_0 );
    exit ( 1 );
  }
  pkn_BandmQRDecomposeMatrixf ( nr, nc, aprof, aa, qprof, qa, rprof, ra );

/* Use these factors to solve the least-squares problem. */
  for ( k = 0; k < ncurves; k++ ) {
    memcpy ( aa, &inctlpoints[k*inpitch],
             spdimen*(inlastknot-degree)*sizeof(float)  );
    pkn_multiBandmReflectVectorf ( nc, qprof, qa, spdimen, aa );
    pkn_multiBandmMultInvUTMVectorf ( nc, rprof, ra, spdimen, aa,
                                      &outctlpoints[k*outpitch] );
  }

  pkv_SetScratchMemTop ( stp );
} /*mbs_multiOsloRemoveKnotsLSQf*/

