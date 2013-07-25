
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2007                                  */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

#include <math.h>
#include <string.h>
#include <stdlib.h>

#include "pkvaria.h"
#include "pknum.h"

#include "msgpool.h"

/* //////////////////////////////////////////////// */
/* solving a linear least squares problem with      */
/* a band matrix and constraints                    */

void pkn_multiBandmSolveCRLSQf ( int nrows, int ncols,
                                 const bandm_profile *aprof, const float *a,
                                 int nconstr, int cpitch, const float *c,
                                 int nrsides, int spdimen,
                                 int bpitch, const float *b,
                                 int dpitch, const float *d,
                                 int xpitch, float *x )
{
  void          *sp;
  int           qsize, rsize, i;
  bandm_profile *qprof, *rprof;
  float         *qa, *ra, *ea, *eaa, *da, *y;

  sp = pkv_GetScratchMemTop ();
    /* working area allocation */
  pkn_BandmFindQRMSizes ( ncols, aprof, &qsize, &rsize );
  qprof = (bandm_profile*)pkv_GetScratchMem ( (ncols+1)*sizeof(bandm_profile) );
  rprof = (bandm_profile*)pkv_GetScratchMem ( (ncols+1)*sizeof(bandm_profile) );
  qa = pkv_GetScratchMemf ( qsize );
  ra = pkv_GetScratchMemf ( rsize );
  ea = pkv_GetScratchMemf ( nconstr*ncols );
  eaa = pkv_GetScratchMemf ( 2*nconstr );
  da = pkv_GetScratchMemf ( nconstr*spdimen );
  y = pkv_GetScratchMemf ( nrows*spdimen );
  if ( !qprof || !rprof || !qa || !ra || !ea || !eaa || !da || !y ) {
    pkv_SignalError ( LIB_PKNUM, 16, ERRMSG_0 );
    exit ( 1 );
  }
    /* QR decomposition of the band matrix A */
  pkn_BandmQRDecomposeMatrixf ( nrows, ncols, aprof, a, qprof, qa, rprof, ra );
    /* solve the system R^TE = C^T */
  pkv_TransposeMatrixf ( nconstr, ncols, cpitch, c, nconstr, ea );
  pkn_multiBandmMultInvTrUTMVectorf ( ncols, rprof, ra,
                                      nconstr, ea, ea );
    /* QR decomposition of the full matrix E */
  pkn_QRDecomposeMatrixf ( ncols, nconstr, ea, eaa );

  for ( i = 0; i < nrsides; i++ ) {
    /* solving a regular least squares problem (without constraints) */
    memcpy ( y, &b[i*bpitch], nrows*spdimen*sizeof(float) );
    pkn_multiBandmReflectVectorf ( ncols, qprof, qa, spdimen, y );
    pkn_multiBandmMultInvUTMVectorf ( ncols, rprof, ra, spdimen, y,
                                      &x[i*xpitch] );
    /* solving the system with the Schur matrix S c = d - C e, S = -E^TE */
    pkn_MultMatrixf ( nconstr, ncols, cpitch, c, spdimen, spdimen,
                      &x[i*xpitch], spdimen, da );
    if ( d )
      pkn_SubtractMatrixf ( nconstr, spdimen, spdimen, da, spdimen, &d[i*dpitch],
                            spdimen, da );
    pkn_multiMultInvTrUTVectorf ( nconstr, ea, spdimen, spdimen, da,
                                  spdimen, da );
    pkn_multiMultInvUTVectorf ( nconstr, ea, spdimen, spdimen, da,
                                spdimen, da );
    /* solving the system R^T h = C^T c and R k = h */
    pkn_MultTMatrixf ( nconstr, ncols, ncols, c, spdimen, spdimen, da,
                       spdimen, y );
    pkn_multiBandmMultInvTrUTMVectorf ( ncols, rprof, ra, spdimen, y, y );
    pkn_multiBandmMultInvUTMVectorf ( ncols, rprof, ra, spdimen, y, y );
    /* adding the correction */
    pkn_SubtractMatrixf ( ncols, spdimen, spdimen, &x[i*xpitch], spdimen, y,
                          spdimen, &x[i*xpitch] );
  }

  pkv_SetScratchMemTop ( sp );
} /*pkn_multiBandmSolveCRLSQf*/

