
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2007, 2013                            */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

#include <math.h>
#include <string.h>
#include <stdlib.h>

#include "pkvaria.h"
#include "pknum.h"

/* //////////////////////////////////////////////// */
/* solving a linear least squares problem with      */
/* a band matrix and constraints                    */

boolean pkn_multiBandmSolveCRLSQd ( int nrows, int ncols,
                                    const bandm_profile *aprof, const double *a,
                                    int nconstr, int cpitch, const double *c,
                                    int nrsides, int spdimen,
                                    int bpitch, const double *b,
                                    int dpitch, const double *d,
                                    int xpitch, double *x )
{
  void          *sp;
  int           qsize, rsize, i;
  bandm_profile *qprof, *rprof;
  double        *qa, *ra, *ea, *eaa, *da, *y;

  sp = pkv_GetScratchMemTop ();
    /* working area allocation */
  pkn_BandmFindQRMSizes ( ncols, aprof, &qsize, &rsize );
  qprof = (bandm_profile*)pkv_GetScratchMem ( (ncols+1)*sizeof(bandm_profile) );
  rprof = (bandm_profile*)pkv_GetScratchMem ( (ncols+1)*sizeof(bandm_profile) );
  qa = pkv_GetScratchMemd ( qsize );
  ra = pkv_GetScratchMemd ( rsize );
  ea = pkv_GetScratchMemd ( nconstr*ncols );
  eaa = pkv_GetScratchMemd ( 2*nconstr );
  da = pkv_GetScratchMemd ( nconstr*spdimen );
  y = pkv_GetScratchMemd ( nrows*spdimen );
  if ( !qprof || !rprof || !qa || !ra || !ea || !eaa || !da || !y ) {
    PKV_SIGNALERROR ( LIB_PKNUM, 2, ERRMSG_2 );
    goto failure;
  }
    /* QR decomposition of the band matrix A */
  pkn_BandmQRDecomposeMatrixd ( nrows, ncols, aprof, a, qprof, qa, rprof, ra );
    /* solve the system R^TE = C^T */
  pkv_TransposeMatrixd ( nconstr, ncols, cpitch, c, nconstr, ea );
  pkn_multiBandmMultInvTrUTMVectord ( ncols, rprof, ra,
                                      nconstr, ea, ea );
    /* QR decomposition of the full matrix E */
  pkn_QRDecomposeMatrixd ( ncols, nconstr, ea, eaa );

  for ( i = 0; i < nrsides; i++ ) {
    /* solving a regular least squares problem (without constraints) */
    memcpy ( y, &b[i*bpitch], nrows*spdimen*sizeof(double) );
    pkn_multiBandmReflectVectord ( ncols, qprof, qa, spdimen, y );
    pkn_multiBandmMultInvUTMVectord ( ncols, rprof, ra, spdimen, y,
                                      &x[i*xpitch] );
    /* solving the system with the Schur matrix S c = d - C e, S = -E^TE */
    pkn_MultMatrixd ( nconstr, ncols, cpitch, c, spdimen, spdimen,
                      &x[i*xpitch], spdimen, da );
    if ( d )
      pkn_SubtractMatrixd ( nconstr, spdimen, spdimen, da, spdimen, &d[i*dpitch],
                            spdimen, da );
    pkn_multiMultInvTrUTVectord ( nconstr, ea, spdimen, spdimen, da,
                                  spdimen, da );
    pkn_multiMultInvUTVectord ( nconstr, ea, spdimen, spdimen, da,
                                spdimen, da );
    /* solving the system R^T h = C^T c and R k = h */
    pkn_MultTMatrixd ( nconstr, ncols, ncols, c, spdimen, spdimen, da,
                       spdimen, y );
    pkn_multiBandmMultInvTrUTMVectord ( ncols, rprof, ra, spdimen, y, y );
    pkn_multiBandmMultInvUTMVectord ( ncols, rprof, ra, spdimen, y, y );
    /* adding the correction */
    pkn_SubtractMatrixd ( ncols, spdimen, spdimen, &x[i*xpitch], spdimen, y,
                          spdimen, &x[i*xpitch] );
  }

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*pkn_multiBandmSolveCRLSQd*/

