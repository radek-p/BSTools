
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2005, 2013                            */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdio.h>  /* for testing */

#include "pkvaria.h"
#include "pkgeom.h"
#include "multibs.h"

/* /////////////////////////////////////////// */
/* constructing spline curves of interpolation, solving the Hermite */
/* interpolation problem with knots of interpolation at the end points */
/* of the domain, with arbitrary multiplicities. */

void mbs_multiInterp2knHermiteBSd ( int ncurves, int spdimen, int degree,
                                    int lastknot, const double *knots,
                                    int nlbc, int lbcpitch, const double *lbc,
                                    int nrbc, int rbcpitch, const double *rbc,
                                    int pitch, double *ctlpoints )
{
  int i, j;

  if ( nlbc+nrbc != lastknot-degree ) {
    PKV_SIGNALERROR ( LIB_MULTIBS, ERRCODE_2, ERRMSG_2 );
    exit ( 1 );
  }

        /* get the interpolation conditions */
  pkv_Selectd ( ncurves, spdimen*nlbc, lbcpitch, pitch, lbc, ctlpoints );
          /* at the right end reverse */
  for ( i = 0; i < ncurves; i++ )
    pkv_Selectd ( nrbc, spdimen, spdimen, -spdimen, &rbc[i*rbcpitch],
                  &ctlpoints[i*pitch+(lastknot-degree-1)*spdimen] );

        /* now the proper computation */
          /* at the left end */
  for ( j = nlbc-2; j >= 0; j-- )
    for ( i = 0; i <= nlbc-2-j; i++ )
      pkn_AddMatrixMd ( ncurves, spdimen, pitch, &ctlpoints[(j+i)*spdimen],
                        pitch, &ctlpoints[(j+i+1)*spdimen],
                        (knots[i+degree+1]-knots[i+j+1])/(double)(degree-j),
                        pitch, &ctlpoints[(j+i+1)*spdimen] );
          /* at the right end */
  for ( j = nrbc-2; j >= 0; j-- )
    for ( i = lastknot-degree-j-1; i >= lastknot-degree-nrbc+1; i-- )
      pkn_AddMatrixMd ( ncurves, spdimen, pitch, &ctlpoints[i*spdimen],
                        pitch, &ctlpoints[(i-1)*spdimen],
                        -(knots[i+degree]-knots[i+j])/(double)(degree-j),
                        pitch, &ctlpoints[(i-1)*spdimen] );
} /*mbs_multiInterp2knHermiteBSd*/

