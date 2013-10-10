
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2005, 2013                            */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

#include <math.h>
#include <string.h>
#include <stdio.h>  /* for testing */

#include "pkvaria.h"
#include "pkgeom.h"
#include "multibs.h"


/* /////////////////////////////////////////// */
/* Horner scheme for Bezier curves and patches */

boolean mbs_multiBCHornerDer3d ( int degree, int ncurves, int spdimen, int pitch,
                                 const double *ctlpoints, double t,
                                 double *p, double *d1, double *d2, double *d3 )
{
  void   *sp;
  int    i;
  double *a, s;

  sp = pkv_GetScratchMemTop ();
  if ( degree <= 2 ) {
    if ( !mbs_multiBCHornerDer2d ( degree, ncurves, spdimen, pitch, ctlpoints,
                                   t, p, d1, d2 ) )
      goto failure;
    memset ( d3, 0, ncurves*spdimen*sizeof(double) );
  }
  else if ( (a = pkv_GetScratchMemd ( 4*spdimen )) ) {
    s = 1.0-t;
    for ( i = 0;
          i < ncurves;
          i++, ctlpoints += pitch,
          p += spdimen, d1 += spdimen, d2 += spdimen, d3 += spdimen ) {
      if ( !mbs_multiBCHornerd ( degree-3, 4, spdimen, spdimen, ctlpoints, t, a ) )
        goto failure;
      pkn_SubtractMatrixd ( 1, spdimen, 0, &a[3*spdimen], 0, a, 0, d3 );
      pkn_SubtractMatrixd ( 1, spdimen, 0, &a[spdimen], 0, &a[2*spdimen], 0, d2 );
      pkn_AddMatrixMd ( 1, spdimen, 0, d3, 0, d2, 3.0, 0, d3 );
      pkn_MultMatrixNumd ( 1, spdimen, 0, d3,
                           (double)(degree*(degree-1)*(degree-2)), 0, d3 );
      pkn_MatrixLinCombd ( 1, spdimen, 0, a, s, 0, &a[spdimen], t, 0, a );
      pkn_MatrixLinCombd ( 1, spdimen, 0, &a[spdimen], s,
                           0, &a[2*spdimen], t, 0, &a[spdimen] );
      pkn_MatrixLinCombd ( 1, spdimen, 0, &a[2*spdimen], s,
                           0, &a[3*spdimen], t, 0, &a[2*spdimen] );

      pkn_AddMatrixd ( 1, spdimen, 0, a, 0, &a[2*spdimen], 0, d2 );
      pkn_AddMatrixMd ( 1, spdimen, 0, d2, 0, &a[spdimen], -2.0, 0, d2 );
      pkn_MultMatrixNumd ( 1, spdimen, 0, d2,
                          (double)(degree*(degree-1)), 0, d2 );
      pkn_MatrixLinCombd ( 1, spdimen, 0, a, s, 0, &a[spdimen], t, 0, a );
      pkn_MatrixLinCombd ( 1, spdimen, 0, &a[spdimen], s,
                           0, &a[2*spdimen], t, 0, &a[spdimen] );
      pkn_MatrixMDifferenced ( 1, spdimen, 0, &a[spdimen], 0, a,
                               (double)degree, 0, d1 );
      pkn_MatrixLinCombd ( 1, spdimen, 0, a, s, 0, &a[spdimen], t, 0, p );
    }
  }
  else
    goto failure;
  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*mbs_multiBCHornerDer3d*/

