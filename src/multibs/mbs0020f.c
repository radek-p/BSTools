
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2005, 2015                            */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

#include <math.h>
#include <string.h>
#include <stdio.h>  /* for testing */

#include "pkvaria.h"
#include "pkgeom.h"
#include "multibs.h"

/* workspace length must be at least 3*spdimen floats */
boolean _mbs_multiBCHornerDer2f ( int degree, int ncurves, int spdimen, int pitch,
                                  const float *ctlpoints,
                                  float t, float *p, float *d1, float *d2,
                                  float *workspace )
{
  int   i;
  float s;

  if ( degree <= 1 ) {
    if ( _mbs_multiBCHornerDerf ( degree, ncurves, spdimen, pitch, ctlpoints,
                                  t, p, d1, workspace ) ) {
      memset ( d2, 0, ncurves*spdimen*sizeof(float) );
      return true;
    }
    else
      return false;
  }
  else {
    s = (float)(1.0-t);
    for ( i = 0;
          i < ncurves;
          i++, ctlpoints += pitch,
          p += spdimen, d1 += spdimen, d2 += spdimen ) {
      if ( !mbs_multiBCHornerf ( degree-2, 3, spdimen, spdimen, ctlpoints, t,
                                 workspace ) )
        return false;
      pkn_AddMatrixf ( 1, spdimen, 0, workspace, 0, &workspace[2*spdimen],
                       0, d2 );
      pkn_AddMatrixMf ( 1, spdimen, 0, d2, 0, &workspace[spdimen], -2.0,
                        0, d2 );
      pkn_MultMatrixNumf ( 1, spdimen, 0, d2,
                          (float)(degree*(degree-1)), 0, d2 );
      pkn_MatrixLinCombf ( 1, spdimen, 0, workspace, s, 0, &workspace[spdimen],
                           t, 0, workspace );
      pkn_MatrixLinCombf ( 1, spdimen, 0, &workspace[spdimen], s,
                           0, &workspace[2*spdimen], t,
                           0, &workspace[spdimen] );
      pkn_MatrixMDifferencef ( 1, spdimen, 0, &workspace[spdimen], 0, workspace,
                               (float)degree, 0, d1 );
      pkn_MatrixLinCombf ( 1, spdimen, 0, workspace, s,
                           0, &workspace[spdimen], t, 0, p );
    }
    return true;
  }
} /*_mbs_multiBCHornerDer2f*/

