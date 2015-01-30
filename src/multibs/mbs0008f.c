
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2005, 2015                            */
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

boolean _mbs_multiBCHornerDerf ( int degree, int ncurves, int spdimen, int pitch,
                                 const float *ctlpoints,
                                 float t, float *p, float *d,
                                 float *workspace )
{
  int i;

  if ( degree == 0 ) {
    pkv_Selectf ( ncurves, spdimen, pitch, spdimen, ctlpoints, p );
    memset ( d, 0, ncurves*spdimen*sizeof(float) );
    return true;
  }
  else {
    for ( i = 0;
          i < ncurves;
          i++, ctlpoints += pitch, p += spdimen, d += spdimen ) {
      if ( !mbs_multiBCHornerf ( degree-1, 2, spdimen, spdimen, ctlpoints, t,
                                 workspace ) )
        return false;
      pkn_MatrixMDifferencef ( 1, spdimen, 0, &workspace[spdimen],
                               0, &workspace[0], (float)degree, 0, d );
      pkn_MatrixLinCombf ( 1, spdimen, 0, &workspace[0], 1.0-t, 0,
                           &workspace[spdimen], t, 0, p );
    }
  }
  return true;
} /*_mbs_multiBCHornerDerf*/
