
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

boolean _mbs_multiBCHornerDerd ( int degree, int ncurves, int spdimen, int pitch,
                                 const double *ctlpoints,
                                 double t, double *p, double *d,
                                 double *workspace )
{
  int i;

  if ( degree == 0 ) {
    pkv_Selectd ( ncurves, spdimen, pitch, spdimen, ctlpoints, p );
    memset ( d, 0, ncurves*spdimen*sizeof(double) );
    return true;
  }
  else {
    for ( i = 0;
          i < ncurves;
          i++, ctlpoints += pitch, p += spdimen, d += spdimen ) {
      if ( !mbs_multiBCHornerd ( degree-1, 2, spdimen, spdimen, ctlpoints, t,
                                 workspace ) )
        return false;
      pkn_MatrixMDifferenced ( 1, spdimen, 0, &workspace[spdimen],
                               0, &workspace[0], (double)degree, 0, d );
      pkn_MatrixLinCombd ( 1, spdimen, 0, &workspace[0], 1.0-t, 0,
                           &workspace[spdimen], t, 0, p );
    }
  }
  return true;
} /*_mbs_multiBCHornerDerd*/

