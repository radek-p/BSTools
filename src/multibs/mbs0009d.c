
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

boolean mbs_multiBCHornerDerd ( int degree, int ncurves, int spdimen, int pitch,
                                const double *ctlpoints,
                                double t, double *p, double *d )
{
  void    *sp;
  double  *aux;
  boolean result;

  if ( degree == 0 )
    return _mbs_multiBCHornerDerd ( degree, ncurves, spdimen, pitch,
                                    ctlpoints, t, p, d, NULL );
  else {
    sp = pkv_GetScratchMemTop ();
    aux = pkv_GetScratchMemd ( 2*spdimen );
    if ( aux )
      result = _mbs_multiBCHornerDerd ( degree, ncurves, spdimen, pitch,
                                        ctlpoints, t, p, d, aux );
    else
      result = false;
    pkv_SetScratchMemTop ( sp );
    return result;
  }
} /*mbs_multiBCHornerDerd*/

