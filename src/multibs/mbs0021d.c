
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

boolean mbs_multiBCHornerDer2d ( int degree, int ncurves, int spdimen, int pitch,
                                 const double *ctlpoints,
                                 double t, double *p, double *d1, double *d2 )
{
  void    *sp;
  double  *aux;
  boolean result;

  sp = aux = pkv_GetScratchMemd ( 3*spdimen );
  if ( aux ) {
    result = _mbs_multiBCHornerDer2d ( degree, ncurves, spdimen, pitch,
                                       ctlpoints, t, p, d1, d2, aux );
    pkv_SetScratchMemTop ( sp );
    return result;
  }
  else
    return false;
} /*mbs_multiBCHornerDer2d*/

