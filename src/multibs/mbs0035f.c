
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


boolean mbs_multiBCHornerDer3f ( int degree, int ncurves, int spdimen, int pitch,
                                 const float *ctlpoints, float t,
                                 float *p, float *d1, float *d2, float *d3 )
{
  void    *sp;
  float   *workspace;
  boolean result;

  sp = workspace = pkv_GetScratchMemf ( 4*spdimen );
  if ( workspace ) {
    result = _mbs_multiBCHornerDer3f ( degree, ncurves, spdimen, pitch,
                                       ctlpoints, t, p, d1, d2, d3, workspace );
    pkv_SetScratchMemTop ( sp );
    return result;
  }
  else
    return false;
} /*mbs_multiBCHornerDer3f*/

