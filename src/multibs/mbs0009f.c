
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

boolean mbs_multiBCHornerDerf ( int degree, int ncurves, int spdimen, int pitch,
                                const float *ctlpoints,
                                float t, float *p, float *d )
{
  void    *sp;
  float   *aux;
  boolean result;

  if ( degree == 0 )
    return _mbs_multiBCHornerDerf ( degree, ncurves, spdimen, pitch,
                                    ctlpoints, t, p, d, NULL );
  else {
    sp = pkv_GetScratchMemTop ();
    aux = pkv_GetScratchMemf ( 2*spdimen );
    if ( aux )
      result = _mbs_multiBCHornerDerf ( degree, ncurves, spdimen, pitch,
                                        ctlpoints, t, p, d, aux );
    else
      result = false;
    pkv_SetScratchMemTop ( sp );
    return result;
  }
} /*mbs_multiBCHornerDerf*/

