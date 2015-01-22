
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

boolean mbs_BCHornerDer2Pf ( int degreeu, int degreev, int spdimen,   
                             const float *ctlpoints,   
                             float u, float v,   
                             float *p, float *du, float *dv,
                             float *duu, float *duv, float *dvv )
{
  void    *sp;
  float   *aux;
  int     size;
  boolean result;

  size = (36+3*(degreev+1))*spdimen;
  sp = aux = pkv_GetScratchMemf ( size );
  if ( aux ) {
    result = _mbs_BCHornerDer2Pf ( degreeu, degreev, spdimen, ctlpoints,
                                   u, v, p, du, dv, duu, duv, dvv, aux );
    pkv_SetScratchMemTop ( sp );
    return result;
  }
  else
    return false;
} /*mbs_BCHornerDer2Pf*/

