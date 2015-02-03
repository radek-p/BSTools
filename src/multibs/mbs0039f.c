
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

#undef CONST_
#define CONST_

#include "multibs.h"

boolean mbs_BCHornerDer3Pf ( int degreeu, int degreev, int spdimen,
                             CONST_ float *ctlpoints,
                             float u, float v,
                             float *p, float *pu, float *pv,
                             float *puu, float *puv, float *pvv,
                             float *puuu, float *puuv, float *puvv, float *pvvv )
{
  void    *sp;
  float   *workspace;
  boolean result;

  sp = pkv_GetScratchMemTop ();
  workspace = pkv_GetScratchMemf ( (40+4*(degreev+1))*spdimen );
  if ( !workspace ) {
    PKV_SIGNALERROR ( LIB_MULTIBS, ERRCODE_2, ERRMSG_2 );
    return false;
  }
  result = _mbs_BCHornerDer3Pf ( degreeu, degreev, spdimen, ctlpoints,
                                 u, v, p, pu, pv, puu, puv, pvv,
                                 puuu, puuv, puvv, pvvv, workspace );
  pkv_SetScratchMemTop ( sp );
  return result;
} /*mbs_BCHornerDer3Pf*/

