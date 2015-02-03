
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

boolean mbs_BCHornerDer3Pd ( int degreeu, int degreev, int spdimen,
                             CONST_ double *ctlpoints,
                             double u, double v,
                             double *p, double *pu, double *pv,
                             double *puu, double *puv, double *pvv,
                             double *puuu, double *puuv, double *puvv, double *pvvv )
{
  void    *sp;
  double  *workspace;
  boolean result;

  sp = pkv_GetScratchMemTop ();
  workspace = pkv_GetScratchMemd ( (40+4*(degreev+1))*spdimen );
  if ( !workspace ) {
    PKV_SIGNALERROR ( LIB_MULTIBS, ERRCODE_2, ERRMSG_2 );
    return false;
  }
  result = _mbs_BCHornerDer3Pd ( degreeu, degreev, spdimen, ctlpoints,
                                 u, v, p, pu, pv, puu, puv, pvv,
                                 puuu, puuv, puvv, pvvv, workspace );
  pkv_SetScratchMemTop ( sp );
  return result;
} /*mbs_BCHornerDer3Pd*/

