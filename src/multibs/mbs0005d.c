
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

boolean mbs_BCHornerPd ( int degreeu, int degreev, int spdimen,
                         const double *ctlpoints,
                         double u, double v, double *ppoint )
{
  void    *sp;
  double  *workspace;
  boolean result;

  sp = pkv_GetScratchMemTop ();
  workspace = pkv_GetScratchMemd ( (degreev+1)*spdimen );
  if ( workspace )
    result = _mbs_BCHornerPd ( degreeu, degreev, spdimen, ctlpoints,
                               u, v, ppoint, workspace );
  else
    result = false;
  pkv_SetScratchMemTop ( sp );
  return result;
} /*mbs_BCHornerPd*/

