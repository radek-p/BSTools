
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2005, 2015                            */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include "pkvaria.h"
#include "pknum.h"
#include "pkgeom.h"
#include "multibs.h"

boolean mbs_BezC2CoonsFindCornersd ( int spdimen,
                                     int degc00, const double *c00,
                                     int degc01, const double *c01,
                                     int degc02, const double *c02,
                                     int degc10, const double *c10,
                                     int degc11, const double *c11,
                                     int degc12, const double *c12,
                                     double *pcorners )
{
  void    *sp;
  double  *workspace;
  boolean result;

  sp = workspace = pkv_GetScratchMemd ( 3*spdimen );
  if ( !workspace )
    return false;
  result = _mbs_BezC2CoonsFindCornersd ( spdimen,
                   degc00, c00, degc01, c01, degc02, c02,
                   degc10, c10, degc11, c11, degc12, c12,
                   pcorners, workspace );
  pkv_SetScratchMemTop ( sp );
  return result;
} /*mbs_BezC2CoonsFindCornersd*/

