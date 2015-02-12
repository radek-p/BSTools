
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

boolean mbs_BezC1CoonsFindCornersf ( int spdimen,
                                     int degc00, const float *c00,
                                     int degc01, const float *c01,
                                     int degc10, const float *c10,
                                     int degc11, const float *c11,
                                     float *pcorners )
{
  void    *sp;
  float   *workspace;
  boolean result;
  
  sp = workspace = pkv_GetScratchMemf ( 2*spdimen );
  if ( !workspace )
    return false;
  result = _mbs_BezC1CoonsFindCornersf ( spdimen,
                   degc00, c00, degc01, c01, degc10, c10, degc11, c11,
                   pcorners, workspace );
  pkv_SetScratchMemTop ( sp );
  return result;
} /*mbs_BezC1CoonsFindCornersf*/

