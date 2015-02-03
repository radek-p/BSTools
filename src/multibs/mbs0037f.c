
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

boolean mbs_FindBezPatchDiagFormf ( int degreeu, int degreev, int spdimen,
                                    CONST_ float *cpoints,
                                    int k, int l, float u, float v,
                                    float *dfcp )
{
  void    *sp;
  float   *workspace;
  int     size;
  boolean result;

  sp = pkv_GetScratchMemTop ();
  if ( k < degreeu ) {
    size = (k+1)*(degreev+1)*spdimen;
    workspace = pkv_GetScratchMemf ( size );
    if ( !workspace ) {
      PKV_SIGNALERROR ( LIB_MULTIBS, ERRCODE_2, ERRMSG_2 );
      return false;
    }
  }
  else
    workspace = NULL;
  result = _mbs_FindBezPatchDiagFormf ( degreeu, degreev, spdimen, cpoints,
                                        k, l, u, v, dfcp, workspace );
  pkv_SetScratchMemTop ( sp );
  return result;
} /*mbs_FindBezPatchDiagFormf*/

