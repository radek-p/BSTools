
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

boolean _mbs_FindBezPatchDiagFormf ( int degreeu, int degreev, int spdimen,
                                     CONST_ float *cpoints,
                                     int k, int l, float u, float v,
                                     float *dfcp, float *workspace )
{
  int i, pitch, ptch;

  if ( k > degreeu || l > degreev ) {
    PKV_SIGNALERROR ( LIB_MULTIBS, ERRCODE_5, ERRMSG_5 );
    return false;
  }
  pitch = (degreev+1)*spdimen;
  if ( k < degreeu ) {
    if ( !mbs_multiBCHornerf ( degreeu-k, k+1, pitch, pitch, cpoints, u,
                               workspace ) )
      return false;
  }
  else
    workspace = cpoints;
  if ( l < degreev ) {
    ptch = (l+1)*spdimen;
    for ( i = 0; i <= k; i++ )
      if ( !mbs_multiBCHornerf ( degreev-l, l+1, spdimen, spdimen,
                                 &workspace[i*pitch], v, &dfcp[i*ptch] ) )
        return false;
  }
  else
    memcpy ( dfcp, workspace, (k+1)*(l+1)*sizeof(float) );
  return true;
} /*_mbs_FindBezPatchDiagFormf*/

