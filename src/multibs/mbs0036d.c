
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

boolean _mbs_FindBezPatchDiagFormd ( int degreeu, int degreev, int spdimen,
                                     CONST_ double *cpoints,
                                     int k, int l, double u, double v,
                                     double *dfcp, double *workspace )
{
  int i, pitch, ptch;

  if ( k > degreeu || l > degreev ) {
    PKV_SIGNALERROR ( LIB_MULTIBS, ERRCODE_5, ERRMSG_5 );
    return false;
  }
  pitch = (degreev+1)*spdimen;
  if ( k < degreeu ) {
    if ( !mbs_multiBCHornerd ( degreeu-k, k+1, pitch, pitch, cpoints, u,
                               workspace ) )
      return false;
  }
  else
    workspace = cpoints;
  if ( l < degreev ) {
    ptch = (l+1)*spdimen;
    for ( i = 0; i <= k; i++ )
      if ( !mbs_multiBCHornerd ( degreev-l, l+1, spdimen, spdimen,
                                 &workspace[i*pitch], v, &dfcp[i*ptch] ) )
        return false;
  }
  else
    memcpy ( dfcp, workspace, (k+1)*(l+1)*sizeof(double) );
  return true;
} /*_mbs_FindBezPatchDiagFormd*/

