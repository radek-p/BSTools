
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2005, 2013                            */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdio.h>

#include "pkvaria.h"
#include "pknum.h"
#include "pkgeom.h"
#include "multibs.h"

/* ////////////////////////////////////////// */
/* degree reduction of Bezier patches         */

boolean mbs_BCDegRedPf ( int spdimen,
                         int indegreeu, int indegreev, const float *inctlp,
                         int deltadegu, int deltadegv,
                         int *outdegreeu, int *outdegreev,
                         float *outctlp )
{
  void  *stp;
  float *acp;

  if ( deltadegu < 0 || deltadegu > indegreeu ||
       deltadegv < 0 || deltadegv > indegreev )
    return false;
  stp = pkv_GetScratchMemTop ();
  acp = pkv_GetScratchMemf ( (indegreeu-deltadegu+1)*(indegreev+1)*spdimen );
  if ( !acp ) {
    PKV_SIGNALERROR ( LIB_MULTIBS, ERRCODE_2, ERRMSG_2 );
    goto failure;
  }
  mbs_multiBCDegRedf ( 1, spdimen*(indegreev+1), 0, indegreeu,
                       (float*)inctlp, deltadegu, 0, outdegreeu, acp );
  mbs_multiBCDegRedf ( indegreeu-deltadegu+1, spdimen,
                       spdimen*(indegreev+1), indegreev, acp, deltadegv,
                       spdimen*(indegreev-deltadegv+1), outdegreev,
                       (float*)outctlp );
  pkv_SetScratchMemTop ( stp );
  return true;

failure:
  pkv_SetScratchMemTop ( stp );
  return false;
} /*mbs_BCDegRedPf*/

