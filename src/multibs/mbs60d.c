
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

boolean mbs_BCDegRedPd ( int spdimen,
                         int indegreeu, int indegreev, const double *inctlp,
                         int deltadegu, int deltadegv,
                         int *outdegreeu, int *outdegreev,
                         double *outctlp )
{
  void   *stp;
  double *acp;

  if ( deltadegu < 0 || deltadegu > indegreeu ||
       deltadegv < 0 || deltadegv > indegreev )
    return false;
  stp = pkv_GetScratchMemTop ();
  acp = pkv_GetScratchMemd ( (indegreeu-deltadegu+1)*(indegreev+1)*spdimen );
  if ( !acp ) {
    PKV_SIGNALERROR ( LIB_MULTIBS, ERRCODE_2, ERRMSG_2 );
    goto failure;
  }
  mbs_multiBCDegRedd ( 1, spdimen*(indegreev+1), 0, indegreeu,
                       (double*)inctlp, deltadegu, 0, outdegreeu, acp );
  mbs_multiBCDegRedd ( indegreeu-deltadegu+1, spdimen,
                       spdimen*(indegreev+1), indegreev, acp, deltadegv,
                       spdimen*(indegreev-deltadegv+1), outdegreev,
                       (double*)outctlp );
  pkv_SetScratchMemTop ( stp );
  return true;

failure:
  pkv_SetScratchMemTop ( stp );
  return false;
} /*mbs_BCDegRedPd*/

