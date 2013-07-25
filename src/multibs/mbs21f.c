
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2005, 2009                            */
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

#include "msgpool.h"

/* ////////////////////////////////////////// */
/* degree elevation of Bezier patches         */

void mbs_BCDegElevPf ( int spdimen,
                       int indegreeu, int indegreev, const float *inctlp,
                       int deltadegu, int deltadegv,
                       int *outdegreeu, int *outdegreev,
                       float *outctlp )
{
  void  *stp;
  float *acp;

  stp = pkv_GetScratchMemTop ();
  acp = pkv_GetScratchMemf ( (indegreeu+deltadegu+1)*(indegreev+1)*spdimen );
  if ( !acp ) {
    pkv_SignalError ( LIB_MULTIBS, 21, ERRMSG_0 );
    exit ( 1 );
  }
  mbs_multiBCDegElevf ( 1, spdimen*(indegreev+1), 0, indegreeu, (float*)inctlp,
                        deltadegu, 0, outdegreeu, acp );
  mbs_multiBCDegElevf ( indegreeu+deltadegu+1, spdimen, spdimen*(indegreev+1),
                        indegreev, acp, deltadegv,
                        spdimen*(indegreev+deltadegv+1), outdegreev,
                        (float*)outctlp );
  pkv_SetScratchMemTop ( stp );
} /*mbs_BCDegElevPf*/

