
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2011, 2012                            */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

#include <limits.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

#include "pkvaria.h"
#include "pknum.h"
#include "pkgeom.h"
#include "multibs.h"
#include "bsmesh.h"
#include "egholed.h"
#include "g2blendingd.h"
#include "g2mblendingd.h"

#include "g2blprivated.h"
#include "g2mblprivated.h"
#include "g2mblmlprivated.h"

/* ///////////////////////////////////////////////////////////////////////// */
boolean g2mbl_MLSSuggestNLevels ( int nv, BSMvertex *mv, int *mvhei,
                                  int nhe, BSMhalfedge *mhe,
                                  int nfac, BSMfacet *mfac, int *mfhei,
                                  byte *mkcp,
                                  int *minlev, int *maxlev )
{
  int nvcp, mi, ma;

  nvcp = g2mbl_GetNvcp ( nv, mv, mvhei, nhe, mhe, nfac, mfac, mfhei, mkcp );
  if ( nvcp <= 1 )
    return false;
  mi = ma = 1;
  while ( nvcp > MAX_NVCPS ) {
    mi ++;
    ma ++;
    nvcp = (SUBBLOCK_STEP*nvcp)/(2*SUBBLOCK_STEP-1);
                  /* the same formula as in _g2mbl_MLDivideBlock */
  }
  while ( nvcp > MIN_NVCP ) {
    ma ++;
    nvcp = (SUBBLOCK_STEP*nvcp)/(2*SUBBLOCK_STEP-1);
                  /* the same formula as above */
  }
  *minlev = min ( mi, G2MBL_MAX_LEVELS );
  *maxlev = min ( ma, G2MBL_MAX_LEVELS );
  return true;
} /*g2mbl_MLSSuggestNLevels*/

