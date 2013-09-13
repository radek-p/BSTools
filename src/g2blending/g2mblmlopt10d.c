
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2011, 2013                            */
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
#include "msgpool.h"

/* ///////////////////////////////////////////////////////////////////////// */
boolean _g2mbl_MLSetupBlockHessiansd ( mesh_ml_optdata *d )
{
  int          i, nblocks;
  mlblock_desc *bd;

  nblocks = d->nblocks;
  for ( i = 0; 2*i+1 < nblocks; i++ ) {
    bd = &d->bd[i];
    if ( bd->nvcp <= MAX_NVCP ) {
      if ( !_g2mbl_MLSetupBlockCholHessiand ( d, i ) ) {
        PKV_SIGNALERROR ( LIB_G2BLENDING, ERRCODE_12, ERRMSG_12 );
        return false;
      }
    }
    else {
      if ( !_g2mbl_MLSetupBlockCGHessiand ( d, i ) ) {
        PKV_SIGNALERROR ( LIB_G2BLENDING, ERRCODE_13, ERRMSG_13 );
        return false;
      }
    }
  }
  for ( ; i < nblocks; i++ )
    if ( !_g2mbl_MLSetupBlockCholHessiand ( d, i ) ) {
      PKV_SIGNALERROR ( LIB_G2BLENDING, ERRCODE_12, ERRMSG_12 );
      return false;
    }
  return true;
} /*_g2mbl_MLSetupBlockHessiansd*/

