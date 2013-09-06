
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

#include <sys/times.h>
#include <unistd.h>

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
#ifdef _DEBUG
static clock_t tic, tt;
#endif

boolean g2mbl_MLSOptIterd ( void *data, boolean *finished )
{
  void            *sp;
  mesh_ml_optdata *d;
  int             bl;

#ifdef _DEBUG
pkv_Tic ( &tic );
#endif
  sp = pkv_GetScratchMemTop ();
  d = (mesh_ml_optdata*)data;
  bl = d->currentblock;
  if ( bl < 0 || bl >= d->nblocks )
    goto failure;
#ifdef _DEBUG
printf ( "%3d:", bl );
#endif
  d->bd[bl].fghflag &= ~FLAG_H;
  if ( d->bd[bl].iHbl ) { /* a "big" block */
    if ( !g2mbl_MLSOptBlockAd ( data, bl ) )
      goto failure;
  }
  else {                  /* a "small" block */
    if ( !g2mbl_MLSOptBlockBd ( data, bl ) )
      goto failure;
  }
  *finished = d->currentblock == -1;

  pkv_SetScratchMemTop ( sp );
#ifdef _DEBUG
tt = pkv_Toc ( &tic );
printf ( " time = %5.2f\n", (double)pkv_Seconds ( tt ) );
#endif
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*g2mbl_MLSOptIterd*/

