
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
void g2mbl_MLOptDeallocated ( void ** data )
{
  mesh_ml_optdata *d;
  mlblock_desc    *bd;
  int             i, nbl;

  if ( *data ) {
    d = (mesh_ml_optdata*)*data;
    nbl = d->nblocks;
    if ( d->bd ) {
      bd = d->bd;
      for ( i = 0; i < nbl; i++ ) {
        if ( bd[i].vncpi )    PKV_FREE ( bd[i].vncpi );
        if ( bd[i].domelind ) PKV_FREE ( bd[i].domelind );
        if ( bd[i].iHbl )     PKV_FREE ( bd[i].iHbl );
        if ( bd[i].cHbl )     PKV_FREE ( bd[i].cHbl );
        if ( bd[i].hprof )    PKV_FREE ( bd[i].hprof );
        if ( bd[i].hrows )    PKV_FREE ( bd[i].hrows );
        if ( bd[i].wncpi )    PKV_FREE ( bd[i].wncpi );
        if ( bd[i].rmnzi )    PKV_FREE ( bd[i].rmnzi );
      }
      PKV_FREE ( d->bd );
    }
        /* this one is not NULL only for the shape only optimization algorithm */
    if ( d->mvcpn )   PKV_FREE ( d->mvcpn );

    if ( d->aqcoeff ) PKV_FREE ( d->aqcoeff );
    if ( d->domelem ) PKV_FREE ( d->domelem );
    if ( d->nncpi )   PKV_FREE ( d->nncpi );
    if ( d->Hbl )     PKV_FREE ( d->Hbl );
    if ( d->ftab1 )   PKV_FREE ( d->ftab1 );
    if ( d->mvtag )   PKV_FREE ( d->mvtag );

    PKV_FREE ( *data );
  }
} /*g2mbl_MLOptDeallocated*/

/* ///////////////////////////////////////////////////////////////////////// */
void g2mbl_MLSetNextBlock ( void *data, int nbl )
{
  mesh_ml_optdata *d;

  d = (mesh_ml_optdata*)data;
  nbl = max ( nbl, 0 );
  nbl = min ( nbl, d->nblocks-1 );
  d->currentblock = nbl;
  d->nextlevel = d->nlevels-1;
  d->dirtyblock = 0;
} /*g2mbl_MLSetNextBlock*/

void g2mbl_MLSetLogLeveld ( void *data, short level )
{
  ((mesh_ml_optdata*)data)->log_level = level;
} /*g2mbl_MLSetLogLeveld*/

short g2mbl_MLGetInfod ( void *data )
{
  return ((mesh_ml_optdata*)data)->error_code;
} /*g2mbl_MLGetInfod*/

/* ///////////////////////////////////////////////////////////////////////// */
void _g2mbl_MLInvalSmallBlocks ( mesh_ml_optdata *d, int bl )
{
  if ( bl < d->nblocks ) {
    d->bd[bl].fghflag = 0;
    _g2mbl_MLInvalSmallBlocks ( d, 2*bl+1 );
    _g2mbl_MLInvalSmallBlocks ( d, 2*bl+2 );
  }
} /*_g2mbl_MLInvalSmallBlocks*/

