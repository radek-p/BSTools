
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2010                                  */
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

int g2mbl_GetBLMBlockNumd ( void *data, int *lastblock )
{
  mesh_lmt_optdata *d;

  d = data;
  if ( lastblock ) {
    if ( d->last_ntn )
      *lastblock = -1;
    else
      *lastblock = d->nextblock;
  }
  return d->nbl;
} /*g2mbl_GetBLMBlockNumd*/

void g2mbl_GetBLMTBlockInfod ( void *data,
                               int bln, int *nv, int *nvcp, int **nncpi,
                               int *c0, int *bnvcp,
                               int **vncpi, int **bvncpi, int **vpermut )
{
  mesh_lmt_optdata *d;
  block_desc       *bd;

  d = data;
  *nncpi = d->nncpi;
  *nv = d->nv;
  *nvcp = d->nvcp;
  *vncpi = d->vncpi;
  if ( bln >= 0 && bln < d->nbl ) {
    bd = &d->block[bln];
    *bnvcp = bd->ncp;
    *c0 = bd->n0 / 3;
    *bvncpi = bd->vncpi;
    *vpermut = bd->vpermut;
  }
  else {
    *bnvcp = d->nvcp;
    *c0 = 0;
    *bvncpi = *vpermut = NULL;
  }
} /*g2mbl_GetBLMTBlockInfod*/

