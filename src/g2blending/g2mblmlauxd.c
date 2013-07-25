
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2011                                  */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

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
int g2mbl_MLGetLastBlockd ( void *data )
{
  return ((mesh_ml_optdata*)data)->lastblock;
} /*g2mbl_MLGetLastBlockd*/

boolean g2mbl_MLGetBlockVCPNumbersd ( void *data, int bl,
                                      int *nvcp, int **vncpi, int *seed )
{
  mesh_ml_optdata *d;
  mlblock_desc    *bd;

  d = (mesh_ml_optdata*)data;
  if ( bl < 0 || bl >= d->nblocks )
    return false;
  else {
    bd = &d->bd[bl];
    *nvcp = bd->nvcp;
    *vncpi = bd->vncpi;
    if ( seed )
      *seed = bd->seed;
    return true;
  } 
} /*g2mbl_MLGetBlockVCPNumbersd*/

