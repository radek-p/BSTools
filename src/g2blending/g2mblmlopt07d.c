
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2011                                  */
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
boolean _g2mbl_MLFindBlockElementsd ( mesh_ml_optdata *d )
{
  void         *sp;
  int          nv, nvcp, ndomelems, ndomelems1, nblocks;
  int          i, j, j1, k, k0, k1, l;
  meshdom_elem *domelem;
  int          *domelind, *domelind1, *domelcpind, *vncpi;
  boolean      *mkel, *mkv;

  sp = pkv_GetScratchMemTop ();
  nv         = d->nv;
  ndomelems  = d->ndomelems;
  domelem    = d->domelem;
  domelcpind = d->domelcpind;
  nblocks    = d->nblocks;
  mkel = pkv_GetScratchMem ( ndomelems*sizeof(boolean) );
  mkv  = pkv_GetScratchMem ( nv*sizeof(boolean) );
  if ( !mkel || !mkv )
    goto failure;
        /* for the topmost block all elements are relevant */
  PKV_MALLOC ( d->bd[0].domelind, ndomelems*sizeof(int) );
  if ( !d->bd[0].domelind )
    goto failure;
  d->bd[0].ndomel = ndomelems;
  domelind = d->bd[0].domelind;
  for ( i = 0; i < ndomelems; i++ )
    domelind[i] = i;

  for ( j = 1; j < nblocks; j++ ) {
    j1 = (j-1) / 2;      /* number of the block above */
    ndomelems1 = d->bd[j1].ndomel;
    domelind1 = d->bd[j1].domelind;
    nvcp  = d->bd[j].nvcp;
    vncpi = d->bd[j].vncpi;
    ndomelems = 0;
    for ( i = 0; i < ndomelems1; i++ )
      mkel[domelind1[i]] = false;
    memset ( mkv, 0, nv*sizeof(boolean) );
    for ( i = 0; i < nvcp; i++ )
      mkv[vncpi[i]] = true;
    for ( i = l = 0;  i < ndomelems1;  i++ ) {
      k0 = domelem[domelind1[i]].firstcpi;
      k1 = k0 + domelem[domelind1[i]].ncp;
      for ( k = k0; k < k1; k++ )
        if ( mkv[domelcpind[k]] ) {
          mkel[domelind1[i]] = true;
          l++;
          break;
        }
    }
    d->bd[j].ndomel = l;
    PKV_MALLOC ( d->bd[j].domelind, l*sizeof(int) );
    if ( !d->bd[j].domelind )
      goto failure;
    domelind = d->bd[j].domelind;
    for ( i = l = 0; i < ndomelems1; i++ )
      if ( mkel[domelind1[i]] )
        domelind[l++] = domelind1[i];
  }

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*_g2mbl_MLFindBlockElementsd*/

