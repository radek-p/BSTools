
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
#include "msgpool.h"

/* ///////////////////////////////////////////////////////////////////////// */
boolean _g2mbl_MLSetupBlockCholHessiand ( mesh_ml_optdata *d, int bn )
{
  void         *sp;
  meshdom_elem *domelem;
  mlblock_desc *bd;
  int          nv, nvcp, nvars, *vncpi, ndomelems,
               *domelcpind, ndomel, *domelind, *nncpi;
  int          hsize;
  int          nzcdsize;
  byte         *nzcdistr;
  int          i, j, k, vj, vk, fcp, ncp, el;
  vertex_desc  *vertd;

  sp = pkv_GetScratchMemTop ();
  ndomelems  = d->ndomelems;
  domelem    = d->domelem;
  domelcpind = d->domelcpind;
  nv         = d->nv;
  bd         = &d->bd[bn];
  nvcp       = bd->nvcp;
  nvars      = 3*nvcp;  /* == 3*bd->nvars */
  vncpi      = bd->vncpi;
  ndomel     = bd->ndomel;
  domelind   = bd->domelind;
  nzcdsize   = pkn_TMBSize ( nvcp );
  nzcdistr   = pkv_GetScratchMem ( nzcdsize );
  nncpi      = pkv_GetScratchMemi ( nv );
  vertd      = pkv_GetScratchMem ( nvcp*sizeof(vertex_desc) );
  if ( !nzcdistr || !nncpi || !vertd ) {
printf ( "%s\n", ERRMSG_0 );
    goto failure;
  }
        /* find the number of variable vertices in the block */
        /* for each vertex, -1 for the vertices not from the block */
  memset ( nncpi, 0xFF, nv*sizeof(int) );
  for ( i = 0; i < nvcp; i++ )
    nncpi[vncpi[i]] = i;
        /* find the distribution of the nonzero Hessian 3x3 blocks */
        /* and the neighbours for each variable vertex */
  for ( i = 0; i < nvcp; i++ ) {
    vertd[i].vnum = vncpi[i];
    vertd[i].firstsel = nvcp;
    vertd[i].nneigh = 0;
  }
  memset ( nzcdistr, 0, nzcdsize );
  for ( i = 0; i < nvcp; i++ )
    pkn_TMBElemSet ( nzcdistr, i, i );
  for ( i = 0; i < ndomel; i++ ) {
    el = domelind[i];
    fcp = domelem[el].firstcpi;
    ncp = domelem[el].ncp;
    for ( j = 0; j < ncp; j++ ) {
      vj = nncpi[domelcpind[fcp+j]];
      if ( vj >= 0 ) {
        for ( k = 0; k < ncp; k++ ) {
          vk = nncpi[domelcpind[fcp+k]];
          if ( vk > vj )
            if ( !pkn_TMBTestAndSet ( nzcdistr, vj, vk ) ) {
              vertd[vj].nneigh ++;
              vertd[vk].nneigh ++;
            }
        }
      }
    }
  }
  if ( g2mbl_outputnzdistr )
    g2mbl_outputnzdistr ( d->nblocks, bn, false, d->nvcp, nvcp, nzcdistr );

  PKV_MALLOC ( bd->hprof, 3*nvcp*sizeof(int) );
  if ( !bd->hprof ) {
printf ( "%s\n", ERRMSG_1 );
    goto failure;
  }
  if ( !_g2mbl_OrderCPoints ( nv, nvcp, 0, nzcdsize, nzcdistr, nncpi, vncpi,
                              vertd, ndomelems, domelem, domelcpind,
                              NULL, 3, &bd->hsize, bd->hprof, NULL ) ) {
printf ( "%s\n", ERRMSG_11 );
    goto failure;
  }
  hsize = bd->hsize;

  if ( g2mbl_outputnzdistr ) {
    memset ( nzcdistr, 0, nzcdsize );
    for ( i = 0; i < nvcp; i++ )
      pkn_TMBElemSet ( nzcdistr, i, i );
    for ( i = 0; i < ndomel; i++ ) {
      el = domelind[i];
      fcp = domelem[el].firstcpi;
      ncp = domelem[el].ncp;
      for ( j = 0; j < ncp; j++ ) {
        vj = nncpi[domelcpind[fcp+j]];
        if ( vj >= 0 )
          for ( k = 0; k < ncp; k++ ) {
            vk = nncpi[domelcpind[fcp+k]];
            if ( vk > vj )
              pkn_TMBElemSet ( nzcdistr, vk, vj );
          }
      }
    }
    g2mbl_outputnzdistr ( d->nblocks, bn, true, d->nvcp, nvcp, nzcdistr );
  }
  PKV_MALLOC ( bd->hrows, 2*(nvars*sizeof(double*)+hsize*sizeof(double)) );
  if ( !bd->hrows ) {
printf ( "%s\n", ERRMSG_1 );
    goto failure;
  }
  bd->lhrows = &bd->hrows[nvars];
  pkn_NRBFindRowsd ( nvars, bd->hprof, (double*)&bd->lhrows[nvars], bd->hrows );
  pkn_NRBFindRowsd ( nvars, bd->hprof, &bd->hrows[0][hsize], bd->lhrows );

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*_g2mbl_MLSetupBlockCholHessiand*/

