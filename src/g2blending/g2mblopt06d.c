
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2010, 2012                            */
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

#define _DEBUG

/* ///////////////////////////////////////////////////////////////////////// */
boolean _g2mbl_SetupHbl3x3d ( mesh_lmt_optdata *d )
{
  void          *sp;
  meshdom_elem  *domelem;
  int           nvcp, *nncpi, ndomelems, *domelcpind;
  int           Hblsize;
  nzHbl_rowdesc *iHbl;
  int           *cHbl;
  double        *Hbl;
  int           nzcdsize;
  byte          *nzcdistr;
  int           i, j, k, vi, vb, ei, eb;

  sp = pkv_GetScratchMemTop ();
  nvcp = d->nvcp;
  nncpi = d->nncpi;
  ndomelems = d->ndomelems;
  domelem = d->domelem;
  domelcpind = d->domelcpind;
  PKV_MALLOC ( d->iHbl, nvcp*sizeof(nzHbl_rowdesc) );
  if ( !d->iHbl )
    goto failure;
  iHbl = d->iHbl;
  nzcdsize = pkn_TMBSize ( nvcp );
  nzcdistr = pkv_GetScratchMem ( nzcdsize );
  if ( !nzcdistr )
    goto failure;

        /* count the preceding neighbours of each variable control point */
  memset ( iHbl, 0, nvcp*sizeof(nzHbl_rowdesc) );
  memset ( nzcdistr, 0, nzcdsize );
  for ( i = 0; i < ndomelems; i++ ) {
    vi = domelem[i].ncp;
    vb = domelem[i].firstcpi;
    for ( j = 0; j < vi; j++ ) {
      ei = nncpi[domelcpind[vb+j]];
      if ( ei >= 0 ) {
        for ( k = 0; k < vi; k++ ) {
          eb = nncpi[domelcpind[vb+k]];
          if ( eb >= ei ) {
            if ( !pkn_TMBTestAndSet ( nzcdistr, eb, ei ) )
              iHbl[eb].nhbl ++;
          }
        }
      }
    }
  }
  for ( i = 0; i < nvcp-1; i++ )
    iHbl[i+1].firsthbl = iHbl[i].firsthbl+iHbl[i].nhbl;
  Hblsize = d->Hblsize = iHbl[nvcp-1].firsthbl + iHbl[nvcp-1].nhbl;
#ifdef _DEBUG
printf ( "%d nonzero blocks out of %d\n", Hblsize, ((nvcp+1)*nvcp)/2 );
#endif
  PKV_MALLOC ( d->Hbl, Hblsize*(9*sizeof(double)+sizeof(int)) );
  if ( !d->Hbl )
    goto failure;
  Hbl = d->Hbl;
  cHbl = d->cHbl = (int*)&Hbl[Hblsize*9];
  memset ( cHbl, 0, Hblsize*sizeof(int) );
        /* setup the lists of the preceding neighbours */
  for ( i = 0; i < nvcp; i++ )
    iHbl[i].nhbl = 0;
  memset ( nzcdistr, 0, nzcdsize );
  for ( i = 0; i < ndomelems; i++ ) {
    vi = domelem[i].ncp;
    vb = domelem[i].firstcpi;
    for ( j = 0; j < vi; j++ ) {
      ei = nncpi[domelcpind[vb+j]];
      if ( ei >= 0 ) {
        for ( k = 0; k < vi; k++ ) {
          eb = nncpi[domelcpind[vb+k]];
          if ( eb >= ei ) {
            if ( !pkn_TMBTestAndSet ( nzcdistr, eb, ei ) ) {
              cHbl[iHbl[eb].firsthbl+iHbl[eb].nhbl] = ei;
              iHbl[eb].nhbl ++;
            }
          }
        }
      }
    }
  }
        /* sort the lists of the preceding neighbours */
  for ( i = 0; i < nvcp; i++ )
    if ( iHbl[i].nhbl > 1 ) {
      if ( pkv_SortFast ( sizeof(int), ID_UNSIGNED, sizeof(int),
                          0, iHbl[i].nhbl, &cHbl[iHbl[i].firsthbl] ) != SORT_OK )
        goto failure;
    }

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*_g2mbl_SetupHbl3x3d*/

boolean _g2mbl_SetupBlockHessianProfiled ( mesh_lmt_optdata *d, int bnum )
{
  void         *sp;
  meshdom_elem *domelem;
  block_desc   *bd;
  int          *nncpi, *vncpi, *domelcpind,
               *bvncpi, ndomelems, *domelind;
  int          *vpermut;
  vertex_desc  *nvcpi;
  int          nzcdsize;
  byte         *nzcdistr;
  int          i, j, k, fcp, ncp, vi, vb, ei, eb;

  sp = pkv_GetScratchMemTop ();
  nncpi = d->nncpi;
  vncpi = d->vncpi;
  domelem = d->domelem;
  domelcpind = d->domelcpind;
  bd = &d->block[bnum];
  fcp = bd->n0 / 3;
  ncp = bd->ncp;
  bvncpi = bd->vncpi;
  nvcpi = pkv_GetScratchMem ( ncp*sizeof(vertex_desc) );
  nzcdsize = pkn_TMBSize ( ncp );
  nzcdistr = pkv_GetScratchMem ( nzcdsize );
  if ( !nvcpi || !nzcdistr )
    goto failure;
  ndomelems = bd->ndomelems;
  domelind = bd->domelind;

        /* find the distribution of nonzero coefficients in the block */
        /* and count the in-block neighbours of each variable control point */
  for ( i = 0; i < ncp; i++ ) {
    nvcpi[i].vnum = vncpi[fcp+i];
    nvcpi[i].firstsel = ncp;
    nvcpi[i].nneigh = 0;
  }
  memset ( nzcdistr, 0, nzcdsize );
  for ( i = 0; i < ncp; i++ )
    pkn_TMBElemSet ( nzcdistr, i, i );
  for ( i = 0; i < ndomelems; i++ ) {
    vi = domelem[domelind[i]].ncp;
    vb = domelem[domelind[i]].firstcpi;
    for ( j = 0; j < vi; j++ ) {
      ei = nncpi[domelcpind[vb+j]] - fcp;
      if ( ei >= 0 && ei < ncp )
        for ( k = 0; k < vi; k++ ) {
          eb = nncpi[domelcpind[vb+k]] - fcp;
          if ( eb >= ei && eb < ncp ) {
            if ( !pkn_TMBTestAndSet ( nzcdistr, eb, ei ) ) {
              nvcpi[eb].nneigh ++;
              nvcpi[ei].nneigh ++;
            }
          }
        }
    }
  }
  if ( g2mbl_outputnzdistr )
    g2mbl_outputnzdistr ( d->nbl, bnum, false, d->nvcp, ncp, nzcdistr );

  if ( !_g2mbl_OrderCPoints ( d->nv, bd->ncp, fcp, nzcdsize, nzcdistr,
                              nncpi, bvncpi, nvcpi, ndomelems, domelem, domelcpind,
                              domelind, 3, &bd->hsize, bd->hprof, bd->vpermut ) )
    goto failure;
  vpermut = bd->vpermut;
  if ( g2mbl_outputnzdistr ) {
    memset ( nzcdistr, 0, nzcdsize );
    for ( i = 0; i < ncp; i++ )
        pkn_TMBElemSet ( nzcdistr, i, i );
    for ( i = 0; i < ndomelems; i++ ) {
      vi = domelem[domelind[i]].ncp;
      vb = domelem[domelind[i]].firstcpi;
      for ( j = 0; j < vi; j++ ) {
        ei = nncpi[domelcpind[vb+j]] - fcp;
        if ( ei >= 0 && ei < ncp )
          for ( k = 0; k < vi; k++ ) {
            eb = nncpi[domelcpind[vb+k]] - fcp;
            if ( eb >= ei && eb < ncp )
              pkn_TMBElemSet ( nzcdistr, vpermut[eb], vpermut[ei] );
          }
      }
    }
    g2mbl_outputnzdistr ( d->nbl, bnum, true, d->nvcp, ncp, nzcdistr );
  }

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*_g2mbl_SetupBlockHessianProfiled*/

boolean _g2mbl_SetupBlockRowsd ( mesh_lmt_optdata *d )
{
  void       *sp;
  int        nbl, thsize, i;
  double     *bh;
  block_desc *bd;
#ifdef DEBUG
FILE *f;
int  j;
#endif

  sp = pkv_GetScratchMemTop ();
  nbl = d->nbl;
  thsize = d->block[0].hsize;
  for ( i = 1; i < nbl; i++ )
    thsize += d->block[i].hsize;
  PKV_MALLOC ( d->BHessian, 2*thsize*sizeof(double) );
  if ( !d->BHessian )
    goto failure;
  bh = d->BHessian;
  for ( i = 0; i < nbl; i++ ) {
    bd = &d->block[i];
    pkn_NRBFindRowsd ( bd->nvars, bd->hprof, bh, bd->hrows );
    bh = &bh[bd->hsize];
    pkn_NRBFindRowsd ( bd->nvars, bd->hprof, bh, bd->lhrows );
    bh = &bh[bd->hsize];
  }
#ifdef DEBUG
f = fopen ( "tab.txt", "w+" );
fprintf ( f, "nncpi:\n" );
for ( i = 0; i < d->nv; i++ )
  fprintf ( f, "%4d,", d->nncpi[i] );
fprintf ( f, "\n" );
fprintf ( f, "vncpi:\n" );
for ( i = 0; i < d->nvcp; i++ )
  fprintf ( f, "%4d,", d->vncpi[i] );
fprintf ( f, "\n" );
for ( j = 0; j < d->nbl; j++ ) {
  fprintf ( f, "block %d:\n", j );
  fprintf ( f, "vpermut:\n" );
  for ( i = 0; i < d->block[j].ncp; i++ )
    fprintf ( f, "%4d,", d->block[j].vpermut[i] );
  fprintf ( f, "\n" );
  fprintf ( f, "vncpi:\n" );
  for ( i = 0; i < d->block[j].ncp; i++ )
    fprintf ( f, "%4d,", d->block[j].vncpi[i] );
  fprintf ( f, "\n" );
}
fclose ( f );
#endif

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*_g2mbl_SetupBlockRowsd*/

