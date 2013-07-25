
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2010, 2011                            */
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
#define _g2mbl_FindDistances bsm_FindVertexDistances2

boolean _g2mbl_SetupAltBlocks ( int nv, BSMvertex *mv, int *mvhei,
                                int nhe, BSMhalfedge *mhe,
                                int nfac, BSMfacet *mfac, int *mfhei,
                                int *nncpi, char nbl, int *bltag, int *blseed )
{
  void *sp;
  int       *sdist[G2MBL_MAX_BLOCKS];
  int       bn[G2MBL_MAX_BLOCKS];
  int       i, j, k, s, extv, nvcp;
  int       d, dmax;
  int       deg, fhe, fdeg, ffhe, f, t, v;
  int       m;
  pkv_queue *q;

  sp = pkv_GetScratchMemTop ();
  q = NULL;
  if ( nbl < 2 )
    goto failure;
        /* count the variable control points */
  for ( i = nvcp = 0;  i < nv;  i++ )
    if ( nncpi[i] >= 0 )
      nvcp ++;

  sdist[0] = pkv_GetScratchMemi ( (nbl+1)*nv );
  if ( !sdist[0] )
    goto failure;
  for ( i = 1; i < nbl; i++ )
    sdist[i] = &sdist[i-1][nv];
  blseed[0] = 0;
        /* choose seeds and find the distances from these seeds in the mesh */
  for ( i = 0; i < nbl; i++ ) {
    if ( !_g2mbl_FindDistances ( nv, mv, mvhei, nhe, mhe, nfac, mfac, mfhei,
                                 blseed[i], sdist[i] ) )
      goto failure;
    if ( i < nbl-1 ) {
      s = 0;
      dmax = 0;
      for ( j = 0; j < nv; j++ )
        if ( nncpi[j] >= 0 ) {
          d = INT_MAX;
          for ( k = 0; k < i; k++ )
            d = min ( d, sdist[k][j] );
          if ( d > dmax ) {
            dmax = d;
            s = j;
          }
        }
      blseed[i+1] = s;
    }
  }
  for ( i = 0; i < nbl; i++ ) {
    s = 0;
    dmax = 0;
    for ( j = 0; j < nv; j++ )
      if ( nncpi[j] >= 0 ) {
        d = INT_MAX;
        for ( k = 0; k < nbl; k++ )
          if ( k != i )
            d = min ( d, sdist[k][j] );
        if ( d > dmax ) {
          dmax = d;
          s = j;
        }
      }
    blseed[i] = s;
    if ( !_g2mbl_FindDistances ( nv, mv, mvhei, nhe, mhe, nfac, mfac, mfhei,
                                 blseed[i], sdist[i] ) )
      goto failure;
  }
        /* assign the vertices to the blocks - Voronoi diagram */
  memset ( bltag, 0, nv*sizeof(int) );
  memset ( bn, 0, nbl*sizeof(int) );
  for ( i = 0; i < nv; i++ )
    if ( nncpi[i] >= 0 ) {
      for ( k = 0, d = sdist[0][i], j = 1;  j < nbl; j++ )
        if ( sdist[j][i] < d ) {
          d = sdist[j][i];
          k = j;
        }
      bltag[i] |= 0x0001 << k;
      bn[k] ++;
    }
        /* extend the blocks - to be fine tuned */
  if ( nbl == 2 )
    extv = nvcp / 4;
  else
    extv = (3*nvcp)/(4*nbl); /* nvcp / nbl ??? */
  q = pkv_InitQueue ( nv, sizeof(int) );
  if ( !q )
    goto failure;
  for ( k = 0; k < nbl; k++ ) {
    m = 0x0001 << k;
    for ( i = 0; i < nv; i++ )
      if ( bltag[i] & m )
        pkv_QueueInsert ( q, &i );
    for ( j = 0;  j < extv && !pkv_QueueEmpty ( q ); ) {
      pkv_QueueRemoveFirst ( q, &i );
      deg = mv[i].degree;
      fhe = mv[i].firsthalfedge;
      for ( s = 0; s < deg; s++ ) {
        f = mhe[mvhei[fhe+s]].facetnum;
        fdeg = mfac[f].degree;
        ffhe = mfac[f].firsthalfedge;
        for ( t = 0; t < fdeg; t++ ) {
          v = mhe[mfhei[ffhe+t]].v0;
          if ( nncpi[v] >= 0 && !(bltag[v] & m) ) {
            bltag[v] |= m;
            pkv_QueueInsert ( q, &v );
            j ++;
          }
        }
      }
    }
    pkv_ResetQueue ( q );
  }

  PKV_FREE ( q );
  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  if ( q ) PKV_FREE ( q );
  pkv_SetScratchMemTop ( sp );
  return false;
} /*_g2mbl_SetupAltBlocks*/

boolean _g2mbl_SetupAltBlockHessianProfiled ( mesh_lmt_optdata *d, int bnum )
{
  void         *sp;
  meshdom_elem *domelem;
  block_desc   *bd;
  int          *nncpi, *vncpi, *domelcpind, ndomelems;
  vertex_desc  *nvcpi;
  int          nzcdsize;
  byte         *nzcdistr;
  int          i, j, k, nv, ncp, vi, vb, ei, eb;
  int          *bltag, mask;

  sp = pkv_GetScratchMemTop ();
  nv = d->nv;
  domelem = d->domelem;
  domelcpind = d->domelcpind;
  bltag = d->bltag;
  bd = &d->block[bnum];
  ncp = bd->ncp;
  nncpi = bd->nncpi;
  vncpi = bd->vncpi;
  ndomelems = d->ndomelems;
  nvcpi = pkv_GetScratchMem ( ncp*sizeof(vertex_desc) );
  nzcdsize = pkn_TMBSize ( ncp );
  nzcdistr = pkv_GetScratchMem ( nzcdsize );
  if ( !nvcpi || !nzcdistr )
    goto failure;

        /* find the distribution of nonzero coefficients in the block */
        /* and count the in-block neighbours of each variable control point */
  mask = 0x0001 << bnum;
  for ( i = j = 0;  i < nv;  i++ )
    if ( bltag[i] & mask ) {
      nvcpi[j].vnum = i;
      nvcpi[j].firstsel = ncp;
      nvcpi[j].nneigh = 0;
      nncpi[i] = j ++;
    }
    else
      nncpi[i] = -1;
  memset ( nzcdistr, 0, nzcdsize );
  for ( i = 0; i < ncp; i++ )
    pkn_TMBElemSet ( nzcdistr, i, i );
  for ( i = 0; i < ndomelems; i++ ) {
    vi = domelem[i].ncp;
    vb = domelem[i].firstcpi;
    for ( j = 0; j < vi; j++ ) {
      ei = nncpi[domelcpind[vb+j]];
      if ( ei >= 0 )
        for ( k = 0; k < vi; k++ ) {
          eb = nncpi[domelcpind[vb+k]];
          if ( eb > ei ) {
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

  if ( !_g2mbl_OrderCPoints ( d->nv, bd->ncp, 0, nzcdsize, nzcdistr,
                              nncpi, vncpi, nvcpi, ndomelems, domelem, domelcpind,
                              NULL, 3, &bd->hsize, bd->hprof, NULL ) )
    goto failure;
  if ( g2mbl_outputnzdistr ) {
    memset ( nzcdistr, 0, nzcdsize );
    for ( i = 0; i < ncp; i++ )
      pkn_TMBElemSet ( nzcdistr, i, i );
    for ( i = 0; i < ndomelems; i++ ) {
      vi = domelem[i].ncp;
      vb = domelem[i].firstcpi;
      for ( j = 0; j < vi; j++ ) {
        ei = nncpi[domelcpind[vb+j]];
        if ( ei >= 0 )
          for ( k = 0; k < vi; k++ ) {
            eb = nncpi[domelcpind[vb+k]];
            if ( eb > ei )
              pkn_TMBElemSet ( nzcdistr, eb, ei );
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
} /*_g2mbl_SetupAltBlockHessianProfiled*/

boolean _g2mbl_SetupAltBlockDescription ( mesh_lmt_optdata *d, int nbl )
{
  void         *sp;
  int          nv, ndomelems, bndomelems, nbcpi;
  int          size;
  block_desc   *bd;
  meshdom_elem *domelem;
  int          *domelcpind;
  int          i, j, k, l, fcp, ncp;
  int          *bltag, mask;
  int          *belind, *bnncpi, *bvcpi, *bhprof;
  double       **bhrows, **blhrows;

  sp = pkv_GetScratchMemTop ();
  d->nbl = nbl;
  nv = d->nv;
  bd = d->block;
  ndomelems = d->ndomelems;
  domelem = d->domelem;
  domelcpind = d->domelcpind;
        /* count variable control points in each block */
  bltag = d->bltag;
  for ( i = nbcpi = 0; i < nbl; i++ ) {
    memset ( &bd[i], 0, sizeof(block_desc) );
    mask = 0x0001 << i;
    for ( j = k = 0;  j < nv;  j++ )
      if ( bltag[j] & mask )
        k ++;
    bd[i].ncp = k;
    nbcpi += k;
    bd[i].nvars = 3*k;
  }
        /* count domain elements associated with each block */
  bndomelems = 0;
  for ( j = 0; j < ndomelems; j++ )
    domelem[j].blmask = 0;
  for ( i = 0; i < nbl; i++ ) {
    mask = 0x0001 << i;
    for ( j = k = 0;  j < ndomelems;  j++ ) {
      fcp = domelem[j].firstcpi;
      ncp = domelem[j].ncp;
      for ( l = 0; l < ncp; l++ )
        if ( bltag[domelcpind[fcp+l]] & mask ) {
          domelem[j].blmask |= mask;
          k ++;
          break;
        }
    }
    bd[i].ndomelems = k;
    bndomelems += k;
  }
  size = bndomelems*sizeof(int) + nbcpi*(4*sizeof(int)+6*sizeof(double*)) +
         nv*nbl*sizeof(int);
  PKV_MALLOC ( d->belind, size );
  if ( !d->belind )
    goto failure;
  belind = d->belind;
  bnncpi = &belind[bndomelems];
  bvcpi = &bnncpi[nbl*nv];
  bhprof = &bvcpi[nbcpi];
  bhrows = (double**)&bhprof[3*nbcpi];
  blhrows = &bhrows[3*nbcpi];
  for ( i = j = k = 0; i < nbl; i++ ) {
    bd[i].domelind = &belind[j];
    j += bd[i].ndomelems;
    bd[i].nncpi = &bnncpi[i*nv];
    bd[i].vncpi = &bvcpi[k];
    bd[i].hprof = &bhprof[3*k];
    bd[i].hrows = &bhrows[3*k];
    bd[i].lhrows = &blhrows[3*k];
    k += bd[i].ncp;
  }
        /* find elements associated with each block */
  for ( j = 0; j < nbl; j++ )
    bd[j].ndomelems = 0;
  for ( i = 0; i < ndomelems; i++ ) {
    for ( j = 0; j < nbl; j++ ) {
      mask = 0x0001 << j;
      if ( domelem[i].blmask & mask ) {
        k = bd[j].ndomelems;
        bd[j].domelind[k] = i;
        bd[j].ndomelems = k+1;
      }
    }
  }
        /* setup profiles for all blocks */
  for ( i = 0; i < nbl; i++ ) {
    if ( !_g2mbl_SetupAltBlockHessianProfiled ( d, i ) )
      goto failure;
    bd[i].func = MYINFINITY;
    bd[i].nu = -1.0;
    bd[i].accurate = d->nkn1 == d->nkn2;
    bd[i].hblstat = 0;
  }

#ifdef _DEBUG
for ( i = 0; i < nbl; i++ )
  printf ( "block %d: ncp = %d, nvars = %d, ndomelems = %d, hsize = %d\n",
           i, bd[i].ncp, bd[i].nvars, bd[i].ndomelems, bd[i].hsize );
#endif

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*_g2mbl_SetupAltBlockDescription*/

boolean g2mbl_InitBlSurfaceOptAltBLMTd ( int nv, BSMvertex *mv, int *mvhei,
                                         point3d *mvcp, int nhe, BSMhalfedge *mhe,
                                         int nfac, BSMfacet *mfac, int *mfhei,
                                         byte *mkcp,
                                         double C, double dO, double dM,
                                         int nkn1, int nkn2, int nbl, void **data )
{
  void             *sp;
  mesh_lmt_optdata *d;
  int              blseed[G2MBL_MAX_BLOCKS];

  sp = pkv_GetScratchMemTop ();
  *data = d = NULL;
  PKV_MALLOC ( *data, sizeof(mesh_lmt_optdata) );
  d = *data;
  if ( !d )
    goto failure;
  memset ( d, 0, sizeof(mesh_lmt_optdata) );  /* clear all pointer types */
                                              /* and the eltypes array */
  if ( !_g2mbl_AssignMeshd ( d, nv, mv, mvhei, mvcp, nhe, mhe,
                             nfac, mfac, mfhei, mkcp ) )
    goto failure;
  if ( !_g2mbl_SetupHessianProfiled ( d, true ) )
    goto failure;
  if ( g2mbl_outputnzdistr )
    _g2mbl_OutputNZDistr ( d );
  if ( !_g2mbl_SetupHbl3x3d ( d ) )
    goto failure;
  if ( !_g2mbl_AllocBFArraysd ( d, nkn1, nkn2 ) )
    goto failure;
  if ( !_g2mbl_SetupElemConstd ( d, dM, dO, C ) )
    goto failure;
  PKV_MALLOC ( d->bltag, nv*sizeof(int) );
  if ( !d->bltag )
    goto failure;
  if ( !_g2mbl_SetupAltBlocks ( nv, mv, mvhei, nhe, mhe, nfac, mfac, mfhei,
                                d->nncpi, nbl, d->bltag, blseed ) )
    goto failure;
  if ( !_g2mbl_SetupAltBlockDescription ( d, nbl ) )
    goto failure;
  if ( !_g2mbl_SetupBlockRowsd ( d ) )
    goto failure;
  d->newpoint = true;
  if ( (d->accurate = nkn1 == nkn2) )
    d->ibl = 0;
  else
    d->ibl = nbl;
  d->itn = 0;
  d->next_entire = nbl;
  d->func = MYINFINITY;
  d->nu[0] = d->nu[1] = -1.0;
  d->nextblock = 0;
  d->last_ntn = false;

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  g2mbl_OptLMTDeallocated ( data );
  pkv_SetScratchMemTop ( sp );
  return false;
} /*g2mbl_InitBlSurfaceOptAltBLMTd*/

