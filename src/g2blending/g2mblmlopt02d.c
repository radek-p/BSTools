
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

#define NSEEDC 24

/* ///////////////////////////////////////////////////////////////////////// */
static void _g2mbl_MLMarkBlockVert ( int nv, BSMvertex *mv, int *mvhei,
                               int nhe, BSMhalfedge *mhe,
                               int nfac, BSMfacet *mfac, int *mfhei,
                               int nvcp, int *vncpi,
                               byte *mk, int *dist, pkv_queue *q )
{
  int i, j, k, d, v, vdeg, vfhe, f, fdeg, ffhe;

        /* mark the subset of vertices to define the distance */
  pkv_ResetQueue ( q );
  memset ( mk, 0, nv*sizeof(byte) );
  for ( i = 0; i < nv; i++ )
    dist[i] = nv;
  for ( i = 0; i < nvcp; i++ ) {
    j = vncpi[i];
    mk[j] = 0x03;
    dist[j] = 0;
    pkv_QueueInsert ( q, &j );
  }
  for ( ; !pkv_QueueEmpty ( q ); ) {
    pkv_QueueRemoveFirst ( q, &i );
    d = dist[i]+1;
    if ( d > 1 )
      break;
    vdeg = mv[i].degree;
    vfhe = mv[i].firsthalfedge;
    for ( j = 0; j < vdeg; j++ ) {
      f = mhe[mvhei[vfhe+j]].facetnum;
      fdeg = mfac[f].degree;
      ffhe = mfac[f].firsthalfedge;
      for ( k = 0; k < fdeg; k++ ) {
        v = mhe[mfhei[ffhe+k]].v0;
        if ( !(mk[v] & 0x02) ) {
          mk[v] |= 0x02;
          dist[v] = d;
          pkv_QueueInsert ( q, &v );
        }
      }
    }
  }
} /*_g2mbl_MLMarkBlockVert*/

static void _g2mbl_FindDist ( int nv, BSMvertex *mv, int *mvhei, byte *mk,
                              int nhe, BSMhalfedge *mhe,
                              int nfac, BSMfacet *mfac, int *mfhei,
                              int v, int *dist, boolean init, pkv_queue *q )
{
/* this procedure is similar to bsm_FindVertexDistance2, but the set of */
/* vertices is restricted to the ones with (mk[v] & 0x02) == 0x02 */
  int i, j, k, f, fhe, deg, ffhe, fdeg;
  int d;

  if ( init ) {
    for ( i = 0; i < nv; i++ )
      dist[i] = nv;
  }
  dist[v] = 0;
  pkv_ResetQueue ( q );
  pkv_QueueInsert ( q, &v );
  do {
    pkv_QueueRemoveFirst ( q, &i );
    d = dist[i]+1;
    deg = mv[i].degree;
    fhe = mv[i].firsthalfedge;
    for ( j = 0; j < deg; j++ ) {
      f = mhe[mvhei[fhe+j]].facetnum;
      fdeg = mfac[f].degree;
      ffhe = mfac[f].firsthalfedge;
      for ( k = 0; k < fdeg; k++ ) {
        v = mhe[mfhei[ffhe+k]].v0;
        if ( (mk[v] & 0x02) && dist[v] >= nv ) {
          dist[v] = d;
          pkv_QueueInsert ( q, &v );
        }
      }
    }
  } while ( !pkv_QueueEmpty ( q ) );
} /*_g2mbl_FindDist*/

static void _g2mbl_MLCountShr ( int nvcp, int *vncpi, int *dist1, int *dist2,
                                int sh, int *vc1, int *vc2 )
{
  int i, j, _vc1, _vc2;

  _vc1 = _vc2 = 0;
  for ( i = 0; i < nvcp; i++ ) {
    j = vncpi[i];
    if ( dist1[j]+sh < dist2[j] ) _vc1 ++;
    if ( dist2[j]-sh < dist1[j] ) _vc2 ++;
  }
  *vc1 = _vc1;
  *vc2 = _vc2;
} /*_g2mbl_MLCountShr*/

static int _g2mbl_MLCountShrb ( int nvcp, int *vncpi, int *dist1, int *dist2,
                                int sh1, int sh2 )
{
  int i, j, cnt;

  cnt = 0;
  for ( i = 0; i < nvcp; i++ ) {
    j = vncpi[i];
    if ( dist1[j]+sh1 >= dist2[j] &&
         dist2[j]-sh2 >= dist1[j] )
      cnt += 2;
    else if ( dist1[j]+sh1-1 == dist2[j] ||
              dist2[j]-sh2-1 == dist1[j] )
      cnt ++;
  }
  return cnt;
} /*_g2mbl_MLCountShrb*/

static int _g2mbl_MLFindShr ( int nvcp, int *vncpi, int *dist1, int *dist2,
                              int *sh1, int *sh2 )
{
  int _sh1, _sh2, _vc1, _vc2;

  _sh1 = 0;
  _g2mbl_MLCountShr ( nvcp, vncpi, dist1, dist2, _sh1, &_vc1, &_vc2 );
  while ( _vc2 > _vc1 && 2*_vc2 > nvcp ) {
    _sh1 --;
    _g2mbl_MLCountShr ( nvcp, vncpi, dist1, dist2, _sh1, &_vc1, &_vc2 );
  }
  while ( _vc1 > _vc2 && 2*_vc1 > nvcp ) {
    _sh1 ++;
    _g2mbl_MLCountShr ( nvcp, vncpi, dist1, dist2, _sh1, &_vc1, &_vc2 );
  }
  _sh2 = _sh1;
  while ( _vc2 > _vc1 && 2*_vc2 > nvcp ) {
    _sh2 --;
    _g2mbl_MLCountShr ( nvcp, vncpi, dist1, dist2, _sh2, &_vc1, &_vc2 );
  }
  *sh1 = _sh1;  *sh2 = _sh2;
  return _g2mbl_MLCountShrb ( nvcp, vncpi, dist1, dist2, _sh1, _sh2 );
} /*_g2mbl_MLFindShr*/

static void _g2mbl_MLFindBDist ( int nv, BSMvertex *mv, int *mvhei,
                                 int nhe, BSMhalfedge *mhe,
                                 int nfac, BSMfacet *mfac, int *mfhei,
                                 byte *mk, byte mask1, byte mask2,
                                 int nvbcp, int *nvbcpi, int *dist,
                                 pkv_queue *q )
{
  int  i, j, v, vdeg, vfhe, f, fdeg, ffhe, e, ee, d;

  pkv_ResetQueue ( q );
  for ( i = 0; i < nv; i++ )
    dist[i] = nv;
  for ( i = 0; i < nvbcp; i++ ) {
    j = nvbcpi[i];
    if ( (mk[j] & mask1) && !(mk[j] & mask2) ) {
      dist[j] = 0;
      pkv_QueueInsert ( q, &j );
    }
  }
  for ( ; !pkv_QueueEmpty ( q ); ) {
    pkv_QueueRemoveFirst ( q, &v );
    vdeg = mv[v].degree;
    vfhe = mv[v].firsthalfedge;
    d = dist[v]+1;
    for ( e = 0; e < vdeg; e++ ) {
      f = mhe[mvhei[vfhe+e]].facetnum;
      fdeg = mfac[f].degree;
      ffhe = mfac[f].firsthalfedge;
      for ( ee = 0; ee < fdeg; ee ++ ) {
        j = mhe[mfhei[ffhe+ee]].v0;
        if ( (mk[j] & 0x02) && dist[j] > d ) {
          dist[j] = d;
          pkv_QueueInsert ( q, &j );
        }
      }
    }
  }
} /*_g2mbl_MLFindBDist*/

static int _g2mbl_MLFindMDist ( int nv, BSMvertex *mv, int *mvhei,
                                int nhe, BSMhalfedge *mhe,
                                int nfac, BSMfacet *mfac, int *mfhei,
                                byte *mk, byte mask1, byte mask2,
                                int nvbcp, int *nvbcpi, int *dist,
                                pkv_queue *q )
{
  int  i, j, k;
  long ad;

  _g2mbl_MLFindBDist ( nv, mv, mvhei, nhe, mhe, nfac, mfac, mfhei,
                       mk, mask1, mask2, nvbcp, nvbcpi, dist, q );
  for ( i = k = 0, ad = 0;  i < nvbcp;  i++ ) {
    j = nvbcpi[i];
    if ( (mk[j] & mask1) && !(mk[j] & mask2) ) {
      ad += dist[j];
      k ++;
    }
  }
  return ad / k;
} /*_g2mbl_MLFindMDist*/

static boolean _g2mbl_MLMakeOverlap ( int nv, BSMvertex *mv, int *mvhei,
                                      int nhe, BSMhalfedge *mhe,
                                      int nfac, BSMfacet *mfac, int *mfhei,
                                      byte *mk, byte mask,
                                      int nvbcp, int vc, int *vnbcpi,
                                      pkv_queue *q )
{
  int i, j, v, vdeg, vfhe, f, fdeg, ffhe, e, ee;

  pkv_ResetQueue ( q );
  for ( i = 0; i < nv; i++ )
    mk[i] &= ~mask;
        /* store the vertices in the queue in this order */
  for ( i = 0; i < vc; i++ ) {
    j = vnbcpi[i];
    pkv_QueueInsert ( q, &j );
    mk[j] |= mask;
  }
  while ( vc < nvbcp && !pkv_QueueEmpty ( q ) ) {
    pkv_QueueRemoveFirst ( q, &v );
    vdeg = mv[v].degree;
    vfhe = mv[v].firsthalfedge;
    for ( e = 0; e < vdeg; e++ ) {
      f = mhe[mvhei[vfhe+e]].facetnum;
      fdeg = mfac[f].degree;
      ffhe = mfac[f].firsthalfedge;
      for ( ee = 0; ee < fdeg; ee++ ) {
        j = mhe[mfhei[ffhe+ee]].v0;
        if ( (mk[j] & 0x02) && !(mk[j] & mask) ) {
          pkv_QueueInsert ( q, &j );
          mk[j] |= mask;
          if ( mk[j] & 0x01 ) {
            vnbcpi[vc++] = j;
            if ( vc == nvbcp )
              return true;
          }
        }
      }
    }
  }
  return false;
} /*_g2mbl_MLMakeOverlap*/

boolean _g2mbl_MLDivideBlock ( int nv, BSMvertex *mv, int *mvhei,
                               int nhe, BSMhalfedge *mhe,
                               int nfac, BSMfacet *mfac, int *mfhei,
                               int nvcp, int *vncpi,
                               int *nvbcp1, int **vnbcp1, int *seed1,
                               int *nvbcp2, int **vnbcp2, int *seed2,
                               int step )
{
  void      *sp;
  int       _nvbcp1, *_vnbcp1, _nvbcp2, *_vnbcp2;
  int       seed[NSEEDC], _seed1, _seed2, *dist, *dist0, *dist1, *dist2;
  int       d, sh1, sh2;
  int       i, j, k, e, ee, f, v, vdeg, vfhe, fdeg, ffhe;
  int       vc1, vc2;
  pkv_queue *q, *q2;
  byte      *mk;

  sp = pkv_GetScratchMemTop ();
  q = q2 = NULL;  _vnbcp1 = _vnbcp2 = NULL;
  if ( nvcp < max(5,NSEEDC) )
    goto failure;
  q = pkv_InitQueue ( nv, sizeof(int) );
  *nvbcp1 = *nvbcp2 = _nvbcp1 = _nvbcp2 = (step*nvcp)/(2*step-1);
  PKV_MALLOC ( *vnbcp1, _nvbcp1*sizeof(int) );
  PKV_MALLOC ( *vnbcp2, _nvbcp2*sizeof(int) );
  if ( !q || !*vnbcp1 || !*vnbcp2 )
    goto failure;
  _vnbcp1 = *vnbcp1;
  _vnbcp2 = *vnbcp2;
  dist = pkv_GetScratchMemi ( 4*nv );
  mk = pkv_GetScratchMem ( nv*sizeof(byte) );
  if ( !dist || !mk )
    goto failure;
  dist0 = &dist[nv];
  dist1 = &dist0[nv];
  dist2 = &dist1[nv];

  _g2mbl_MLMarkBlockVert ( nv, mv, mvhei, nhe, mhe, nfac, mfac, mfhei,
                           nvcp, vncpi, mk, dist, q );

        /* choose the first candidate for a seed */
  seed[0] = vncpi[0];
  _g2mbl_FindDist ( nv, mv, mvhei, mk, nhe, mhe, nfac, mfac, mfhei,
                    seed[0], dist, true, q );
        /* test if the block is connected */
  for ( i = 0; i < nvcp; i++ )
    if ( dist[vncpi[i]] >= nv ) {
      printf ( "not connected\n" );
      goto failure;
    }
        /* choose more candidates for seeds */
  v = seed[0];
  for ( i = 1; i < NSEEDC; i++ ) {
    d = 0;
    for ( j = 0; j < nvcp; j++ ) {
      if ( dist[vncpi[j]] > d ) {
        v = vncpi[j];
        d = dist[v];
      }
    }
    seed[i] = v;
    _g2mbl_FindDist ( nv, mv, mvhei, mk, nhe, mhe, nfac, mfac, mfhei,
                      v, dist, false, q );
  }
        /* choose the seeds - the two candidates, whose regions */
        /* have the smallest common boundary */
  sh1 = sh2 = 0;
  _seed1 = seed[0];
  _seed2 = seed[1];
  d = 3*nv;
  for ( i = 0; i < NSEEDC; i++ ) {
    _g2mbl_FindDist ( nv, mv, mvhei, mk, nhe, mhe, nfac, mfac, mfhei,
                      seed[i], dist, true, q );
    for ( j = i+1; j < NSEEDC; j++ ) {
      _g2mbl_FindDist ( nv, mv, mvhei, mk, nhe, mhe, nfac, mfac, mfhei,
                        seed[j], dist0, true, q );
      k = _g2mbl_MLFindShr ( nvcp, vncpi, dist, dist0, &vc1, &vc2 );
      if ( k < d ) {
        d = k;
        _seed1 = seed[i];  sh1 = vc1;
        memcpy ( dist1, dist, nv*sizeof(int) );
        _seed2 = seed[j];  sh2 = vc2;
        memcpy ( dist2, dist0, nv*sizeof(int) );
      }
    }
  }
  if ( d >= 3*nv )
    goto failure;
  *seed1 = _seed1;
  *seed2 = _seed2;
/*
printf ( "seed1 = %d, seed2 = %d\n", _seed1, _seed2 );
*/
        /* assign the vertices of the Voronoi regions to the blocks */
  q2 = pkv_InitQueue ( nv, sizeof(int) );
  if ( !q2 )
    goto failure;
  pkv_ResetQueue ( q );
  for ( i = vc1 = vc2 = 0;  i < nvcp;  i++ ) {
    j = vncpi[i];
    if ( dist1[j]+sh1 < dist2[j] ) {
      pkv_QueueInsert ( q, &j );
      mk[j] |= 0x04;
      _vnbcp1[vc1++] = j;
    }
    else if ( dist2[j]-sh2 < dist1[j] ) {
      pkv_QueueInsert ( q2, &j );
      mk[j] |= 0x08;
      _vnbcp2[vc2++] = j;
    }
  }
        /* assign the vertices not assigned yet */
  while ( vc1 < _nvbcp1 || vc2 < _nvbcp2 ) {
    if ( vc1 < vc2 && !pkv_QueueEmpty ( q ) ) {
      pkv_QueueRemoveFirst ( q, &v );
      vdeg = mv[v].degree;
      vfhe = mv[v].firsthalfedge;
      for ( e = 0; e < vdeg; e++ ) {
        f = mhe[mvhei[vfhe+e]].facetnum;
        fdeg = mfac[f].degree;
        ffhe = mfac[f].firsthalfedge;
        for ( ee = 0; ee < fdeg; ee++ ) {
          j = mhe[mfhei[ffhe+ee]].v0;
          if ( (mk[j] & 0x02) && !(mk[j] & 0x0C) ) {
            pkv_QueueInsert ( q, &j );
            mk[j] |= 0x04;
            if ( mk[j] & 0x01 ) {
              _vnbcp1[vc1++] = j;
              if ( vc1 == _nvbcp1 )
                goto qq1;
            }
          }
        }
      }
    }
    else if ( vc2 <= vc1 && !pkv_QueueEmpty ( q2 ) ) {
      pkv_QueueRemoveFirst ( q2, &v );
      vdeg = mv[v].degree;
      vfhe = mv[v].firsthalfedge;
      for ( e = 0; e < vdeg; e++ ) {
        f = mhe[mvhei[vfhe+e]].facetnum;
        fdeg = mfac[f].degree;
        ffhe = mfac[f].firsthalfedge;
        for ( ee = 0; ee < fdeg; ee++ ) {
          j = mhe[mfhei[ffhe+ee]].v0;
          if ( (mk[j] & 0x02) && !(mk[j] & 0x0C) ) {
            pkv_QueueInsert ( q2, &j );
            mk[j] |= 0x08;
            if ( mk[j] & 0x01 ) {
              _vnbcp2[vc2++] = j;
              if ( vc2 == _nvbcp2 )
                goto qq1;
            }
          }
        }
      }
    }
    else break;
qq1:
    ;
  }
        /* extend the blocks to create overlaps */
  if ( vc1 < _nvbcp1 ) {
    if ( !_g2mbl_MLMakeOverlap ( nv, mv, mvhei, nhe, mhe, nfac, mfac, mfhei,
                                 mk, 0x04, _nvbcp1, vc1, _vnbcp1, q ) )
      goto failure;
    vc1 = _nvbcp1;
  }
  if ( vc2 < _nvbcp2 ) {
    if ( !_g2mbl_MLMakeOverlap ( nv, mv, mvhei, nhe, mhe, nfac, mfac, mfhei,
                                 mk, 0x08, _nvbcp2, vc2, _vnbcp2, q ) )
      goto failure;
    vc2 = _nvbcp2;
  }
        /* trim the overlaps and then extend them */
  d = _g2mbl_MLFindMDist ( nv, mv, mvhei, nhe, mhe, nfac, mfac, mfhei,
                           mk, 0x04, 0x08, _nvbcp1, _vnbcp1, dist1, q );
  d = (step*d)/(step-1);
  for ( i = 0; i < vc1; ) {
    if ( dist1[_vnbcp1[i]] > d )
      _vnbcp1[i] = _vnbcp1[--vc1];
    else
      i++;
  }
  if ( vc1 < _nvbcp1 ) {
    if ( !_g2mbl_MLMakeOverlap ( nv, mv, mvhei, nhe, mhe, nfac, mfac, mfhei,
                                 mk, 0x04, _nvbcp1, vc1, _vnbcp1, q ) )
      goto failure;
    vc1 = _nvbcp1;
  }
  d = _g2mbl_MLFindMDist ( nv, mv, mvhei, nhe, mhe, nfac, mfac, mfhei,
                           mk, 0x08, 0x04, _nvbcp2, _vnbcp2, dist2, q );
  d = (step*d)/(step-1);
  for ( i = 0; i < vc2; ) {
    if ( dist2[_vnbcp2[i]] > d )
      _vnbcp2[i] = _vnbcp2[--vc2];
    else
      i++;
  }
  if ( vc2 < _nvbcp2 ) {
    if ( !_g2mbl_MLMakeOverlap ( nv, mv, mvhei, nhe, mhe, nfac, mfac, mfhei,
                                 mk, 0x08, _nvbcp2, vc2, _vnbcp2, q ) )
      goto failure;
    vc2 = _nvbcp2;
  }

  *nvbcp1 = vc1;
  *nvbcp2 = vc2;

        /* verification - do the subblocks cover the block? */
/*
  for ( i = 0; i < nvcp; i++ )
    if ( !(mk[vncpi[i]] & 0x0C) )
      goto failure;
*/
  PKV_FREE ( q );
  PKV_FREE ( q2 );
  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  if ( q ) PKV_FREE ( q );
  if ( q2 ) PKV_FREE ( q2 );
  if ( *vnbcp1 ) PKV_FREE ( *vnbcp1 );
  if ( *vnbcp2 ) PKV_FREE ( *vnbcp2 );
  pkv_SetScratchMemTop ( sp );
  return false;
} /*_g2mbl_MLDivideBlock*/

