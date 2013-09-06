
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2012, 2013                            */
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
boolean _g2mbl_CMPFindCoarseMeshBlock (
                     int fnv, int cnv, int rmnnz, index2 *rmnzi,
                     int nvbcp, int *vnbcp,
                     int *nwbcp, int **wnbcp )
{
  void *sp;
  char *fmkv, *cmkv;
  int  i, nwb;
  int  *_wnbcp;

  sp = pkv_GetScratchMemTop ();
        /* allocate arrays */
  fmkv = pkv_GetScratchMem ( fnv + cnv );
  if ( !fmkv )
    goto failure;
  cmkv = &fmkv[fnv];
        /* determine the vertices of the coarse mesh, which influence */
        /* only the vertices of the fine mesh, which belong to the */
        /* indicated block */
  memset ( fmkv, 0, fnv+cnv );
  for ( i = 0; i < nvbcp; i++ )
    fmkv[vnbcp[i]] = 1;
        /* searching the coarse mesh vertices, which influence only the */
        /* vertices of the indicated block of the fine mesh */
  for ( i = 0;  i < rmnnz;  i++ )
    if ( !fmkv[rmnzi[i].i] )
      cmkv[rmnzi[i].j] = 1;
  for ( i = nwb = 0; i < cnv; i++ )
    if ( !cmkv[i] )
      nwb ++;
  if ( !nwb )
    goto failure;
  *nwbcp = nwb;
  PKV_MALLOC ( _wnbcp, nwb*sizeof(int) );
  if ( !_wnbcp )
    goto failure;
  *wnbcp = _wnbcp;
  for ( i = nwb = 0;  i < cnv;  i++ )
    if ( !cmkv[i] )
      _wnbcp[nwb++] = i;
  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*_g2mbl_CMPFindCoarseMeshBlock*/

boolean _g2mbl_CMPFindBlockCGPmatrix (
                             int fnv, int cnv, int rmnnz, index2 *rmnzi,
                             int nvbcp, int *vnbcp, int nwbcp, int *wnbcp,
                             int *smnnz, index3 **smi )
{
  void         *sp;
  unsigned int *permut;
  int          *cols, *rows;
  int          i, j, k, _smnnz;
  index3       *_smi;

  sp = pkv_GetScratchMemTop ();
        /* allocate the arrays */
  permut = (unsigned int*)pkv_GetScratchMemi ( rmnnz+fnv+cnv+1 );
  if ( !permut )
    goto failure;
  rows = (int*)&permut[rmnnz];
  cols = &rows[fnv];

  if ( !pkn_SPMFindCols ( fnv, cnv, rmnnz, rmnzi, permut, false, cols ) )
    goto failure;
  for ( i = 0; i < fnv; i++ )
    rows[i] = -1;
  for ( i = 0; i < nvbcp; i++ )
    rows[vnbcp[i]] = i;
  for ( j = _smnnz = 0;  j < nwbcp;  j++ ) {
    k = wnbcp[j];
    for ( i = cols[k]; i < cols[k+1]; i++ )
      if ( rows[rmnzi[permut[i]].i] >= 0 )
        _smnnz ++;
  }
  *smnnz = _smnnz;
  PKV_MALLOC ( _smi, _smnnz*sizeof(index3) );
  if ( !_smi )
    goto failure;
  *smi = _smi;
  for ( j = _smnnz = 0;  j < nwbcp; j++ ) {
    k = wnbcp[j];
    for ( i = cols[k]; i < cols[k+1]; i++ )
      if ( rows[rmnzi[permut[i]].i] >= 0 ) {
        _smi[_smnnz].i = rows[rmnzi[permut[i]].i];
        _smi[_smnnz].j = j;
        _smi[_smnnz++].k = permut[i];
      }
  }

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*_g2mbl_CMPFindBlockCGPmatrix*/

boolean _g2mbl_CMPOrderCoarsePoints ( int nwcp, int nnz, index2 *nzci,
                                      int brmnnz, index3 *brmnzi, int blsize,
                                      int *hsize, int *hprof )
{
  void         *sp;
  int          nzcdsize;
  byte         *nzcdistr;
  unsigned int *rthrpermut, *permut, *ipermut;
  int          *rthrrows;
  vertex_desc  *vertd;
  int          i, j, k, l, ei, eb, vi;

  sp = pkv_GetScratchMemTop ();
  rthrpermut = (unsigned int*)pkv_GetScratchMemi ( nnz+nwcp+1 );
  if ( !rthrpermut )
    goto failure;
  rthrrows = (int*)&rthrpermut[nnz];
  if ( !pkn_SPMFindRows ( nwcp, nwcp, nnz, nzci, rthrpermut, false, rthrrows ) )
    goto failure;
  vertd = pkv_GetScratchMem ( nwcp*sizeof(vertex_desc) );
  if ( !vertd )
    goto failure;
  for ( i = 0; i < nwcp; i++ ) {
    vertd[i].vnum = i;
    vertd[i].firstsel = nwcp;
    vertd[i].fneigh = rthrrows[i];
    vertd[i].nneigh = vertd[i].nncn = rthrrows[i+1]-rthrrows[i];
  }

  nzcdsize = pkn_TMBSize ( nwcp );
  nzcdistr = pkv_GetScratchMem ( nzcdsize );
  if ( !nzcdistr )
    goto failure;
  memset ( nzcdistr, 0, nzcdsize );
  for ( i = 0; i < nnz; i++ )
    if ( nzci[i].i >= nzci[i].j )
      pkn_TMBElemSet ( nzcdistr, nzci[i].i, nzci[i].j );

#ifdef DEBUG
if ( g2mbl_outputnzdistr )
  g2mbl_outputnzdistr ( 0, 0, false, nwcp, nwcp, nzcdistr );
#endif

  permut = (unsigned int*)pkv_GetScratchMemi ( 2*nwcp );
  if ( !permut )
    goto failure;
  ipermut = &permut[nwcp];
  for ( i = 0; i < nwcp; i++ )
    permut[i] = i;
  j = 0;
  k = vertd[0].nncn;
  for ( i = 1; i < nwcp; i++ )
    if ( vertd[i].nncn < k ) { j = i;  k = vertd[i].nncn; }
  if ( j )
    { permut[0] = j;  permut[j] = 0; }
  vertd[j].firstsel = 0;
  for ( i = 0; i < nwcp-1; i++ ) {
    k = permut[i];
    ei = vertd[k].fneigh;
    eb = vertd[k].nneigh;
    for ( j = 0; j < eb; j++ ) {

if ( k != nzci[rthrpermut[ei+j]].i )
  exit ( 1 );

      l = nzci[rthrpermut[ei+j]].j;
      if ( pkn_TMBTestAndClear ( nzcdistr, k, l ) ) {
        vertd[l].firstsel = min ( vertd[l].firstsel, i );
        vertd[l].nncn --;
      }
    }
        /* selection */
    k = i+1;
    vi = permut[k];
    ei = vertd[vi].nncn;
    eb = vertd[vi].firstsel;
    for ( j = i+2; j < nwcp; j++ ) {
      l = permut[j];
      if ( vertd[l].firstsel < eb ||
           (vertd[l].firstsel == eb && vertd[l].nncn < ei) ) {
        k = j;
        ei = vertd[l].nncn;
        eb = vertd[l].firstsel;
      }
    }
    if ( k != i+1 )
      { l = permut[i+1];  permut[i+1] = permut[k];  permut[k] = l; }
    if ( vertd[permut[i+1]].firstsel > i+1 )
      vertd[permut[i+1]].firstsel = i+1;
  }
        /* setup the coarse Hessian matrix profile */
  for ( i = k = 0;  i < nwcp;  i++, k += blsize ) {
    j = permut[i];
    ipermut[j] = i;
    hprof[k] = blsize*vertd[j].firstsel;
    for ( l = 1; l < blsize; l++ )
      hprof[k+l] = hprof[k];
  }
#ifdef DEBUG
  for ( i = 0; i < nnz; i++ ) {
    nzci[i].i = ipermut[nzci[i].i];
    nzci[i].j = ipermut[nzci[i].j];
  }
#endif
  *hsize = pkn_NRBArraySize ( blsize*nwcp, hprof );
        /* reorder the refinement submatrix columns */
  for ( i = 0; i < brmnnz; i++ )
    brmnzi[i].j = ipermut[brmnzi[i].j];

#ifdef DEBUG
memset ( nzcdistr, 0, nzcdsize );
for ( i = 0; i < nnz; i++ )
  if ( ipermut[nzci[rthrpermut[i]].i] >= ipermut[nzci[rthrpermut[i]].j] )
    pkn_TMBElemSet ( nzcdistr, ipermut[nzci[rthrpermut[i]].i],
                               ipermut[nzci[rthrpermut[i]].j] );
if ( g2mbl_outputnzdistr )
  g2mbl_outputnzdistr ( 0, 0, false, nwcp, nwcp, nzcdistr );
#endif

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*_g2mbl_CMPOrderCoarsePoints*/

