
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
#define AVE_VERTN 64  /* upper estimate for now */

typedef struct {
    int k;
    int left, right;
  } tree_el;

static int _g2mbl_MLTreeCountd ( tree_el *tree, int root )
{
  int l, r;

  if ( root >= 0 ) {
    l = _g2mbl_MLTreeCountd ( tree, tree[root].left );
    r = _g2mbl_MLTreeCountd ( tree, tree[root].right );
    return 1+l+r;
  }
  else
    return 0;
} /*_g2mbl_MLTreeCountd*/

static void _g2mbl_MLTreeOutputAscending ( tree_el *tree, int root,
                                           int **destp )
{
  if ( root >= 0 ) {
    _g2mbl_MLTreeOutputAscending ( tree, tree[root].left, destp );
    **destp = tree[root].k;
    (*destp) ++;
    _g2mbl_MLTreeOutputAscending ( tree, tree[root].right, destp );
  }
} /*_g2mbl_MLTreeOutputAscending*/

boolean _g2mbl_MLFindVCPNeighboursd ( mesh_ml_optdata *d,
                                      nzHbl_rowdesc **vn, int **vni )
{
  void          *sp;
  int           nv, nvcp;
  nzHbl_rowdesc *_vn;
  int           *_vni;
  meshdom_elem  *domelem;
  int           ndomel, *domelcpind, *domelind;
  int           *vncpi, *nncpi;
  int           i, el, fcp, ncp, j, k, jv, kv, vnilgt;
  tree_el       *tree;
  int           *vroot, ntreeel, ltreeel, *destp;

#define TREE_INSERT(v1,v2) { \
  int *vp; \
  vp = &vroot[v1];          /* select the tree */ \
  while ( *vp >= 0 ) {      /* search in the tree */ \
    if ( tree[*vp].k == v2 ) \
      goto out;             /* already present - no copies */ \
    else if ( tree[*vp].k > v2 ) \
      vp = &tree[*vp].left; \
    else \
      vp = &tree[*vp].right; \
  } \
  if ( ltreeel >= ntreeel ) /* not enough storage space allocated */ \
    goto failure; \
  *vp = ltreeel++;          /* allocate tree element and insert */ \
  tree[*vp].k = v2; \
  tree[*vp].left = tree[*vp].right = -1; \
out: ; \
}

  sp = pkv_GetScratchMemTop ();
  nv         = d->nv;
  nvcp       = d->nvcp;
  PKV_MALLOC ( *vn, nvcp*sizeof(nzHbl_rowdesc) );
  if ( !*vn )
    goto failure;
  _vn = *vn;

  domelem    = d->domelem;
  domelcpind = d->domelcpind;
  ndomel     = d->bd[0].ndomel;
  domelind   = d->bd[0].domelind;
  vncpi      = d->bd[0].vncpi;

  nncpi = pkv_GetScratchMemi ( nv );
  vroot = pkv_GetScratchMemi ( nvcp );
  if ( !nncpi || !vroot )
    goto failure;
  for ( i = 0; i < nv; i++ )
    nncpi[i] = -1;
  for ( i = 0; i < nvcp; i++ ) {
    nncpi[vncpi[i]] = i;
    vroot[i] = -1;      /* all binary search trees initially empty */
  }
  ntreeel = AVE_VERTN*nvcp;
  tree = pkv_GetScratchMem ( ntreeel*sizeof(tree_el) );
  if ( !tree )
    goto failure;
  ltreeel = 0;
        /* now lookup the elements and insert neighbour indices to the trees */
  for ( i = 0; i < ndomel; i++ ) {
    el = domelind[i];
    fcp = domelem[el].firstcpi;
    ncp = domelem[el].ncp;
    for ( j = 0; j < ncp; j++ ) {
      jv = nncpi[domelcpind[fcp+j]];
      if ( jv >= 0 ) {
        for ( k = 0; k < ncp; k++ ) {
          kv = nncpi[domelcpind[fcp+k]];
          if ( kv > jv )
            TREE_INSERT ( kv, jv )
        }
      }
    }
  }
        /* count the neighbours of each vertex, i.e. the elements in its tree */
  for ( i = 0; i < nvcp; i++ )
    _vn[i].nhbl = _g2mbl_MLTreeCountd ( tree, vroot[i] ) + 1;
  _vn[0].firsthbl = 0;
  for ( i = 1; i < nvcp; i++ )
    _vn[i].firsthbl = _vn[i-1].firsthbl+_vn[i-1].nhbl;
  vnilgt = _vn[nvcp-1].firsthbl+_vn[nvcp-1].nhbl;
  PKV_MALLOC ( *vni, 2*vnilgt*sizeof(int) );
  if ( !*vni )
    goto failure;
  _vni = *vni;
        /* transcript the neighbours from the trees to the *vni array */
  for ( i = 0; i < nvcp; i++ ) {
    destp = &_vni[_vn[i].firsthbl];
    _g2mbl_MLTreeOutputAscending ( tree, vroot[i], &destp );
    *destp = i;
  }
  d->Hblsize = vnilgt;
  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*_g2mbl_MLFindVCPNeighboursd*/

boolean _g2mbl_MLSetupBlockCGHessiand ( mesh_ml_optdata *d, int bn )
{
  void          *sp;
  int           Hblsize, i, j, fnz, nnz, fnz1, nnz1;
  mlblock_desc  *bd, *bd1;
  int           nbcp1, *vncpi1, *ind1, *cHbl1, *tHbl1;
  int           nbcp, *vncpi, *ind, *cHbl, *tHbl;
  nzHbl_rowdesc *iHbl, *iHbl1;

  sp = pkv_GetScratchMemTop ();
  bd = &d->bd[bn];
  if ( bn == 0 ) {
    if ( !_g2mbl_MLFindVCPNeighboursd ( d, &bd->iHbl, &bd->cHbl ) ) {
printf ( "%s\n", ERRMSG_10 );
      goto failure;
    }
    Hblsize = d->Hblsize;
    PKV_MALLOC ( d->Hbl, Hblsize*9*sizeof(double) );
    if ( !d->Hbl ) {
printf ( "%s\n", ERRMSG_1 );
      goto failure;
    }
    bd->nHbl = Hblsize;
    tHbl = bd->tHbl = &bd->cHbl[Hblsize];
    for ( i = 0; i < Hblsize; i++ )
      tHbl[i] = 9*i;
  }
  else {
    bd1 = &d->bd[(bn-1)/2];  /* number of the block 1 level up in the hierarchy */
    nbcp   = bd->nvcp;
    vncpi  = bd->vncpi;
    nbcp1  = bd1->nvcp;
    vncpi1 = bd1->vncpi;
    iHbl1  = bd1->iHbl;
    cHbl1  = bd1->cHbl;
    tHbl1  = bd1->tHbl;
    ind  = pkv_GetScratchMemi ( nbcp );
    ind1 = pkv_GetScratchMemi ( nbcp1 );
    if ( !ind || !ind1 ) {
printf ( "%s\n", ERRMSG_0 );
      goto failure;
    }
        /* vncpi1 contains an increasing sequence of integers (indices of */
        /* control points of the block above) and vncpi contains a subsequence */
        /* the first task is to find the indices cross-referencing the arrays */
    for ( i = j = 0;  i < nbcp1;  i++ )
      if ( vncpi[j] == vncpi1[i] ) {
        ind[j] = i;
        ind1[i] = j ++;
      }
      else
        ind1[i] = -1;
        /* with that the job is to select the rows and vertices of the Hessian */
        /* first, count the nonzero 3x3 blocks in the selected rows and columns */
    Hblsize = 0;
    for ( i = 0; i < nbcp1; i++ )
      if ( ind1[i] >= 0 ) {
        fnz1 = iHbl1[i].firsthbl;
        nnz1 = iHbl1[i].nhbl;
        for ( j = 0; j < nnz1; j++ )
          if ( ind1[cHbl1[fnz1+j]] >= 0 )
            Hblsize ++;
      }
    PKV_MALLOC ( bd->iHbl, nbcp*sizeof(nzHbl_rowdesc) );
    PKV_MALLOC ( bd->cHbl, 2*Hblsize*sizeof(int) );
    if ( !bd->iHbl || !bd->cHbl ) {
printf ( "%s\n", ERRMSG_1 );
      goto failure;
    }
    iHbl = bd->iHbl;
    cHbl = bd->cHbl;
    tHbl = bd->tHbl = &cHbl[Hblsize];
        /* now select the blocks, i.e. store their indices in the arrays */
    for ( i = fnz = 0;  i < nbcp1;  i++ )
      if ( ind1[i] >= 0 ) {
        fnz1 = iHbl1[i].firsthbl;
        nnz1 = iHbl1[i].nhbl;
        for ( j = nnz = 0;  j < nnz1;  j++ )
          if ( ind1[cHbl1[fnz1+j]] >= 0 ) {
            cHbl[fnz+nnz] = ind1[cHbl1[fnz1+j]];
            tHbl[fnz+nnz] = tHbl1[fnz1+j];
            nnz ++;
          }
        iHbl[ind1[i]].firsthbl = fnz;
        iHbl[ind1[i]].nhbl = nnz;
        fnz += nnz;
      }
  }
  bd->nHbl = Hblsize;

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*_g2mbl_MLSetupBlockCGHessiand*/

