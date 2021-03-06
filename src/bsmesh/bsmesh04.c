
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2010, 2012                            */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "pkvaria.h"
#include "pknum.h"
#include "bsmesh.h"
#include "bsmprivate.h"


int _bsm_AveragingFVN ( int inv, BSMvertex *imv, int *imvhei,
                        int inhe, BSMhalfedge *imhe,
                        int infac, BSMfacet *imfac, int *imfhei,
                        int fn, char *vtag )
{
/* this procedure computes the number of vertices corresponding to */
/* a facet generated by the averaing operation; actually this is the number */
/* of sequences of inner facet vertices separated by boundary vertices */
/* bsm_TagMesh must be called before using this procedure */
  int  i, j, k, fhe, d, v0, v1;
  char b0, b1;

  d = imfac[fn].degree;
  fhe = imfac[fn].firsthalfedge;
  for ( i = 0; i < d; i++ ) {
    v0 = imhe[imfhei[fhe+i]].v0;
    if ( vtag[v0] )
      goto stage2;
  }
  return 1;  /* all vertices of this facet are inner */

stage2:
  b0 = vtag[v0];
  for ( j = k = 0;  j < d;  j++ ) {
    v1 = imhe[imfhei[fhe+i]].v1;
    b1 = vtag[v1];
    if ( b0 && !b1 )
      k ++;
    v0 = v1;  b0 = b1;
    i = i >= d-1 ? 0 : i+1;  /* increase i modulo d */
  }
  return k;
} /*_bsm_AveragingFVN*/

boolean bsm_AveragingNum ( int inv, BSMvertex *imv, int *imvhei,
                           int inhe, BSMhalfedge *imhe,
                           int infac, BSMfacet *imfac, int *imfhei,
                           int *onv, int *onhe, int *onfac )
{
  void *sp;
  int  vi, vb, ei, eb;
  int  _onv, _onhe;
  int  i;
  char *vtag, *ftag;

  sp = pkv_GetScratchMemTop ();
  vtag = pkv_GetScratchMem ( inv );
  ftag = pkv_GetScratchMem ( infac );
  if ( !vtag || !ftag )
    goto failure;
  bsm_TagMesh ( inv, imv, imvhei, inhe, imhe, infac, imfac, imfhei,
                vtag, ftag, &vi, &vb, &ei, &eb );
  if ( vb == 0 ) {  /* all vertices are inner */
    *onv = infac;
    *onhe = inhe;
    *onfac = inv;
  }
  else {
        /* number of output vertices */
    for ( i = _onv = 0;  i < infac;  i++ )
      _onv += _bsm_AveragingFVN ( inv, imv, imvhei, inhe, imhe,
                                  infac, imfac, imfhei, i, vtag );
    *onv = _onv;
        /* number of output facets */
    *onfac = vi;
        /* number of output halfedges */
    for ( i = _onhe = 0;  i < inhe;  i++ )
      if ( imhe[i].otherhalf >= 0 && !vtag[imhe[i].v1] )
        _onhe ++;
    *onhe = _onhe;
  }
  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*bsm_AveragingNum*/

int bsm_AveragingMatSize ( int inv, BSMvertex *imv, int *imvhei,
                           int inhe, BSMhalfedge *imhe,
                           int infac, BSMfacet *imfac, int *imfhei )
{
  void *sp;
  int  vi, vb, ei, eb, i, _namat;
  char *vtag, *ftag;

  sp = pkv_GetScratchMemTop ();
  vtag = pkv_GetScratchMem ( inv );
  ftag = pkv_GetScratchMem ( infac );
  if ( vtag && ftag ) {
    bsm_TagMesh ( inv, imv, imvhei, inhe, imhe, infac, imfac, imfhei,
                  vtag, ftag, &vi, &vb, &ei, &eb );
    _namat = 0;
    if ( vb == 0 ) { /* there are no boundary vertices and edges */
      for ( i = 0; i < infac; i++ )
        _namat += imfac[i].degree;
    }
    else {
      for ( i = 0; i < infac; i++ )
        _namat += imfac[i].degree *
                  _bsm_AveragingFVN ( inv, imv, imvhei, inhe, imhe,
                                      infac, imfac, imfhei, i, vtag );
    }
  }
  else
    _namat = 0;
  pkv_SetScratchMemTop ( sp );
  return _namat;
} /*bsm_AveragingMatSize*/

