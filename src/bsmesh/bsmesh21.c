
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

/* ////////////////////////////////////////////////////////////////////////// */
static void _bsm_rTagBZVert ( int nv, BSMvertex *mv, int *mvhei,
                              int nhe, BSMhalfedge *mhe, int v, char d,
                              char *vtag )
{
  char l;
  int  j, deg, fhe, u;

  l = vtag[v]+1;
  if ( l < d ) {
    deg = mv[v].degree;
    fhe = mv[v].firsthalfedge;
    for ( j = 0; j < deg; j++ ) {
      u = mhe[mvhei[fhe+j]].v1;
      if ( vtag[u] > l ) {
        vtag[u] = l;
        _bsm_rTagBZVert ( nv, mv, mvhei, nhe, mhe, u, d, vtag );
      }
    }
  }
} /*_bsm_rTagBZVert*/

void bsm_TagBoundaryZoneVertices ( int nv, BSMvertex *mv, int *mvhei,
                                   int nhe, BSMhalfedge *mhe,
                                   char d, char *vtag )
{
  int i, deg, fhe;

  for ( i = 0; i < nv; i++ )
    vtag[i] = d;
  for ( i = 0; i < nv; i++ ) {
    deg = mv[i].degree;
    fhe = mv[i].firsthalfedge;
    if ( mhe[mvhei[fhe+deg-1]].otherhalf < 0 ) {
      vtag[i] = 0;
      _bsm_rTagBZVert ( nv, mv, mvhei, nhe, mhe, i, d, vtag );
    }
  }
} /*bsm_TagBoundaryZoneVertices*/

