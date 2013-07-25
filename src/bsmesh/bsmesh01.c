
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

void bsm_TagMesh ( int nv, BSMvertex *mv, int *mvhei,
                   int nhe, BSMhalfedge *mhe,
                   int nfac, BSMfacet *mfac, int *mfhei,
                   char *vtag, char *ftag,
                   int *vi, int *vb, int *ei, int *eb )
{
  int _vi, _vb, _ei, _eb;
  int i, j, k, d;

  _vi = _vb = _ei = _eb = 0;
  for ( i = 0; i < nv; i++ ) {
    j = mv[i].firsthalfedge;
    d = mv[i].degree;
    if ( mhe[mvhei[j+d-1]].otherhalf < 0 ) {
      vtag[i] = 1;
      _vb ++;
    }
    else {
      vtag[i] = 0;
      _vi ++;
    }
  }
  for ( i = 0; i < nhe; i++ ) {
    if ( mhe[i].otherhalf == -1 )
      _eb ++;
    else if ( mhe[i].otherhalf > i )
      _ei ++;
  }
  for ( i = 0; i < nfac; i++ ) {
    j = mfac[i].firsthalfedge;
    d = mfac[i].degree;
    ftag[i] = 0;
    for ( k = j; k < j+d; k++ )
      if ( mhe[mfhei[k]].otherhalf < 0 )
        break;
    if ( k < j+d )
      ftag[i] = 1;
    else
      ftag[i] = 0;
  }
  *vi = _vi;
  *vb = _vb;
  *ei = _ei;
  *eb = _eb;
} /*bsm_TagMesh*/

