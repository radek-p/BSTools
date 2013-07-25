
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2010                                  */
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

boolean bsm_DoublingNum ( int inv, BSMvertex *imv, int *imvhei,
                          int inhe, BSMhalfedge *imhe,
                          int infac, BSMfacet *imfac, int *imfhei,
                          int *onv, int *onhe, int *onfac )
{
  void *sp;
  int  ei, eb, vi, vb, vout;
  char *vtag, *ftag;
  int  i, d;

  sp = pkv_GetScratchMemTop ();
  vtag = pkv_GetScratchMem ( inv );
  ftag = pkv_GetScratchMem ( infac );
  if ( !vtag || !ftag )
    goto failure;
  bsm_TagMesh ( inv, imv, imvhei, inhe, imhe, infac, imfac, imfhei,
                vtag, ftag, &vi, &vb, &ei, &eb );
        /* the number of output facets */
  *onfac = infac + ei + eb + inv;
        /* the number of output halfedges */
  *onhe = 8*ei + 6*eb + 2*vb;
        /* the number of output verticess */
  vout = 0;
  for ( i = 0; i < inv; i++ ) {
    d = imv[i].degree;
    if ( vtag[i] )
      vout += d + 2;
    else
      vout += d;
  }
  *onv = vout;
  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*bsm_DoublingNum*/

int bsm_DoublingMatSize ( int inv, BSMvertex *imv, int *imvhei,
                          int inhe, BSMhalfedge *imhe,
                          int infac, BSMfacet *imfac, int *imfhei )
{
  void *sp;
  int  vi, vb, ei, eb, ndmat, i;
  char *vtag, *ftag;

  sp = pkv_GetScratchMemTop ();
  vtag = pkv_GetScratchMem ( inv );
  ftag = pkv_GetScratchMem ( infac );
  if ( vtag && ftag ) {
    bsm_TagMesh ( inv, imv, imvhei, inhe, imhe, infac, imfac, imfhei,
                  vtag, ftag, &vi, &vb, &ei, &eb );
    ndmat = 0;
    for ( i = 0; i < inv; i++ ) {
      if ( vtag[i] )
        ndmat += imv[i].degree + 2;
      else
        ndmat += imv[i].degree;
    }
  }
  else
    ndmat = 0;
  pkv_SetScratchMemTop ( sp );
  return ndmat;
} /*bsm_DoublingMatSize*/

