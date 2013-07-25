
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

/* ///////////////////////////////////////////////////////////////////////// */
/* count the variable vertices in the mesh */
int g2mbl_GetNvcp ( int nv, BSMvertex *mv, int *mvhei,
                    int nhe, BSMhalfedge *mhe,
                    int nfac, BSMfacet *mfac, int *mfhei,
                    byte *mkcp )
{
  void *sp;
  int  i, vi, vb, ei, eb;
  char *vtag, *ftag;

  sp = pkv_GetScratchMemTop ();
  vtag = pkv_GetScratchMem ( nv );
  ftag = pkv_GetScratchMem ( nfac );
  if ( !vtag || !ftag )
    goto failure;
  bsm_TagMesh ( nv, mv, mvhei, nhe, mhe, nfac, mfac, mfhei,
                vtag, ftag, &vi, &vb, &ei, &eb );
  for ( i = 0; i < nv; i++ )
    if ( vtag[i] == 1 ) {
      if ( mv[i].degree >= 4 )
        goto failure;
      else
        _g2mbl_TagBoundaryCondVert ( nv, mv, mvhei, nhe, mhe, i, vtag );
    }
    else if ( mv[i].degree < 3 || mv[i].degree > GH_MAX_K )
      goto failure;
  if ( mkcp ) {
    for ( i = 0; i < nv; i++ )
      if ( mkcp[i] )
        vtag[i] = 1;
  }
  for ( i = vi = 0; i < nv; i++ )
    if ( !vtag[i] )
      vi ++;
  pkv_SetScratchMemTop ( sp );
  return vi;

failure:
  pkv_SetScratchMemTop ( sp );
  return -1;
} /*g2mbl_GetNvcp*/

