
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
#include <math.h>
#include <malloc.h>
#include <string.h>
#include <ctype.h>

#include "pkvaria.h"
#include "pknum.h"
#include "pkgeom.h"
#include "camerad.h"
#include "multibs.h"
#include "bsmesh.h"
#include "bsfile.h"

boolean bsf_ReadBSMesh4d ( int maxnv, int maxnhe, int maxnfac,
                           int *degree,
                           int *nv, BSMvertex *mv, int *mvhei, point4d *vc,              
                           int *nhe, BSMhalfedge *mhe,      
                           int *nfac, BSMfacet *mfac, int *mfhei,
                           int *spdimen, boolean *rational,
                           byte *mkv, char *name )
{
  boolean _name, deg, vertices, halfedges, facets, dimen, _mkv;
  int     nitems, nmkv;
  int     i, k, d;

  if ( bsf_nextsymbol != BSF_SYMB_BSMESH )
    goto failure;
  bsf_GetNextSymbol ();
  if ( bsf_nextsymbol != BSF_SYMB_LBRACE )
    goto failure;
  bsf_GetNextSymbol ();

        /* nothing has been read in yet */
  _name = deg = vertices = halfedges = facets = dimen = *rational = _mkv = false;
  if ( name )
    *name = 0;
  *degree = 0;
  nmkv = -1;
  for (;;) {
    switch ( bsf_nextsymbol ) {
case BSF_SYMB_NAME:
      if ( _name )
        goto failure;
      bsf_GetNextSymbol ();
      if ( bsf_nextsymbol == BSF_SYMB_STRING ) {
        strcpy ( name, bsf_namebuffer );
        bsf_GetNextSymbol ();
      }
      else
        goto failure;
      _name = true;
      break;

case BSF_SYMB_DIM:
      if ( dimen )
        goto failure;
      bsf_GetNextSymbol ();
      if ( !bsf_ReadSpaceDim ( 4, spdimen ) )
        goto failure;
      dimen = true;
      break;

case BSF_SYMB_DEGREE:
      if ( deg )
        goto failure;
      bsf_GetNextSymbol ();
      if ( !bsf_ReadCurveDegree ( 255, degree ) )
        goto failure;
      deg = true;
      break;

case BSF_SYMB_VERTICES:
      if ( vertices )
        goto failure;
      bsf_GetNextSymbol ();
      if ( bsf_nextsymbol != BSF_SYMB_INTEGER )
        goto failure;
      nitems = bsf_nextint;
      if ( nitems < 0 || nitems > maxnv )
        goto failure;
      bsf_GetNextSymbol ();
      if ( bsf_nextsymbol != BSF_SYMB_LBRACE )
        goto failure;
      bsf_GetNextSymbol ();
          /* read in the vertices */
      k = 0;
      for ( i = 0; i < nitems; i++ ) {
        if ( bsf_nextsymbol != BSF_SYMB_LBRACE )
          goto failure;
        d = 0;
        mv[i].firsthalfedge = k;
        do {
          bsf_GetNextSymbol ();
          if ( k >= maxnhe || bsf_nextsymbol != BSF_SYMB_INTEGER )
            goto failure;
          mvhei[k++] = bsf_nextint;
          d ++;
          bsf_GetNextSymbol ();
        } while ( bsf_nextsymbol == BSF_SYMB_COMMA );
        mv[i].degree = d;
        if ( !bsf_ReadPointd ( 4, &vc[i].x, &d ) )
          goto failure;
        if ( bsf_nextsymbol != BSF_SYMB_RBRACE )
          goto failure;
        bsf_GetNextSymbol ();
        if ( i < nitems-1 ) {
          if ( bsf_nextsymbol != BSF_SYMB_COMMA )
            goto failure;
          bsf_GetNextSymbol ();
        }
      }
      if ( bsf_nextsymbol != BSF_SYMB_RBRACE )
        goto failure;
      bsf_GetNextSymbol ();
      *nv = nitems;
      vertices = true;
      break;

case BSF_SYMB_HALFEDGES:
      if ( halfedges )
        goto failure;
      bsf_GetNextSymbol ();
      if ( bsf_nextsymbol != BSF_SYMB_INTEGER )
        goto failure;
      nitems = bsf_nextint;
      if ( nitems < 0 || nitems > maxnhe )
        goto failure;
      bsf_GetNextSymbol ();
      if ( bsf_nextsymbol != BSF_SYMB_LBRACE )
        goto failure;
      bsf_GetNextSymbol ();
        /* read in the halfedges */
      for ( i = 0; i < nitems; i++ ) {
        if ( bsf_nextsymbol != BSF_SYMB_LBRACE )
          goto failure;
        bsf_GetNextSymbol ();
        if ( bsf_nextsymbol != BSF_SYMB_INTEGER )
          goto failure;
        mhe[i].v0 = bsf_nextint;
        bsf_GetNextSymbol ();
        if ( bsf_nextsymbol != BSF_SYMB_COMMA )
          goto failure;
        bsf_GetNextSymbol ();
        if ( bsf_nextsymbol != BSF_SYMB_INTEGER )
          goto failure;
        mhe[i].v1 = bsf_nextint;
        bsf_GetNextSymbol ();
        if ( bsf_nextsymbol != BSF_SYMB_COMMA )
          goto failure;
        bsf_GetNextSymbol ();
        if ( bsf_nextsymbol != BSF_SYMB_INTEGER )
          goto failure;
        mhe[i].facetnum = bsf_nextint;
        bsf_GetNextSymbol ();
        if ( bsf_nextsymbol != BSF_SYMB_COMMA )
          goto failure;
        bsf_GetNextSymbol ();
        if ( !bsf_ReadIntNumber ( &mhe[i].otherhalf ) )  /* this may be -1 or >= 0 */
          goto failure;
        if ( mhe[i].otherhalf < -1 )
          goto failure;
        if ( bsf_nextsymbol != BSF_SYMB_RBRACE )
          goto failure;
        bsf_GetNextSymbol ();
        if ( i < nitems-1 ) {
          if ( bsf_nextsymbol != BSF_SYMB_COMMA )
            goto failure;
          bsf_GetNextSymbol ();
        }
      }
      if ( bsf_nextsymbol != BSF_SYMB_RBRACE )
        goto failure;
      bsf_GetNextSymbol ();
      *nhe = nitems;
      halfedges = true;
      break;

case BSF_SYMB_FACETS:
      if ( facets )
        goto failure;
      bsf_GetNextSymbol ();
      if ( bsf_nextsymbol != BSF_SYMB_INTEGER )
        goto failure;
      nitems = bsf_nextint;
      if ( nitems < 0 || nitems > maxnfac )
        goto failure;
      bsf_GetNextSymbol ();
      if ( bsf_nextsymbol != BSF_SYMB_LBRACE )
        goto failure;
      bsf_GetNextSymbol ();
          /* read in the facets */
      k = 0;
      for ( i = 0; i < nitems; i++ ) {
        if ( bsf_nextsymbol != BSF_SYMB_LBRACE )
          goto failure;
        d = 0;
        mfac[i].firsthalfedge = k;
        do {
          bsf_GetNextSymbol ();
          if ( k >= maxnhe || bsf_nextsymbol != BSF_SYMB_INTEGER )
            goto failure;
          mfhei[k++] = bsf_nextint;
          d ++;
          bsf_GetNextSymbol ();
        } while ( bsf_nextsymbol == BSF_SYMB_COMMA );
        mfac[i].degree = d;
        if ( bsf_nextsymbol != BSF_SYMB_RBRACE )
          goto failure;
        bsf_GetNextSymbol ();
        if ( i < nitems-1 ) {
          if ( bsf_nextsymbol != BSF_SYMB_COMMA )
            goto failure;
          bsf_GetNextSymbol ();
        }
      }
      if ( bsf_nextsymbol != BSF_SYMB_RBRACE )
        goto failure;
      bsf_GetNextSymbol ();
      *nfac = nitems;
      facets = true;
      break;

case BSF_SYMB_CPOINTSMK:
      if ( _mkv )
        goto failure;
      bsf_GetNextSymbol ();
      nmkv = bsf_ReadPointsMK ( maxnv, mkv );
      _mkv = true;
      break;

case BSF_SYMB_RBRACE:
      bsf_GetNextSymbol ();
      if ( nmkv >= 0 && nmkv != *nv )
        goto failure;
      if ( vertices && halfedges && facets ) {
        if ( !_mkv && mkv )
          memset ( mkv, 0, *nv );
        return true;
      }
      else
        goto failure;

case BSF_SYMB_RATIONAL:
      *rational = true;
      bsf_GetNextSymbol ();
      break;

case BSF_SYMB_COLOR:
case BSF_SYMB_COLOUR:
      if ( !_bsf_ReadColour ( bsf_current_readers ) )
        goto failure;
      break;

default:
      goto failure;  /* anything else is forbidden */
    }
  }

failure:
  return false;
} /*bsf_ReadBSMesh4d*/

