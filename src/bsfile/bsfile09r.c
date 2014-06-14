
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2009, 2014                            */
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
#include "bsfile.h"

boolean bsf_ReadBSplineHoled ( int maxk, int *hole_k, double *knots,
                               point2d *domain_cp, point4d *hole_cp,
                               int *spdimen, boolean *rational,
                               char *name, int *ident,
                               bsf_UserReaders *readers )
{
  boolean _name, sides, dimen, domcp, surfcp, _id;
  int     ns, nk, ncp, lkn, dim, cpdimen;

  if ( bsf_nextsymbol != BSF_SYMB_BSHOLE )
    goto failure;
  bsf_GetNextSymbol ();
  if ( bsf_nextsymbol != BSF_SYMB_LBRACE )
    goto failure;
  bsf_GetNextSymbol ();

        /* nothing has ben read yet */
  _name = sides = dimen = domcp = surfcp = *rational = _id = false;
  if ( name )
    *name = 0;
  ns = nk = 0;
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

case BSF_SYMB_IDENT:
      if ( _id )
        goto failure;
      if ( !bsf_ReadIdent ( ident ) )
        goto failure;
      _id = true;
      break;

case BSF_SYMB_DIM:
      if ( dimen )
        goto failure;
      bsf_GetNextSymbol ();
      if ( !bsf_ReadSpaceDim ( 4, spdimen ) )
        goto failure;
      dimen = true;
      break;

case BSF_SYMB_SIDES:
      if ( sides )
        goto failure;
      bsf_GetNextSymbol ();
      if ( bsf_nextsymbol != BSF_SYMB_INTEGER )
        goto failure;
      if ( bsf_nextint < 3 || bsf_nextint > maxk )
        goto failure;
      *hole_k = ns = bsf_nextint;
      bsf_GetNextSymbol ();
      sides = true;
      break;

case BSF_SYMB_KNOTS:
      if ( !sides )
        goto failure;
      if ( nk >= ns )
        goto failure;
      bsf_GetNextSymbol ();
      if ( !bsf_ReadKnotSequenced ( 10, &lkn, &knots[nk*11], NULL ) )
        goto failure;
      if ( lkn != 10 )
        goto failure;
      nk ++;
      break;

case BSF_SYMB_DOMAIN:
      if ( domcp || !sides )
        goto failure;
      bsf_GetNextSymbol ();
      if ( bsf_nextsymbol != BSF_SYMB_CPOINTS )
        goto failure;
      bsf_GetNextSymbol ();
      ncp = bsf_ReadPointsd ( 2, 12*ns+1, &domain_cp[0].x, &dim );
      if ( ncp != 12*ns+1 )
        goto failure;
      domcp = true;
      break;

case BSF_SYMB_CPOINTS:
      if ( surfcp || !sides )
        goto failure;
      bsf_GetNextSymbol ();
      ncp = bsf_ReadPointsd ( 4, 12*ns+1, &hole_cp[0].x, &cpdimen );
      if ( ncp != 12*ns+1 )
        goto failure;
      surfcp = true;
      break;

case BSF_SYMB_RBRACE:
      bsf_GetNextSymbol ();
      if ( (*rational && *spdimen < cpdimen-1) ||
           (!*rational && *spdimen < cpdimen) )
        goto failure;  
      if ( sides && domcp && surfcp && nk == ns )
        return true;
      else
        goto failure;

case BSF_SYMB_RATIONAL:
      bsf_GetNextSymbol ();
      *rational = true;
      break;

        /* optional attributes */
case BSF_SYMB_CPOINTSMK:
      if ( !_bsf_ReadCPMark ( readers, 12*maxk+1 ) )
        goto failure;
      break;

case BSF_SYMB_COLOR:
case BSF_SYMB_COLOUR:
      if ( !_bsf_ReadColour ( readers ) )
        goto failure;
      break;

default:
      goto failure;  /* anything else is forbidden */
    }
  }

failure:
  return false;
} /*bsf_ReadBSplineHoled*/

