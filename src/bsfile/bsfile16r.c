
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2014                                  */
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

#include "bsfprivate.h"

boolean bsf_ReadPolyline4d ( int maxvert, int *nvert, point4d *vert,
                             int *spdimen, boolean *rational, boolean *closed,
                             char *name, int *ident,
                             bsf_UserReaders *readers )
{
  boolean _name, _dimen, _id, _vert;
  int     nv, dim, cpdimen;

  if ( bsf_nextsymbol != BSF_SYMB_POLYLINE )
    goto failure;
  bsf_GetNextSymbol ();
  if ( bsf_nextsymbol != BSF_SYMB_LBRACE )
    goto failure;
  bsf_GetNextSymbol ();

        /* nothing has been read in yet */
  _name = _dimen = _id = _vert = false;
  nv = 0;
  if ( name )   *name = 0;
  if ( *ident ) *ident = -1;
  *rational = false;
  for (;;) {
    switch ( bsf_nextsymbol ) {
case BSF_SYMB_NAME:
      if ( _name )
        goto failure;
      bsf_GetNextSymbol ();
      if ( bsf_nextsymbol == BSF_SYMB_STRING ) {
        strcpy ( name, bsf_name );
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
      if ( _dimen )
        goto failure;
      bsf_GetNextSymbol ();
      if ( !bsf_ReadSpaceDim ( 4, &dim ) )
        goto failure;
      _dimen = true;
      break;

case BSF_SYMB_POINTS:
      if ( _vert )
        goto failure;
      bsf_GetNextSymbol ();
      nv = bsf_ReadPointsd ( 4, maxvert, &vert[0].x, &cpdimen );
      if ( !nv )
        goto failure;
      *nvert = nv;
      _vert = true;
      break;

case BSF_SYMB_RBRACE:
      bsf_GetNextSymbol ();
      if ( !_dimen || !_vert )
        goto failure;
      if ( (*rational && dim < cpdimen-1) || (!*rational && dim < cpdimen) )
        goto failure;
      *spdimen = dim;
      return true;

case BSF_SYMB_RATIONAL:
      bsf_GetNextSymbol ();
      *rational = true;
      break;

        /* optional attributes */
case BSF_SYMB_CPOINTSMK:
      if ( !_bsf_ReadCPMark ( readers, maxvert ) )
        goto failure;
      break;

case BSF_SYMB_COLOR:
case BSF_SYMB_COLOUR:
      if ( !_bsf_ReadColour ( readers ) )
        goto failure;
      break;

default:
      goto failure;
    }
  }

failure:
  return false;
} /*bsf_ReadPolyline4d*/

