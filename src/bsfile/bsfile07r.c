
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2009, 2013                            */
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

boolean bsf_ReadBezierPatch4d ( int maxdeg, int *udeg, int *vdeg,
                                int *pitch, point4d *cpoints,
                                int *spdimen, boolean *rational,
                                byte *mk, char *name )
{
  boolean _name, deg, c_points, dimen, _mk;
  int     ncpoints, nmk, dim, cpdimen;

  if ( bsf_nextsymbol != BSF_SYMB_BPATCH )
    goto failure;
  bsf_GetNextSymbol ();
  if ( bsf_nextsymbol != BSF_SYMB_LBRACE )
    goto failure;
  bsf_GetNextSymbol ();

        /* nothing has been read in yet */
  _name = deg = c_points = dimen = *rational = _mk = false;
  ncpoints = 0;
  if ( name )
    *name = 0;
  nmk = -1;
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

case BSF_SYMB_DEGREE:
      if ( deg )
        goto failure;
      bsf_GetNextSymbol ();
      if ( !bsf_ReadPatchDegree ( maxdeg, udeg, vdeg ) )
        goto failure;
      deg = true;
      break;

case BSF_SYMB_DIM:
      if ( dimen )
        goto failure;
      bsf_GetNextSymbol ();
      if ( !bsf_ReadSpaceDim ( 4, &dim ) )
        goto failure;
      dimen = true;
      break;

case BSF_SYMB_CPOINTS:
      if ( c_points )
        goto failure;
      bsf_GetNextSymbol ();
      ncpoints = bsf_ReadPointsd ( 4, (maxdeg+1)*(maxdeg+1),
                                   &cpoints[0].x, &cpdimen );
      c_points = true;
      if ( !ncpoints )
        goto failure;
      break;

case BSF_SYMB_CPOINTSMK:
      if ( _mk )
        goto failure;
      bsf_GetNextSymbol ();
      nmk = bsf_ReadPointsMK ( (maxdeg+1)*(maxdeg+1), mk );
      _mk = true;
      break;

case BSF_SYMB_RBRACE:
      bsf_GetNextSymbol ();
      if ( nmk >= 0 && nmk != ncpoints )
        goto failure;
      if ( (*rational && dim < cpdimen+1) || (!*rational && dim < cpdimen) )
        goto failure;  
      *spdimen = dim;
      if ( deg && c_points && ncpoints == (*udeg+1)*(*vdeg+1) ) {
        *pitch = 4*(*vdeg+1);
        if ( !_mk && mk )
          memset ( mk, 0, ncpoints );
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
      goto failure;
    }
  }

failure:
  return false;
} /*bsf_ReadBezierPatch4d*/

