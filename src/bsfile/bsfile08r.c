
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

#include "bsfprivate.h"

boolean bsf_ReadBSplinePatch4d ( int maxdeg, int maxlastknot, int maxncpoints,
                                 int *udeg, int *lastknotu, double *knotsu,
                                 int *vdeg, int *lastknotv, double *knotsv,
                                 boolean *closed_u, boolean *closed_v,
                                 int *pitch, point4d *cpoints,
                                 int *spdimen, boolean *rational,
                                 char *name, int *ident,
                                 bsf_UserReaders *readers )
{
  boolean _name, deg, knots_u, knots_v, c_points, dimen, _id, dep;
  int     ncpoints, dim, cpdimen;

  if ( bsf_nextsymbol != BSF_SYMB_BSPATCH )
    goto failure;
  bsf_GetNextSymbol ();
  if ( bsf_nextsymbol != BSF_SYMB_LBRACE )
    goto failure;
  bsf_GetNextSymbol ();

        /* nothing has been read in yet */
  _name = deg = knots_u = knots_v = c_points = *rational =
  dimen = _id = dep = false;
  ncpoints = 0;
  if ( name )
    *name = 0;
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

case BSF_SYMB_KNOTS_U:
      if ( knots_u )
        goto failure;
      bsf_GetNextSymbol ();
      if ( !bsf_ReadKnotSequenced ( maxlastknot, lastknotu, knotsu, closed_u ) )
        goto failure;
      knots_u = true;
      break;

case BSF_SYMB_KNOTS_V:
      if ( knots_v )
        goto failure;
      bsf_GetNextSymbol ();
      if ( !bsf_ReadKnotSequenced ( maxlastknot, lastknotv, knotsv, closed_v ) )
        goto failure;
      knots_v = true;
      break;

case BSF_SYMB_CPOINTS:
      if ( c_points )
        goto failure;
      bsf_GetNextSymbol ();
      ncpoints = bsf_ReadPointsd ( 4, maxncpoints, &cpoints[0].x, &cpdimen );
      c_points = true;
      if ( !ncpoints )
        goto failure;
      break;

case BSF_SYMB_RBRACE:
      bsf_GetNextSymbol ();
      if ( (*rational && dim < cpdimen-11) || (!*rational && dim < cpdimen) )
        goto failure;
      *spdimen = dim;
      if ( deg && knots_u && knots_v && c_points &&
           ncpoints == (*lastknotu-*udeg)*(*lastknotv-*vdeg) ) {
        *pitch = 4*(*lastknotv-*vdeg);
        return true;
      }
      else
        goto failure;

case BSF_SYMB_RATIONAL:
      *rational = true;
      bsf_GetNextSymbol ();
      break;

        /* trimmed patch domain */
case BSF_SYMB_TRIMMED:
      if ( !bsf_ReadTrimmedDomaind ( readers ) )
        goto failure;
      break;

        /* optional attributes */
case BSF_SYMB_CPOINTSMK:
      if ( !_bsf_ReadCPMark ( readers, maxncpoints ) )
        goto failure;
      break;

case BSF_SYMB_COLOR:
case BSF_SYMB_COLOUR:
      if ( !_bsf_ReadColour ( readers ) )
        goto failure;
      break;

        /* dependencies */
case BSF_SYMB_SPHERICAL_PRODUCT:
      if ( dep )
        goto failure;
      if ( !bsf_ReadDependencies ( readers ) )
        goto failure;
      dep = true;
      break;

default:
      goto failure;  /* anything else is forbidden */
    }
  }

failure:
  return false;
} /*bsf_ReadBSplinePatch4d*/

