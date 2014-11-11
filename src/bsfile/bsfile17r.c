
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

boolean bsf_ReadTrimmedDomaind ( bsf_UserReaders *readers )
{
  void    *sp;
  boolean one_contour;
  double  *knots;
  point4d *points;
  int     maxp, maxlkn, maxdeg, deg, dim, lkn, nvert;
  boolean closed, rational, in_contour;
  mbs_polycurved elem;

  sp = pkv_GetScratchMemTop ();
  if ( bsf_nextsymbol != BSF_SYMB_TRIMMED )
    goto failure;
  bsf_GetNextSymbol ();
  if ( bsf_nextsymbol != BSF_SYMB_DOMAIN )
    goto failure;
  bsf_GetNextSymbol ();
  if ( bsf_nextsymbol != BSF_SYMB_LBRACE )
    goto failure;
  bsf_GetNextSymbol ();
  if ( bsf_nextsymbol == BSF_SYMB_LBRACE ) {
    one_contour = false;
    bsf_GetNextSymbol ();
  }
  else
    one_contour = true;
  if ( readers->BeginReader )
    readers->BeginReader ( readers->userData, BSF_TRIMMED_DOMAIN );
  in_contour = true;
  for (;;) {
    elem.ident = -1;
    switch ( bsf_nextsymbol ) {
case BSF_SYMB_BCURVE:
      maxdeg = readers->bc_maxdeg;
      maxp = maxdeg+1;
      points = pkv_GetScratchMem ( maxp*sizeof(point4d) );
      if ( !points )
        goto failure;
      elem.points = &points[0].x;
      if ( !bsf_ReadBezierCurve4d ( maxdeg, &deg, points, &dim,
                                    &rational, NULL, &elem.ident, readers ) )
        goto failure;
      if ( readers->TrimmedReader ) {
        elem.closing = bsf_nextsymbol == BSF_SYMB_RBRACE;
        if ( dim != 2 )
          goto failure;
        if ( rational )
          dim ++;
        elem.cpdimen = dim;
        pkv_Rearranged ( deg+1, dim, 4, dim, points );
        elem.degree = deg;
        elem.knots = NULL;
        elem.lastknot = -1;
        readers->TrimmedReader ( readers->userData, &elem );
      }
      pkv_SetScratchMemTop ( sp );
      break;
case BSF_SYMB_BSCURVE:
      maxdeg = readers->bsc_maxdeg;
      maxp = maxlkn = readers->bsc_maxlkn;
      points = pkv_GetScratchMem ( maxp*sizeof(point4d) );
      knots = elem.knots = pkv_GetScratchMemd ( maxlkn+1 );
      if ( !knots || !points )
        goto failure;
      elem.points = &points[0].x;
      if ( !bsf_ReadBSplineCurve4d ( maxdeg, maxlkn, maxp,
                                     &deg, &lkn, knots,
                                     &closed, points, &dim,
                                     &rational, NULL, &elem.ident, readers ) )
        goto failure;
      if ( readers->TrimmedReader ) {
        elem.closing = bsf_nextsymbol == BSF_SYMB_RBRACE;
        if ( dim != 2 )
          goto failure;
        if ( rational )
          dim ++;
        elem.cpdimen = dim;
        pkv_Rearranged ( lkn-deg, dim, 4, dim, points );
        elem.degree = deg;
        elem.lastknot = lkn;
        readers->TrimmedReader ( readers->userData, &elem );
      }
      pkv_SetScratchMemTop ( sp );
      break;
case BSF_SYMB_POLYLINE:
      maxp = readers->poly_maxvert;
      points = pkv_GetScratchMem ( maxp*sizeof(point4d) );
      if ( !points )
        goto failure;
      elem.points = &points[0].x;
      if ( !bsf_ReadPolyline4d ( maxp, &nvert, points, &dim,
                                 &rational, &closed, NULL, &elem.ident, readers) )
        goto failure;
      if ( readers->TrimmedReader ) {
        elem.closing = bsf_nextsymbol == BSF_SYMB_RBRACE;
        if ( dim != 2 )
          goto failure;
        if ( rational )
          dim ++;
        elem.cpdimen = dim;
        pkv_Rearranged ( nvert, dim, 4, dim, points );
        elem.knots = NULL;
        elem.degree = 1;
        elem.lastknot = nvert-1;
        readers->TrimmedReader ( readers->userData, &elem );
      }
      pkv_SetScratchMemTop ( sp );
      break;
case BSF_SYMB_LBRACE:
      bsf_GetNextSymbol ();
      if ( in_contour )     /* nesting too deep */
        goto failure;
      in_contour = true;    /* begin a next contour */
      break;
case BSF_SYMB_RBRACE:
      bsf_GetNextSymbol ();
      if ( one_contour || !in_contour ) {  /* end of the trimmed domain description */
        pkv_SetScratchMemTop ( sp );
        if ( readers->EndReader )
          readers->EndReader ( readers->userData, BSF_TRIMMED_DOMAIN, true );
        return true;
      }
      else {                /* end of a contour */
        in_contour = false;
        break;
      }
default:
      goto failure;
    }
  }

failure:
  pkv_SetScratchMemTop ( sp );
  if ( readers->EndReader )
    readers->EndReader ( readers->userData, BSF_TRIMMED_DOMAIN, false );
  return false;
} /*bsf_ReadTrimmedDomaind*/

