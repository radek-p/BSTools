
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2013, 2014                            */
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

#include "bsfprivate.h"

boolean bsf_ReadCamera ( CameraRecd *Camera, int *ident )
{
  boolean frame, pos, orient, depth, projdata, parallel, _id;
  int     width, height, xmin, ymin, dim;
  point3d position, angles;
  point2d depthrange;
  double  f, xi0, eta0, wdt, hgh, diag;
  int     dim_case;

  if ( bsf_nextsymbol != BSF_SYMB_CAMERA )
    goto failure;
  bsf_GetNextSymbol ();
  if ( bsf_nextsymbol != BSF_SYMB_LBRACE )
    goto failure;
  bsf_GetNextSymbol ();

        /* nothing has been read yet */
  frame = pos = orient = depth = projdata = _id = false;
  parallel = false;
  for (;;) {
    switch ( bsf_nextsymbol ) {
case BSF_SYMB_IDENT:
      if ( _id )
        goto failure;
      if ( !bsf_ReadIdent ( ident ) )
        goto failure;
      _id = true;
      break;

case BSF_SYMB_FRAME: /* pixel dimensions of the camera frame */
           /* two, three or four integers in braces are expected */
      if ( frame )
        goto failure;
      bsf_GetNextSymbol ();
      if ( bsf_nextsymbol != BSF_SYMB_LBRACE )
        goto failure;
      bsf_GetNextSymbol ();
      xmin = ymin = 0;  /* default value if absent in the file */
      if ( bsf_nextsymbol != BSF_SYMB_INTEGER )
        goto failure;
      width = bsf_nextint;
      bsf_GetNextSymbol ();
      if ( bsf_nextsymbol != BSF_SYMB_COMMA )
        goto failure;
      bsf_GetNextSymbol ();
      if ( bsf_nextsymbol != BSF_SYMB_INTEGER )
        goto failure;
      height = bsf_nextint;
      bsf_GetNextSymbol ();
      if ( bsf_nextsymbol == BSF_SYMB_COMMA ) {
        bsf_GetNextSymbol ();
        if ( bsf_nextsymbol == BSF_SYMB_INTEGER ) {
          xmin = bsf_nextint;
          bsf_GetNextSymbol ();
        }
        else
          goto failure;
      }
      if ( bsf_nextsymbol == BSF_SYMB_COMMA ) {
        bsf_GetNextSymbol ();
        if ( bsf_nextsymbol == BSF_SYMB_INTEGER ) {
          ymin = bsf_nextint;
          bsf_GetNextSymbol ();
        }
        else
          goto failure;
      }
      if ( bsf_nextsymbol != BSF_SYMB_RBRACE )
        goto failure;
      bsf_GetNextSymbol ();
      if ( width <= 0 || height <= 0 )  /* is it correct? */
        goto failure;
      frame = true;
      break;

case BSF_SYMB_POSITION:  /* viewer position */
      if ( pos )
        goto failure;
      bsf_GetNextSymbol ();
      if ( !bsf_ReadPointd ( 3, &position.x, &dim ) )
        goto failure;
      pos = true;
      break;

case BSF_SYMB_EULERANGLES:  /* orientation represented by Euler angles */
      if ( orient )
        goto failure;
      bsf_GetNextSymbol ();
        /* these are 3 real numbers in braces, */
        /* so a point reading procedure is used */
      if ( !bsf_ReadPointd ( 3, &angles.x, &dim ) )
        goto failure;
      orient = true;
      break;

case BSF_SYMB_DEPTH:  /* depth range, optional */
      if ( depth )
        goto failure;
      bsf_GetNextSymbol ();
      if ( !bsf_ReadPointd ( 2, &depthrange.x, &dim ) )
        goto failure;
      if ( depthrange.x <= 0.0 || depthrange.y <= depthrange.x ) /* is it correct? */
        goto failure;
      depth = true;
      break;

case BSF_SYMB_PARALLEL:  /* data for a parallel projection */
      if ( projdata )
        goto failure;
      bsf_GetNextSymbol ();
      if ( bsf_nextsymbol != BSF_SYMB_LBRACE )
        goto failure;
      bsf_GetNextSymbol ();
      if ( !bsf_ReadDoubleNumber ( &wdt ) )
        goto failure;
      if ( bsf_nextsymbol != BSF_SYMB_COMMA )
        goto failure;
      bsf_GetNextSymbol ();
      if ( !bsf_ReadDoubleNumber ( &hgh ) )
        goto failure;
      if ( bsf_nextsymbol != BSF_SYMB_COMMA )
        goto failure;
      bsf_GetNextSymbol ();
      if ( !bsf_ReadDoubleNumber ( &diag ) )
        goto failure;
      if ( bsf_nextsymbol != BSF_SYMB_COMMA )
        goto failure;
      bsf_GetNextSymbol ();
      if ( bsf_nextsymbol != BSF_SYMB_INTEGER )
        goto failure;
      dim_case = bsf_nextint;
      bsf_GetNextSymbol ();
      if ( bsf_nextsymbol != BSF_SYMB_RBRACE )
        goto failure;
      bsf_GetNextSymbol ();
      if ( wdt <= 0.0 || hgh <= 0.0 || diag <= 0.0 ||  /* is it correct? */
           dim_case < 0 || dim_case > 2 )
        goto failure;
      projdata = parallel = true;
      break;

case BSF_SYMB_PERSPECTIVE:  /* data for a perspective projection */
      if ( projdata )
        goto failure;
      bsf_GetNextSymbol ();
      if ( bsf_nextsymbol != BSF_SYMB_LBRACE )
        goto failure;
      bsf_GetNextSymbol ();
      if ( !bsf_ReadDoubleNumber ( &f ) )
        goto failure;
      if ( bsf_nextsymbol != BSF_SYMB_COMMA )
        goto failure;
      bsf_GetNextSymbol ();
      if ( !bsf_ReadDoubleNumber ( &xi0 ) )
        goto failure;
      if ( bsf_nextsymbol != BSF_SYMB_COMMA )
        goto failure;
      bsf_GetNextSymbol ();
      if ( !bsf_ReadDoubleNumber ( &eta0 ) )
        goto failure;
      if ( bsf_nextsymbol != BSF_SYMB_RBRACE )
        goto failure;
      bsf_GetNextSymbol ();
      if ( f <= 0.0 )
        goto failure;
      projdata = true;
      parallel = false;
      break;

case BSF_SYMB_RBRACE:  /* if everything has been read in, */
                       /* setup the camera and return */
      bsf_GetNextSymbol ();
      if ( !frame || !pos || !orient || !projdata )
        goto failure;
                       /* default values for the other parameters are */
                       /* used, i.e. aspect is 1.0, no upside down, */
                       /* magnification is 1 etc. An application will */
                       /* probably change the camera to suit its needs, */
                       /* but it is supposed to obtain a ready to use one */
      CameraInitFramed ( Camera, parallel, false,
                         width, height, xmin, ymin, 1.0,
                         depth ? 6 : 4 );
      if ( depth )
        CameraSetDepthRanged ( Camera, depthrange.x, depthrange.y );
      CameraInitPosd ( Camera );
      Camera->position = position;
      Camera->psi   = angles.x;
      Camera->theta = angles.y;
      Camera->phi   = angles.z;
      if ( parallel ) {
        Camera->vd.para.wdt = wdt;
        Camera->vd.para.wdt = hgh;
        Camera->vd.para.wdt = diag;
        Camera->vd.para.dim_case = dim_case;
      }
      else {
        Camera->vd.persp.f = f;
        Camera->vd.persp.xi0 = xi0;
        Camera->vd.persp.eta0 = eta0;
        Camera->vd.persp.dxi0 = Camera->vd.persp.deta0 = 0.0;
      }
      if ( CameraSetMappingd ( Camera ) )
        goto success;
      else
        goto failure;

default:
      goto failure;
    }
  }
success:
  return true;

failure:
  return false;
} /*bsf_ReadCamera*/

