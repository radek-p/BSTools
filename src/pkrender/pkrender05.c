
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2015                                  */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#include "pkvaria.h"
#include "pkvthreads.h"
#include "pknum.h"
#include "pkgeom.h"
#include "multibs.h"
#include "camera.h"
#include "raybez.h"
#include "pkrender.h"

#include "pkrenderprivate.h"

boolean RendReset ( pkRenderer *rend )
{
  int i;

  rend->RenderingIsOn = false;
  DestroyTNodeTree ( rend->root );
  rend->root = NULL;
  for ( i = 0; i < rend->nobjects; i++ )
    switch ( rend->obj_tab[i].type ) {
  case obj_TRIANGLE:
      PKV_FREE ( rend->obj_tab[i].triang.trdata );
      break;
  case obj_BSPATCH:
      rbez_DestroyBezPatchTreed ( rend->obj_tab[i].bsp.ptree );
      break;
  case obj_RBSPATCH:
      rbez_DestroyRBezPatchTreed ( rend->obj_tab[i].rbsp.ptree );
      break;
  case obj_BEZCURVE:
      rbez_DestroyBezCurveTreed ( rend->obj_tab[i].bezc.ctree );
      break;
  case obj_RBEZCURVE:
      rbez_DestroyRBezCurveTreed ( rend->obj_tab[i].rbezc.ctree );
      break;
  default:
      break;
    }
  PKV_FREE ( rend->obj_tab );
  rend->obj_tab_length = 100;
  PKV_MALLOC ( rend->obj_tab, rend->obj_tab_length*sizeof(renderobj) );
  if ( !rend->obj_tab )
    rend->obj_tab_length = 0;
  rend->nobjects = 0;
  return rend->obj_tab_length > 0;
} /*RendReset*/

