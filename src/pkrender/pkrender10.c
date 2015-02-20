
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

boolean RendBegin ( pkRenderer *rend, boolean antialias, short nthr )
{
  if ( rend->root )
    DestroyTNodeTree ( rend->root );
  if ( !BuildObjectTree ( rend ) )
    return false;
  rend->nthr = nthr;
  if ( nthr > 1 )
    raybez_EnablePThreads ();
  else
    raybez_DisablePThreads ();
  FindMinMaxShapeFunc ( rend, 10 );
  pkv_Tic ( &rend->tic );
  rend->y = rend->CPos.ymin;
  rend->RenderingIsOn = true;
  if ( (rend->swAntialias = antialias) )
    InitRenderingAA ( rend );
  return true;
} /*RendBegin*/

boolean RendRestart ( pkRenderer *rend )
{
  if ( rend->RenderingIsOn ) {
    pkv_Tic ( &rend->tic );
    rend->y = rend->CPos.ymin;
    if ( rend->swAntialias )
      InitRenderingAA ( rend );
    return true;
  }
  else
    return false;
} /*RendRestart*/

int RenderLine ( pkRenderer *rend )
{
  if ( rend->swAntialias )
    return RenderLineAA ( rend );
  else
    return RenderLineA ( rend );
} /*RenderLine*/

