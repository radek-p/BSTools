
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

boolean RendInit ( pkRenderer *rend, void *private_data, short nthr,
                   short maxwidth, short maxheight,
                   void (*SetPixel)( void *private_data,
                                     short x, short y, byte r, byte g, byte b ) )
{
  memset ( rend, 0, sizeof(pkRenderer) );
  rend->maxwidth = maxwidth;
  rend->maxheight = maxheight;
  rend->private_data = private_data;
  rend->nthr = nthr;
  rend->SetPixel = SetPixel;
  rend->c_shape_func = rend->d_shape_func = shapefunc_NONE;
  PKV_MALLOC ( rend->aabuf, 7*3*(3*maxwidth+4) );
  if ( !rend->aabuf )
    goto failure;
  rend->obj_tab_length = 100;
  PKV_MALLOC ( rend->obj_tab, rend->obj_tab_length*sizeof(renderobj) );
  if ( !rend->obj_tab )
    goto failure;
  rend->nobjects = 0;
  return true;

failure:
  if ( rend->obj_tab ) PKV_FREE ( rend->obj_tab );
  if ( rend->aabuf )   PKV_FREE ( rend->aabuf );
  return false;
} /*RendInit*/

void RendDestroy ( pkRenderer *rend )
{
  RendReset ( rend );
  PKV_FREE ( rend->obj_tab );
  PKV_FREE ( rend->aabuf );
} /*RendDestroy*/

