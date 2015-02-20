
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

static boolean DoLineA ( void *usrdata, int4 *jobnum )
{
  pkRenderer   *rend;
  int          nthr, xmax, x;
  double       y;
  ray3d        ray;
  unsigned int r, g, b;

  rend = (pkRenderer*)usrdata;
  nthr = rend->nthr;
  xmax = rend->CPos.xmin+rend->CPos.width;
  y = (double)rend->y;
  for ( x = rend->CPos.xmin + jobnum->x;
        x < xmax;
        x += nthr ) {
      /* get the ray */
    CameraRayOfPixeld ( &rend->CPos, (double)x, y, &ray );
      /* compute the pixel colour */
    GetPixelColour ( rend, &ray, &r, &g, &b );
    rend->SetPixel ( rend->private_data, x, rend->y, r, g, b );
  }
  return true;
} /*DoLineA*/

int RenderLineA ( pkRenderer *rend )
{
  short int     nthr;
  int4          jobsize;
  boolean       success;

  if ( !rend->RenderingIsOn )
    return rend->y;

  nthr = rend->nthr;
  if ( nthr == 1 ) {
    jobsize.x = 0;
    DoLineA ( (void*)rend, &jobsize );
  }
  else {
    jobsize.x = nthr;
    pkv_SetPThreadsToWork ( 1, &jobsize, nthr, PKRENDER_STACK, PKRENDER_SCRATCHMEM,
                            (void*)rend, DoLineA, NULL, NULL, &success );
  }
  rend->y ++;
  if ( rend->y >= rend->CPos.ymin+rend->CPos.height ) {
    if ( nthr > 1 )
      raybez_DisablePThreads ();
    rend->ticks2 = pkv_Toc ( &rend->tic );
    rend->RenderingIsOn = false;
  }
  return rend->y;
} /*RenderLineA*/

