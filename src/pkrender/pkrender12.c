
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

boolean OverThreshold1 ( byte *a, byte *b )
{
#define THR 5
  int  i;

  for ( i = 0; i < 3; i++ ) {
    if ( abs ( a[i]-b[i] ) > THR )
      return true;
  }
  return false;
#undef THR
} /*OverThreshold1*/

boolean OverThreshold2 ( byte *a, byte *b, byte *c, byte *d )
{
#define THR 5
  byte cmin, cmax;
  int  i;

  for ( i = 0; i < 3; i++ ) {
    if ( a[i] < b[i] ) { cmin = a[i];  cmax = b[i]; }
                  else { cmin = b[i];  cmax = a[i]; }
    if ( c[i] < d[i] ) { cmin = min ( cmin, c[i] );  cmax = max ( cmax, d[i] ); }
                  else { cmin = min ( cmin, d[i] );  cmax = max ( cmax, c[i] ); }
    if ( cmax-cmin > THR )
      return true;
  }
  return false;
#undef THR
} /*OverThreshold2*/

void Interpol2_3 ( byte *a, byte *b, byte *c )
{
  int i;

  for ( i = 0; i < 3; i++ )
    c[i] = (2*a[i]+b[i])/3;
} /*Interpol2_3*/

void GetSample ( pkRenderer *rend, double x, double y, byte *rgb )
{
  ray3d ray;
  unsigned int r, g, b;

  CameraRayOfPixeld ( &rend->CPos, x, y, &ray );
  GetPixelColour ( rend, &ray, &r, &g, &b );
  rgb[0] = (byte)r;
  rgb[1] = (byte)g;
  rgb[2] = (byte)b;
} /*GetSample*/

typedef struct {
    pkRenderer *rend;
    double     y;
    byte       *aaline;
  } aa_data;

static boolean DoSubpixelLine1 ( void *usrdata, int4 *jobnum )
{
  aa_data    *aadata;
  pkRenderer *rend;
  int        nthr, nthr9, xmin, xmax, x, i;
  double     y;
  byte       *aaline;

  aadata = (aa_data*)usrdata;
  rend = aadata->rend;
  nthr = rend->nthr;
  nthr9 = 9*nthr;
  xmin = rend->CPos.xmin-1+jobnum->x;
  xmax = rend->CPos.xmin+rend->CPos.width;
  y = aadata->y;
  aaline = aadata->aaline;
  for ( x = xmin, i = 9*jobnum->x;  x <= xmax;  x += nthr, i += nthr9 )
    GetSample ( rend, (double)x, y, &aaline[i] );
  return true;
} /*DoSubpixelLine1*/

static boolean DoSubpixelLine2 ( void *usrdata, int4 *jobnum )
{
  aa_data    *aadata;
  pkRenderer *rend;
  int        nthr, nthr9, xmin, xmax, x, i;
  double     y;
  byte       *aaline;

  aadata = (aa_data*)usrdata;
  rend = aadata->rend;
  nthr = rend->nthr;
  nthr9 = 9*nthr;
  xmin = rend->CPos.xmin-1+jobnum->x;
  xmax = rend->CPos.xmin+rend->CPos.width;
  y = aadata->y;
  aaline = aadata->aaline;
  for ( x = xmin, i = 9*jobnum->x;  x < xmax;  x += nthr, i += nthr9 ) {
    if ( OverThreshold1 ( &aaline[i], &aaline[i+9] ) ) {
      GetSample ( rend, (double)x+(double)(1.0/3.0), y, &aaline[i+3] );
      GetSample ( rend, (double)x+(double)(2.0/3.0), y, &aaline[i+6] );
    }
    else {
      Interpol2_3 ( &aaline[i], &aaline[i+9], &aaline[i+3] );
      Interpol2_3 ( &aaline[i+9], &aaline[i], &aaline[i+6] );
    }
  }
  return true;
} /*DoSubpixelLine2*/

void RenderSubpixelLine ( pkRenderer *rend, double y, byte *aaline )
{
  short int nthr;
  aa_data   aadata;
  int4      jobsize;
  boolean   success;

  if ( !rend->RenderingIsOn )
    return;

  nthr = rend->nthr;
  aadata.rend = rend;
  aadata.y = y;
  aadata.aaline = aaline;
  if ( nthr == 1 ) {
    jobsize.x = 0;
    DoSubpixelLine1 ( (void*)&aadata, &jobsize );
    DoSubpixelLine2 ( (void*)&aadata, &jobsize );
  }
  else {
    jobsize.x = nthr;
    pkv_SetPThreadsToWork ( 1, &jobsize, nthr, PKRENDER_STACK, PKRENDER_SCRATCHMEM,
                            (void*)&aadata, DoSubpixelLine1, NULL, NULL, &success );
    pkv_SetPThreadsToWork ( 1, &jobsize, nthr, PKRENDER_STACK, PKRENDER_SCRATCHMEM,
                            (void*)&aadata, DoSubpixelLine2, NULL, NULL, &success );
  }
} /*RenderSubpixelLine*/

static boolean DoSubpixelLine3 ( void *usrdata, int4 *jobnum )
{
  aa_data    *aadata;
  pkRenderer *rend;
  int        nthr, nthr9, xmin, xmax, x, i;
  double     y;
  byte       **aaline;

  aadata = (aa_data*)usrdata;
  rend = aadata->rend;
  nthr = rend->nthr;
  nthr9 = 9*nthr;
  xmin = rend->CPos.xmin-1+jobnum->x;
  xmax = rend->CPos.xmin+rend->CPos.width;
  y = aadata->y;
  aaline = rend->aaline;
  for ( x = xmin, i = 9*jobnum->x;  x <= xmax;  x += nthr, i += nthr9 ) {
    if ( OverThreshold1 ( &aaline[3][i], &aaline[6][i] ) ) {
      GetSample ( rend, (double)x, y+(double)(1.0/3.0), &aaline[4][i] );
      GetSample ( rend, (double)x, y+(double)(2.0/3.0), &aaline[5][i] );
    }
    else {
      Interpol2_3 ( &aaline[3][i], &aaline[6][i], &aaline[4][i] );
      Interpol2_3 ( &aaline[6][i], &aaline[3][i], &aaline[5][i] );
    }
  }
  return true;
} /*DoSubpixelLine3*/

static boolean DoSubpixelLine4 ( void *usrdata, int4 *jobnum )
{
  aa_data    *aadata;
  pkRenderer *rend;
  int        nthr, nthr9, xmin, xmax, x, i;
  double     y, xx;
  byte       **aaline;

  aadata = (aa_data*)usrdata;
  rend = aadata->rend;
  nthr = rend->nthr;
  nthr9 = 9*nthr;
  xmin = rend->CPos.xmin-1+jobnum->x;
  xmax = rend->CPos.xmin+rend->CPos.width;
  y = aadata->y;
  aaline = rend->aaline;
  for ( x = xmin, i = 9*jobnum->x;  x < xmax;  x += nthr, i += nthr9 ) {
    xx = (double)x;
    if ( OverThreshold2 ( &aaline[3][i], &aaline[3][i+9],
                          &aaline[6][i], &aaline[6][i+9] ) ) {
      GetSample ( rend, xx+(double)(1.0/3.0), y+(double)(1.0/3.0), &aaline[4][i+3] );
      GetSample ( rend, xx+(double)(2.0/3.0), y+(double)(1.0/3.0), &aaline[4][i+6] );
      GetSample ( rend, xx+(double)(1.0/3.0), y+(double)(2.0/3.0), &aaline[5][i+3] );
      GetSample ( rend, xx+(double)(2.0/3.0), y+(double)(2.0/3.0), &aaline[5][i+6] );
    }
    else {
      Interpol2_3 ( &aaline[4][i], &aaline[4][i+9], &aaline[4][i+3] );
      Interpol2_3 ( &aaline[4][i+0], &aaline[4][i], &aaline[4][i+6] );
      Interpol2_3 ( &aaline[5][i], &aaline[5][i+9], &aaline[5][i+3] );
      Interpol2_3 ( &aaline[5][i+0], &aaline[5][i], &aaline[5][i+6] );
    }
  }
  return true;
} /*DoSubpixelLine4*/

void GetSupersamples ( pkRenderer *rend, double y )
{
  short int nthr;
  aa_data   aadata;
  int4      jobsize;
  boolean   success;

  if ( !rend->RenderingIsOn )
    return;

  nthr = rend->nthr;
  aadata.rend = rend;
  aadata.y = y;
  if ( nthr == 1 ) {
    jobsize.x = 0;
    DoSubpixelLine3 ( (void*)&aadata, &jobsize );
    DoSubpixelLine4 ( (void*)&aadata, &jobsize );
  }
  else {
    jobsize.x = nthr;
    pkv_SetPThreadsToWork ( 1, &jobsize, nthr, PKRENDER_STACK, PKRENDER_SCRATCHMEM,
                            (void*)&aadata, DoSubpixelLine3, NULL, NULL, &success );
    pkv_SetPThreadsToWork ( 1, &jobsize, nthr, PKRENDER_STACK, PKRENDER_SCRATCHMEM,
                            (void*)&aadata, DoSubpixelLine4, NULL, NULL, &success );
  }
} /*GetSupersamples*/

void InitRenderingAA ( pkRenderer *rend )
{
  int i, lgt;

  lgt = 3*(3*rend->CPos.width+4);
  for ( i = 0; i < 7; i++ )
    rend->aaline[i] = &rend->aabuf[i*lgt];
  RenderSubpixelLine ( rend, (double)(rend->y-1), rend->aaline[3] );
  RenderSubpixelLine ( rend, (double)rend->y, rend->aaline[6] );
  GetSupersamples ( rend, (double)(rend->y-1) );
} /*InitRenderingAA*/

static boolean FilterSupersamples ( void *usrdata, int4 *jobnum )
{
      /* Gaussian filter coefficients for convolution */
  const static unsigned int fc[7] = {4,22,60,84,60,22,4};
  pkRenderer   *rend;
  unsigned int b[8][3];
  int          xmin, xmax, nthr, nthr9, i, j, k, l, x, y;

  rend = (pkRenderer*)usrdata;
  nthr = rend->nthr;
  nthr9 = 9*nthr;
  xmin = rend->CPos.xmin-1+jobnum->x;
  xmax = rend->CPos.xmin+rend->CPos.width;
  y = rend->y;
  for ( x = xmin, i = 9*jobnum->x;  x < xmax;  x += nthr, i += nthr9 ) {
    memset ( b, 0, 24*sizeof(unsigned int) );
    for ( j = 0; j < 7; j++ )
      for ( k = 0; k < 7; k++ )
        for ( l = 0; l < 3; l++ )
          b[j][l] += fc[k]*rend->aaline[j][i+3*k+l];
    for ( j = 0; j < 7; j++ )
      for ( l = 0; l < 3; l++ )
        b[7][l] += fc[j]*b[j][l];
    rend->SetPixel ( rend->private_data, x, y,
                     b[7][0] >> 16, b[7][1] >> 16, b[7][2] >> 16 );
  }
  return true;
} /*FilterSupersamples*/

int RenderLineAA ( pkRenderer *rend )
{
  byte    *a;
  int4    jobsize;
  boolean success;

  if ( rend->RenderingIsOn ) {
        /* move previously computed samples up */
    a = rend->aaline[0];  rend->aaline[0] = rend->aaline[3];
    rend->aaline[3] = rend->aaline[6];  rend->aaline[6] = rend->aaline[2];
    rend->aaline[2] = rend->aaline[5];  rend->aaline[5] = rend->aaline[1];
    rend->aaline[1] = rend->aaline[4];  rend->aaline[4] = a;
        /* get the samples for the current and next line */
    RenderSubpixelLine ( rend, (double)(rend->y+1), rend->aaline[6] );
    GetSupersamples ( rend, (double)rend->y );
        /* filtering and writing the pixels to the image */
    if ( rend->nthr == 1 ) {
      jobsize.x = 0;
      FilterSupersamples ( (void*)rend, &jobsize );
    }
    else {
      jobsize.x = rend->nthr;
      pkv_SetPThreadsToWork ( 1, &jobsize, rend->nthr,
                              PKRENDER_STACK, PKRENDER_SCRATCHMEM,
                              (void*)rend, FilterSupersamples, NULL, NULL, &success );
    }
        /* get to the next line or finish */
    rend->y ++;
    if ( rend->y >= rend->CPos.ymin+rend->CPos.height ) {
      rend->ticks2 = pkv_Toc ( &rend->tic );
      rend->RenderingIsOn = false;
    }
  }
  return rend->y;
} /*RenderLineAA*/

