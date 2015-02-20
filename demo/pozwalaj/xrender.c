
/* //////////////////////////////////////////////////// */
/* This file is a part of the BSTools procedure package */
/* written by Przemyslaw Kiciak                         */
/* //////////////////////////////////////////////////// */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <X11/Xlib.h>
#include <X11/Xutil.h>

#include "pkvaria.h"
#include "pkvthreads.h"
#include "pknum.h"
#include "pkgeom.h"
#include "multibs.h"
#include "bsmesh.h"
#include "camera.h"
#include "raybez.h"
#include "xgedit.h"
#include "pkrender.h"

extern int    rendering_npthreads;

pkRenderer    rend;
XImage        *rendimage = NULL;
static char   *imagedata = NULL;

boolean       swGaussian_c = false, swMean_c = false, swLambert_c = false,
              swReflection_c = false, swHighlight_c = false, swSections_c = false;
              boolean swGaussian_d = false, swMean_d = false, swLambert_d = false,
              swReflection_d = false, swHighlight_d = false, swSections_d = false;
boolean       swAntialias = true, swShadows = false;
double        render_cfrange[2] = {0.0,1.0};
double        render_dfsf = 0.5;

static void SetPixel ( void *private_data, short x, short y,
                       byte r, byte g, byte b )
{
  XPutPixel ( rendimage, x, y, xge_PixelColour ( r, g, b ) );
} /*SetPixel*/

boolean InitXRenderer ( void )
{
  int nplanes;

  nplanes = XDisplayPlanes ( xgedisplay, xgescreen );
  rendimage = XCreateImage ( xgedisplay, xgevisual, nplanes, ZPixmap, 0,
                             NULL, xge_MAX_WIDTH, xge_MAX_HEIGHT, 8, 0 );
  imagedata = (char*)malloc ( xge_MAX_HEIGHT*rendimage->bytes_per_line );
  if ( rendimage && imagedata )
    rendimage->data = imagedata;
  else
    return false;
  return RendInit ( &rend, (void*)rendimage, rendering_npthreads,
                    xge_MAX_WIDTH, xge_MAX_HEIGHT, SetPixel );
} /*InitXRenderer*/

void DestroyXRenderer ( void )
{
  RendDestroy ( &rend );
  XDestroyImage ( rendimage );
  raybez_DisablePThreads ();
} /*DestroyXRenderer*/

void SetRendererSwitches ( void )
{
  rend.swShadows = swShadows;
  rend.swAntialias = swAntialias;
  if ( swGaussian_c )         rend.c_shape_func = shapefunc_GAUSSIAN;
  else if ( swMean_c )        rend.c_shape_func = shapefunc_MEAN;
  else if ( swLambert_c )     rend.c_shape_func = shapefunc_LAMBERTIAN;
  else if ( swReflection_c )  rend.c_shape_func = shapefunc_REFLECTION;
  else if ( swHighlight_c )   rend.c_shape_func = shapefunc_HIGHLIGHT;
  else if ( swSections_c )    rend.c_shape_func = shapefunc_CROSSECTIONS;
  else                        rend.c_shape_func = shapefunc_NONE;
  if ( swGaussian_d )         rend.d_shape_func = shapefunc_GAUSSIAN;
  else if ( swMean_d )        rend.d_shape_func = shapefunc_MEAN;
  else if ( swLambert_d )     rend.d_shape_func = shapefunc_LAMBERTIAN;
  else if ( swReflection_d )  rend.d_shape_func = shapefunc_REFLECTION;
  else if ( swHighlight_d )   rend.d_shape_func = shapefunc_HIGHLIGHT;
  else if ( swSections_d )    rend.d_shape_func = shapefunc_CROSSECTIONS;
  else                        rend.d_shape_func = shapefunc_NONE;
  rend.cfrange[0] = render_cfrange[0];
  rend.cfrange[1] = render_cfrange[1];
  rend.dfsf = xge_LogSlidebarValued ( R_MINDFSF, R_MAXDFSF, render_dfsf );
} /*SetRendererSwitches*/

