
/* //////////////////////////////////////////////////// */
/* This file is a part of the BSTools procedure package */
/* written by Przemyslaw Kiciak                         */
/* //////////////////////////////////////////////////// */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <malloc.h>
#include <string.h>
#include <unistd.h>
#include <sys/times.h>

#include <X11/Xlib.h>
#include <X11/Xutil.h>

#include "pkvaria.h"
#include "pknum.h"
#include "pkgeom.h"
#include "camerad.h"
#include "multibs.h"
#include "eg1holed.h"
#include "eg2holed.h"
#include "xgedit.h"

#include "render.h"
#include "edpolicz.h"
#include "edwidgets.h"
#include "splhole.h"

static clock_t tick;
static long    ticksps;

/* ///////////////////////////////////////////////////////////////////////// */
void StartRendering ( void )
{
  struct tms tt;

  ticksps = sysconf ( _SC_CLK_TCK );
  tick    = times ( &tt );
        /* send the patches and camera to the renderer */
  RendReset ();
  SendPatchesToRenderer ();
  RendEnterCamerad ( &swind.CPos[3], swind.cwin[3] );
  RendEnterLightsd ( R_NLIGHTS, light_dir, light_int );
  if ( RendBegin ( swShadows, swAntialias ) ) {
        /* prepare the background */
    swind_picture = false;
    xge_SetClipping ( swind.cwin[3] );
    swind.cwin[3]->redraw ( swind.cwin[3], false );
    XGetSubImage ( xgedisplay, xgepixmap, swind.cwin[3]->x, swind.cwin[3]->y,
                   swind.cwin[3]->w, swind.cwin[3]->h, 0xFFFFFFFF, ZPixmap,
                   rendimage, swind.cwin[3]->x, swind.cwin[3]->y  );
        /* launch rendering of the first picture line */
    swind_picture = true;
    xge_PostIdleCommand ( IDLE_COMMAND_RENDER, 0, 0 );
  }
} /*StartRendering*/

void ContRendering ( void )
{
  clock_t    tick1;
  struct tms tt;
/*  int        y; */
  xge_widget *er;

/*  y =*/ RenderLine ();
/*  printf ( "%d\n", y ); */
  tick1 = times ( &tt );
  if ( (tick1-tick) >= ticksps || !RenderingIsOn ) {
        /* show the picture on the screen once per second */
    er = swind.cwin[3];
    xge_SetClipping ( er );
    er->redraw ( er, false );
    xge_RedrawPopups ();
    xge_SetClipping ( er );
    xgeCopyRectOnScreen ( er->w, er->h, er->x, er->y );
    tick = tick1;
  }
  if ( RenderingIsOn )
    xge_PostIdleCommand ( IDLE_COMMAND_RENDER, 0, 0 );
  else {
    BreakRendering ( true );
  }
} /*ContRendering*/

void BreakRendering ( boolean redraw )
{
  if ( RenderingIsOn )
    RendReset ();
  renderbtn0->data0 = renderbtn1->data0 = txtRender;
  if ( redraw ) {
    if ( menu1->data1 == menu01elist ) {
      xge_SetClipping ( renderbtn0 );
      renderbtn0->redraw ( renderbtn0, true );
    }
    else if ( menu1->data1 == menu01flist ) {
      xge_SetClipping ( renderbtn1 );
      renderbtn1->redraw ( renderbtn1, true );
    }
  }
} /*BreakRendering*/

