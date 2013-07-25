
/* //////////////////////////////////////////////////// */
/* This file is a part of the BSTools procedure package */
/* written by Przemyslaw Kiciak                         */
/* //////////////////////////////////////////////////// */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <malloc.h>
#include <string.h>

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


void ResizeWinStatus ( int win )
{
  short w, h;

  xge_SetWindow ( win );
  if ( win == win0 ) {
    if ( status_0_sw ) {
      status0->data0 = statustext0;
      w = xge_current_width;
      h = xge_current_height-36;
    }
    else {
      status0->data0 = txtNULL;
      w = 110;
      h = xge_current_height-20;
    }
    swind.fww.er->msgproc ( swind.fww.er, xgemsg_RESIZE, 0,
                            (short)(xge_current_width-110), h );
    menu0->msgproc ( menu0, xgemsg_RESIZE, 0, xge_current_width, 20 );
    scroll_constr_sw.er->msgproc ( scroll_constr_sw.er, xgemsg_RESIZE, 0,
                                   110, xge_current_height-286 );
    ConfigureConstraintWidgets ( false );
    menu1->msgproc ( menu1, xgemsg_RESIZE, 0,
                     110, (short)(xge_current_height-36) );
    status0->w = (short)(w-110);
    menu2->y = (short)(xge_current_height-16);
    menu2->msgproc ( menu2, xgemsg_RESIZE, 0, w, 16 );
  }
  else {
    if ( status_1_sw ) {
      status1->data0 = statustext1;
      w = xge_current_width;
      h = (short)(xge_current_height-36);
      knotwind->y = (short)(xge_current_height-16-knotwind->h);
    }
    else {
      status1->data0 = txtNULL;
      w = 110;
      h = (short)(xge_current_height-20);
      knotwind->y = (short)(xge_current_height-knotwind->h);
    }
    domwind.er->msgproc ( domwind.er, xgemsg_RESIZE, 0,
                          (short)(xge_current_width-110),
                          (short)(h-knotwind->h-4) );
    knotwind->msgproc ( knotwind, xgemsg_RESIZE, 0,
                      (short)(xge_current_width-110), knotwind->h );
    menu3->msgproc ( menu3, xgemsg_RESIZE, 0, xge_current_width, 20 );
    menu4->msgproc ( menu4, xgemsg_RESIZE, 0,
                     110, (short)(xge_current_height-36) );
    status1->w = (short)(w-110);
    menu5->y = (short)(xge_current_height-16);
    menu5->msgproc ( menu5, xgemsg_RESIZE, 0, w, 16 );
  }
} /*ResizeWinStatus*/

void StatusLineOnOff ( int win )
{
  ResizeWinStatus ( win );
  xge_Redraw ();
} /*StatusLineOnOff*/

void SetStatusLineText ( int win, const char *text, boolean onscreen )
{
  char       *st;
  xge_widget *er;
  int        l;
  boolean    show;

        /* select the window */
  if ( win == win0 ) {
    er = status0;
    st = statustext0;
    show = status_0_sw;
  }
  else {
    er = status1;
    st = statustext1;
    show = status_1_sw;
  }
        /* copy the text */
  l = strlen ( text );
  l = min ( l, STATUS_LINE_LENGTH );
  strncpy ( st, text, l );
  st[l] = 0;
        /* redisplay */
  if ( onscreen && show ) {
    xge_SetWindow ( win );
    xge_SetClipping ( er );
    er->redraw ( er, true );
  }
} /*SetStatusLineText*/

/* ////////////////////////////////////////////////////////////////////////// */
void Notify2DTrans ( int win, xge_2Dwind *_2Dwin )
{
  char str[64];

  switch ( _2Dwin->current_tool ) {
case xge_2DWIN_MOVING_TOOL:
    sprintf ( str, "%f,%f", _2Dwin->trans_params.x, _2Dwin->trans_params.y );
    break;
case xge_2DWIN_SCALING_TOOL:
    sprintf ( str, "%f,%f", _2Dwin->trans_params.x, _2Dwin->trans_params.y );
    break;
case xge_2DWIN_ROTATING_TOOL:
    sprintf ( str, "%f", _2Dwin->trans_params.x );
    break;
default:
    str[0] = 0;
    break;
  }
  SetStatusLineText ( win, str, true );
} /*Notify2DTrans*/

void Notify2DTransChange ( int win, xge_2Dwind *_2Dwin )
{
  char str[64];

  switch ( _2Dwin->current_tool ) {
case xge_2DWIN_MOVING_TOOL:
    str[0] = 0;
    break;
case xge_2DWIN_SCALING_TOOL:
    sprintf ( str, "%f,%f", _2Dwin->scaling_centre.x, _2Dwin->scaling_centre.y );
    break;
case xge_2DWIN_ROTATING_TOOL:
    sprintf ( str, "%f,%f", _2Dwin->rotating_centre.x, _2Dwin->rotating_centre.y );
    break;
default:
    str[0] = 0;
    break;
  }
  SetStatusLineText ( win, str, true );
} /*Notify2DTransChange*/

void Notify3DTrans ( int win, xge_3Dwind *_3Dwin )
{
  char str[64];

    str[0] = 0;
  switch ( _3Dwin->current_tool ) {
case xge_3DWIN_MOVING_TOOL:
    sprintf ( str, "%f,%f,%f", _3Dwin->trans_params.x, _3Dwin->trans_params.y,
              _3Dwin->trans_params.z );
    break;
case xge_3DWIN_SCALING_TOOL:
    sprintf ( str, "%f,%f,%f", _3Dwin->trans_params.x, _3Dwin->trans_params.y,
              _3Dwin->trans_params.z );
    break;
case xge_3DWIN_ROTATING_TOOL:
    sprintf ( str, "%f", _3Dwin->trans_params.x );
    break;
default:
    str[0] = 0;
    break;
  }
  SetStatusLineText ( win, str, true );
} /*Notify3DTrans*/

void Notify3DTransChange ( int win, xge_3Dwind *_3Dwin )
{
  char str[64];

  switch ( _3Dwin->current_tool ) {
case xge_3DWIN_MOVING_TOOL:
    str[0] = 0;
    break;
case xge_3DWIN_SCALING_TOOL:
    sprintf ( str, "%f,%f,%f", _3Dwin->scaling_centre.x,
              _3Dwin->scaling_centre.y, _3Dwin->scaling_centre.z );
    break;
case xge_3DWIN_ROTATING_TOOL:
    sprintf ( str, "%f,%f,%f", _3Dwin->rotating_centre.x,
              _3Dwin->rotating_centre.y, _3Dwin->rotating_centre.z );
    break;
default:
    str[0] = 0;
    break;
  }
  SetStatusLineText ( win, str, true );
} /*Notify3DTransChange*/

void NotifyDoubleNumber ( int win, double x, boolean onscreen )
{
  char str[40];

  sprintf ( str, "%4.2f", x );
  SetStatusLineText ( win, str, onscreen );
} /*NotifyDoubleNumber*/

