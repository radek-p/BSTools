
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2007, 2012                            */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <malloc.h>
#include <string.h>

#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <X11/cursorfont.h>

#include <GL/gl.h>
#include <GL/glx.h>

#include "pkgeom.h"
#include "multibs.h"

#include "xgedit.h"
#include "xgeprivate.h"

/* ///////////////////////////////////////////////////////////////////////// */
boolean _xge_background_msg ( xge_widget *er,
                             int msg, int key, short x, short y )
{
  switch ( msg ) {
case xgemsg_ENTERING:
    xge_SetCurrentWindowCursor ( xgeCURSOR_DEFAULT );
    return true;
default:
    break;
  }
  return false;
} /*_xge_background_msg*/

/* ///////////////////////////////////////////////////////////////////////// */
boolean xge_IntersectXRectangles ( XRectangle *r1, XRectangle *r2 )
{
  short x1, y1, x2, y2;

  x1 = r1->x+r1->width;  y1 = r1->y+r1->height;
  x2 = r2->x+r2->width;  y2 = r2->y+r2->height;
  x1 = min ( x1, x2 );   y1 = min ( y1, y2 );
  if ( r2->x > r1->x ) r1->x = r2->x;
  if ( r2->y > r1->y ) r1->y = r2->y;
  x2 = x1-r1->x;   if ( x2 < 0 ) x2 = 0;
  y2 = y1-r1->y;   if ( y2 < 0 ) y2 = 0;
  r1->width = x2;  r1->height = y2;
  return x2 > 0 && y2 > 0;
} /*xge_IntersectXRectangles*/

boolean xge_SetClipping ( xge_widget *er )
{
  XRectangle rect1, rect2;

  if ( er->window_num != -1 && er->window_num != xge_current_win )
    xge_SetWindow ( er->window_num );
  rect1.width = er->w;  rect1.height = er->h;
  rect1.x = er->x;  rect1.y = er->y;
  for ( er = er->up; er; er = er->up ) {
    rect2.width = er->w;  rect2.height = er->h;
    rect2.x = er->x;  rect2.y = er->y;
    xge_IntersectXRectangles ( &rect1, &rect2 );
  }
  XSetClipRectangles ( xgedisplay, xgegc, 0, 0, &rect1, 1, Unsorted );
  return rect1.width > 0 && rect1.height > 0;
} /*xge_SetClipping*/

void xge_ResetClipping ( void )
{
  XSetClipMask ( xgedisplay, xgegc, None );
} /*xge_ResetClipping*/

void xge_RedrawPopups ( void )
{
  xge_widget *rp;

  for ( rp = xge_windesc[xge_current_win].popup0;  rp;  rp = rp->next ) {
    xge_SetClipping ( rp );
    rp->redraw ( rp, false );
  }
} /*xge_RedrawPopups*/

void xge_Redraw ( void )
{
  xge_widget *rp;

  xge_ResetClipping ();
  xgeSetForeground ( xgec_Black );
  xgeFillRectangle ( xge_current_width, xge_current_height, 0, 0 );
  for ( rp = xge_windesc[xge_current_win].win_edr0;  rp;  rp = rp->next ) {
    xge_SetClipping ( rp );
    rp->redraw ( rp, false );
  }
  xge_RedrawPopups ();
  xge_ResetClipping ();
  xgeCopyRectOnScreen ( xge_current_width, xge_current_height, 0, 0 );
} /*xge_Redraw*/

void xge_RedrawAll ( void )
{
  int win, i;

  win = xge_current_win;
  for ( i = 0; i < xge_winnum; i++ ) {
    xge_SetWindow ( i );
    xge_Redraw ();
  }
  xge_SetWindow ( win );
} /*xge_RedrawAll*/

/* ///////////////////////////////////////////////////////////////////////// */
void xge_BoundPoint ( xge_widget *er, short *x, short *y )
{
  int a;

  a = er->x + 2;
  *x = (short)(*x > a ? *x : a);
  a = er->x + er->w - 3;
  *x = (short)(*x < a ? *x : a);
  a = er->y + 2;
  *y = (short)(*y > a ? *y : a);
  a = er->y + er->h - 3;
  *y = (short)(*y < a ? *y : a);
} /*xge_BoundPoint*/

/* ///////////////////////////////////////////////////////////////////////// */
boolean xge_PointInRect ( xge_widget *edr, short x, short y )
{
  return (boolean)((x >= edr->x) && (x < edr->x+edr->w) &&
                   (y >= edr->y) && (y < edr->y+edr->h));
} /*xge_PointInRect*/

boolean xge_RectanglesIntersect ( short wa, short ha, short xa, short ya,
                                  short wb, short hb, short xb, short yb )
{
  if ( xa >= xb+wb || xb >= xa+wa || ya >= yb+hb || yb >= ya+ha )
    return false;
  else
    return true;
} /*xge_RectanglesIntersect*/

/* ///////////////////////////////////////////////////////////////////////// */
boolean xge_CallMsgProc ( xge_widget **last, xge_widget *er,
                          int msg, int key, short x, short y )
{

  switch ( msg ) {
case xgemsg_NULL:
case xgemsg_KEY:
case xgemsg_SPECIAL_KEY:
case xgemsg_MMOVE:
case xgemsg_MCLICK:
case xgemsg_OTHEREVENT:
    if ( er != *last ) {
      if ( *last )
        (*last)->msgproc ( *last, xgemsg_EXITING, 0, x, y );
      if ( er ) {
        er->msgproc ( er, xgemsg_ENTERING, 0, x, y );
        *last = er;
      }
      else
        *last = NULL;
    }
    break;

default:
    break;
  }
  if ( er )
    return er->msgproc ( er, msg, key, x, y );
  else
    return false;
} /*xge_CallMsgProc*/

void xge_dispatch_message ( unsigned int msg, unsigned int key, short x, short y )
{
  xge_widget *rp;
  int        fsp;

  if ( msg == xgemsg_NULL )
    return;

  fsp = (int)xge_windesc[xge_current_win].fsp;
  if ( fsp > 0 ) {
/* printf ( "fsp = %d,", fsp ); */
    rp = xge_windesc[xge_current_win].focusstack[fsp-1];
/* printf ( "  rp = %x,", (unsigned int)rp ); */
    if ( (xge_notinfocus = (boolean)(rp->window_num != xge_current_win)) ) {
/* printf ( "  rp->window_num = %d", rp->window_num ); */
      xge_SetWindow ( rp->window_num );
      x = (short)-x;  /* this is to prevent a wrong interpretation of the input */
      y = (short)-y;
    }
/* printf ( "\n" ); */
    if ( xge_CallMsgProc ( &xge_lastwin, rp, msg, key, x, y ) ) {
      if ( rp->window_num != xge_current_win )
        xge_SetWindow ( xge_current_win );
      if ( !xge_windesc[xge_current_win].fsp ) {
        msg = xgemsg_MMOVE;
        key = xge_mouse_buttons;
        x = (short)xge_mouse_x;
        y = (short)xge_mouse_y;
        goto nonfocus;
      }
    }
    else
      xge_callback ( NULL, msg, key, x, y );
  }
  else {
nonfocus:
    for ( rp = xge_windesc[xge_current_win].popup1;  rp;  rp = rp->prev )
      if ( xge_PointInRect ( rp, x, y ) ) {
        if ( xge_CallMsgProc ( &xge_lastwin, rp, msg, key, x, y ) )
          goto store_pos;
      }
    for ( rp = xge_windesc[xge_current_win].win_edr1;  rp;  rp = rp->prev )
      if ( xge_PointInRect ( rp, x, y ) ) {
        if ( xge_CallMsgProc ( &xge_lastwin, rp, msg, key, x, y ) )
          goto store_pos;
      }
    xge_callback ( NULL, msg, key, x, y );
  }
store_pos:
  if ( msg == xgemsg_MMOVE || msg == xgemsg_MCLICK ) {
    xge_prevx = (short)xge_mouse_x;
    xge_prevy = (short)xge_mouse_y;
  }
} /*xge_dispatch_message*/

void xge_PostIdleCommand ( unsigned int key, short x, short y )
{
  XEvent event;

  event.type = ClientMessage;
  event.xclient.window = xgewindow;
  event.xclient.message_type = xgemsg_IDLE_COMMAND;
  event.xclient.format = 32;
  event.xclient.data.l[0] = (int)key;
  event.xclient.data.l[1] = (int)x;
  event.xclient.data.l[2] = (int)y;
  XSendEvent ( xgedisplay, xgewindow, false, 0, &event );
} /*xge_PostIdleCommand*/

/* ///////////////////////////////////////////////////////////////////////// */
int xge_NewWindow ( char *title )
{
  XSetWindowAttributes swa;
  Window win, root;
  Pixmap pix;
  int    x, y;
  unsigned int border_width, depth;

  if ( xge_winnum < xge_MAX_WINDOWS ) {
    swa.colormap = xgecolormap;
    xge_background = swa.background_pixel = swa.border_pixel =
      BlackPixel ( xgedisplay, xgescreen );
    swa.event_mask = ButtonPressMask | ExposureMask | ButtonReleaseMask |
                     PointerMotionMask | KeyPressMask | StructureNotifyMask ;
    xge_foreground = WhitePixel ( xgedisplay, xgescreen );
    win = XCreateWindow ( xgedisplay, xgeroot, 0, 0, xge_WIDTH, xge_HEIGHT, 0,
                          xge_nplanes, InputOutput, xgevisual,
                          CWBackPixel | CWBorderPixel | CWColormap | CWEventMask,
                          &swa );

    XSetStandardProperties ( xgedisplay, win, title, title,
                             None, _xge_argv, _xge_argc, &xgehints );
    XSelectInput ( xgedisplay, win,
                   ButtonPressMask | ExposureMask | ButtonReleaseMask |
                   PointerMotionMask | KeyPressMask | StructureNotifyMask );
    XMapWindow ( xgedisplay, win );
    pix = XCreatePixmap ( xgedisplay, win,
                          xge_MAX_WIDTH, xge_MAX_HEIGHT, xge_nplanes );

    xge_windesc[xge_winnum].thewindow = win;
    xge_windesc[xge_winnum].thepixmap = pix;
    xge_windesc[xge_winnum].win_rect_num = 0;
    xge_windesc[xge_winnum].win_edr0 = xge_windesc[xge_winnum].win_edr1 =
    xge_windesc[xge_winnum].popup0 = xge_windesc[xge_winnum].popup1 = NULL;
    xge_windesc[xge_winnum].fsp = 0;
    XGetGeometry ( xgedisplay, win, &root, &x, &y,
                   &xge_current_width,
                   &xge_current_height,
                   &border_width, &depth );
    xge_windesc[xge_winnum].thewinrect.width  = xge_current_width;
    xge_windesc[xge_winnum].thewinrect.height = xge_current_height;
    win = xge_winnum;
    xge_winnum++;
    return win;
  }
  else return -1;
} /*xge_NewWindow*/

boolean xge_SetWindow ( int win )
{
  if ( win >= 0 && win < xge_winnum ) {
    if ( win != xge_current_win ) {
      if ( xge_current_win >= 0 && xge_current_win < xge_winnum ) {
        xge_windesc[xge_current_win].thewindow         = xgewindow;
        xge_windesc[xge_current_win].thepixmap         = xgepixmap;
        xge_windesc[xge_current_win].thewinrect.width  = xge_current_width;
        xge_windesc[xge_current_win].thewinrect.height = xge_current_height;
      }
      xgewindow          = xge_windesc[win].thewindow;
      xgepixmap          = xge_windesc[win].thepixmap;
      xge_current_width  = xge_windesc[win].thewinrect.width;
      xge_current_height = xge_windesc[win].thewinrect.height;
      xge_current_win    = win;
    }
    return true;
  }
  else
    return false;
} /*xge_SetWindow*/

int xge_CurrentWindow ( void )
{
  return xge_current_win;
} /*xge_CurrentWindow*/

static int xge_FindWindowNum ( Window thewin )
{
  int i;

  for ( i = 0; i < xge_winnum; i++ )
    if ( thewin == xge_windesc[i].thewindow )
      return i;
  return -1;
} /*xge_FindWindowNum*/

void  xge_SetWinEdRect ( xge_widget *edr )
{
  xge_windesc[xge_current_win].win_edr1 = edr;
  if ( edr ) {  /* find the list end */
    while ( edr->prev )
      edr = edr->prev;
    xge_windesc[xge_current_win].win_edr0 = edr;
  }
  else
    xge_windesc[xge_current_win].win_edr0 = NULL;
} /*xge_SetWinEdRect*/

/* ///////////////////////////////////////////////////////////////////////// */
static boolean xge_ProcessSpecialKey ( char key, unsigned int keysym )
{
  int        dx, dy, fsp;
  xge_widget *rp;

  if ( (fsp = (int)xge_windesc[xge_current_win].fsp) > 0 ) {
    rp = xge_windesc[xge_current_win].focusstack[fsp-1];
    if ( key ) {
      if ( xge_CallMsgProc ( &xge_lastwin, rp, xgemsg_KEY, key, 0, 0 ) )
        return true;
    }
    else {
      if ( xge_CallMsgProc ( &xge_lastwin, rp,
                             xgemsg_SPECIAL_KEY, keysym, 0, 0 ) )
        return true;
    }
  }
  if ( !key ) {
    switch ( keysym ) {
case 0xFF51: /* key left */
      dx = -1;  dy = 0;
      goto cont;
case 0xFF52: /* key up */
      dx = 0;  dy = -1;
      goto cont; 
case 0xFF53: /* key right */
      dx = +1;  dy = 0;
      goto cont;
case 0xFF54: /* key down */
      dx = 0;  dy = +1;
cont:
      XWarpPointer ( xgedisplay, None, None, 0, 0, 0, 0, dx, dy );
      xge_dispatch_message ( xgemsg_MMOVE, xge_mouse_buttons, 
                             (short)(xge_mouse_x+dx), (short)(xge_mouse_y+dy) );
      return true;

default:
      xge_dispatch_message ( xgemsg_SPECIAL_KEY, keysym,
                             xge_mouse_x, xge_mouse_y );
      break;
    }
  }
  return false;
} /*xge_ProcessSpecialKey*/

void xge_get_message ( unsigned int *msg, unsigned int *key, short *x, short *y )
{
  Window         r_w, c_w;
  int            win;
  int            x_r, y_r;
  short          w, h;
  unsigned int   xbutton_mask;
  unsigned int   newmask = 0;

  XNextEvent ( xgedisplay, &xgeevent );
  win = xge_FindWindowNum ( xgeevent.xany.window );
  if ( win != xge_current_win )
    xge_SetWindow ( win );

  switch ( xgeevent.type )
  {
case Expose:
    if ( xgeevent.xexpose.count == 0 ) {
      XSetClipRectangles ( xgedisplay, xgegc, 0, 0,
                           &xge_windesc[xge_current_win].thewinrect, 1, Unsorted );
      xgeCopyRectOnScreen ( xge_windesc[xge_current_win].thewinrect.width,
                            xge_windesc[xge_current_win].thewinrect.height, 0, 0 );
    }
    *msg = xgemsg_NULL;
    break;

case ButtonPress:
    xge_mouse_x = xgeevent.xbutton.x;
    xge_mouse_y = xgeevent.xbutton.y;
    switch ( xgeevent.xbutton.button )
    {
  case Button1:
      newmask = xge_mouse_buttons | xgemouse_LBUTTON_DOWN;
      *key = newmask | xgemouse_LBUTTON_CHANGE;
      break;
  case Button2:
      newmask = xge_mouse_buttons | xgemouse_MBUTTON_DOWN;
      *key = newmask | xgemouse_MBUTTON_CHANGE;
      break;
  case Button3:
      newmask = xge_mouse_buttons | xgemouse_RBUTTON_DOWN;
      *key = newmask | xgemouse_RBUTTON_CHANGE;
      break;
  case Button4:
      newmask = xge_mouse_buttons | xgemouse_WHEELFW_DOWN;
      *key = newmask | xgemouse_WHEELFW_CHANGE;
      break;
  case Button5:
      newmask = xge_mouse_buttons | xgemouse_WHEELBK_DOWN;
      *key = newmask | xgemouse_WHEELBK_CHANGE;
      break;
    }
    xge_mouse_buttons = newmask;
    *msg = xgemsg_MCLICK;
    *x = (short)xge_mouse_x;
    *y = (short)xge_mouse_y;
    break;

case ButtonRelease:
    xge_mouse_x = xgeevent.xbutton.x;
    xge_mouse_y = xgeevent.xbutton.y;
    switch ( xgeevent.xbutton.button )
    {
  case Button1: 
      newmask = xge_mouse_buttons & (~xgemouse_LBUTTON_DOWN);
      *key = newmask | xgemouse_LBUTTON_CHANGE;
      break;
  case Button2:
      newmask = xge_mouse_buttons & (~xgemouse_MBUTTON_DOWN);
      *key = newmask | xgemouse_MBUTTON_CHANGE;
      break;
  case Button3:
      newmask = xge_mouse_buttons & (~xgemouse_RBUTTON_DOWN);
      *key = newmask | xgemouse_RBUTTON_CHANGE;
      break;
  case Button4:
      newmask = xge_mouse_buttons & (~xgemouse_WHEELFW_DOWN);
      *key = newmask | xgemouse_WHEELFW_CHANGE;
      break;
  case Button5:
      newmask = xge_mouse_buttons & (~xgemouse_WHEELBK_DOWN);
      *key = newmask | xgemouse_WHEELBK_CHANGE;
      break;
    }
    xge_mouse_buttons = newmask;
    *msg = xgemsg_MCLICK;
    *x = (short)xge_mouse_x;
    *y = (short)xge_mouse_y;
    break;

case MotionNotify:
        /* Getting the pointer position by calling XQueryPointer  */
        /* instead of taking it from the XMotion is supposed to   */
        /* result in the position at the moment of processing the */
        /* event, not at the event generation. In this way        */
        /* a number of motion events may be combined into one,    */
        /* which results in a faster program execution.           */
    XQueryPointer ( xgedisplay, xgewindow,
                    &r_w, &c_w, &x_r, &y_r, &xge_mouse_x, &xge_mouse_y, &xbutton_mask );
    if ( xge_mouse_x != xge_prevx || xge_mouse_y != xge_prevy ) {
      *x = (short)xge_mouse_x;
      *y = (short)xge_mouse_y;
      *key = xge_mouse_buttons;
      *msg = xgemsg_MMOVE;
    }
    else *msg = xgemsg_NULL;
    break;

case ConfigureNotify:
    w = (short)max ( xgeevent.xconfigure.width, xge_WIDTH );
    w = (short)min ( w, xge_MAX_WIDTH );
    h = (short)max ( xgeevent.xconfigure.height, xge_HEIGHT );
    h = (short)min ( h, xge_MAX_HEIGHT );
    if ( (w != xge_windesc[win].thewinrect.width ||
          h != xge_windesc[win].thewinrect.height) ) {
      if ( win == xge_errmsg_win && xge_errmsg_edr->state == xgestate_MESSAGE ) {
        xge_RemoveErrorMessage ();
        xge_current_width  = w;
        xge_current_height = h;
        xge_callback ( NULL, xgemsg_RESIZE, 1, w, h );
        if ( xge_errmsg_msgtext )
          _xge_DisplayErrorMessage ( xge_errmsg_msgtext, xge_msgbkcolour,
                                     xge_usr_msg_key );
        else if ( xge_info_msgtext )
          xge_DisplayInfoMessage ( xge_info_msgtext, xge_usr_msg_key );
      }
      else {
        xge_current_width = w;
        xge_current_height = h;
        xge_callback ( NULL, xgemsg_RESIZE, 1, w, h );
      }
      xge_windesc[win].thewinrect.width  = w;
      xge_windesc[win].thewinrect.height = h;
    }
    *msg = xgemsg_NULL;
    break;

case KeyPress:
    *x = (short)(xge_mouse_x = xgeevent.xkey.x);
    *y = (short)(xge_mouse_y = xgeevent.xkey.y);
    *key = 0;
    XLookupString ( &xgeevent.xkey, (char*)key, 1, &xgekeysym, NULL );
    if ( xge_ProcessSpecialKey ( *((char*)key), xgekeysym ) )
      *msg = xgemsg_NULL;
    else
      *msg = xgemsg_KEY;
    break;

case ClientMessage:
    xge_callback ( NULL, xgeevent.xclient.message_type,
                   xgeevent.xclient.data.l[0],
                   (short)xgeevent.xclient.data.l[1],
                   (short)xgeevent.xclient.data.l[2] );
    *msg = xgemsg_NULL;
    break;

default:
    xge_callback ( NULL, xgemsg_OTHEREVENT, 0, 0, 0 );
    *msg = xgemsg_NULL;
    break;
  }
} /*xge_get_message*/

void xge_GetWindowSize ( void )
{
  Window root;
  int x, y;
  unsigned int w, h, border_width, depth;

  XGetGeometry ( xgedisplay, xgewindow, &root, &x, &y,
                 &w, &h, &border_width, &depth );
  xge_current_width = (unsigned short)min ( w, xge_MAX_WIDTH );
  xge_current_height = (unsigned short)min ( h, xge_MAX_HEIGHT );
} /*xge_GetWindowSize*/

void _xge_FindAspect ( void )
{
#ifdef XGE_AUTO_ASPECT
  float dw, dh, dwmm, dhmm;
#endif

#ifdef XGE_AUTO_ASPECT
  /* this may not work properly if the virtual screen has a different size */
  /* than the physical screen; this is to be fixed when I find out, how to */
  /* obtain a necessary information */
  dw = (float)DisplayWidth ( xgedisplay, xgescreen );
  dh = (float)DisplayHeight ( xgedisplay, xgescreen );
  dwmm = (float)DisplayWidthMM ( xgedisplay, xgescreen );
  dhmm = (float)DisplayHeightMM ( xgedisplay, xgescreen );
  xge_aspect = (dwmm/dw)/(dhmm/dh);
#else
  /* if the above does not work properly, comment out the #definition */
  /* XGE_AUTO_ASPECT in xgeprivate.h */
  xge_aspect = XGE_DEFAULT_ASPECT;
#endif
} /*_xge_FindAspect*/

void xge_MessageLoop ( void )
{
  unsigned int msg, key;
  short        x, y;

  while ( !xge_done ) {
    xge_get_message ( &msg, &key, &x, &y );
    xge_dispatch_message ( msg, key, x, y );
  }
} /*xge_MessageLoop*/

