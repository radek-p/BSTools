
/* //////////////////////////////////////////////////// */
/* This file is a part of the BSTools procedure package */
/* written by Przemyslaw Kiciak.                        */
/* //////////////////////////////////////////////////// */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <malloc.h>
#include <string.h>

#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <X11/cursorfont.h>

#include "pkgeom.h"
#include "multibs.h"

#include "oldxgedit.h"

int  _argc;
char **_argv;

typedef struct {
    Window        thewindow;
    int           win_rect_num;
    ed_rect       *win_edr;
    unsigned int  current_width, current_height;
  } WinDesc;

int     winnum = 0;            /* current number of windows */
int     current_win = -1;      /* window currently active */
int     focus_win = -1;        /* window with the input focus */
WinDesc windesc[MAX_WINDOWS];  /* window descriptors */

Display    *thedisplay;
Window     thewindow;
Pixmap     thepixmap;
int        thescreen;
GC         thegc;
Visual     *thevisual;
XSizeHints thehints;
Cursor     thecursor0, thecursor1;
KeySym     thekeysym;

unsigned int current_width, current_height;  /* window dimensions */
unsigned int mouse_buttons = 0;
int  mouse_x, mouse_y, xx, yy;

char       *p_name;
char       done = 0;
int        prevx = -1, prevy = -1;
/*char       prevl, prevp;*/
int        slidebarid;

Colormap   thecolormap;
int        ncolors, nplanes;
unsigned long foreground, background;
unsigned long palette[256];   /* my private palette */
boolean    notinfocus;

ed_rect *focus = NULL;
int stan = STATE_NOTHING;
static char **info_msgtext = NULL;
static int InfoNLines, InfoMaxLength;


/* ///////////////////////////////////////////////////////////////////////// */
XEvent  theevent;

static int     win_rect_num = 0;
static ed_rect *win_edr = NULL;

static ed_rect *lastwin = NULL;
static ed_rect *lastmenu = NULL;

static ed_rect errmsg_edr;
static int     errmsg_win;
static char    *errmsg_msgtext;

/* ///////////////////////////////////////////////////////////////////////// */
void redraw ()
{
  int i;

  for ( i = 0; i < win_rect_num; i++ )
    win_edr[i].redraw ( &win_edr[i] );
  if ( stan == STATE_MESSAGE && current_win == errmsg_win )
    errmsg_edr.redraw ( &errmsg_edr );
} /*redraw*/

void redraw_all ()
{
  int win, i;

  win = current_win;
  for ( i = 0; i < winnum; i++ ) {
    SetWindow ( i );
    redraw ();
  }
  SetWindow ( win );
} /*redraw_all*/

/* ///////////////////////////////////////////////////////////////////////// */
void DrawEmpty ( ed_rect *er )
{
} /*DrawEmpty*/

void EmptyMsg ( ed_rect *er, int msg, int key, int x, int y )
{
} /*EmptyMsg*/

void DrawMenu ( ed_rect *er )
{
  int     i, menu_rect_num;
  ed_rect *menu;

  XSetForeground ( thedisplay, thegc, c_black );
  XFillRectangle ( thedisplay, thepixmap, thegc,
                   er->x, er->y, er->w, er->h );

  GetMenuData ( er, &menu_rect_num, &menu );
  if ( !menu )
    return;
  for ( i = 0; i < menu_rect_num; i++ )
    menu[i].redraw ( &menu[i] );

  XCopyArea ( thedisplay, thepixmap, thewindow, thegc,
              er->x, er->y, er->w, er->h,
              er->x, er->y );
} /*DrawMenu*/

void MenuMsg ( ed_rect *er, int msg, int key, int x, int y )
{
  int     i, menu_rect_num;
  ed_rect *menu;

  GetMenuData ( er, &menu_rect_num, &menu );
  if ( !menu )
    return;

  if ( msg == msg_EXITING ) {
    if ( lastmenu )
      lastmenu->msgproc ( lastmenu, msg, key, x, y );
    lastmenu = NULL; 
    return;
  }

  for ( i = 0;  i < menu_rect_num;  i++ ) {
    if ( (x >= menu[i].x) && (x < menu[i].x+menu[i].w) &&
         (y >= menu[i].y) && (y < menu[i].y+menu[i].h) ) {
      if ( &menu[i] != lastmenu ) {
        if ( lastmenu )
          lastmenu->msgproc ( lastmenu, msg_EXITING, key, x, y );
        menu[i].msgproc ( &menu[i], msg_ENTERING, key, x, y );
        lastmenu = &menu[i];
      }
      menu[i].msgproc ( &menu[i], msg, key, x, y );
      return;
    }
  }
} /*MenuMsg*/

void DrawSwitch ( ed_rect *er )
{
  char    *title;
  boolean state, *switchvar;
  void    (*switchproc)(ed_rect *er);

  GetSwitchData ( er, &title, &switchvar, &switchproc );
  state = *switchvar;
  XSetForeground ( thedisplay, thegc, c_black );
  XFillRectangle ( thedisplay, thepixmap, thegc,
                   er->x, er->y, er->w, er->h );
  XSetForeground ( thedisplay, thegc, c_blue );
  XFillRectangle ( thedisplay, thepixmap, thegc,
             er->x+1, er->y+1, er->h-2, er->h-2 );
  XSetForeground ( thedisplay, thegc, c_white );
  XDrawRectangle ( thedisplay, thepixmap, thegc,
             er->x, er->y, er->h-1, er->h-1 );
  if ( state )
    XFillRectangle ( thedisplay, thepixmap, thegc,
                     er->x+4, er->y+4, er->h-8, er->h-8 );
  XDrawString ( thedisplay, thepixmap, thegc,
                     er->x+er->h+2, er->y+er->h-4,
                     title, strlen(title) );
  XCopyArea ( thedisplay, thepixmap, thewindow, thegc,
              er->x, er->y, er->w, er->h,
              er->x, er->y );
} /*DrawSwitch*/

void SwitchMsg ( ed_rect *er, int msg, int key, int x, int y )
{
  char    *title;
  boolean state, *switchvar;
  void    (*switchproc)(ed_rect *er);

  if ( (msg == msg_MCLICK && key & mouse_LBUTTON_DOWN) ||
       (msg == msg_KEY && key == 13) ) {
    GetSwitchData ( er, &title, &switchvar, &switchproc );
    state = *switchvar;
    *switchvar = !state;
    (*switchproc)( er );
  }
} /*SwitchMsg*/

void DrawButton ( ed_rect *er )
{
  char *title;
  void (*buttonproc)(ed_rect *er);

  GetButtonData ( er, &title, &buttonproc );
  XSetForeground ( thedisplay, thegc, c_dk_red );
  XFillRectangle ( thedisplay, thepixmap, thegc,
             er->x+1, er->y+1, er->w-2, er->h-2 );
  XSetForeground ( thedisplay, thegc, c_white );
  XDrawRectangle ( thedisplay, thepixmap, thegc,
             er->x, er->y, er->w-1, er->h-1 );
  XDrawString ( thedisplay, thepixmap, thegc,
                     er->x+2, er->y+er->h-5,
                     title, strlen(title) );
  XCopyArea ( thedisplay, thepixmap, thewindow, thegc,
              er->x, er->y, er->w, er->h,
              er->x, er->y );
} /*DrawButton*/

void ButtonMsg ( ed_rect *er, int msg, int key, int x, int y )
{
  char *title;
  void (*buttonproc)(ed_rect *er);

  if ( (msg == msg_MCLICK && key & mouse_LBUTTON_DOWN) ||
       (msg == msg_KEY && key == 13) ) {
    GetButtonData ( er, &title, &buttonproc );
    (*buttonproc)( er );
  }
} /*ButtonMsg*/

void BoundPoint ( ed_rect *er, int *x, int *y )
{
  int a;

  a = er->x + 2;
  *x = *x > a ? *x : a;
  a = er->x + er->w - 3;
  *x = *x < a ? *x : a;
  a = er->y + 2;
  *y = *y > a ? *y : a;
  a = er->y + er->h - 3;
  *y = *y < a ? *y : a;
} /*BoundPoint*/

void DrawSlidebar ( ed_rect *er )
{
  int   x;
  float *slipos;
  void  (*slideproc)(ed_rect *er);

  GetSlidebarData ( er, &slipos, &slideproc );
  XSetForeground ( thedisplay, thegc, c_blue );
  XFillRectangle ( thedisplay, thepixmap, thegc,
                   er->x+1, er->y+1, er->w-2, er->h-2 );
  XSetForeground ( thedisplay, thegc, c_white );
  XDrawRectangle ( thedisplay, thepixmap, thegc,
                   er->x, er->y, er->w-1, er->h-1 );
  x = er->x + 2 + (int)((*slipos)*(float)(er->w - 10));
  XFillRectangle ( thedisplay, thepixmap, thegc,
                   x, er->y+2, 6, 6 );
  XCopyArea ( thedisplay, thepixmap, thewindow, thegc,
              er->x, er->y, er->w, er->h,
              er->x, er->y );
} /*DrawSlidebar*/

void SlidebarMsg ( ed_rect *er, int msg, int key, int x, int y )
{
  float z;
  float *slipos;
  void  (*slideproc)(ed_rect *er);

  GetSlidebarData ( er, &slipos, &slideproc );
  if ( stan == STATE_NOTHING ) {
    if ( msg == msg_MCLICK && key & mouse_LBUTTON_DOWN ) {
      if ( x < er->x+5 ) x = er->x+5;
      else if ( x > er->x+er->w-5 ) x = er->x+er->w-5;
      stan = STATE_MOVINGSLIDE;
      focus = er;
      slidebarid = er->id;
      goto update;
    }
  }
  else if ( stan == STATE_MOVINGSLIDE && er->id == slidebarid ) {
    if ( msg == msg_MMOVE || msg == msg_MCLICK ) {
      if ( key & mouse_LBUTTON_DOWN ) {
        if ( x < er->x+5 ) x = er->x+5;
        else if ( x > er->x+er->w-5 ) x = er->x+er->w-5;
        if ( x != xx ) {
  update:
          xx = x;
          z = (float)(x-er->x-5)/(float)(er->w-10);
          if ( *slipos != z ) {
            *slipos = z;
            (*slideproc)( er );
            DrawSlidebar ( er );
          }
        }
      }
      else {
        stan = STATE_NOTHING;
        focus = NULL;
      }
    }
  }
} /*SlidebarMsg*/

float LinSlidebarValue ( float xmin, float xmax, float t )
{
  return xmin + t*(xmax-xmin);
} /*LinSlidebarValue*/

float LogSlidebarValue ( float xmin, float xmax, float t )
{
  return xmin + (xmax-xmin)*(exp(t)-1.0)/(exp(1.0)-1.0);
} /*LogSlidebarValue*/

void DrawText ( ed_rect *er )
{
  char    *title;

  GetTextTitle ( er, &title );
  XSetForeground ( thedisplay, thegc, c_white );
  XDrawString ( thedisplay, thepixmap, thegc, er->x, er->y+er->h-4,
                title, strlen(title) );
  XCopyArea ( thedisplay, thepixmap, thewindow, thegc,
              er->x, er->y, er->w, er->h,
              er->x, er->y );
} /*DrawText*/

/* ///////////////////////////////////////////////////////////////////////// */
static void RedrawErrorMessage ( ed_rect *er )
{
  int sl;

  XSetForeground ( thedisplay, thegc, c_red );
  XFillRectangle ( thedisplay, thepixmap, thegc,
                   er->x, er->y, er->w, er->h );
  XSetForeground ( thedisplay, thegc, c_white );
  XDrawRectangle ( thedisplay, thepixmap, thegc,
                   er->x, er->y, er->w-1, er->h-1 );
  sl = strlen ( errmsg_msgtext );
  XDrawString ( thedisplay, thepixmap, thegc,
                er->x+((er->w-6*sl)/2), er->y+16,
                errmsg_msgtext, sl );
  XDrawString ( thedisplay, thepixmap, thegc,
                er->x+(er->w/2)-7, er->y+er->h-6, "OK", 2 );
  XCopyArea ( thedisplay, thepixmap, thewindow, thegc,
              er->x, er->y, er->w, er->h, er->x, er->y );
} /*RedrawErrorMessage*/

static void RemoveErrorMessage ()
{
  if ( stan == STATE_MESSAGE ) {
    SetWindow ( errmsg_win );
    errmsg_win = -1;
    stan = STATE_NOTHING;
    focus = NULL;
    redraw ();
  }
} /*RemoveErrorMessage*/

static void ErrorMessageMsg ( ed_rect *er, int msg, int key, int x, int y )
{
  if ( current_win == errmsg_win && stan == STATE_MESSAGE && !notinfocus ) {
    if ( (msg == msg_MCLICK && PointInRect ( er, x, y )) ||
         (msg == msg_KEY && key == 0x0d) )
      RemoveErrorMessage ();
  }
} /*ErrorMessageMsg*/

void DisplayErrorMessage ( char *message )
{
  focus_win = errmsg_win = current_win;
  errmsg_msgtext = message;
  info_msgtext   = NULL;

  SetupEdRect ( &errmsg_edr, 0, 0, WIDTH-20, 60,
                (current_width-WIDTH+20)/2, (2*current_height)/3-30,
                &ErrorMessageMsg, &RedrawErrorMessage );
  stan = STATE_MESSAGE;
  focus = &errmsg_edr;
  if ( theevent.type != ConfigureNotify )
    XRaiseWindow ( thedisplay, thewindow );
  RedrawErrorMessage ( &errmsg_edr );
} /*DisplayErrorMessage*/

static void RedrawInfoMessage ( ed_rect *er )
{
  int  i, sl;
  char *msgstr;

  XSetForeground ( thedisplay, thegc, c_dk_grey );
  XFillRectangle ( thedisplay, thepixmap, thegc,
                   er->x, er->y, er->w, er->h );
  XSetForeground ( thedisplay, thegc, c_white );
  XDrawRectangle ( thedisplay, thepixmap, thegc,
                   er->x, er->y, er->w-1, er->h-1 );
  for ( i = 0; i < InfoNLines; i++ ) {
    msgstr = info_msgtext[i];
    if ( (sl = strlen ( msgstr )) )
      XDrawString ( thedisplay, thepixmap, thegc,
                    er->x+((er->w-6*InfoMaxLength)/2), er->y+16*(i+1),
                    msgstr, sl );
  }

  XDrawString ( thedisplay, thepixmap, thegc,
                er->x+(er->w/2)-7, er->y+er->h-6, "OK", 2 );
  XCopyArea ( thedisplay, thepixmap, thewindow, thegc,
              er->x, er->y, er->w, er->h, er->x, er->y );
} /*RedrawInfoMessage*/

void DisplayInfoMessage ( char **msglines )
{
  int height;

  focus_win = errmsg_win = current_win;
  errmsg_msgtext = NULL;
  info_msgtext   = msglines;
  for ( InfoNLines = InfoMaxLength = 0;  *msglines;  msglines++, InfoNLines++ )
    InfoMaxLength = max ( InfoMaxLength, strlen(*msglines) );

  height = 44 + 16*InfoNLines;
  SetupEdRect ( &errmsg_edr, 0, 0, WIDTH-20, height,
                (current_width-WIDTH+20)/2, (current_height-height)/2,
                &ErrorMessageMsg, &RedrawInfoMessage );
  stan = STATE_MESSAGE;
  focus = &errmsg_edr;
  if ( theevent.type != ConfigureNotify )
    XRaiseWindow ( thedisplay, thewindow );
  RedrawInfoMessage ( &errmsg_edr );
} /*DisplayInfoMessage*/

/* ///////////////////////////////////////////////////////////////////////// */
void SetupEdRect ( ed_rect *edr, int en, int id, int w, int h, int x, int y,
                   message_proc msgproc, redraw_proc redraw )
{
  edr[en].id = id;
  edr[en].w = w;
  edr[en].h = h;
  edr[en].x = x;
  edr[en].y = y;
  edr[en].msgproc = msgproc;
  edr[en].redraw = redraw;
} /*SetupEdRect*/

void SetWinEdRect ( int rect_num, ed_rect *edr )
{
  windesc[current_win].win_rect_num = win_rect_num = rect_num;
  windesc[current_win].win_edr = win_edr = edr;
} /*SetWinEdRect*/

/* ///////////////////////////////////////////////////////////////////////// */
boolean PointInRect ( ed_rect *edr, int x, int y )
{
  return (x >= edr->x) && (x < edr->x+edr->w) &&
         (y >= edr->y) && (y < edr->y+edr->h);
} /*PointInRect*/

boolean RectanglesIntersect ( int wa, int ha, int xa, int ya,
                              int wb, int hb, int xb, int yb )
{
  if ( xa >= xb+wb || xb >= xa+wa || ya >= yb+wb || yb >= ya+wa )
    return false;
  else
    return true;
} /*RectanglesIntersect*/

void dispatch_message ( unsigned int msg, unsigned int key, int x, int y )
{
  int i;

  if ( msg == msg_NULL )
    return;

  if ( msg == msg_KEY && stan != STATE_MESSAGE )
    process_key ( &key );

  if ( focus ) {
    if ( (notinfocus = focus_win != current_win) )
      SetWindow ( focus_win );
    focus->msgproc ( focus, msg, key, x, y );
    if ( focus_win != current_win )
      SetWindow ( current_win );
  }
  else {
    for ( i = 0;  i < win_rect_num;  i++ ) {
      if ( PointInRect ( &win_edr[i], x, y ) ) {
        if ( &win_edr[i] != lastwin ) {
          if ( lastwin )
            lastwin->msgproc ( lastwin, msg_EXITING, key, x, y );
          win_edr[i].msgproc ( &win_edr[i], msg_ENTERING, key, x, y );
          lastwin = &win_edr[i];
        }
        win_edr[i].msgproc ( &win_edr[i], msg, key, x, y );
        if ( focus )
          focus_win = current_win;
        return;
      }
    }
  }
} /*dispatch_message*/

/* ///////////////////////////////////////////////////////////////////////// */
int NewWindow ( char *p_name )
{
  Window win, root;
  int x, y;
  unsigned int border_width, depth;

  if ( winnum < MAX_WINDOWS ) {
    background = BlackPixel ( thedisplay, thescreen );
    foreground = WhitePixel ( thedisplay, thescreen );
    win = XCreateSimpleWindow ( thedisplay, DefaultRootWindow (thedisplay),
                                thehints.x, thehints.y,
                                thehints.width, thehints.height,
                                BORDER_WIDTH, foreground, background );
    XSetStandardProperties ( thedisplay, win, p_name, p_name,
                             None, _argv, _argc, &thehints );
    XSelectInput ( thedisplay, win,
                   ButtonPressMask | ExposureMask | ButtonReleaseMask |
                   PointerMotionMask | KeyPressMask | StructureNotifyMask );
    XMapWindow ( thedisplay, win );

    windesc[winnum].thewindow = win;
    windesc[winnum].win_rect_num = 0;
    windesc[winnum].win_edr = NULL;
    XGetGeometry ( thedisplay, win, &root, &x, &y,
                   &windesc[winnum].current_width,
                   &windesc[winnum].current_height,
                   &border_width, &depth );

    win = winnum;
    winnum++;
    return win;
  }
  else return -1;
} /*NewWindow*/

boolean SetWindow ( int win )
{
  if ( win >= 0 && win < winnum ) {
    if ( win != current_win ) {
      if ( current_win >= 0 && current_win < winnum ) {
        windesc[current_win].thewindow      = thewindow;
        windesc[current_win].win_rect_num   = win_rect_num;
        windesc[current_win].win_edr        = win_edr;
        windesc[current_win].current_width  = current_width;
        windesc[current_win].current_height = current_height;
      }
      thewindow      = windesc[win].thewindow;
      win_rect_num   = windesc[win].win_rect_num;
      win_edr        = windesc[win].win_edr;
      current_width  = windesc[win].current_width;
      current_height = windesc[win].current_height;
      current_win    = win;
    }
    return true;
  }
  else
    return false;
} /*SetWindow*/

int CurrentWindow ()
{
  return current_win;
} /*CurrentWindow*/

int FindWindowNum ( Window thewin )
{
  int i;

  for ( i = 0; i < winnum; i++ )
    if ( thewin == windesc[i].thewindow )
      return i;
  return -1;
} /*FindWindowNum*/

/* ///////////////////////////////////////////////////////////////////////// */
void get_message ( unsigned int *msg, unsigned int *key, int *x, int *y )
{
  Window         r_w, c_w;
  int            win;
  int            x_r, y_r;
  unsigned int   xbutton_mask;
  unsigned int   newmask = 0;
  XComposeStatus status;

  XNextEvent ( thedisplay, &theevent );
  win = FindWindowNum ( theevent.xany.window );
  if ( win != current_win )
    SetWindow ( win );

  switch ( theevent.type )
  {
case Expose:
    if ( theevent.xexpose.count == 0 )
      redraw ();
    *msg = msg_NULL;
    break;

case ButtonPress:
    switch ( theevent.xbutton.button )
    {
  case Button1:
      newmask = mouse_buttons | mouse_LBUTTON_DOWN;
      *key = newmask | mouse_LBUTTON_CHANGE;
      break;
  case Button2:
      newmask = mouse_buttons | mouse_MBUTTON_DOWN;
      *key = newmask | mouse_MBUTTON_CHANGE;
      break;
  case Button3:
      newmask = mouse_buttons | mouse_RBUTTON_DOWN;
      *key = newmask | mouse_RBUTTON_CHANGE;
      break;
    }
    mouse_buttons = newmask;
    *msg = msg_MCLICK;
    *x = mouse_x;
    *y = mouse_y;
    break;

case ButtonRelease:
    switch ( theevent.xbutton.button )
    {
  case Button1: 
      newmask = mouse_buttons & (~mouse_LBUTTON_DOWN);
      *key = newmask | mouse_LBUTTON_CHANGE;
      break;
  case Button2:
      newmask = mouse_buttons & (~mouse_MBUTTON_DOWN);
      *key = newmask | mouse_MBUTTON_CHANGE;
      break;
  case Button3:
      newmask = mouse_buttons & (~mouse_RBUTTON_DOWN);
      *key = newmask | mouse_RBUTTON_CHANGE;
      break;
    }
    mouse_buttons = newmask;
    *msg = msg_MCLICK;
    *x = mouse_x;
    *y = mouse_y;
    break;

case MotionNotify:
    XQueryPointer ( thedisplay, thewindow,
                    &r_w, &c_w, &x_r, &y_r, &mouse_x, &mouse_y, &xbutton_mask );
    if ( mouse_x != prevx || mouse_y != prevy ) {
      *x = prevx = mouse_x;
      *y = prevy = mouse_y;
      *key = mouse_buttons;
      *msg = msg_MMOVE;
    }
    else *msg = msg_NULL;
    break;

case ConfigureNotify:
    if ( (theevent.xconfigure.width != windesc[win].current_width ||
          theevent.xconfigure.height != windesc[win].current_height) ) {
      if ( win == errmsg_win && stan == STATE_MESSAGE ) {
        RemoveErrorMessage ();
        resize_edwin ();
        if ( errmsg_msgtext )
          DisplayErrorMessage ( errmsg_msgtext );
        else if ( info_msgtext )
          DisplayInfoMessage ( info_msgtext );
      }
      else
        resize_edwin ();
    }
    *msg = msg_NULL;
    break;

case KeyPress:
    *msg = msg_KEY;
    *key = 0;
    XLookupString ( &theevent.xkey, (char*)key, 1, &thekeysym, &status );
    break;

default:
    ProcessOtherEvent ();
    *msg = msg_NULL;
    break;
  }
} /*get_message*/

void MakeTCPalette ()
{
  XColor        cdefs[32];
  int           i, r, g, b;
  unsigned long mask;
  unsigned long red_mask, green_mask, blue_mask;
  unsigned long red_step = 0, green_step = 0, blue_step = 0;
  int           red_cols, green_cols, blue_cols;

  for ( i = 0, mask = 0x01; i < nplanes; i++, mask = mask << 1 )
    cdefs[i].pixel = mask;

  nplanes = nplanes < 32 ? nplanes : 32;
  XQueryColors ( thedisplay, thecolormap, cdefs, nplanes );

  red_mask   = 0;
  green_mask = 0;
  blue_mask  = 0;
  red_cols   = 1;
  green_cols = 1;
  blue_cols  = 1;
  for ( i = 0, mask = 0x01; i < nplanes; i++, mask = mask << 1 )
  {
/*
    printf ( "R: %5d, G: %5d, B: %5d\n", cdefs[i].red, cdefs[i].green, cdefs[i].blue );
*/
    if (cdefs[i].red)
    {
      if ( !red_mask ) red_step = mask;
      red_mask   |= mask;
      red_cols   *= 2;
    }
    if (cdefs[i].green)
    {
      if ( !green_mask ) green_step = mask;
      green_mask |= mask;
      green_cols *= 2;
    }
    if (cdefs[i].blue)
    {
      if ( !blue_mask ) blue_step = mask;
      blue_mask  |= mask;
      blue_cols  *= 2;
    }
  }
  for ( r = 0; r < 8; r++ )
    for ( g = 0; g < 8; g++ )
      for ( b = 0; b < 4; b++ )
      {
        i = r + (g << 3) + (b << 6);
        palette[i] = (r*(red_cols-1)/7)*red_step +
                     (g*(green_cols-1)/7)*green_step +
                     (b*(blue_cols-1)/3)*blue_step;
      }
} /*MakeTCPalette*/

void GetWindowSize ()
{
  Window root;
  int x, y;
  unsigned int border_width, depth;

  XGetGeometry ( thedisplay, thewindow, &root, &x, &y,
                 &current_width, &current_height, &border_width, &depth );
  if ( current_width > MAX_WIDTH )
    current_width = MAX_WIDTH;
  if ( current_height > MAX_HEIGHT )
    current_height = MAX_HEIGHT;
} /*GetWindowSize*/

void init ( int argc, char *argv[] )
{
  _argc = argc;
  _argv = argv;

  p_name = argv[0];
  thedisplay = XOpenDisplay ( "" );
  thescreen = DefaultScreen ( thedisplay );

  thehints.flags = PPosition | PSize | PMinSize;
  thehints.height = HEIGHT;
  thehints.width = WIDTH;
  thehints.x = 0;
  thehints.y = 100;
  thehints.min_height = HEIGHT;
  thehints.min_width = WIDTH;

  NewWindow ( p_name );
  SetWindow ( 0 );

  thegc = XCreateGC ( thedisplay, thewindow, 0, 0 );

  thevisual = XDefaultVisual ( thedisplay, thescreen );
  thecolormap = XDefaultColormap ( thedisplay, thescreen );
  nplanes = XDisplayPlanes ( thedisplay, thescreen );
  ncolors = XDisplayCells ( thedisplay, thescreen );
  thepixmap = XCreatePixmap ( thedisplay, thewindow,
                              MAX_WIDTH, MAX_HEIGHT, nplanes );
  thecursor0 = XCreateFontCursor ( thedisplay, XC_crosshair );
  thecursor1 = XCreateFontCursor ( thedisplay, XC_watch );
  XDefineCursor ( thedisplay, thewindow, None );
  GetWindowSize ();

  switch ( thevisual->class )
  {
case PseudoColor :
    break;
case DirectColor :
    break;
case TrueColor :
    MakeTCPalette ();
    break;
default :
    ;
  }
} /*init*/

void cleanup ()
{
  XFreePixmap ( thedisplay, thepixmap );
  XFreeGC ( thedisplay, thegc );
  XUndefineCursor ( thedisplay, thewindow );
  XDestroyWindow ( thedisplay, thewindow );
  XCloseDisplay ( thedisplay );
} /*cleanup*/

int main ( int argc, char *argv[] )
{
  unsigned int msg, key;
  int x, y;

  init ( argc, argv );
  init_edwin ( argc, argv );
  while ( !done ) {
    get_message ( &msg, &key, &x, &y );
    dispatch_message ( msg, key, x, y );
  }
  destroy_edwin ();
  cleanup ();
  exit ( 0 );
} /*main*/

