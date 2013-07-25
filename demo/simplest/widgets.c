
/* //////////////////////////////////////////////////// */
/* This file is a part of the BSTools procedure package */
/* written by Przemyslaw Kiciak.                        */
/* //////////////////////////////////////////////////// */

#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <malloc.h>
#include <string.h>

#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <X11/cursorfont.h>

#include "pkgeom.h"
#include "xgedit.h"

#include "widgets.h"

xge_fourww       fourww;
xge_string_ed    stred;
xge_int_widget   intwidg0, intwidg1;
quaterniond      q;
trans3d          qtr;
xge_quatrotballd qball;
char             text[21] = "cosik";
float            dialpos1, dialpos2;
int              intvalue0 = 1, intvalue1 = -1;

xge_listbox    filelist, dirlist;
xge_widget     *mainmenu, *fourwwedr, *popupmenu, *colourpopup, *colourwdg;

float sl2f[2] = {0.0,1.0};

float colour[3] = {1.0,1.0,1.0};

Cursor cursor[XC_num_glyphs/2];

char filter[] = "*-*.?z";

/* ///////////////////////////////////////////////////////////////////////// */
boolean ColourWMsg ( xge_widget *er, int msg, int key, short x, short y )
{
  return false;
} /*ColourWMsg*/

void ColourWRedraw ( xge_widget *er, boolean onscreen )
{
  float *colour;

  colour = er->data0;
  xgeSetForeground ( xge_PixelColourf ( colour[0], colour[1], colour[2] ) );
  xgeFillRectangle ( er->w-2, er->h-2, er->x+1, er->y+1 );
  xgeSetForeground ( xgec_White );
  xgeDrawRectangle ( er->w-1, er->h-1, er->x, er->y );
  if ( onscreen )
    xgeCopyRectOnScreen ( er->w, er->h, er->x, er->y );
} /*ColourWRedraw*/

/* ///////////////////////////////////////////////////////////////////////// */
void RysujOkno ( xge_widget *er, boolean onscreen )
{
  xgeSetForeground ( xgec_Blue6 );
  xgeFillRectangle ( er->w-2, er->h-2, er->x+1, er->y+1 );
  xgeSetForeground ( xgec_CornflowerBlue );
  xgeDrawRectangle ( er->w-1, er->h-1, er->x, er->y );
  if ( onscreen )
    xgeCopyRectOnScreen ( er->w, er->h, er->x, er->y );
} /*RysujOkno*/

boolean OknoMsg ( xge_widget *er, int msg, int key, short x, short y )
{
/*  printf ( "%d, %d, %d, %d, %d\n", er->id, msg, key, x, y ); */
  switch ( msg ) {
case xgemsg_ENTERING:
    xge_SetCurrentWindowCursor ( xgeCURSOR_ARROW );
    break;
case xgemsg_EXITING:
    xge_SetCurrentWindowCursor ( xgeCURSOR_DEFAULT );
    break;
default:
    break;
  }
  return true;
} /*OknoMsg*/

boolean MyMenuMsg ( xge_widget *er, int msg, int key, short x, short y )
{
  switch ( msg ) {
case xgemsg_ENTERING:
    if ( intvalue1 >= 0 )
      xge_SetCurrentWindowCursor ( cursor[intvalue1] );
    break;
case xgemsg_EXITING:
    if ( intvalue1 >= 0 )
      xge_SetCurrentWindowCursor ( None );
    break;
default:
    break;    
  }
  return xge_MenuMsg ( er, msg, key, x, y );
} /*MyMenuMsg*/

/* ///////////////////////////////////////////////////////////////////////// */
int CallBack ( xge_widget *er, int msg, int key, short x, short y )
{
  if ( er ) {
/*    printf ( "%d, %d, %d, %d, %d\n", er->id, msg, key, x, y ); */
    switch ( msg ) {
case xgemsg_BUTTON_COMMAND:
      switch ( er->id ) {
  case 0:
        xge_done = 1;
        break;

  case 6:
        xge_SetupFileList ( &filelist, ".", filter );
        xge_SetupDirList ( &dirlist, ".", NULL, NULL );
        xge_AddPopup ( popupmenu );
        xge_GrabFocus ( popupmenu, true );
        break;

  case 7:
        xge_RemovePopup ( true );
        xge_ClearListBox ( &filelist );
        xge_ClearListBox ( &dirlist );
        break;

  case 8:
        xge_AddPopup ( colourpopup );
        xge_GrabFocus ( colourpopup, true );
        break;

  case 16:
        xge_RemovePopup ( true );
        break;

  default:
        break;
      }
      break;

case xgemsg_SLIDEBAR_COMMAND:
/* RGB slidebars are supposed to have identifiers 3n, 3n+1 and 3n+2 */
      switch ( er->id ) {
  case 12:
  case 13:
  case 14:
        xge_SetClipping ( colourwdg );
        colourwdg->redraw ( colourwdg, true );
        break;
  default:
        break;
      }
      break;

case xgemsg_INT_WIDGET_COMMAND:
      switch ( er->id ) {
  case 4:
        intvalue0 = key;
        printf ( "%d\n", intvalue0 );
        break;

  case 5:
        intvalue1 = key;
        if ( key < 0 )
          xge_SetCurrentWindowCursor ( None );
        else
          xge_SetCurrentWindowCursor ( cursor[key] );
        break;

  default:
        break;
      }
      break;

case xgemsg_LISTBOX_ITEM_PICK:
      switch ( er->id ) {
  case 8:
        printf ( "dir: %s\n",
                 &dirlist.itemstr[dirlist.itemind[dirlist.currentitem]] );
        if ( !chdir ( &dirlist.itemstr[dirlist.itemind[dirlist.currentitem]] ) ) {
          xge_SetupFileList ( &filelist, ".", filter );
          xge_SetupDirList ( &dirlist, ".", NULL, NULL );
          dirlist.er->redraw ( dirlist.er, true );
          filelist.er->redraw ( filelist.er, true );
        }
        break;

  case 9:
        printf ( "file: %s\n",
                 &filelist.itemstr[filelist.itemind[filelist.currentitem]] );
        break;

  default:
        break;
      }
      break;

case xgemsg_TEXT_EDIT_VERIFY:
      return 1;

case xgemsg_TEXT_EDIT_ENTER:
      return 1;

case xgemsg_QUATROTBALL_COMMAND:
      return 1;

default:
      break;
    }
  }
  else {
/*    printf ( "%d, %d, %d, %d\n", msg, key, x, y ); */
    switch ( msg ) {
case xgemsg_KEY:
      switch ( key ) {
  case 'Q': case 'q':
        xge_done = 1;
        break;
  default:
        break;
      }
      break;

case xgemsg_RESIZE:
      fourww.er->msgproc ( fourww.er, xgemsg_RESIZE, 0,
                           (short)(xge_current_width-110), xge_current_height );
      mainmenu->msgproc ( mainmenu, xgemsg_RESIZE, 0,
                          110, xge_current_height );
      xge_Redraw ();
      break;

default:
      break;
    }
  }
  return 0;
} /*CallBack*/

/* ///////////////////////////////////////////////////////////////////////// */
static char b0[] = "Quit";
static char dialtext1[] = "aside";
static char dialtext2[] = "below";
static char intwtext0[] = "int widget";
static char intwtext1[] = "cursor";
static char b6[] = "File";
static char b7[] = "Cancel";
static char b8[] = "Colour";
static char b14[] = "OK";

void PrepareCursors ( void )
{
  int i;

  for ( i = 0; i < XC_num_glyphs; i += 2 )
    cursor[i/2] = XCreateFontCursor ( xgedisplay, i );
} /*PrepareCursors*/

void init_edwin ( void  )
{
  xge_widget *rp;

  pkv_InitScratchMem ( 655360 );
  PrepareCursors ();
        /* the main menu */
  rp = xge_NewButton ( 0, NULL, 0, 60, 18, 0, 0, b0 );
  rp = xge_NewStringEd ( 0, rp, 1, 100, 18, 0, 20, 20, text, &stred );
  rp = xge_NewDialf ( 0, rp, 2, 110, 40, 0, 40, dialtext1, &dialpos1 );
  rp = xge_NewDialf ( 0, rp, 3, 40, 60, 0, 82, dialtext2, &dialpos2 );
  rp = xge_NewIntWidget ( 0, rp, 4, 110, 18, 0, 144, 0, 20,
                          &intwidg0, intwtext0, &intvalue0 );
  rp = xge_NewIntWidget ( 0, rp, 5, 110, 18, 0, 164, -1, (XC_num_glyphs/2)-1,
                          &intwidg1, intwtext1, &intvalue1 );
  rp = xge_NewButton ( 0, rp, 6, 60, 18, 0, 184, b6 );
  rp = xge_NewVSlidebar2f ( 0, rp, 7, 10, 109, 0, 204, sl2f );
  rp = xge_NewButton ( 0, rp, 8, 60, 18, 0, 320, b8 );
  rp = xge_NewQuatRotBalld ( 0, rp, 9, 80, 80, 14, 204, 80,
                             &qball, &q, &qtr, NULL );
  mainmenu = xge_NewMenu ( 0, NULL, 1, 110, xge_HEIGHT, 0, 0, rp );
  mainmenu->msgproc = MyMenuMsg;

        /* the four window widget */
  rp = xge_NewWidget ( 0, NULL, 0, xge_WIDTH-110, xge_HEIGHT, 110, 0,
                       NULL, NULL, OknoMsg, RysujOkno );
  rp = xge_NewWidget ( 0, rp, 1, xge_WIDTH-110, xge_HEIGHT, 110, 0,
                       NULL, NULL, OknoMsg, RysujOkno );
  rp = xge_NewWidget ( 0, rp, 2, xge_WIDTH-110, xge_HEIGHT, 110, 0,
                       NULL, NULL, OknoMsg, RysujOkno );
  rp = xge_NewWidget ( 0, rp, 3, xge_WIDTH-110, xge_HEIGHT, 110, 0,
                       NULL, NULL, OknoMsg, RysujOkno );
  fourwwedr = xge_NewFourWW ( 0, mainmenu, 2, xge_WIDTH-110, xge_HEIGHT, 110, 0,
                              rp, &fourww );

        /* the popup menu */
  rp = xge_NewButton ( 0, NULL, 7, 60, 18,
                       20+(xge_WIDTH-40-60)/2, 40+160-24, b7 );
  rp = xge_NewListBox ( 0, rp, 8, 180, 115, 20+20, 40+10, &dirlist );
  rp = xge_NewListBox ( 0, rp, 9, 180, 115, 20+236, 40+10, &filelist );
  popupmenu = xge_NewFMenu ( 0, NULL, 2, xge_WIDTH-40, 160, 20, 40, rp );

        /* the colour popup */
/* RGB slidebars are supposed to have identifiers 3n, 3n+1 and 3n+2 */
  rp = xge_NewSlidebarfRGB ( 0, NULL, 12, 109, 10, 50, 40, &colour[0] );
  rp = xge_NewSlidebarfRGB ( 0, rp, 13, 109, 10, 50, 54, &colour[1] );
  rp = xge_NewSlidebarfRGB ( 0, rp, 14, 109, 10, 50, 68, &colour[2] );
  rp = colourwdg = xge_NewWidget ( 0, rp, 15, 50, 50, 180, 40, colour, NULL,
                       ColourWMsg, ColourWRedraw );
  rp = xge_NewButton ( 0, rp, 16, 60, 18, 70, 88, b14 );
  colourpopup = xge_NewFMenu ( 0, NULL, 3, 240, 100, 30, 30, rp );

  xge_SetWinEdRect ( fourwwedr );
  xge_Redraw ();
} /*init_edwin*/

void destroy_edwin ( void )
{
  pkv_DestroyScratchMem ();
} /*destroy_edwin*/

int main ( int argc, char *argv[] )
{
  setvbuf ( stdout, NULL, _IONBF, 0 );
  xge_Init ( argc, argv, CallBack, NULL );
  init_edwin ();
  xge_MessageLoop ();
  destroy_edwin ();  
  xge_Cleanup ();  
  exit ( 0 );
} /*main*/   

