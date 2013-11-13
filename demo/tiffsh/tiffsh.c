
/* //////////////////////////////////////////////////// */
/* This file is a part of the BSTools procedure package */
/* written by Przemyslaw Kiciak.                        */
/* //////////////////////////////////////////////////// */
 
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <tiffio.h>
#include <sys/times.h>
#include <X11/Xlib.h>
#include <X11/Xutil.h>

#include "pkvaria.h"
#include "pknum.h"
#include "pkgeom.h"
#include "multibs.h"
#include "camera.h"
#include "xgedit.h"

#include "tiffsh.h"

xge_widget *imagewin;
XImage     *xim1, *xim2;
char       *imagedata1, *imagedata2;

const char slopt_file[] = "tiffsh.pr";
const char chapter_file[] = "tiffsh.ch";
int        nchapters;
int        ch[MAXCHAPTERS];
int        nsl, speed;
char       slopt[MAXSLIDES];

const char  file_filter[] = "*.tif";
xge_widget  *flist, *dlist, *popup;
xge_listbox dlbox, flbox;
char current_directory[MAX_PATH_LGT+1];
char current_dir[MAX_PATH_SHRT+1];
char filename[MAX_FILENAME_LGT+1];
boolean popup_is_on = false, cursor_is_on = false;
boolean showhidden = false;

char txtLoad[] = "Load";
char txtExit[] = "Exit";
char txtHidden[] = "hidden";

int     animadir = 0, animapos = 0;
double  acp[4] = {0.0,0.0,1.0,1.0};
clock_t time_cnt = 0;
double  time_start;


void ReadSlOptFile ( void )
{
  FILE *f;
  int  r, i;
  char c;

  memset ( slopt, 1, MAXSLIDES );
  f = fopen ( slopt_file, "r+" );
  if ( f ) {
    for (;;) {
      r = fscanf ( f, "%d", &i );
      if ( r == 1 && i >= 0 && i < MAXSLIDES )
        slopt[i] = 0;
      else if ( r == EOF )
        break;
      else {
        do
          r = fscanf ( f, "%c", &c );
        while ( r == 1 && c != ':' );
        if ( r != 1 )
          break;
        r = fscanf ( f, "%d", &i );
        if ( r == 1 && i >= 0 && i < MAXSLIDES )
          slopt[i] = 3;
        else
          break;
      }
    }
    fclose ( f );
  }
} /*ReadSlOptFile*/

void ReadChapterFile ( void )
{
  FILE *f;
  int  r;

  f = fopen ( chapter_file, "r+" );
  nchapters = 0;
  if ( f ) {
    for (;;) {
      r = fscanf ( f, "%d", &ch[nchapters] );
      if ( r == 1 && nchapters < MAXCHAPTERS )
        nchapters ++;
      else
        break;
    }
    fclose ( f );
  }
} /*ReadChapterFile*/

void MakeFileList ( void )
{
  if ( !xge_SetupFileList ( &flbox, current_directory, file_filter, showhidden ) ) {
    printf ( "cannot create the file list\n" );
    exit ( 1 );
  }
/*  printf ( "nitems = %d\n", flbox.nitems );*/
  xge_GetCurrentListBoxString ( &flbox, filename );
/*  printf ( "%s\n", filename ); */
} /*MakeFileList*/

void ClearImage ( void )
{
  memset ( xim1->data, 0, IMHEIGHT*xim1->bytes_per_line );
} /*ClearImage*/

boolean ReadImage ( char *fn )
{
  TIFF   *tif;
  int    w, h, p, npix, x, y, ww, hh;
  uint32 *data;
  byte   *b;

  data = NULL;
  if ( fn[0] ) {
    memset ( xim2->data, 0, IMHEIGHT*xim2->bytes_per_line );
    tif = TIFFOpen ( fn, "r" );
    if ( tif ) {
      TIFFGetField ( tif, TIFFTAG_IMAGEWIDTH, &w );
      TIFFGetField ( tif, TIFFTAG_IMAGELENGTH, &h );
      npix = w*h;
      data = (uint32*)_TIFFmalloc ( npix*sizeof(uint32) );
      if ( !data )
        goto failure;
      memset ( data, 0, npix*sizeof(uint32) );
      if ( !TIFFReadRGBAImage ( tif, w, h, data, 0 ) )
        goto failure;
      TIFFClose ( tif );
      ww = min ( w, IMWIDTH );
      hh = min ( h, IMHEIGHT );
      for ( y = p = 0;  y < hh;  y++ )
        for ( x = 0;  x < ww;  x++, p++ ) {
          b = (byte*)&data[p];
          XPutPixel ( xim2, x, h-1-y, xge_PixelColour ( *b, *(b+1), *(b+2) ) );
        }
      _TIFFfree ( data );
      return true;
    }
    else
      return false;

  failure:
    if ( data )
      _TIFFfree ( data );
    TIFFClose ( tif );
    return false;
  }
  else {
    ClearImage ();
    return true;
  }
} /*ReadImage*/

void DrawImageWin ( xge_widget *er, boolean onscreen )
{
  switch ( animadir ) {
case 0:
    XPutImage ( xgedisplay, xgepixmap, xgegc, xim1,
                0, 0, 0, 0, IMWIDTH, IMHEIGHT );
    break;
case 1:
    XPutImage ( xgedisplay, xgepixmap, xgegc, xim1,
                0, animapos, 0, 0, IMWIDTH, IMHEIGHT-animapos );
    XPutImage ( xgedisplay, xgepixmap, xgegc, xim2,
                0, 0, 0, IMHEIGHT-animapos, IMWIDTH, animapos );
    break;
case 2:
    XPutImage ( xgedisplay, xgepixmap, xgegc, xim1,
                0, 0, 0, animapos, IMWIDTH, IMHEIGHT-animapos );
    XPutImage ( xgedisplay, xgepixmap, xgegc, xim2,
                0, IMHEIGHT-animapos, 0, 0, IMWIDTH, animapos );
    break;
case 3:
    XPutImage ( xgedisplay, xgepixmap, xgegc, xim1,
                animapos, 0, 0, 0, IMWIDTH-animapos, IMHEIGHT );
    XPutImage ( xgedisplay, xgepixmap, xgegc, xim2,
                0, 0, IMWIDTH-animapos, 0, animapos, IMHEIGHT );
    break;
case 4:
    XPutImage ( xgedisplay, xgepixmap, xgegc, xim1,
                0, 0, animapos, 0, IMWIDTH-animapos, IMHEIGHT );
    XPutImage ( xgedisplay, xgepixmap, xgegc, xim2,
                IMWIDTH-animapos, 0, 0, 0, animapos, IMHEIGHT );
    break;
  }
  if ( onscreen )
    xgeCopyRectOnScreen ( er->w, er->h, er->x, er->y );
} /*DrawImageWin*/

void InitAnimation ( int dir, int _speed )
{
  struct tms tt;

  ReadImage ( filename );
  animadir = dir;
  speed = _speed;
  if ( dir == 0 ) {
    animapos = 0;
    memcpy ( imagedata1, imagedata2, IMHEIGHT*xim1->bytes_per_line );
    xge_Redraw ();
  }
  else {
    if ( !animapos )
      time_cnt = times ( &tt );
    xge_PostIdleCommand ( 0, 0, 0 );
  }
} /*InitAnimation*/

void ContinueAnimation ( void )
{
  struct tms tt;
  clock_t    ct;
  double     time_var, pos;

  ct = times ( &tt );
  if ( speed ) {
    time_var = ((double)(ct-time_cnt))/(speed*ATIME*(double)sysconf(_SC_CLK_TCK));
    time_var = max ( time_var, 0.0 );
    time_var = min ( time_var, 1.0 );
  }
  else
    time_var = 1.0;
  mbs_BCHornerC1d ( 3, acp, time_var, &pos );
/*  animapos += ANIMASTEP; */
  switch ( animadir ) {
case 1:
case 2:
    animapos = (int)((double)IMHEIGHT*pos+0.5);
    xge_Redraw ();
    if ( animapos >= IMHEIGHT ) {
      memcpy ( imagedata1, imagedata2, IMHEIGHT*xim1->bytes_per_line );
      animapos = 0;
      animadir = 0;
    }
    else
      xge_PostIdleCommand ( 0, 0, 0 );
    return;
case 3:
case 4:
    animapos = (int)((double)IMWIDTH*pos+0.5);
    xge_Redraw ();
    if ( animapos >= IMWIDTH ) {
      memcpy ( imagedata1, imagedata2, IMHEIGHT*xim1->bytes_per_line );
      animapos = 0;
      animadir = 0;
    }
    else
      xge_PostIdleCommand ( 0, 0, 0 );
    return;
  }
} /*ContinueAnimation*/

void ShowNextPicture ( boolean anim )
{
  int i;

  i = flbox.currentitem;
  if ( i < flbox.nitems-1 ) {
    xge_MoveInListBox ( &flbox, 1 );
    xge_GetCurrentListBoxString ( &flbox, filename );
    ReadImage ( filename );
    if ( anim )
      InitAnimation ( 1, slopt[i] );
    else
      InitAnimation ( 0, 0 );
  }
} /*ShowNextPicture*/

void ShowPreviousPicture ( boolean anim )
{
  int i;

  if ( flbox.currentitem > 0 ) {
    xge_MoveInListBox ( &flbox, -1 );
    i = flbox.currentitem;
    xge_GetCurrentListBoxString ( &flbox, filename );
    ReadImage ( filename );
    if ( anim )
      InitAnimation ( 2, slopt[i] );
    else
      InitAnimation ( 0, 0 );
  }
} /*ShowPreviousPicture*/

void ShowNextChapter ( void )
{
  int i;

  if ( nchapters == 0 )
    ShowNextPicture ( true );
  else {
    for ( i = 0; i < nchapters; i++ )
      if ( ch[i] > flbox.currentitem )
        break;
    if ( ch[i] < flbox.nitems ) {
      xge_MoveInListBox ( &flbox, ch[i]-flbox.currentitem );
      xge_GetCurrentListBoxString ( &flbox, filename );
      ReadImage ( filename );
      InitAnimation ( 3, 1 );
    }
  }
} /*ShowNextChapter*/

void ShowPreviousChapter ( void )
{
  int i;

  if ( nchapters == 0 )
    ShowPreviousPicture ( true );
  else {
    if ( flbox.currentitem == 0 )
      i = nchapters-1;
    else {
      for ( i = nchapters-1; i >= 0; i-- )
        if ( ch[i] < flbox.currentitem )
          break;
    }
    if ( ch[i] >= 0 ) {
      xge_MoveInListBox ( &flbox, ch[i]-flbox.currentitem );
      xge_GetCurrentListBoxString ( &flbox, filename );
      ReadImage ( filename );
      InitAnimation ( 4, 1 );
    }
  }
} /*ShowPreviousChapter*/

void DrawBackground ( void )
{
  int i;

  xge_SetWinEdRect ( imagewin );
  xge_ResetClipping ();
  xgeSetForeground ( xgec_Black );
  xgeFillRectangle ( IMWIDTH, IMHEIGHT, 0, 0 );
  xgeSetForeground ( xgec_RoyalBlue1 );
  for ( i = -IMHEIGHT+32; i < IMWIDTH+IMHEIGHT; i += 64 ) {
    xgeDrawLine ( i, 0, i+IMHEIGHT, IMHEIGHT );
    xgeDrawLine ( i, 0, i-IMHEIGHT, IMHEIGHT );
  }
  XGetSubImage ( xgedisplay, xgepixmap, 0, 0, IMWIDTH, IMHEIGHT, 0xFFFFFFFF,
                 ZPixmap, xim1, 0, 0 );
} /*DrawBackground*/

void SetCurrentDir ( void )
{
  int l;

  getcwd ( current_directory, MAX_PATH_LGT+1 );
  l = strlen ( current_directory );
  if ( l < MAX_PATH_SHRT )
    strcpy ( current_dir, current_directory );
  else {
    strcpy ( current_dir, "..." );
    strcpy ( &current_dir[3], &current_directory[l-MAX_PATH_SHRT+3] );
  }
} /*SetCurrentDir*/

void SetupTheLists ( void )
{
  xge_SetupFileList ( &flbox, ".", file_filter, showhidden );
  xge_SetupDirList ( &dlbox, ".", NULL, showhidden, current_directory );
  SetCurrentDir ();
  DrawBackground ();
  xge_Redraw ();
} /*SetupTheLists*/

void ChangeDir ( void )
{
  if ( !chdir ( &dlbox.itemstr[dlbox.itemind[dlbox.currentitem]] ) )
    SetupTheLists ();
} /*ChangeDir*/

void AltChangeDir ( short x )
{
  short l1, l2;

  l1 = strlen ( current_directory );
  l2 = strlen ( current_dir );
  x += l2 - l1;
  while ( x < l1 && current_directory[x] != '/' )
    x ++;
  if ( x < l1 ) {
    current_directory[x] = 0;
    if ( chdir ( current_directory ) )
      current_directory[x] = '/';
    else
      SetupTheLists ();
  }
} /*AltChangeDir*/

boolean ProcessKey ( int key )
{
  switch ( key ) {
case 'c':  case 'C':
    cursor_is_on = (!cursor_is_on) | popup_is_on;
    if ( cursor_is_on )
      xge_SetCurrentWindowCursor ( xgeCURSOR_DEFAULT );
    else
      xge_SetCurrentWindowCursor ( xgeCURSOR_INVISIBLE );
    return 1;
case 'q':  case 'Q':
    xge_done = 1;
    return 1;
case ' ':  /* space */
    ShowNextPicture ( true );
    return 1;
case 8:    /* backspace */
    ShowPreviousPicture ( true );
    return 1;
case 'f':  case 'F':
    if ( !popup_is_on ) {
      xge_AddPopup ( popup );
      popup_is_on = cursor_is_on = true;
      xge_SetCurrentWindowCursor ( xgeCURSOR_DEFAULT );
    }
    else {
      xge_RemovePopup ( true );
      popup_is_on = cursor_is_on = false;
      xge_SetCurrentWindowCursor ( xgeCURSOR_INVISIBLE );
    }
    xge_Redraw ();
    return 1;
case 46: /* laser */
case 9:  /* tab */
    ShowNextChapter ();
    return 1;
case 0:    /* my laser device or special key */
    switch ( xgeevent.xkey.keycode ) {
  case 71:  /* laser */
  case 23:  /* shift+tab */
      ShowPreviousChapter ();
      return 1;
  case 105:
  case 117:
      ShowNextPicture ( true );
      return 1;
  case 99:
  case 112:
      ShowPreviousPicture ( true );
      return 1;
  default:
/*printf ( "special key = %d\n", xgeevent.xkey.keycode );*/
      return 0;
    }
default:
    return 0;
  }
} /*ProcessKey*/

boolean ImageWinMsg ( xge_widget *er, int msg, int key, short x, short y )
{
  switch ( msg ) {
case xgemsg_KEY:
    return ProcessKey ( key );

case xgemsg_MCLICK:
    if ( (key & xgemouse_RBUTTON_DOWN) && (key & xgemouse_RBUTTON_CHANGE) ) {
      if ( !popup_is_on ) {
        xge_AddPopup ( popup );
        popup_is_on = cursor_is_on = true;
        xge_SetCurrentWindowCursor ( xgeCURSOR_DEFAULT );
      }
      else {
        xge_GetCurrentListBoxString ( &flbox, filename );
        ReadImage ( filename );
        xge_RemovePopup ( true );
        popup_is_on = cursor_is_on = false;
        xge_SetCurrentWindowCursor ( xgeCURSOR_INVISIBLE );
        InitAnimation ( 0, 0 );
      }
      return 1;
    }
    else if ( (key & xgemouse_WHEELFW_DOWN) && (key && xgemouse_WHEELFW_CHANGE) ) {
      ShowPreviousPicture ( false );
      return 1;
    }
    else if ( (key & xgemouse_WHEELBK_DOWN) && (key && xgemouse_WHEELBK_CHANGE) ) {
      ShowNextPicture ( false );
      return 1;
    }
    else
      return 0;

default:
    return false;
  }
} /*ImageWinMsg*/

int CallBack ( xge_widget *er, int msg, int key, short x, short y )
{
  switch ( msg ) {
case xgemsg_BUTTON_COMMAND:
    switch ( er->id ) {
  case btnP01LOAD:
      MakeFileList ();
      ReadSlOptFile ();
      ReadChapterFile ();
      xge_GetCurrentListBoxString ( &flbox, filename );
      ReadImage ( filename );
      xge_RemovePopup ( true );
      popup_is_on = cursor_is_on = false;
      xge_SetCurrentWindowCursor ( xgeCURSOR_INVISIBLE );
      InitAnimation ( 0, 0 );
      return 1;
  case btnP01EXIT:
      xge_done = 1;
      return 1;
  default:
      return 0;
    }

case xgemsg_SWITCH_COMMAND:
    switch ( er->id ) {
  case sw01HIDDEN:
      SetupTheLists ();
      return 1;
  default:
      return 0;
    }

case xgemsg_TEXT_WIDGET_CLICK:
    switch ( er->id ) {
  case txtDIRSTR:
      AltChangeDir ( x );
      return 1;
  default:
      return 0;
    }

case xgemsg_LISTBOX_ITEM_PICK:
    switch ( er->id ) {
  case lb01DIRLIST:
      ChangeDir ();
      return 1;
  case lb01FILELIST:
      xge_GetCurrentListBoxString ( &flbox, filename );
      ReadImage ( filename );
      xge_RemovePopup ( true );
      popup_is_on = cursor_is_on = false;
      xge_SetCurrentWindowCursor ( xgeCURSOR_INVISIBLE );
      InitAnimation ( 0, 0 );
      return 1;
  default:
      return 0;
    }

case xgemsg_KEY:
    return ProcessKey ( key );

case xgemsg_IDLE_COMMAND:
    ContinueAnimation ();
    return 1;

default:
    return 0;
  }
} /*CallBack*/

void init_win ( void )
{
  int        nplanes;
  xge_widget *w;

  pkv_InitScratchMem ( 8388608 ); /* 8 MB */
  XSetWindowBorderWidth ( xgedisplay, xgewindow, 0 );
  XMoveResizeWindow ( xgedisplay, xgewindow, 0, 0, IMWIDTH, IMHEIGHT );
  imagewin = xge_NewWidget ( 0, NULL, 0, IMWIDTH, IMHEIGHT, 0, 0,
                             NULL, NULL, ImageWinMsg, DrawImageWin );
  nplanes = XDisplayPlanes ( xgedisplay, xgescreen );
  xim1 = XCreateImage ( xgedisplay, xgevisual, nplanes, ZPixmap, 0,
                       NULL, IMWIDTH, IMHEIGHT, 8, 0 );
  xim2 = XCreateImage ( xgedisplay, xgevisual, nplanes, ZPixmap, 0,
                       NULL, IMWIDTH, IMHEIGHT, 8, 0 );
  if ( !xim1 || !xim2 )
    exit ( 1 );
  imagedata1 = (char*)malloc ( IMHEIGHT*xim1->bytes_per_line );
  imagedata2 = (char*)malloc ( IMHEIGHT*xim2->bytes_per_line );
  if ( !imagedata1 || !imagedata2 )
    exit ( 1 );
  xim1->data = imagedata1;
  xim2->data = imagedata2;
  ClearImage ();
  DrawBackground ();

    /* setup the popup menu */
  w = xge_NewTextWidget ( 0, NULL, txtDIRSTR, 380, 16, 20+10, 40+10, current_dir );
  w = xge_NewSwitch ( 0, w, sw01HIDDEN, 60, 16, 20+10, 40+30, txtHidden, &showhidden );
  w = xge_NewButton ( 0, w, btnP01LOAD, 58, 19, 20+91, 40+360-30, txtLoad );
  w = xge_NewButton ( 0, w, btnP01EXIT, 58, 19, 20+251, 40+360-30, txtExit );
  dlist = xge_NewListBox ( 0, w, lb01DIRLIST, 180, 259, 20+10, 40+60, &dlbox );
  flist = xge_NewListBox ( 0, dlist, lb01FILELIST, 180, 259, 220+10, 40+60, &flbox );
  popup = xge_NewFMenu ( 0, NULL, POPUP01, 400, 360, 20, 40, flist );
  popup->msgproc = xge_PopupMenuMsg;
  getcwd ( current_directory, MAX_PATH_LGT+1 );
  strncpy ( current_dir, current_directory, MAX_PATH_SHRT );
  xge_SetupDirList ( &dlbox, ".", NULL, false, current_directory );
  MakeFileList ();
  ReadImage ( filename );
  xge_AddPopup ( popup );
  popup_is_on = true;
} /*init_win*/

void destroy_win ( void )
{
  pkv_DestroyScratchMem ();
  XDestroyImage ( xim1 );
  XDestroyImage ( xim2 );
} /*destroy_win*/

int main ( int argc, char *argv[] )
{
  xge_Init ( argc, argv, CallBack, "tiffsh" );
  init_win ();
  xge_MessageLoop ();
  destroy_win ();
  xge_Cleanup ();
  exit ( 0 );
} /*main*/

