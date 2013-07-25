
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

#include "pkgeom.h"
#include "xgedit.h"


#define MAX_PATH_LGT     1024
#define MAX_PATH_SHRT      63
#define MAX_FILENAME_LGT   64

#define SIDE_MENU_WIDTH    60

#define WDT 8192
#define HGH 8192

#define NPARA 11

/* widget identifiers */
#define MENU0             0x100
#define M0_BTN_FILE       (MENU0+1)
#define M0_SWR0           (MENU0+2)

#define POPUP0            0x200
#define P0_BTN_OPEN       (POPUP0+1)
#define P0_BTN_SAVE       (POPUP0+2)
#define P0_BTN_EXIT       (POPUP0+3)

#define POPUP1            0x300
#define P1_BTN_OPEN       (POPUP1+1)
#define P1_BTN_CANCEL     (POPUP1+2)
#define P1_LB_DIRLIST     (POPUP1+3)
#define P1_LB_FILELIST    (POPUP1+4)

#define POPUP2            0x400
#define P2_BTN_SAVE       (POPUP2+1)
#define P2_BTN_CANCEL     (POPUP2+2)
#define P2_LB_DIRLIST     (POPUP2+3)
/*#define P2_LB_FILELIST    (POPUP2+4)*/
#define P2_TXT_SAVE_AS    (POPUP2+5)
#define P2_TXTED_FILENAME (POPUP2+6)

#define CWIND0            0x500

char txtFile[]   = "File";
char txtOpen[]   = "Open";
char txtSave[]   = "Save";
char txtSaveAs[] = "Save as";
char txtCancel[] = "Cancel";
char txtExit[]   = "Exit";

xge_2Dwind cwind;
xge_widget *menu, *popup0, *popup1, *popup2;

xge_listbox dirlist, filelist;
const char file_filter[] = "*.ifs";
const char file_ext[] = ".ifs";
char current_directory[MAX_PATH_LGT+1] = "";
char current_dir[MAX_PATH_SHRT+1] = "";
char filename[MAX_FILENAME_LGT+1] = "";
xge_string_ed filename_editor;

boolean swr[NPARA] = {true,true,true,true,false,false,false,false,false,false};
point2d rcp[3*(NPARA+1)] =
  {{-0.8,-0.8},{1.6,0.0},{0.0,1.6},
   {0.2,0.1},{0.5,0.0},{0.0,0.5},
   {0.1,0.2},{0.5,0.0},{0.0,0.5},
   {0.3,0.2},{0.5,0.0},{0.0,0.5},
   {0.2,0.3},{0.5,0.0},{0.0,0.5},
   {0.1,0.3},{0.5,0.0},{0.0,0.5},
   {0.3,0.1},{0.5,0.0},{0.0,0.5},
   {0.2,0.4},{0.5,0.0},{0.0,0.5},
   {0.1,0.4},{0.5,0.0},{0.0,0.5},
   {0.4,0.4},{0.5,0.0},{0.0,0.5},
   {0.4,0.1},{0.5,0.0},{0.0,0.5}};

trans2d tr[NPARA], trb[NPARA];

int current_point = -1;

/* ///////////////////////////////////////////////////////////////////////// */
byte *bmap;
int transx = 0, transy = 0;

void SetupBitTranslation ( void )
{
  transx = (WDT-cwind.CPos.width)/2-cwind.CPos.xmin;
  transy = (HGH-cwind.CPos.height)/2-cwind.CPos.ymin;
/*printf ( "transx = %d, transy = %d\n", transx, transy ); */
} /*SetupBitTranslation*/

boolean InitMyBitmap ( void )
{
  PKV_MALLOC ( bmap, (WDT/8)*HGH );
  return bmap != NULL;
} /*InitMyBitmap*/

void SetMyPixel ( short x, short y, char on )
{
  int  l;
  byte mask;

  x += transx;
  y += transy;
  if ( x >= 0 && x < WDT && y >= 0 && y < HGH ) {
    l = (WDT/8)*(int)y + (int)x/8;
    mask = 0x01 << (x % 8);
    if ( on )
      bmap[l] |= mask;
    else
      bmap[l] &= ~mask;
  }
} /*SetMyPixel*/

char GetMyPixel ( short x, short y )
{
  int  l;
  byte mask;

  x += transx;
  y += transy;
  if ( x >= 0 && x < WDT && y >= 0 && y < HGH ) {
    l = (WDT/8)*(int)y + (int)x/8;
    mask = 0x01 << (x % 8);
    return (bmap[l] & mask) != 0;
  }
  else
    return 0;
} /*GetMyPixel*/

/* ///////////////////////////////////////////////////////////////////////// */
void MakeTr ( point2d *p, point2d *q, trans2d *tr, trans2d *trb )
{
  double  a[9], b[9];
  trans2d ctr, ctri, qq;

  a[0] = p[0].x;  a[1] = p[0].y;  a[2] = 1.0;
  a[3] = p[1].x;  a[4] = p[1].y;  a[5] = 0.0;
  a[6] = p[2].x;  a[7] = p[2].y;  a[8] = 0.0;
  b[0] = q[0].x;  b[1] = q[0].y;  b[2] = 1.0;
  b[3] = q[1].x;  b[4] = q[1].y;  b[5] = 0.0;
  b[6] = q[2].x;  b[7] = q[2].y;  b[8] = 0.0;
  pkn_multiGaussSolveLinEqd ( 3, a, 3, 3, b );
  tr->U0.a11 = b[0];  tr->U0.a12 = b[3];  tr->U0.a13 = b[6];
  tr->U0.a21 = b[1];  tr->U0.a22 = b[4];  tr->U0.a23 = b[7];

  ctr.U0.a11 = cwind.CPos.CTr.U0.a11;
  ctr.U0.a12 = cwind.CPos.CTr.U0.a12;
  ctr.U0.a13 = cwind.CPos.CTr.U0.a14;
  ctr.U0.a21 = cwind.CPos.CTr.U0.a21;
  ctr.U0.a22 = cwind.CPos.CTr.U0.a22;
  ctr.U0.a23 = cwind.CPos.CTr.U0.a24;
  ctri = ctr;
  InvertTrans2d ( &ctri );
  CompTrans2d ( &qq, tr, &ctri );
  CompTrans2d ( trb, &ctr, &qq );
} /*MakeTr*/

void IterIFS ( void )
{
  pkv_queue *q;
  xpoint    p0, p1;
  double    a[4], b[2];
  point2d   q0, q1;
  int       i;

  memset ( bmap, 0, (WDT/8)*HGH );
  q = pkv_InitQueue ( WDT*HGH, sizeof(xpoint) );
  for ( i = 1; i < NPARA; i++ )
    if ( swr[i] ) {
      a[0] = 1.0-tr[i].U0.a11;  a[1] = -tr[i].U0.a12;
      a[2] = -tr[i].U0.a21;     a[3] = 1.0-tr[i].U0.a22;
      b[0] = tr[i].U0.a13;
      b[1] = tr[i].U0.a23;
      if ( pkn_multiGaussSolveLinEqd ( 2, a, 1, 1, b ) ) {
        CameraProjectPoint2d ( &cwind.CPos, (point2d*)&b[0], &q0 );
        p0.x = (int)(q0.x+0.5);
        p0.y = (int)(q0.y+0.5);
        if ( p0.x >= -transx && p0.x < WDT-transx &&
             p0.y >= -transy && p0.y < HGH-transy )
          break;
      }
    }
  if ( i >= NPARA )
    goto way_out;

  SetMyPixel ( p0.x, p0.y, 1 );
  pkv_QueueInsert ( q, &p0 );
  for ( ; !pkv_QueueEmpty ( q ); ) {
    pkv_QueueRemoveFirst ( q, &p0 );
    SetPoint2d ( &q0, p0.x, p0.y );
    CameraUnProjectPoint2d ( &cwind.CPos, &q0, &q0 );
    for ( i = 1; i < NPARA; i++ )
      if ( swr[i] ) {
        TransPoint2d ( &tr[i], &q0, &q1 );
        CameraProjectPoint2d ( &cwind.CPos, &q1, &q1 );
        p1.x = (int)(q1.x+0.5);
        p1.y = (int)(q1.y+0.5);
        if ( p1.x >= -transx && p1.x < WDT-transx &&
             p1.y >= -transy && p1.y < HGH-transy ) {
          if ( !GetMyPixel ( p1.x, p1.y ) ) {
            SetMyPixel ( p1.x, p1.y, 1 );
            pkv_QueueInsert ( q, &p1 );
          }
        }
      }
  }

way_out:
  PKV_FREE ( q );
} /*IterIFS*/

/* ///////////////////////////////////////////////////////////////////////// */
void DrawParallelogram ( xge_2Dwind *_2Dwin, point2d *para )
{
  point2d p[4];
  XPoint  q[5];
  int     i;

  CameraProjectPoint2d ( &_2Dwin->CPos, para, &p[0] );
  AddVector2d ( &para[0], &para[1], &p[1] );
  AddVector2d ( &para[0], &para[2], &p[3] );
  AddVector2d ( &p[1], &para[2], &p[2] );
  for ( i = 1; i < 4; i++ )
    CameraProjectPoint2d ( &_2Dwin->CPos, &p[i], &p[i] );
  for ( i = 0; i < 4; i++ ) {
    q[i].x = (short)(p[i].x+0.5);
    q[i].y = (short)(p[i].y+0.5);
  }
  q[4] = q[0];
  xgeSetForeground ( xgec_Green );
  xgeDrawLines ( 5, q );
  xgeSetForeground ( xgec_Red );
  xgeFillRectangle ( 3, 3, q[0].x-1, q[0].y-1 );
  xgeFillRectangle ( 3, 3, q[1].x-1, q[1].y-1 );
  xgeSetForeground ( xgec_Yellow );
  xgeFillRectangle ( 3, 3, q[3].x-1, q[3].y-1 );
} /*DrawParallelogram*/

void RysujOkno ( xge_widget *er, boolean onscreen )
{
  xge_2Dwind *_2Dwin;
  int        i, j;

  _2Dwin = er->data0;
  xge_DrawGeomWinBackground ( er );
  if ( _2Dwin->inside && _2Dwin->display_coord )
    xge_2DwindDrawCursorPos ( _2Dwin, xge_xx, xge_yy );

  _pkv_OutputPixels = xge_OutPixels;
  xgeSetForeground ( xgec_White );
  for ( j = er->y; j < er->y+er->h; j++ )
    for ( i = er->x; i < er->x+er->w; i++ )
      if ( GetMyPixel ( i, j ) )
        PKV_SETPIXEL ( i, j );
  PKV_FLUSH;

  for ( i = 0; i < NPARA; i++ )
    if ( swr[0] && swr[i] )
      DrawParallelogram ( _2Dwin, &rcp[3*i] );

  xge_DrawGeomWinSelectionRect ( er, &_2Dwin->selection_rect );
  xge_2DwindDrawGeomWidgets ( er );
  xge_DrawGeomWinFrame ( er, onscreen );
} /*RysujOkno*/

boolean FindNearestPoint ( xge_widget *er, short x, short y, short mindist )
{
  xge_2Dwind *_2Dwin;
  int        i, j;
  int        d, e;
  point2d    p[3];

  _2Dwin = er->data0;
  current_point = -1;
  d = mindist+1;
  for ( i = 1; i < NPARA; i++ )
    if ( swr[i] ) {
      CameraProjectPoint2d ( &_2Dwin->CPos, &rcp[3*i], &p[0] );
      AddVector2d ( &rcp[3*i], &rcp[3*i+1], &p[1] );
      AddVector2d ( &rcp[3*i], &rcp[3*i+2], &p[2] );
      for ( j = 1; j < 3; j++ )
        CameraProjectPoint2d ( &_2Dwin->CPos, &p[j], &p[j] );
      for ( j = 0; j < 3; j++ ) {
        e = (int)(fabs(x-p[j].x)+fabs(y-p[j].y));
        if ( e < d ) { d = e;  current_point = 3*i+j; }
      }
    }
  return d <= mindist;
} /*FindNearestPoint*/

void SetCPoint ( xge_widget *er, short x, short y )
{
  xge_2Dwind *_2Dwin;
  point3d    p;
  int        i;

  if ( current_point >= 0 && current_point < 3*NPARA ) {
    _2Dwin = er->data0;
    SetPoint3d ( &p, x, y, 1.0 );
    CameraUnProjectPoint3d ( &_2Dwin->CPos, &p, &p );
    if ( current_point % 3 ) {
      i = (current_point / 3)*3;
      rcp[current_point].x = p.x - rcp[i].x;
      rcp[current_point].y = p.y - rcp[i].y;
    }
    else {
      rcp[current_point].x = p.x;
      rcp[current_point].y = p.y;
    }
  }
} /*SetCPoint*/

/* ///////////////////////////////////////////////////////////////////////// */
void SetCurrentDir ( void )
{
  int l;

  getcwd ( current_directory, MAX_PATH_LGT+1 );
  l = strlen ( current_directory );
  if ( l <= MAX_PATH_SHRT )
    strcpy ( current_dir, current_directory );
  else {
    strcpy ( current_dir, "..." );
    strcpy ( &current_dir[3], &current_directory[l-MAX_PATH_SHRT+3] );
  }
} /*SetCurrentDir*/

boolean OpenTheFile ( void )
{
  FILE *f;
  int  i, j;

  if ( !xge_GetCurrentListBoxString ( &filelist, filename ) )
    return false;
  f = fopen ( filename, "r+" );
  if ( !f )
    return false;
  for ( i = 0; i < NPARA; i++ )
    if ( fscanf ( f, "%lf %lf %lf %lf %lf %lf %d",
           &rcp[3*i].x, &rcp[3*i].y, &rcp[3*i+1].x, &rcp[3*i+1].y,
           &rcp[3*i+2].x, &rcp[3*i+2].y, &j ) == 7 )
      swr[i] = j;
    else
      break;
  fclose ( f );
  for ( i = 1; i < NPARA; i++ )
    MakeTr ( rcp, &rcp[3*i], &tr[i], &trb[i] );
  return true;
} /*OpenTheFile*/

boolean SaveTheFile ( void )
{
  FILE *f;
  int  i;

  f = fopen ( filename, "w+" );
  if ( !f )
    return false;
  for ( i = 0; i < NPARA; i++ )
    fprintf ( f, "%f %f %f %f %f %f %d\n",
              rcp[3*i].x, rcp[3*i].y, rcp[3*i+1].x, rcp[3*i+1].y,
              rcp[3*i+2].x, rcp[3*i+2].y, swr[i] );
  fclose ( f );
  return true;
} /*SaveTheFile*/

boolean ChangeDir ( xge_widget *popup,
                    xge_listbox *dirlist, xge_listbox *filelist )
{
  int l;

  if ( !chdir ( &dirlist->itemstr[dirlist->itemind[dirlist->currentitem]] ) ) {
    if ( filelist )
      xge_SetupFileList ( filelist, ".", file_filter );
    xge_SetupDirList ( dirlist, ".", NULL, current_directory );
    getcwd ( current_directory, MAX_PATH_LGT+1 );
    l = strlen ( current_directory );
    if ( l <= MAX_PATH_SHRT )
      strcpy ( current_dir, current_directory );
    else {
      strcpy ( current_dir, "..." );
      strcpy ( &current_dir[3], &current_directory[l-MAX_PATH_SHRT+3] );
    }
    xge_SetClipping ( popup );
    popup->redraw ( popup, true );
    return true;
  }
  else
    return false;
} /*ChangeDir*/

/* ///////////////////////////////////////////////////////////////////////// */
int CallBack ( xge_widget *er, int msg, int key, short x, short y )
{
  int i;

  if ( er ) {
/*    printf ( "%d, %d, %d, %d, %d\n", er->id, msg, key, x, y );*/
    switch ( msg ) {
case xgemsg_BUTTON_COMMAND:
      switch ( er->id ) {
  case M0_BTN_FILE:
        xge_AddPopup ( popup0 );
        xge_GrabFocus ( popup0, true );
        return true;

  case P0_BTN_OPEN:
        xge_RemovePopup ( true );
        SetCurrentDir ();
        xge_SetupFileList ( &filelist, ".", file_filter );
        xge_SetupDirList ( &dirlist, ".", NULL, NULL );
        xge_AddPopup ( popup1 );
        xge_GrabFocus ( popup1, true );
        return true;

  case P0_BTN_SAVE:
        xge_RemovePopup ( true );
        SetCurrentDir ();
        xge_SetupDirList ( &dirlist, ".", NULL, NULL );
        xge_AddPopup ( popup2 );
        xge_GrabFocus ( popup2, true );
        return true;

  case P0_BTN_EXIT:
        xge_done = 1;
        return true;

  case P1_BTN_OPEN:
        if ( OpenTheFile () ) {
          IterIFS ();
          xge_RemovePopup ( true );
          xge_ClearListBox ( &filelist );
          xge_ClearListBox ( &dirlist ); 
        }
        return true;

  case P1_BTN_CANCEL:  /* cancel opening file */
        xge_RemovePopup ( true );
        xge_ClearListBox ( &filelist );
        xge_ClearListBox ( &dirlist ); 
        return true;

  case P2_BTN_SAVE:
        if ( SaveTheFile () ) {
          xge_RemovePopup ( true );
          xge_ClearListBox ( &dirlist ); 
        }
        return true;

  case P2_BTN_CANCEL: /* cancel saving file */
        xge_RemovePopup ( true );
        xge_ClearListBox ( &dirlist );
        return true;

  default:
        break;
      }
      break;

case xgemsg_SWITCH_COMMAND:
      if ( er->id >= M0_SWR0 && er->id < M0_SWR0+NPARA ) {
        IterIFS ();
        xge_Redraw ();
      }
      break;

case xgemsg_LISTBOX_ITEM_PICK:
      switch ( er->id ) {
  case P1_LB_DIRLIST:
        ChangeDir ( popup1, &dirlist, &filelist );
        break;
  case P1_LB_FILELIST:
        if ( OpenTheFile () ) {
          IterIFS ();
          xge_RemovePopup ( true );
          xge_ClearListBox ( &filelist );
        }
        break;
  case P2_LB_DIRLIST:
        ChangeDir ( popup2, &dirlist, NULL );
        break;
      }
      return true;

case xgemsg_TEXT_EDIT_VERIFY:
      switch ( er->id ) {
  case P2_TXTED_FILENAME:
        return filename[0] != 0;
  default:
        break;
      }
      break;

case xgemsg_2DWIN_PROJCHANGE:
      return true;

case xgemsg_2DWIN_PICK_POINT:
      return FindNearestPoint ( er, x, y, xge_MINDIST );

case xgemsg_2DWIN_MOVE_POINT:
      SetCPoint ( er, x, y );
      MakeTr ( rcp, &rcp[3*(current_point/3)],
               &tr[current_point/3], &trb[current_point/3] );
      IterIFS ();
      xge_Redraw ();
      return true;

default:
      break;
    }
  }
  else {
/*    printf ( "%d, %d, %d, %d\n", msg, key, x, y );*/
    switch ( msg ) {
case xgemsg_RESIZE:
      cwind.er->w = (short)(x-SIDE_MENU_WIDTH);  cwind.er->h = (short)y;
      menu->msgproc ( menu, xgemsg_RESIZE, 0, SIDE_MENU_WIDTH, (short)y );
      xge_2DwindInitProjection ( &cwind, cwind.RefBBox.x0, cwind.RefBBox.x1,
                                 cwind.RefBBox.y0, cwind.RefBBox.y1 );
      SetupBitTranslation ();
      for ( i = 1; i < NPARA; i++ )
        MakeTr ( rcp, &rcp[3*i], &tr[i], &trb[i] );
      IterIFS ();
      xge_Redraw ();
      return true;

case xgemsg_KEY:
      switch ( key ) {
  case 'Q': case 'q':
        xge_done = 1;
        break;
  default:
        break;
      }
      break;

default:
      break;
    }
  }
  return 0;
} /*CallBack*/

/* ///////////////////////////////////////////////////////////////////////// */
void init_edwin ( void )
{
  xge_widget *rp, *edr;
  int        i;

  pkv_InitScratchMem ( 64*1048576 );

    /* setup the first popup */
  rp = xge_NewButton ( 0, NULL, P0_BTN_OPEN, 60, 18, 2, 21, txtOpen );
  rp = xge_NewButton ( 0, rp, P0_BTN_SAVE, 60, 18, 2, 41, txtSave );
  rp = xge_NewButton ( 0, rp, P0_BTN_EXIT, 60, 18, 2, 61, txtExit );
  popup0 = xge_NewFMenu ( 0, NULL, POPUP0, 64, 62, 0, 19, rp );
  popup0->msgproc = xge_PopupMenuMsg;

    /* setup the Open file menu */
  rp = xge_NewTextWidget ( 0, NULL, 0, 380, 16, 20+10, 40+10, current_dir ); 
  rp = xge_NewButton ( 0, rp, P1_BTN_OPEN, 76, 19, 20+82, 40+180-30, txtOpen );
  rp = xge_NewButton ( 0, rp, P1_BTN_CANCEL, 76, 19, 20+242, 40+180-30, txtCancel );
  rp = xge_NewListBox ( 0, rp, P1_LB_DIRLIST, 180, 99, 20+10, 40+40, &dirlist );   
  rp = xge_NewListBox ( 0, rp, P1_LB_FILELIST, 180, 99, 220+10, 40+40, &filelist );
  popup1 = xge_NewFMenu ( 0, NULL, POPUP1, 400, 180, 20, 40, rp );

    /* setup the Save as file menu */
  rp = xge_NewTextWidget ( 0, NULL, 0, 380, 16, 20+10, 40+10, current_dir ); 
  rp = xge_NewButton ( 0, rp, P2_BTN_SAVE, 76, 19, 20+82, 40+180-30, txtSave );
  rp = xge_NewButton ( 0, rp, P2_BTN_CANCEL, 76, 19, 20+242, 40+180-30, txtCancel );
  rp = xge_NewListBox ( 0, rp, P2_LB_DIRLIST, 180, 99, 20+10, 40+40, &dirlist );
  rp = xge_NewTextWidget ( 0, rp, P2_TXT_SAVE_AS, 120, 16, 220+10, 40+38, txtSaveAs );
  rp = xge_NewStringEd ( 0, rp, P2_TXTED_FILENAME, 180, 19, 220+10, 40+56,
                         MAX_FILENAME_LGT, filename, &filename_editor );
  popup2 = xge_NewFMenu ( 0, NULL, POPUP2, 400, 180, 20, 40, rp );

    /* setup the main menu and fractal editing window */
  rp = xge_NewButton ( 0, NULL, M0_BTN_FILE, 60, 18, 0, 0, txtFile );
  for ( i = 0; i < NPARA; i++ )
    rp = xge_NewSwitch ( 0, rp, M0_SWR0+i, 16, 16, 0, 20+20*i, NULL, &swr[i] );
  menu = xge_NewMenu ( 0, NULL, 1, SIDE_MENU_WIDTH, xge_HEIGHT, 0, 0, rp );
  cwind.er = edr = xge_New2Dwind
                       ( 0, menu, CWIND0, xge_WIDTH-SIDE_MENU_WIDTH,
                         xge_HEIGHT, SIDE_MENU_WIDTH, 0,
                         &cwind, RysujOkno );
  xge_SetWinEdRect ( edr );

    /* prepare the canvas */
  SetupBitTranslation ();
  _pkv_InitPixelBuffer ();
  if ( !InitMyBitmap () )
    exit ( 1 );
  for ( i = 1; i < NPARA; i++ )
    MakeTr ( rcp, &rcp[3*i], &tr[i], &trb[i] );
  IterIFS ();

  xge_Redraw ();
} /*init_edwin*/

void destroy_edwin ( void )
{
  pkv_DestroyScratchMem ();
} /*destroy_edwin*/

int main ( int argc, char *argv[] )
{
  xge_Init ( argc, argv, CallBack, NULL );
  init_edwin ();
  xge_MessageLoop ();
  destroy_edwin ();  
  xge_Cleanup ();  
  exit ( 0 );
} /*main*/   

