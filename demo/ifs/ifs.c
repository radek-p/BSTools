
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

#include "ifs.h"

char txtFile[]      = "File";
char txtSelect[]    = "Select";
char txtTranslate[] = "Translate";
char txtScale[]     = "Scale";
char txtRotate[]    = "Rotate";
char txtShear[]     = "Shear";
char txtPanZoom[]   = "Pan/Zoom";
char txtCoord[]     = "Coord.";
char txtAA[]        = "aa";
char txtOpen[]      = "Open";
char txtSave[]      = "Save";
char txtSaveAs[]    = "Save as";
char txtCancel[]    = "Cancel";
char txtExit[]      = "Exit";

xge_2Dwind cwind;
xge_widget *menu, *popup0, *popup1, *popup2;

xge_listbox dirlist, filelist;
const char file_filter[] = "*.ifs";
const char file_ext[] = ".ifs";
char current_directory[MAX_PATH_LGT+1] = "";
char current_dir[MAX_PATH_SHRT+1] = "";
char filename[MAX_FILENAME_LGT+1] = "";
xge_string_ed filename_editor;

boolean antialias = false;

boolean swr[NPARA] = {true,true,true,true,false,false,false,false,false,false};
point2d rcp[3*NPARA] =
  {{-0.8,-0.8},{1.6,0.0},{0.0,1.6},
   {-0.8,-0.8},{0.8,0.0},{0.0,0.8},
   {-0.8,0.0},{0.0,0.8},{0.8,0.0},
   {0.0,-0.8},{0.8,0.0},{0.0,0.8},
   {0.2,0.3},{0.5,0.0},{0.0,0.5},
   {0.1,0.3},{0.5,0.0},{0.0,0.5},
   {0.3,0.1},{0.5,0.0},{0.0,0.5},
   {0.2,0.4},{0.5,0.0},{0.0,0.5},
   {0.1,0.4},{0.5,0.0},{0.0,0.5},
   {0.4,0.4},{0.5,0.0},{0.0,0.5},
   {0.4,0.4},{0.5,0.2},{0.0,0.5},
   {0.4,0.1},{0.5,0.0},{0.0,0.5}};
point2d rcps[3*NPARA];
boolean mkrcp[3*NPARA];

trans2d tr[NPARA], trb[NPARA], traa[NPARA];

int current_point = -1;

/* ///////////////////////////////////////////////////////////////////////// */
byte           *bmap;
unsigned short *bmapaa;

int transx = 0, transy = 0,
    ximin = 0, ximax = WDT, etamin = 0, etamax = HGH;
int transxaa = 0, transyaa = 0,
    ximinaa = 0, ximaxaa = WDT, etaminaa = 0, etamaxaa = HGH;

void SetupBitTranslation ( void )
{
        /* translation and pixel box */
  transx = (WDT-cwind.CPos.width)/2-cwind.CPos.xmin;
  transy = (HGH-cwind.CPos.height)/2-cwind.CPos.ymin;
  ximin = -transx;
  ximax = WDT-transx;
  etamin = -transy;
  etamax = HGH-transy;
        /* the same for the antialiased version */
  transxaa = (WDT-4*cwind.CPos.width)/2-4*cwind.CPos.xmin;
  transyaa = (HGH-4*cwind.CPos.height)/2-4*cwind.CPos.ymin;
  ximinaa = -transxaa;
  ximaxaa = WDT-transxaa;
  etaminaa = -transyaa;
  etamaxaa = HGH-transyaa;
} /*SetupBitTranslation*/

boolean InitMyBitmap ( void )
{
  PKV_MALLOC ( bmap, (WDT/8)*HGH );
  bmapaa = (unsigned short*)bmap;
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

void SetMySubpixelAA ( short x, short y, char on )
{
  unsigned short mask;
  short          xa, xb, ya, yb;
  int            l;

  x += transxaa;
  y += transyaa;
  if ( x >= 0 && x < WDT && y >= 0 && y < HGH ) {
    xa = x >> 2;  xb = x & 0x0003;
    ya = y >> 2;  yb = y & 0x0003;
    l = (WDT/4)*(int)ya + (int)xa;
    mask = 0x0001 << (4*yb+xb);
    if ( on )
      bmapaa[l] |= mask;
    else
      bmapaa[l] &= ~mask;
  }
} /*SetMySubpixelAA*/

char GetMySubpixelAA ( short x, short y )
{
  unsigned short mask;
  short          xa, xb, ya, yb;
  int            l;

  x += transxaa;
  y += transyaa;
  if ( x >= 0 && x < WDT && y >= 0 && y < HGH ) {
    xa = x >> 2;  xb = x & 0x0003;
    ya = y >> 2;  yb = y & 0x0003;
    l = (WDT/4)*(int)ya + (int)xa;
    mask = 0x0001 << (4*yb+xb);
    return (bmapaa[l] & mask) != 0;
  }
  else
    return 0;
} /*GetMySubpixelAA*/

char GetMyPixelAA ( short x, short y )
{
  unsigned short pixel;
  int            l, i;
  char           nbits;

  x += transxaa/4;
  y += transyaa/4;
  if ( x >= 0 && x < WDT/4 && y >= 0 && y < HGH/4 ) {
    l = (WDT/4)*(int)y + (int)x;
    pixel = bmapaa[l];
    for ( i = 0, nbits = 0;  i < 16;  i++, pixel >>= 1 )
      if ( pixel & 0x0001 )
        nbits ++;
    return nbits;
  }
  else
    return 0;
} /*GetMyPixelAA*/

/* ///////////////////////////////////////////////////////////////////////// */
void MakeTr ( point2d *p, point2d *q, trans2d *tr,
              trans2d *trb, trans2d *traa )
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
  ctr.U0.a21 = -cwind.CPos.CTr.U0.a21;
  ctr.U0.a22 = -cwind.CPos.CTr.U0.a22;
  ctr.U0.a23 = 2.0*cwind.CPos.ymin+cwind.CPos.height-cwind.CPos.CTr.U0.a24;
  ctr.U1.detsgn = -cwind.CPos.CTr.U1.detsgn;
  ctri = ctr;
  InvertTrans2d ( &ctri );
  CompTrans2d ( &qq, tr, &ctri );
  CompTrans2d ( trb, &ctr, &qq );
  *traa = *trb;
  traa->U0.a13 *= 4.0;
  traa->U0.a23 *= 4.0;
} /*MakeTr*/

void MakeAllTr ( void )
{
  int i;

  for ( i = 1; i < NPARA; i++ )
    MakeTr ( rcp, &rcp[3*i], &tr[i], &trb[i], &traa[i] );
} /*MakeAllTr*/

void IterIFS ( void )
{
  pkv_queue *q;
  xpoint    p0, p1;
  double    a[4], b[2];
  point2d   q0, q1;
  int       i;

  memset ( bmap, 0, (WDT/8)*HGH );
  q = pkv_InitQueue ( WDT*HGH, sizeof(xpoint) );
  if ( antialias ) {
    for ( i = 1; i < NPARA; i++ )
      if ( swr[i] ) {
        a[0] = 1.0-traa[i].U0.a11;  a[1] = -traa[i].U0.a12;
        a[2] = -traa[i].U0.a21;     a[3] = 1.0-traa[i].U0.a22;
        b[0] = traa[i].U0.a13;
        b[1] = traa[i].U0.a23;
        if ( pkn_multiGaussSolveLinEqd ( 2, a, 1, 1, b ) ) {
          p0.x = (int)(b[0]+0.5);
          p0.y = (int)(b[1]+0.5);
          if ( p0.x >= ximinaa && p0.x < ximaxaa &&
               p0.y >= etaminaa && p0.y < etamaxaa )
            break;
        }
      }
    if ( i >= NPARA )
      goto way_out;

    SetMySubpixelAA ( p0.x, p0.y, 1 );
    pkv_QueueInsert ( q, &p0 );
    for ( ; !pkv_QueueEmpty ( q ); ) {
      pkv_QueueRemoveFirst ( q, &p0 );
      SetPoint2d ( &q0, p0.x, p0.y );
      for ( i = 1; i < NPARA; i++ )
        if ( swr[i] ) {
          TransPoint2d ( &traa[i], &q0, &q1 );
          p1.x = (int)(q1.x+0.5);
          p1.y = (int)(q1.y+0.5);
          if ( p1.x >= ximin && p1.x < ximax &&
               p1.y >= etamin && p1.y < etamax ) {
            if ( !GetMySubpixelAA ( p1.x, p1.y ) ) {
              SetMySubpixelAA ( p1.x, p1.y, 1 );
              pkv_QueueInsert ( q, &p1 );
            }
          }
        }
    }
  }
  else {  /* no antialiasing */
    for ( i = 1; i < NPARA; i++ )
      if ( swr[i] ) {
        a[0] = 1.0-trb[i].U0.a11;  a[1] = -trb[i].U0.a12;
        a[2] = -trb[i].U0.a21;     a[3] = 1.0-trb[i].U0.a22;
        b[0] = trb[i].U0.a13;
        b[1] = trb[i].U0.a23;
        if ( pkn_multiGaussSolveLinEqd ( 2, a, 1, 1, b ) ) {
          p0.x = (int)(b[0]+0.5);
          p0.y = (int)(b[1]+0.5);
          if ( p0.x >= ximin && p0.x < ximax &&
               p0.y >= etamin && p0.y < etamax )
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
      for ( i = 1; i < NPARA; i++ )
        if ( swr[i] ) {
          TransPoint2d ( &trb[i], &q0, &q1 );
          p1.x = (int)(q1.x+0.5);
          p1.y = (int)(q1.y+0.5);
          if ( p1.x >= ximin && p1.x < ximax &&
               p1.y >= etamin && p1.y < etamax ) {
            if ( !GetMyPixel ( p1.x, p1.y ) ) {
              SetMyPixel ( p1.x, p1.y, 1 );
              pkv_QueueInsert ( q, &p1 );
            }
          }
        }
    }
  }

way_out:
  PKV_FREE ( q );
} /*IterIFS*/

/* ///////////////////////////////////////////////////////////////////////// */
void DrawParallelogram ( xge_2Dwind *_2Dwin, point2d *para,
                         xgecolour_int c, boolean *cpmk )
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
  xgeSetForeground ( c );
  xgeDrawLines ( 5, q );
  if ( cpmk[0] )
    xgeSetForeground ( xgec_Red );
  else
    xgeSetForeground ( xgec_DarkOrange4 );
  xgeFillRectangle ( 3, 3, q[0].x-1, q[0].y-1 );
  if ( cpmk[1] )
    xgeSetForeground ( xgec_Red );
  else
    xgeSetForeground ( xgec_DarkOrange4 );
  xgeFillRectangle ( 3, 3, q[1].x-1, q[1].y-1 );
  if ( cpmk[2] )
    xgeSetForeground ( xgec_Yellow );
  else
    xgeSetForeground ( xgec_Goldenrod4 );
  xgeFillRectangle ( 3, 3, q[3].x-1, q[3].y-1 );
} /*DrawParallelogram*/

void DrawBitmap ( xge_widget *er )
{
  int i, j;

  _pkv_OutputPixels = xge_OutPixels;
  xgeSetForeground ( xgec_White );
  for ( j = er->y; j < er->y+er->h; j++ )
    for ( i = er->x; i < er->x+er->w; i++ )
      if ( GetMyPixel ( i, j ) )
        PKV_SETPIXEL ( i, j );
  PKV_FLUSH;
} /*DrawBitmap*/

void DrawBitmapAA ( xge_widget *er )
{
  xgecolour_int cc[17];
  int           i, j;
  byte          r0, g0, b0, r1, g1, b1, c0, c1;

      /* make the palette */
  xge_GetPixelColour ( xgec_Blue5, &r0, &g0, &b0 );
  xge_GetPixelColour ( xgec_White, &r1, &g1, &b1 );
  for ( i = 0; i <= 16; i++ )
    cc[i] = xge_PixelColour (
                (byte)(((16.0-i)*r0 + (double)i*r1)/16.0 + 0.5),
                (byte)(((16.0-i)*g0 + (double)i*g1)/16.0 + 0.5),
                (byte)(((16.0-i)*b0 + (double)i*b1)/16.0 + 0.5) );
      /* draw the bitmap */
  _pkv_OutputPixels = xge_OutPixels;
  c0 = 0;
  for ( j = er->y; j < er->y+er->h; j++ )
    for ( i = er->x; i < er->x+er->w; i++ ) {
      c1 = GetMyPixelAA ( i, j );
      if ( c1 ) {
        if ( c1 != c0 ) {
          PKV_FLUSH;
          xgeSetForeground ( cc[c1] );
          c0 = c1;
        }
        PKV_SETPIXEL ( i, j );
      }
    }
  PKV_FLUSH;
} /*DrawBitmapAA*/

void RysujOkno ( xge_widget *er, boolean onscreen )
{
  xge_2Dwind *_2Dwin;
  int        i;

  _2Dwin = er->data0;
  xge_DrawGeomWinBackground ( er );
  if ( _2Dwin->inside && _2Dwin->display_coord )
    xge_2DwindDrawCursorPos ( _2Dwin, xge_xx, xge_yy );

  if ( antialias )
    DrawBitmapAA ( er );
  else
    DrawBitmap ( er );

  if ( swr[0] )
    DrawParallelogram ( _2Dwin, &rcp[0], xgec_Green4, mkrcp );
  for ( i = 1; i < NPARA; i++ )
    if ( swr[0] && swr[i] )
      DrawParallelogram ( _2Dwin, &rcp[3*i], xgec_Green, &mkrcp[3*i] );

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
void SelectPoints ( xge_widget *er, short x0, short x1, short y0, short y1,
                    boolean select )
{
  xge_2Dwind *_2Dwin;
  int        i, k;
  point2d    p, q;

  _2Dwin = er->data0;
  if ( x1-x0 >= 2 || y1-y0 >= 2 ) {
    for ( i = 1, k = 3;  i < NPARA;  i++, k += 3 ) {
      CameraProjectPoint2d ( &_2Dwin->CPos, &rcp[k], &q );
      if ( q.x >= x0 && q.x <= x1 && q.y >= y0 && q.y <= y1 )
        mkrcp[k] = select;
      AddVector2d ( &rcp[k], &rcp[k+1], &p );
      CameraProjectPoint2d ( &_2Dwin->CPos, &p, &q );
      if ( q.x >= x0 && q.x <= x1 && q.y >= y0 && q.y <= y1 )
        mkrcp[k+1] = select;
      AddVector2d ( &rcp[k], &rcp[k+2], &p );
      CameraProjectPoint2d ( &_2Dwin->CPos, &p, &q );
      if ( q.x >= x0 && q.x <= x1 && q.y >= y0 && q.y <= y1 )
        mkrcp[k+2] = select;
    }
  }
  else {
    if ( FindNearestPoint ( er, x0, y0, xge_MINDIST ) )
      mkrcp[current_point] = select;
  }
} /*SelectPoints*/

void SavePoints ( void )
{
  int i;

  for ( i = 0; i < 3*NPARA; i++ )
    rcps[i] = rcp[i];
/*  memcpy ( &rcps[0].x, &rcp[0].x, 3*NPARA*sizeof(point2d) ); */
} /*SavePoints*/

void TransformPoints ( xge_widget *er )
{
  int        i, k;
  xge_2Dwind *_2Dwin;
  trans2d    *t;

  _2Dwin = er->data0;
  t = &_2Dwin->gwtrans;
  for ( i = 1, k = 3;  i < NPARA;  i++, k += 3 )
    if ( swr[i] ) {
      if ( mkrcp[k] )   TransPoint2d ( t, &rcps[k], &rcp[k] );
      if ( mkrcp[k+1] ) TransVector2d ( t, &rcps[k+1], &rcp[k+1] );
      if ( mkrcp[k+2] ) TransVector2d ( t, &rcps[k+2], &rcp[k+2] );
      if ( mkrcp[k] || mkrcp[k+1] || mkrcp[k+2] )
        MakeTr ( rcp, &rcp[3*i], &tr[i], &trb[i], &traa[i] );
    }
} /*TransformPoints*/

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
  MakeAllTr ();
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
      xge_SetupFileList ( filelist, ".", file_filter, false );
    xge_SetupDirList ( dirlist, ".", NULL, false, current_directory );
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
        xge_SetupFileList ( &filelist, ".", file_filter, false );
        xge_SetupDirList ( &dirlist, ".", NULL, false, NULL );
        xge_AddPopup ( popup1 );
        xge_GrabFocus ( popup1, true );
        return true;

  case P0_BTN_SAVE:
        xge_RemovePopup ( true );
        SetCurrentDir ();
        xge_SetupDirList ( &dirlist, ".", NULL, false, NULL );
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
      switch ( er->id ) {
  case M0_SW_SELECT:
        if ( cwind.selecting_mode )
          xge_2DwindEnableGeomWidget ( &cwind, xge_2DWIN_SELECTING_TOOL );
        else
          xge_2DwindEnableGeomWidget ( &cwind, xge_2DWIN_NO_TOOL );
        xge_Redraw ();
        break;
  case M0_SW_TRANSLATE:
        if ( cwind.moving_tool )
          xge_2DwindEnableGeomWidget ( &cwind, xge_2DWIN_MOVING_TOOL );
        else
          xge_2DwindEnableGeomWidget ( &cwind, xge_2DWIN_NO_TOOL );
        xge_Redraw ();
        break;
  case M0_SW_SCALE:
        if ( cwind.scaling_tool )
          xge_2DwindEnableGeomWidget ( &cwind, xge_2DWIN_SCALING_TOOL );
        else
          xge_2DwindEnableGeomWidget ( &cwind, xge_2DWIN_NO_TOOL );
        xge_Redraw ();
        break;
  case M0_SW_ROTATE:
        if ( cwind.rotating_tool )
          xge_2DwindEnableGeomWidget ( &cwind, xge_2DWIN_ROTATING_TOOL );
        else
          xge_2DwindEnableGeomWidget ( &cwind, xge_2DWIN_NO_TOOL );
        xge_Redraw ();
        break;
  case M0_SW_SHEAR:
        if ( cwind.shear_tool )
          xge_2DwindEnableGeomWidget ( &cwind, xge_2DWIN_SHEAR_TOOL );
        else
          xge_2DwindEnableGeomWidget ( &cwind, xge_2DWIN_NO_TOOL );
        xge_Redraw ();
        break;
  case M0_SW_PANZOOM:
        if ( cwind.panning )
          xge_2DwindEnableGeomWidget ( &cwind, xge_2DWIN_PANNING_TOOL );
        else
          xge_2DwindEnableGeomWidget ( &cwind, xge_2DWIN_NO_TOOL );
        xge_Redraw ();
        break;
  case M0_SW_COORD:
        break;
  case M0_SW_ANTIALIAS:
        IterIFS ();
        xge_Redraw ();
        break;
  default:
        if ( er->id >= M0_SWR0 && er->id < M0_SWR0+NPARA ) {
          IterIFS ();
          xge_Redraw ();
        }
        break;
      }
      return true;

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
      MakeAllTr ();
      IterIFS ();
      xge_Redraw ();
      return true;

case xgemsg_2DWIN_PICK_POINT:
      return FindNearestPoint ( er, x, y, xge_MINDIST );

case xgemsg_2DWIN_MOVE_POINT:
      SetCPoint ( er, x, y );
      MakeTr ( rcp, &rcp[3*(current_point/3)],
               &tr[current_point/3], &trb[current_point/3],
               &traa[current_point/3] );
      IterIFS ();
      xge_Redraw ();
      return true;

case xgemsg_2DWIN_SELECT_POINTS:
      SelectPoints ( er, cwind.selection_rect.x0, cwind.selection_rect.x1,
                     cwind.selection_rect.y0, cwind.selection_rect.y1, true );
      return true;

case xgemsg_2DWIN_UNSELECT_POINTS:
      SelectPoints ( er, cwind.selection_rect.x0, cwind.selection_rect.x1,
                     cwind.selection_rect.y0, cwind.selection_rect.y1, false );
      return true;

case xgemsg_2DWIN_SAVE_POINTS:
      SavePoints ();
      return true;

case xgemsg_2DWIN_TRANSFORM_POINTS:
      TransformPoints ( er );
      IterIFS ();
      return true;

case xgemsg_2DWIN_FIND_REFBBOX:
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
      MakeAllTr ();
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
  int        i, j, k;

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
  rp = xge_NewButton ( 0, NULL, M0_BTN_FILE, 57, 18, 0, 0, txtFile );
  for ( i = k = 0;  k < NPARA;  i++ )
    for ( j = 0; j < 3 && k < NPARA; j++, k++ )
      rp = xge_NewSwitch ( 0, rp, M0_SWR0+k, 16, 16, 20*j, 20+20*i, NULL, &swr[k] );
    /* these 3 switches are to be moved after window resizing */   
  rp = xge_NewSwitch ( 0, rp, M0_SW_SELECT, 109, 16, 0, xge_HEIGHT-156,
                           txtSelect, &cwind.selecting_mode );
  xge_SetWidgetPositioning ( rp, 2, 0, -156 );
  rp = xge_NewSwitch ( 0, rp, M0_SW_TRANSLATE, 109, 16, 0, xge_HEIGHT-136,
                           txtTranslate, &cwind.moving_tool );
  xge_SetWidgetPositioning ( rp, 2, 0, -136 );
  rp = xge_NewSwitch ( 0, rp, M0_SW_SCALE, 109, 16, 0, xge_HEIGHT-116,
                           txtScale, &cwind.scaling_tool );
  xge_SetWidgetPositioning ( rp, 2, 0, -116 );
  rp = xge_NewSwitch ( 0, rp, M0_SW_ROTATE, 109, 16, 0, xge_HEIGHT-96,
                           txtRotate, &cwind.rotating_tool );
  xge_SetWidgetPositioning ( rp, 2, 0, -96 );
  rp = xge_NewSwitch ( 0, rp, M0_SW_SHEAR, 109, 16, 0, xge_HEIGHT-76,
                           txtShear, &cwind.shear_tool );
  xge_SetWidgetPositioning ( rp, 2, 0, -76 );
  rp = xge_NewSwitch ( 0, rp, M0_SW_PANZOOM, 109, 16, 0, xge_HEIGHT-56,
                           txtPanZoom, &cwind.panning );
  xge_SetWidgetPositioning ( rp, 2, 0, -56 );
  rp = xge_NewSwitch ( 0, rp, M0_SW_COORD, 109, 16, 0, xge_HEIGHT-36,
                           txtCoord, &cwind.display_coord );
  xge_SetWidgetPositioning ( rp, 2, 0, -36 );
  rp = xge_NewSwitch ( 0, rp, M0_SW_ANTIALIAS, 109, 16, 0, xge_HEIGHT-16,
                           txtAA, &antialias );
  xge_SetWidgetPositioning ( rp, 2, 0, -16 );
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
  MakeAllTr ();
  IterIFS ();
  memset ( mkrcp, 0, 3*NPARA );

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

