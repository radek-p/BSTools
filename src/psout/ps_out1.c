
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2005, 2012                            */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

#include <math.h>
#include <memory.h>
#include <stdio.h>

#include "pkvaria.h"
#include "pkgeom.h"
#include "psout.h"

#include "psprivate.h"

short ps_dec_digits = 3;


static FILE *_ps_fstack[PS_FILE_STACK_LENGTH];
static int _ps_fsp = 0;

FILE *_ps_f;
unsigned short _ps_dpi;   /* dots per inch */
float _ps_cgray;   /* current grey level */
float _ps_cred, _ps_cgreen, _ps_cblue;  /* current colour components */
float _ps_cwidth;   /* current line width */
long _ps_written;

short _ps_bmp_y;   /* lines to the end of bitmap */
unsigned short _ps_bmp_w;  /* size of the bitmap */
byte _ps_bmp_p;   /* bitmap pixel size */


static float xmin, ymin, xmax, ymax;
static boolean noext;
static char f_NAME[64];

static int bbx1, bby1, bbx2, bby2;
static boolean bbox = false;

void _ps_ExtendRect ( float x, float y )
{
  if ( noext ) {
    xmin = x;
    xmax = x;
    ymin = y;
    ymax = y;
    noext = false;
    return;
  }
  if ( x < xmin )
    xmin = x;
  else if ( x > xmax )
    xmax = x;
  if ( y < ymin )
    ymin = y;
  else if ( y > ymax )
    ymax = y;
} /*_ps_ExtendRect*/

void _ps_OutProc ( ps_proc p )
{
  if ( ((1L << ((long)p)) & _ps_written) != 0 )
    return;
  switch ( p ) {

  case _newpath:
    fprintf ( _ps_f, "/_np { newpath } bdef\n" );
    break;

  case _drawline:
    _ps_OutProc(_newpath);
    fprintf ( _ps_f, "/_dl { _np moveto lineto stroke } bdef\n" );
    break;

  case _dcircle:
    _ps_OutProc(_newpath);
    fprintf ( _ps_f, "/_dc { _np 0 360 arc stroke } bdef\n" );
    break;

  case _fcircle:
    _ps_OutProc(_newpath);
    fprintf ( _ps_f, "/_fc { _np 0 360 arc fill } bdef\n" );
    break;

  case _mcircle:
    _ps_OutProc ( _newpath );
    fprintf ( _ps_f, "/_mc {\n" );
    fprintf ( _ps_f, "  /y exch def\n" );
    fprintf ( _ps_f, "  /x exch def\n" );
    fprintf ( _ps_f, "  currentrgbcolor" );
    fprintf ( _ps_f, "  _np x y %4.1f 0 360 arc fill\n", _ps_dpi * 0.02 );
    fprintf ( _ps_f, "  1 setgray\n" );
    fprintf ( _ps_f, "  _np x y %4.1f 0 360 arc fill\n", _ps_dpi * 0.01333 );
    fprintf ( _ps_f, "  setrgbcolor\n" );
    fprintf ( _ps_f, "} bdef\n" );
    break;

  case _drawrect:
    _ps_OutProc ( _newpath );
    fprintf ( _ps_f, "/_dr {\n" );
    fprintf ( _ps_f, "  /y exch def\n" );
    fprintf ( _ps_f, "  /x exch def\n" );
    fprintf ( _ps_f, "  /h exch def\n" );
    fprintf ( _ps_f, "  /w exch def\n" );
    fprintf ( _ps_f, "  _np\n" );
    fprintf ( _ps_f, "  x y moveto\n" );
    fprintf ( _ps_f, "  x w add y lineto\n" );
    fprintf ( _ps_f, "  x w add y h add lineto\n" );
    fprintf ( _ps_f, "  x y h add lineto closepath stroke\n" );
    fprintf ( _ps_f, "} bdef\n" );
    break;

  case _fillrect:
    _ps_OutProc ( _newpath );
    fprintf ( _ps_f, "/_fr {\n" );
    fprintf ( _ps_f, "  /y exch def\n" );
    fprintf ( _ps_f, "  /x exch def\n" );
    fprintf ( _ps_f, "  /h exch def\n" );
    fprintf ( _ps_f, "  /w exch def\n" );
    fprintf ( _ps_f, "  _np\n" );
    fprintf ( _ps_f, "  x y moveto\n" );
    fprintf ( _ps_f, "  x w add y lineto\n" );
    fprintf ( _ps_f, "  x w add y h add lineto\n" );
    fprintf ( _ps_f, "  x y h add lineto closepath fill\n" );
    fprintf ( _ps_f, "} bdef\n" );
    break;

  case _hatchrect:
    _ps_OutProc ( _newpath );
    fprintf ( _ps_f, "/_hatchrect {\n" );
    fprintf ( _ps_f, "  20 dict begin\n" );
    fprintf ( _ps_f, "  /d exch def /a exch def\n" );
    fprintf ( _ps_f, "  /y exch def /x exch def\n" );
    fprintf ( _ps_f, "  /h exch def /w exch def\n" );
    fprintf ( _ps_f, "  /c a cos def /s a sin def\n" );
    fprintf ( _ps_f, "  /xa x c s mul h mul sub def\n" );
    fprintf ( _ps_f, "  /ya y c c mul h mul add def\n" );
    fprintf ( _ps_f, "  /xab s s mul w mul c s mul h mul add def\n" );
    fprintf ( _ps_f, "  /yab c c mul h mul c s mul w mul add def\n" );
    fprintf ( _ps_f, "  /xda c s mul h mul c c mul w mul add def\n" );
    fprintf ( _ps_f, "  /yda s s mul h mul c s mul w mul add def\n" );
    fprintf ( _ps_f, "  /lab xab dup mul yab dup mul add sqrt def\n" );
    fprintf ( _ps_f, "  /xab xab lab div d mul def\n" );
    fprintf ( _ps_f, "  /yab yab lab div d mul def\n" );
    fprintf ( _ps_f, "  /n lab d div 1 add round cvi def\n" );
    fprintf ( _ps_f, "  gsave\n" );
    fprintf ( _ps_f, "  _np x y moveto w 0 rlineto 0 h rlineto w neg 0 rlineto closepath clip\n" );
    fprintf ( _ps_f, "  n {\n" );
    fprintf ( _ps_f, "    _np xa ya moveto xda yda rlineto stroke\n" );
    fprintf ( _ps_f, "    /xa xa xab add def /ya ya yab sub def\n" );
    fprintf ( _ps_f, "  } repeat\n" );
    fprintf ( _ps_f, "  grestore\n" );
    fprintf ( _ps_f, "  end\n" );
    fprintf ( _ps_f, "} bdef\n" );
    break;

  case _moveto:
    fprintf ( _ps_f, "/_mt { moveto } bdef\n" );
    break;

  case _lineto:
    fprintf ( _ps_f, "/_lt { lineto } bdef\n" );
    break;

  case _shcone:
    _ps_OutProc ( _newpath );
    _ps_OutProc ( _moveto );
    _ps_OutProc ( _lineto );
    fprintf ( _ps_f, "/_shcone {\n");
    fprintf ( _ps_f, "  /y2 exch def /x2 exch def\n" );
    fprintf ( _ps_f, "  /y1 exch def /x1 exch def\n" );
    fprintf ( _ps_f, "  /y exch def  /x exch def\n" );
    fprintf ( _ps_f, "  /gray 1 def\n" );
    fprintf ( _ps_f, "  /dx1 x1 100 div def /dy1 y1 100 div def\n" );
    fprintf ( _ps_f, "  /dx2 x2 100 div def /dy2 y2 100 div def\n" );
    fprintf ( _ps_f, "  1 -0.01 0.8 {\n" );
    fprintf ( _ps_f, "    setgray\n" );
    fprintf ( _ps_f, "    _np\n" );
    fprintf ( _ps_f, "    x  y  _mt\n" );
    fprintf ( _ps_f, "    x x1 add y y1 add _lt\n" );
    fprintf ( _ps_f, "    x x2 add y y2 add _lt\n" );
    fprintf ( _ps_f, "    closepath fill\n" );
    fprintf ( _ps_f, "    /x1 x1 dx1 sub def /y1 y1 dy1 sub def\n" );
    fprintf ( _ps_f, "    /x2 x2 dx2 sub def /y2 y2 dy2 sub def\n" );
    fprintf ( _ps_f, "  } for\n" );
    fprintf ( _ps_f, "} def\n" );
    break;
  }
  _ps_written |= 1L << ((long)p);
} /*_ps_OutProc*/

boolean ps_OpenFile ( const char *filename, unsigned int dpi )
{
  _ps_written = 0x00;
  xmin = 0.0;
  ymin = 0.0;
  xmax = 0.0;
  ymax = 0.0;
  strcpy ( f_NAME, filename );
  if ( _ps_fsp < PS_FILE_STACK_LENGTH ) {
    _ps_f = _ps_fstack[_ps_fsp] = fopen ( f_NAME, "w" );
    if ( !_ps_f ) {
      printf ( "cannot create file %s\n", f_NAME);
      return false;
    }
    _ps_fsp++;
    fprintf ( _ps_f, "%%!PS-Adobe-3.0 EPSF-3.0\n" );
    if ( bbox ) {
      fprintf ( _ps_f, "%%%%BoundingBox: %d %d %d %d\n", bbx1, bby1, bbx2, bby2 );
      bbox = false;
    }
    else
      fprintf ( _ps_f, "%%%%BoundingBox: \n" );
    fprintf ( _ps_f, "%% %s\n", filename );
    fprintf ( _ps_f, "gsave\n" );
    fprintf ( _ps_f, "20 dict begin\n" );
    fprintf ( _ps_f, "/bdef { bind def } def\n" );
    fprintf ( _ps_f, "/resolution %u def\n", dpi );
    fprintf ( _ps_f, "72 resolution div dup scale\n" );
    _ps_written = 0;
    _ps_cgray = _ps_cred = _ps_cgreen = _ps_cblue = -1.0;
    ps_Set_Gray ( 0.0 );
    _ps_cwidth = -1.0;
    ps_Set_Line_Width ( 1.0 );
    _ps_bmp_y = -1;
    noext = true;
    _ps_dpi = (short)dpi;
    _psl_InitPSLib ();
    return true;
  }
  else {
    printf ( "too many open PostScript files.\n" );
    return false;
  }
} /*ps_OpenFile*/

void ps_CloseFile ( void )
{
  if ( _ps_fsp > 0 ) {
    fprintf ( _ps_f, "end\n" );
    fprintf ( _ps_f, "grestore\n" );
    fprintf ( _ps_f, "showpage\n" );
    fprintf ( _ps_f, "%%\n" );
    fclose ( _ps_f );
    _ps_fsp --;
    if ( _ps_fsp > 0 )
      _ps_f = _ps_fstack[_ps_fsp-1];
    else
      _ps_f = NULL;
  }
} /*ps_CloseFile*/

void ps_WriteBBox ( float x1, float y1, float x2, float y2 )
{
  bbx1 = (int)(floor(x1+0.5));
  bby1 = (int)(floor(y1+0.5));
  bbx2 = (int)(floor(x2+0.5));
  bby2 = (int)(floor(y2+0.5));
  bbox = true;
} /*ps_WriteBBox*/

void ps_GetSize ( float *x1, float *y1, float *x2, float *y2 )
{
  *x1 = xmin;
  *y1 = ymin;
  *x2 = xmax;
  *y2 = ymax;
} /*ps_GetSize*/

