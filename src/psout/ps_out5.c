
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2005, 2010                            */
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

void ps_Draw_Circle ( float x, float y, float r )
{
  _ps_OutProc(_dcircle);
  fprintf(_ps_f, "%5.*f %5.*f %5.*f _dc\n",
	  ps_dec_digits, x, ps_dec_digits, y, ps_dec_digits, r);
  _ps_ExtendRect(x + r, y + r);
  _ps_ExtendRect(x - r, y - r);
} /*ps_Draw_Circle*/

void ps_Fill_Circle ( float x, float y, float r )
{
  _ps_OutProc(_fcircle);
  fprintf(_ps_f, "%5.*f %5.*f %5.*f _fc\n",
	  ps_dec_digits, x, ps_dec_digits, y, ps_dec_digits, r);
  _ps_ExtendRect(x + r, y + r);
  _ps_ExtendRect(x - r, y - r);
} /*ps_Fill_Circle*/

void ps_Mark_Circle ( float x, float y )
{
  float r;

  _ps_OutProc ( _mcircle );
  fprintf ( _ps_f, "%5.*f %5.*f _mc\n", ps_dec_digits, x, ps_dec_digits, y );
  r = (float)(_ps_dpi*0.02);
  _ps_ExtendRect ( x + r, y + r );
  _ps_ExtendRect ( x - r, y - r );
} /*ps_Mark_Circle*/

void ps_Draw_Arc ( float x, float y, float r, float a0, float a1 )
{
  fprintf ( _ps_f, "newpath %5.*f %5.*f %5.*f %7.2f %7.2f arc stroke\n",
            ps_dec_digits, x, ps_dec_digits, y, ps_dec_digits, r,
            180.0/PI*a0, 180.0/PI*a1 );
} /*ps_Draw_Arc*/

void ps_Write_Command ( char *command )
{
  fprintf ( _ps_f, "%s\n", command );
} /*ps_Write_Command*/

void ps_Newpath ( void )
{
  _ps_OutProc(_newpath);
  fprintf(_ps_f, "_np\n");
} /*ps_Newpath*/

void ps_MoveTo ( float x, float y )
{
  _ps_OutProc(_moveto);
  _ps_OutProc(_lineto);
  fprintf(_ps_f, "%5.*f %5.*f _mt\n", ps_dec_digits, x, ps_dec_digits, y);
  _ps_ExtendRect(x, y);
} /*ps_MoveTo*/

void ps_LineTo ( float x, float y )
{
  _ps_OutProc(_moveto);
  _ps_OutProc(_lineto);
  fprintf(_ps_f, "%5.*f %5.*f _lt\n", ps_dec_digits, x, ps_dec_digits, y);
  _ps_ExtendRect(x, y);
} /*ps_LineTo*/

void ps_ShCone ( float x, float y, float x1, float y1, float x2, float y2 )
{
  _ps_OutProc(_shcone);
  fprintf(_ps_f, "%5.*f %5.*f %5.*f %5.*f %5.*f %5.*f _shcone\n",
	  ps_dec_digits, x, ps_dec_digits, y, ps_dec_digits, x1,
	  ps_dec_digits, y1, ps_dec_digits, x2, ps_dec_digits, y2);
  _ps_cgray = _ps_cred = _ps_cgreen = _ps_cblue = 0.8;
} /*ps_ShCone*/

void ps_GSave ( void )
{
  fprintf(_ps_f, "gsave\n");
} /*ps_GSave*/

void ps_GRestore ( void )
{
  fprintf(_ps_f, "grestore\n");
  _ps_cgray = _ps_cred = _ps_cgreen = _ps_cblue = -1.0;
  _ps_cwidth = -1.0;
} /*ps_GRestore*/

void ps_BeginDict ( int n )
{
  fprintf ( _ps_f, "%d dict begin\n", n );
} /*ps_BeginDict*/

void ps_EndDict ( void )
{
  fprintf ( _ps_f, "end\n" );
  _ps_written = 0;
} /*ps_EndDict*/

void ps_DenseScreen ( void )
{
  fprintf(_ps_f, "currentscreen\n");
  fprintf(_ps_f, "3 -1 roll\n");
  fprintf(_ps_f, "2 mul\n");
  fprintf(_ps_f, "3 1 roll\n");
  fprintf(_ps_f, "setscreen\n");
} /*ps_DenseScreen*/

