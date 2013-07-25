
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2005, 2009                            */
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

static void _ps_OutputPath2f ( int n, const point2f *p )
{
  int i;

  _ps_OutProc ( _moveto );
  _ps_OutProc ( _lineto );
  fprintf ( _ps_f, "%5.*f %5.*f newpath _mt\n",
	    ps_dec_digits, p[0].x, ps_dec_digits, p[0].y );
  _ps_ExtendRect ( p[0].x, p[0].y );
  for ( i = 1; i < n; i++ ) {
    fprintf ( _ps_f, "%5.*f %5.*f _lt\n",
	      ps_dec_digits, p[i].x, ps_dec_digits, p[i].y );
    _ps_ExtendRect ( p[i].x, p[i].y );
  }
} /*_ps_OutputPath2f*/

static void _ps_OutputPath2d ( int n, const point2d *p )
{
  int i;

  _ps_OutProc ( _moveto );
  _ps_OutProc ( _lineto );
  fprintf ( _ps_f, "%5.*f %5.*f newpath _mt\n",
	    ps_dec_digits, p[0].x, ps_dec_digits, p[0].y );
  _ps_ExtendRect ( (float)p[0].x, (float)p[0].y );
  for ( i = 1; i < n; i++ ) {
    fprintf ( _ps_f, "%5.*f %5.*f _lt\n",
	      ps_dec_digits, p[i].x, ps_dec_digits, p[i].y );
    _ps_ExtendRect ( (float)p[i].x, (float)p[i].y );
  }
} /*_ps_OutputPath2d*/

void ps_Draw_Polyline2f ( int n, const point2f *p )
{
  if ( n < 2 )
    return;
  _ps_OutputPath2f ( n, p );
  fprintf(_ps_f, " stroke\n");
} /*ps_Draw_Polyline2f*/

void ps_Draw_Polyline2d ( int n, const point2d *p )
{
  if ( n < 2 )
    return;
  _ps_OutputPath2d ( n, p );
  fprintf(_ps_f, " stroke\n");
} /*ps_Draw_Polyline2d*/

static void _ps_OutputPath2Rf ( int n, const point3f *p )
{
  int i;

  _ps_OutProc ( _moveto );
  _ps_OutProc ( _lineto );
  fprintf ( _ps_f, "%5.*f %5.*f newpath _mt\n",
	    ps_dec_digits, p[0].x/p[0].z, ps_dec_digits, p[0].y/p[0].z );
  _ps_ExtendRect ( p[0].x/p[0].z, p[0].y/p[0].z );
  for ( i = 1; i < n; i++ ) {
    fprintf ( _ps_f, "%5.*f %5.*f _lt\n",
	      ps_dec_digits, p[i].x/p[i].z, ps_dec_digits, p[i].y/p[i].z );
    _ps_ExtendRect ( p[i].x/p[i].z, p[i].y/p[i].z );
  }
} /*_ps_OutputPath2Rf*/

static void _ps_OutputPath2Rd ( int n, const point3d *p )
{
  int i;

  _ps_OutProc ( _moveto );
  _ps_OutProc ( _lineto );
  fprintf ( _ps_f, "%5.*f %5.*f newpath _mt\n",
	    ps_dec_digits, p[0].x/p[0].z, ps_dec_digits, p[0].y/p[0].z );
  _ps_ExtendRect ( (float)(p[0].x/p[0].z), (float)(p[0].y/p[0].z) );
  for ( i = 1; i < n; i++ ) {
    fprintf ( _ps_f, "%5.*f %5.*f _lt\n",
	      ps_dec_digits, p[i].x/p[i].z, ps_dec_digits, p[i].y/p[i].z );
    _ps_ExtendRect ( (float)(p[i].x/p[i].z), (float)(p[i].y/p[i].z) );
  }
} /*_ps_OutputPatch2Rd*/

void ps_Draw_Polyline2Rf ( int n, const point3f *p )
{
  if ( n < 2 )
    return;
  _ps_OutputPath2Rf ( n, p );
  fprintf(_ps_f, " stroke\n");
} /*ps_Draw_Polyline2Rf*/

void ps_Draw_Polyline2Rd ( int n, const point3d *p )
{
  if ( n < 2 )
    return;
  _ps_OutputPath2Rd ( n, p );
  fprintf(_ps_f, " stroke\n");
} /*ps_Draw_Polyline2Rd*/

void ps_Set_Clip_Polygon2f ( int n, const point2f *p )
{
  if ( n < 2 )
    return;
  _ps_OutputPath2f ( n, p );
  fprintf(_ps_f, " closepath clip\n");
} /*ps_Set_Clip_Polygon2f*/

void ps_Set_Clip_Polygon2d ( int n, const point2d *p )
{
  if ( n < 2 )
    return;
  _ps_OutputPath2d ( n, p );
  fprintf(_ps_f, " closepath clip\n");
} /*ps_Set_Clip_Polygon2d*/

void ps_Set_Clip_Polygon2Rf ( int n, const point3f *p )
{
  if ( n < 2 )
    return;
  _ps_OutputPath2Rf ( n, p );
  fprintf(_ps_f, " closepath clip\n");
} /*ps_Set_Clip_Polygon2Rf*/

void ps_Set_Clip_Polygon2Rd ( int n, const point3d *p )
{
  if ( n < 2 )
    return;
  _ps_OutputPath2Rd ( n, p );
  fprintf(_ps_f, " closepath clip\n");
} /*ps_Set_Clip_Polygon2Rd*/

void ps_Fill_Polygon2f ( int n, const point2f *p )
{
  if ( n < 2 )
    return;
  _ps_OutputPath2f ( n, p );
  fprintf(_ps_f, "closepath fill\n");
} /*ps_Fill_Polygon2f*/

void ps_Fill_Polygon2d ( int n, const point2d *p )
{
  if ( n < 2 )
    return;
  _ps_OutputPath2d ( n, p );
  fprintf(_ps_f, "closepath fill\n");
} /*ps_Fill_Polygon2d*/

void ps_Fill_Polygon2Rf ( int n, const point3f *p )
{
  if ( n < 2 )
    return;
  _ps_OutputPath2Rf ( n, p );
  fprintf(_ps_f, "closepath fill\n");
} /*ps_Fill_Polygon2Rf*/

void ps_Fill_Polygon2Rd ( int n, const point3d *p )
{
  if ( n < 2 )
    return;
  _ps_OutputPath2Rd ( n, p );
  fprintf(_ps_f, "closepath fill\n");
} /*ps_Fill_Polygon2Rd*/

