
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2005, 2009                            */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

#include <math.h>
#include <string.h>
#include <stdio.h>

#include "pkvaria.h"
#include "pkgeom.h"
#include "psout.h"


static float xx1, yy1, xx2, yy2, tt1, tt2, lgt, vx, vy;

void psl_SetLine ( float x1, float y1, float x2, float y2, float t1, float t2 )
{
  float dx, dy;

  xx1 = x1;
  yy1 = y1;
  xx2 = x2;
  yy2 = y2;
  tt1 = t1;
  tt2 = t2;
  dx = x2 - x1;
  dy = y2 - y1;
  lgt = (float)sqrt(dx * dx + dy * dy);
  vx = dx / lgt;
  vy = dy / lgt;
} /*psl_SetLine*/

void psl_GetPointf (float t, float *x, float *y)
{
  *x = xx1 + (xx2 - xx1) * (t - tt1) / (tt2 - tt1);
  *y = yy1 + (yy2 - yy1) * (t - tt1) / (tt2 - tt1);
}  /*psl_GetPointf*/

void psl_GetPointd (double t, double *x, double *y)
{
  *x = xx1 + (xx2 - xx1) * (t - tt1) / (tt2 - tt1);
  *y = yy1 + (yy2 - yy1) * (t - tt1) / (tt2 - tt1);
}  /*psl_GetPointd*/

float psl_GetDParam(float dl)
{
  return (dl / lgt * (tt2 - tt1));
}  /*psl_GetDParam*/

void psl_GoAlong(float s, float *x, float *y)
{
  *x += s * vx;
  *y += s * vy;
}  /*psl_GoAlong*/

void psl_GoPerp(float s, float *x, float *y)
{
  *x -= s * vy;
  *y += s * vx;
}  /*pls_GoPerp*/

void psl_Tick(float t)
{
  float x, y, x1, y1;

  psl_GetPointf ( t, &x, &y );
  x1 = x;
  y1 = y;
  psl_GoPerp ( tickl, &x, &y );
  psl_GoPerp ( -tickl, &x1, &y1 );
  ps_Set_Line_Width ( tickw );
  ps_Draw_Line ( x, y, x1, y1 );
}  /*psl_Tick*/

void psl_BTick ( float t )
{
  float x, y, x1, y1;

  psl_GetPointf ( t, &x, &y );
  x1 = x;
  y1 = y;
  psl_GoPerp ( (float)(tickl + 2.0), &x, &y );
  psl_GoPerp ( (float)(-tickl - 2.0), &x1, &y1 );
  ps_Set_Line_Width ( (float)(tickw + 4.0) );
  ps_Draw_Line ( x, y, x1, y1 );
}  /*psl_BTick*/

void psl_HTick(float t, boolean left)
{
  float x, y, x1, y1;

  psl_GetPointf ( t, &x, &y );
  x1 = x;
  y1 = y;
  if ( left )
    psl_GoPerp ( tickl, &x, &y );
  else
    psl_GoPerp ( -tickl, &x1, &y1 );
  ps_Set_Line_Width ( tickw );
  ps_Draw_Line ( x, y, x1, y1 );
}  /*psl_HTick*/

void psl_Dot(float t)
{
  float x, y;

  psl_GetPointf ( t, &x, &y );
  ps_Fill_Circle ( x, y, dotr );
}  /*psl_Dot*/

void psl_HDot(float t)
{
  float x, y;

  psl_GetPointf ( t, &x, &y );
  ps_Mark_Circle ( x, y );
}  /*psl_HDot*/

void psl_TrMark ( float x, float y )
{
  static boolean first = true;
  char sx[8], sy[8];
  char s[21];

  if ( first ) {
    ps_Write_Command("/trm {");
    ps_Write_Command("  /y exch def");
    ps_Write_Command("  /x exch def");
    ps_Write_Command("  1 setgray");
    ps_Write_Command(
      "  newpath x y moveto 15 -35 rlineto -30 0 rlineto closepath fill");
    ps_Write_Command("  2 setlinewidth");
    ps_Write_Command("  0 setgray");
    ps_Write_Command(
      "  newpath x y moveto 15 -35 rlineto -30 0 rlineto closepath stroke");
    ps_Write_Command("} def");
    first = false;
  }
  sprintf(sx, "%7.2f", x);
  sprintf(sy, "%7.2f", y);
  sprintf(s, "%s %s trm", sx, sy);
  ps_Write_Command(s);
}  /*psl_TrMark*/

void psl_BlackTrMark(float x, float y)
{
  static boolean first = true;
  char sx[8], sy[8];
  char s[23];

  if ( first ) {
    ps_Write_Command("/btrm {");
    ps_Write_Command("  /y exch def");
    ps_Write_Command("  /x exch def");
    ps_Write_Command(
      "  newpath x y moveto 15 -35 rlineto -30 0 rlineto closepath fill");
    ps_Write_Command("} def");
    first = false;
  }
  sprintf(sx, "%7.2f", x);
  sprintf(sy, "%7.2f", y);
  sprintf(s, "%s %s btrm", sx, sy);
  ps_Write_Command(s);
}  /*psl_BlackTrMark*/

void psl_HighTrMark(float x, float y)
{
  static boolean first = true;
  char sx[8], sy[8];
  char s[21];

  if ( first ) {
    ps_Write_Command("/htrm {");
    ps_Write_Command("  /y exch def");
    ps_Write_Command("  /x exch def");
    ps_Write_Command("  1 setgray");
    ps_Write_Command(
      "  newpath x y moveto 15 -70 rlineto -30 0 rlineto closepath fill");
    ps_Write_Command("  2 setlinewidth");
    ps_Write_Command("  0 setgray");
    ps_Write_Command(
      "  newpath x y moveto 15 -70 rlineto -30 0 rlineto closepath stroke");
    ps_Write_Command("} def");
    first = false;
  }
  sprintf(sx, "%7.2f", x);
  sprintf(sy, "%7.2f", y);
  sprintf(s, "%s %s htrm", sx, sy);
  ps_Write_Command(s);
}  /*psl_HighTrMark*/

void psl_BlackHighTrMark(float x, float y)
{
  static boolean first = true;
  char sx[8], sy[8];
  char s[23];

  if ( first ) {
    ps_Write_Command("/bhtrm {");
    ps_Write_Command("  /y exch def");
    ps_Write_Command("  /x exch def");
    ps_Write_Command(
      "  newpath x y moveto 15 -70 rlineto -30 0 rlineto closepath fill");
    ps_Write_Command("} def");
    first = false;
  }
  sprintf(sx, "%7.2f", x);
  sprintf(sy, "%7.2f", y);
  sprintf(s, "%s %s bhtrm", sx, sy);
  ps_Write_Command(s);
}  /*psl_BlackHighTrMark*/

void psl_LTrMark ( float t )
{
  float x, y;

  psl_GetPointf ( t, &x, &y );
  psl_TrMark ( x, y );
}  /*psl_LTrMark*/

void psl_BlackLTrMark ( float t )
{
  float x, y;

  psl_GetPointf ( t, &x, &y );
  psl_BlackTrMark ( x, y );
}  /*psl_BlackLTrMark*/

void psl_HighLTrMark ( float t )
{
  float x, y;

  psl_GetPointf ( t, &x, &y );
  psl_HighTrMark ( x, y );
}  /*psl_HighLTrMark*/

void psl_BlackHighLTrMark ( float t )
{
  float x, y;

  psl_GetPointf ( t, &x, &y );
  psl_BlackHighTrMark ( x, y );
}  /*psl_BlackHighLTrMark*/

void psl_Arrow ( float t, boolean sgn )
{
  point2f pp[3];

  psl_GetPointf ( t, &pp[0].x, &pp[0].y );
  pp[1] = pp[0];
  if (sgn)
    psl_GoAlong ( -arrowl, &pp[1].x, &pp[1].y );
  else
    psl_GoAlong ( arrowl, &pp[1].x, &pp[1].y );
  pp[2] = pp[1];
  psl_GoPerp ( arroww, &pp[1].x, &pp[1].y );
  psl_GoPerp ( -arroww, &pp[2].x, &pp[2].y );
  ps_Fill_Polygon2f ( 3, pp );
}  /*psl_Arrow*/

void psl_BkArrow ( float t, boolean sgn )
{
  point2f pp[4];

  ps_Set_Line_Width(6.0);
  psl_GetPointf ( t, &pp[0].x, &pp[0].y );
  pp[1] = pp[0];
  if (sgn)
    psl_GoAlong ( -arrowl, &pp[1].x, &pp[1].y );
  else
    psl_GoAlong ( arrowl, &pp[1].x, &pp[1].y );
  pp[2] = pp[1];
  psl_GoPerp(arroww, &pp[1].x, &pp[1].y);
  psl_GoPerp(-arroww, &pp[2].x, &pp[2].y);
  pp[3] = pp[0];
  ps_Draw_Polyline2f ( 4, pp );
}  /*psl_BkArrow*/

void psl_Draw ( float ta, float tb, float w )
{
  float x1, y1, x2, y2;

  psl_GetPointf ( ta, &x1, &y1 );
  psl_GetPointf ( tb, &x2, &y2 );
  ps_Set_Line_Width ( w );
  ps_Draw_Line ( x1, y1, x2, y2 );
}  /*psl_Draw*/

void psl_ADraw ( float ta, float tb, float ea, float eb, float w )
{
  float x1, y1, x2, y2;

  psl_GetPointf ( ta, &x1, &y1 );
  psl_GoAlong ( ea, &x1, &y1 );
  psl_GetPointf ( tb, &x2, &y2 );
  psl_GoAlong ( eb, &x2, &y2 );
  ps_Set_Line_Width ( w );
  ps_Draw_Line( x1, y1, x2, y2 );
}  /*psl_ADraw*/

void psl_MapsTo ( float t )
{
#define MAG 0.85
  float x, y, x1, y1, x2, y2, x3, y3;
  char s[80];

  psl_GetPointf ( t, &x, &y );
  ps_Newpath ();
  psl_GoPerp ( (float)(0.4*MAG), &x, &y );
  ps_MoveTo ( x, y );  x1 = x;  y1 = y;
  psl_GoAlong ( (float)(-16.0*MAG), &x1, &y1 );
  psl_GoPerp ( (float)(5.6*MAG), &x1, &y1 );  x2 = x1;  y2 = y1;
  psl_GoAlong ( (float)(-14.0*MAG), &x2, &y2  );
  psl_GoPerp ( (float)(10.0*MAG), &x2, &y2 );  x3 = x2;  y3 = y2;
  psl_GoAlong ( (float)(-4.5*MAG), &x3, &y3 );
  psl_GoPerp ( (float)(7.8*MAG), &x3, &y3 );
  sprintf ( s, "%6.2f %6.2f %6.2f %6.2f %6.2f %6.2f curveto", x1, y1, x2, y2, x3, y3 );
  ps_Write_Command ( s );  x = x3;  y = y3;
  psl_GoAlong ( (float)(-3.4*MAG), &x, &y );
  ps_LineTo ( x, y );  x1 = x;  y1 = y;
  psl_GoAlong ( (float)(3.0*MAG), &x1, &y1 );
  psl_GoPerp ( (float)(-8.0*MAG), &x1, &y1 );  x2 = x1;  y2 = y1;
  psl_GoPerp ( (float)(-8.0*MAG), &x2, &y2 );
  psl_GoAlong ( (float)(8.0*MAG), &x2, &y2 );  x3 = x2;  y3 = y2;
  psl_GoPerp ( (float)(-7.8*MAG), &x3, &y3 );
  psl_GoAlong ( (float)(12.0*MAG), &x3, &y3 );
  sprintf ( s, "%6.2f %6.2f %6.2f %6.2f %6.2f %6.2f curveto", x1, y1, x2, y2, x3, y3 );
  ps_Write_Command ( s );  x1 = x3;  y1 = y3;
  psl_GoAlong ( (float)(-12.0*MAG), &x1, &y1 );
  psl_GoPerp ( (float)(-7.8*MAG), &x1, &y1 );  x2 = x1;  y2 = y1;
  psl_GoAlong ( (float)(-8.0*MAG), &x2, &y2 );
  psl_GoPerp ( (float)(-8.0*MAG), &x2, &y2 );  x3 = x2;  y3 = y2;
  psl_GoAlong ( (float)(-3.0*MAG), &x3, &y3 );
  psl_GoPerp ( (float)(-8.0*MAG), &x3, &y3 );
  sprintf ( s, "%6.2f %6.2f %6.2f %6.2f %6.2f %6.2f curveto", x1, y1, x2, y2, x3, y3 );
  ps_Write_Command ( s );
  psl_GoAlong ( (float)(3.4*MAG), &x3, &y3 );
  ps_LineTo ( x3, y3 );  x1 = x3;  y1 = y3;
  psl_GoAlong ( (float)(4.5*MAG), &x1, &y1 );
  psl_GoPerp ( (float)(7.8*MAG), &x1, &y1 );  x2 = x1;  y2 = y1;
  psl_GoAlong ( (float)(14.0*MAG), &x2, &y2 );
  psl_GoPerp ( (float)(10.0*MAG), &x2, &y2 );  x3 = x2;  y3 = y2;
  psl_GoAlong ( (float)(16.0*MAG), &x3, &y3 );
  psl_GoPerp ( (float)(5.6*MAG), &x3, &y3 );
  sprintf ( s, "%6.2f %6.2f %6.2f %6.2f %6.2f %6.2f curveto", x1, y1, x2, y2, x3, y3 );
  ps_Write_Command ( s );
  ps_Write_Command ( "closepath fill" );
} /*psl_MapsTo*/

#define mx              (float)1.4

static void SetupTrans(trans2f *tr, float mag, float ang, float t, byte cc)
{
  tr->U1.a[0][0] = mx * mag * vx;
  tr->U1.a[0][1] = -mag * vy;
  tr->U1.a[1][0] = mx * mag * vy;
  tr->U1.a[1][1] = mag * vx;
  if (!(cc & 1)) {
    tr->U1.a[0][0] = -tr->U1.a[0][0];
    tr->U1.a[1][0] = -tr->U1.a[1][0];
  }
  if ((cc / 2) & 1) {
    tr->U1.a[0][1] = -tr->U1.a[0][1];
    tr->U1.a[1][1] = -tr->U1.a[1][1];
  }
  tr->U1.a[0][2] = 0.0;
  tr->U1.a[1][2] = 0.0;
  RotTrans2f ( tr, ang );
  psl_GetPointf ( t, &tr->U1.a[0][2], &tr->U1.a[1][2] );
}  /*SetupTrans*/

static void rb(point3f *c, point2f *p)
{
  point3f q[4];
  short i, j, k;
  float t;

  for (k = 0; k <= 30; k++) {
    t = (float)(k / 30.0);
    memmove((char *)q, c, sizeof(point3f) * 4L);
    for (j = 1; j <= 3; j++) {
      for (i = 0; i <= 3 - j; i++)
	InterPoint3f ( &q[i], &q[i+1], t, &q[i] );
    }
    p[k].x = q[0].x / q[0].z;
    p[k].y = q[0].y / q[0].z;
  }
}  /*rb*/

static void CVC(point3f *c, point3f *d, trans2f *tr)
{
  short i;

  for (i = 0; i <= 3; i++) {
    d[i].x = (tr->U1.a[0][0] * c[i].x + tr->U1.a[0][1] * c[i].y + tr->U1.a[0][2]) * c[i].z;
    d[i].y = (tr->U1.a[1][0] * c[i].x + tr->U1.a[1][1] * c[i].y + tr->U1.a[1][2]) * c[i].z;
    d[i].z = c[i].z;
  }
}  /*CVC*/

void psl_DrawEye(float t, byte cc, float mag, float ang)
{
  static point3f c1[4] = {
    { 0.0, -50.0, 1.0 },
    { -50.0, -50.0, 0.333333 },
    { -50.0, 50.0, 0.333333 },
    { 0.0, 50.0, 1.0 }
  };

  static point3f c2[4] = {
    { 0.0, 50.0, 1.0 },
    { 50.0, 50.0, 0.333333 },
    { 50.0, -50.0, 0.333333 },
    { 0.0, -50.0, 1.0 }
  };

  static point3f c3[4] = {
    { 5.0, -25.0, 1.0 },
    { -20.0, -25.0, 0.333333 },
    { -20.0, 25.0, 0.333333 },
    { 5.0, 25.0, 1.0 }
  };

  static point3f c4[4] = {
    { 5.0, 25.0, 1.0 },
    { 30.0, 25.0, 0.333333 },
    { 30.0, -25.0, 0.333333 },
    { 5.0, -25.0, 1.0 }
  };

  static point3f c5[4] = {
    { 7.0, -3.0, 1.0 },
    { 4.0, -3.0, 0.333333 },
    { 4.0, 3.0, 0.333333 },
    { 7.0, 3.0, 1.0 }
  };

  static point3f c6[4] = {
    { 7.0, 3.0, 1.0 },
    { 10.0, 3.0, 0.333333 },
    { 10.0, -3.0, 0.333333 },
    { 7.0, -3.0, 1.0 }
  };

  static point3f c7[4] = {
    { -55.0, -12.0, 1.0 },
    { -20.0, 68.0, 1.0 },
    { 30.0, 61.0, 1.0 },
    { 36.0, 0.0, 1.0 }
  };

  static point3f c8[4] = {
    { 36.0, 0.0, 1.0 },
    { 16.0, -67.0, 1.0 },
    { -32.0, -61.0, 1.0 },
    { -55.0, -12.0, 1.0 }
  };

  point3f d[4];
  point2f p[61];
  point2f *qq;
  trans2f tr;

  SetupTrans(&tr, mag, ang, t, cc);
  ps_Set_Gray(0.0);
  ps_Set_Line_Width(1.0);
  CVC(c1, d, &tr);
  rb(d, p);
  CVC(c2, d, &tr);
  rb(d, (&p[30]));
  qq = &p[3];
  ps_Draw_Polyline2f ( 22, qq );
  qq = &p[36];
  ps_Draw_Polyline2f ( 16, qq );
  CVC(c3, d, &tr);
  rb(d, p);
  CVC(c4, d, &tr);
  rb(d, (&p[30]));
  ps_Fill_Polygon2f ( 61, p );
  ps_Set_Gray(1.0);
  CVC(c5, d, &tr);
  rb(d, p);
  CVC(c6, d, &tr);
  rb(d, (&p[30]));
  ps_Fill_Polygon2f ( 61, p );
  ps_Set_Gray(0.0);
  ps_Set_Line_Width(2.0);
  CVC(c7, d, &tr);
  rb(d, p);
  CVC(c8, d, &tr);
  rb(d, (&p[30]));
  ps_Draw_Polyline2f ( 61, p );
}  /*psl_DrawEye*/

#undef mx


