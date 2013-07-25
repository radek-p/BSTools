
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2012                                  */
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

#include "xgedit.h"
#include "xgeprivate.h"


#define T0 0.5
#define T1 0.53125
#define T2 0.46875
static char xround = 0, yround = 0;

static void _xge_qb_roundpixel ( vector3f *ec, vector3f *v, xpoint *px )
{
  float x, fx, y, fy, t;

  x = ec->x + v->x + 0.5;  fx = floor ( x );
  y = ec->y - v->y + 0.5;  fy = floor ( y );
  switch ( xround ) {
 case 1: t = T1;  break;
 case 2: t = T2;  break;
default: t = T0;  break;
  }
  if ( x-fx < t ) { px->x = (short)fx;      xround = 1; }
             else { px->x = (short)(fx+1);  xround = 2; }
  switch ( yround ) {
 case 1: t = T1;  break;
 case 2: t = T2;  break;
default: t = T0;  break;
  }
  if ( y-fy < t ) { px->y = (short)fy;      yround = 1; }
             else { px->y = (short)(fy+1);  yround = 2; }
} /*_xge_qb_roundpixel*/

static void _xge_r_DrawArc ( vector3f va, vector3f vb, vector3f ec )
{
#define STKSIZE 16

  typedef struct {
      xpoint   px;
      vector3f vx;
      short    level;
    } stackel;

  static float t[9] = {0.7071067812,0.5411961001,0.5097955791,
                       0.5024192862,0.5006029982,0.500150636,
                       0.5000376519,0.5000094125,0.5000023531};
  void     *sp;
  xpoint   last, next;
  float    lastz, nextz;
  short    levela, levelb, levelc;
  vector3f vc;
  stackel  *stack;
  int      stp;

#define NEIGHBOURS(a,b) \
  abs(a.x-b.x) <= 1 && abs(a.y-b.y) <= 1
#define OUTPUT(px,z) \
  { if ( z >= 0.0 ) { \
      if ( _pkv_npix == PKV_BUFSIZE ) { _pkv_npix -= 2;  PKV_FLUSH \
        memcpy ( &_pkv_pixbuf[0], &_pkv_pixbuf[PKV_BUFSIZE-2], 2*sizeof(xpoint) ); \
        _pkv_npix = 2; } \
      if ( _pkv_npix > 1 ) \
        { if ( NEIGHBOURS ( px, _pkv_pixbuf[_pkv_npix-2] ) ) _pkv_npix--; } \
      _pkv_pixbuf[_pkv_npix] = px;  _pkv_npix++; \
    } \
    lastz = z; \
  }
#define PUSH(p,v,l) \
  { if (stp < STKSIZE) { stack[stp].px = p;  stack[stp].vx = (v); \
                         stack[stp].level = l;  stp++; } }
#define READ(p,v,l) \
  { p = stack[stp-1].px;  v = stack[stp-1].vx;  l = stack[stp-1].level; }
#define POP \
  stp--;

  sp = pkv_GetScratchMemTop ();
  stack = pkv_GetScratchMem ( STKSIZE*sizeof(stackel) );
  if ( !stack )
    return;
  stp = 0;  /* stack becomes empty */
  xround = yround = 0;
  _xge_qb_roundpixel ( &ec, &va, &last );
  levela = 0;
  OUTPUT ( last, ec.z+va.z )
  _xge_qb_roundpixel ( &ec, &vb, &next );
  PUSH ( next, vb, 0 )
  do {
    READ ( next, vb, levelb )
    nextz = ec.z+vb.z;
    if ( NEIGHBOURS ( last, next ) && fabs(nextz-lastz) < 0.1 ) {
      POP
      last   = next;
      va     = vb;
      levela = levelb;
      OUTPUT ( last, nextz );
    }
    else {
      levelc = max ( levela, levelb );
      if ( levelc < 9 ) {
        AddVector3f ( &va, &vb, &vc );
        MultVector3f ( t[levelc], &vc, &vc );
      }
      else
        MidPoint3f ( &va, &vb, &vc );
      _xge_qb_roundpixel ( &ec, &vc, &next );
      PUSH ( next, vc, levelc+1 )
    }
  } while ( stp );
  PKV_FLUSH
  pkv_SetScratchMemTop ( sp );
} /*_xge_r_DrawArc*/

static void _xge_QuatRotBallDrawCircle ( short xc, short yc, short r, trans3f *tr,
                                         vector3f *c, vector3f *a, vector3f *b )
{
  vector3f ec, va, vb, uva, uvb;

  TransVector3f ( tr, c, &ec );
  ec.x = (float)xc + (float)r*ec.x;
  ec.y = (float)yc - (float)r*ec.y;
  ec.z *= (float)r;
  TransVector3f ( tr, a, &va );
  MultVector3f ( (float)r, &va, &va );
  SetVector3f ( &uva, -va.x, -va.y, -va.z );
  TransVector3f ( tr, b, &vb );
  MultVector3f ( (float)r, &vb, &vb );
  SetVector3f ( &uvb, -vb.x, -vb.y, -vb.z );
  _xge_r_DrawArc ( va, vb, ec );
  _xge_r_DrawArc ( vb, uva, ec );
  _xge_r_DrawArc ( uva, uvb, ec );
  _xge_r_DrawArc ( uvb, va, ec );
} /*_xge_QuatRotBallDrawCircle*/

void _xge_QuatRotBallDrawCircles ( short xc, short yc, short r, trans3f *tr )
{
  static vector3f z = {0.0, 0.0, 0.0},
    e1 = {1.0, 0.0, 0.0}, e2 = {0.0, 1.0, 0.0}, e3 = {0.0, 0.0, 1.0},
    f1 = {0.5*SQRT2, 0.0, 0.0}, f2 = {0.0, 0.5*SQRT2, 0.0},
    g1 = {0.5*SQRT2, 0.5*SQRT2,0.0}, g2 = {0.5*SQRT2, -0.5*SQRT2, 0.0},
    c1 = {0.0, 0.0, 0.5*SQRT2}, c2 = {0.0, 0.0, -0.5*SQRT2};
  void *sp;

  sp = pkv_GetScratchMemTop ();
  _pkv_InitPixelBuffer ();
  if ( _pkv_pixbuf ) {
    _pkv_OutputPixels = xge_OutPixels;
        /* equator */
    _xge_QuatRotBallDrawCircle ( xc, yc, r, tr, &z, &e1, &e2 );
    _xge_QuatRotBallDrawCircle ( xc, yc, r, tr, &c1, &f1, &f2 );
    _xge_QuatRotBallDrawCircle ( xc, yc, r, tr, &c2, &f1, &f2 );
        /* meridians */
    _xge_QuatRotBallDrawCircle ( xc, yc, r, tr, &z, &e1, &e3 );
    _xge_QuatRotBallDrawCircle ( xc, yc, r, tr, &z, &e2, &e3 );
    _xge_QuatRotBallDrawCircle ( xc, yc, r, tr, &z, &g1, &e3 );
    _xge_QuatRotBallDrawCircle ( xc, yc, r, tr, &z, &g2, &e3 );
  }
  _pkv_DestroyPixelBuffer ();
  pkv_SetScratchMemTop ( sp );
} /*_xge_QuatRotBallDrawCircles*/

