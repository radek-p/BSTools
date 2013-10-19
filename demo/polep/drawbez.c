
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2005                                  */
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

#include "pkvaria.h"
#include "pkgeom.h"
#include "camera.h"
#include "multibs.h" 

#include "xgedit.h"
#include "drawbez.h"


static vector3f pplanes[4];

static CameraRecf *_CPos;


static void SetForeground ( int i, int dd, int colour1, int colour2 )
{
  if ( i == 0 || i == dd )
    xgeSetForeground ( colour1 );
  else if ( i == 1 )
    xgeSetForeground ( colour2 );
} /*SetForeground*/

/* ///////////////////////////////////////////////////////////////////////// */
static boolean DrawArcBC3a ( int degree, const point3f *cpoints )
{
  void    *sp;
  int     i;
  point3f *cp;
  boolean result;

  sp = pkv_GetScratchMemTop ();
  cp = pkv_GetScratchMem ( (degree+1)*sizeof(point3f) );
  if ( !cp )
    exit ( 1 );

  for ( i = 0; i <= degree; i++ ) {
    CameraProjectPoint3f ( _CPos, &cpoints[i], &cp[i] );
    cp[i].x *= cp[i].z;  cp[i].y *= cp[i].z;
  }
  result = xge_DrawBC2Rf ( degree, cp );
  pkv_SetScratchMemTop ( sp );
  return result;
} /*DrawArcBC3a*/

void DrawBezPatch3a ( CameraRecf *CPos,
                      int n, int m, const point3f *p,
                      int dd1, int dd2, int colour1, int colour2 )
{
  void    *sp;
  int     i, k;
  point3f *cp;
  float   t;

  sp = pkv_GetScratchMemTop ( );
  k = max ( n, m ) + 1;
  cp = (point3f*)pkv_GetScratchMem ( k*sizeof(point3f) );
  if ( !cp )
    exit ( 1 );

  _CPos = CPos;
  for ( i = 0; i <= dd1; i++ ) {
    SetForeground ( i, dd1, colour1, colour2 );
    t = (float)i/(float)dd1;
    mbs_multiBCHornerf ( n, 1, 3*(m+1), 0, (float*)p, t, (float*)cp );
    mbs_ClipBC3f ( 4, CPos->cplane, m, cp, DrawArcBC3a );
  }

  for ( i = 0; i <= dd2; i++ ) {
    SetForeground ( i, dd2, colour1, colour2 );
    t = (float)i/(float)dd2;
    mbs_multiBCHornerf ( m, n+1, 3, 3*(m+1), (float*)p, t, (float*)cp );
    mbs_ClipBC3f ( 4, CPos->cplane, n, cp, DrawArcBC3a );
  }

  pkv_SetScratchMemTop ( sp );
} /*DrawBezPatch3a*/

void DrawBezCurve3a ( CameraRecf *CPos,
                      int n, const point3f *p, int colour )
{
  _CPos = CPos;
  xgeSetForeground ( colour );
  mbs_ClipBC3f ( 4, CPos->cplane, n, p, DrawArcBC3a );
} /*DrawBezCurve3a*/

/* ///////////////////////////////////////////////////////////////////////// */
static void GetPClipRect ( CameraRecf *PPos )
{
  SetVector3f ( &pplanes[0], 0.0, -1.0, (float)(PPos->ymin+PPos->height) );
  SetVector3f ( &pplanes[1], 0.0, 1.0, (float)PPos->ymin );
  SetVector3f ( &pplanes[2], -1.0, 0.0, (float)(PPos->xmin+PPos->width) );
  SetVector3f ( &pplanes[3], 1.0, 0.0, (float)PPos->xmin );
} /*GetClipRect*/

void DrawBezPatch3b ( CameraRecf *PPos,
                      int n, int m, const point3f *p,
                      int dd1, int dd2, int colour1, int colour2 )
{
  void    *sp;
  int     i, k;
  point2f *cp, *pp;
  point3f q;
  float   t;

  sp = pkv_GetScratchMemTop ( );
  k = max ( n, m ) + 1;
  cp = (point2f*)pkv_GetScratchMem ( 2*k*sizeof(point2f) );
  pp = (point2f*)pkv_GetScratchMem ( (n+1)*(m+1)*sizeof(point2f) );
  if ( !cp || !pp )
    exit ( 1 );

  GetPClipRect ( PPos );

  for ( i = 0; i < (n+1)*(m+1); i++ ) {
    CameraProjectPoint3f ( PPos, &p[i], &q );
    SetPoint2f ( &pp[i], q.x, q.y );
  }

  for ( i = 0; i <= dd1; i++ ) {
    SetForeground ( i, dd1, colour1, colour2 );
    t = (float)i/(float)dd1;
    mbs_multiBCHornerf ( n, 1, 2*(m+1), 0, (float*)pp, t, (float*)cp );
    mbs_ClipBC2f ( 4, pplanes, m, cp, xge_DrawBC2f );
  }

  for ( i = 0; i <= dd2; i++ ) {
    SetForeground ( i, dd2, colour1, colour2 );
    t = (float)i/(float)dd2;
    mbs_multiBCHornerf ( m, n+1, 2, 2*(m+1), (float*)pp, t, (float*)cp );
    mbs_ClipBC2f ( 4, pplanes, n, cp, xge_DrawBC2f );
  }

  pkv_SetScratchMemTop ( sp );
} /*DrawBezPatch3b*/

void DrawBezCurve3b ( CameraRecf *PPos,
                      int n, const point3f *p, int colour )
{
  int     i;
  point2f *cp;
  point3f q;
  void    *sp;

  sp = pkv_GetScratchMemTop ();
  cp = (point2f*)pkv_GetScratchMem ( (n+1)*sizeof(point2f) );
  if ( !cp )
    exit ( 1 );

  GetPClipRect ( PPos );

  xgeSetForeground ( colour );
  for ( i = 0; i <= n; i++ ) {
    CameraProjectPoint3f ( PPos, &p[i], &q );
    SetPoint2f ( &cp[i], q.x, q.y );
  }
  mbs_ClipBC2f ( 4, pplanes, n, cp, xge_DrawBC2f );

  pkv_SetScratchMemTop ( sp );
} /*DrawBezCurve3b*/

/* ///////////////////////////////////////////////////////////////////////// */
void DrawBezPatch2 ( CameraRecf *PPos,
                     int n, int m, const point2f *p,
                     int dd1, int dd2, int colour1, int colour2 )
{
  void    *sp;
  int     i, k;
  point2f *pp, *cp;
  float   t;

  sp = pkv_GetScratchMemTop ( );
  k = max ( n, m ) + 1;
  cp = (point2f*)pkv_GetScratchMem ( k*sizeof(point2f) );
  pp = (point2f*)pkv_GetScratchMem ( (n+1)*(m+1)*sizeof(point2f) );
  if ( !cp || !pp )
    exit ( 1 );

  GetPClipRect ( PPos );

  for ( i = 0; i < (n+1)*(m+1); i++ )
    CameraProjectPoint2f ( PPos, &p[i], &pp[i] );

  for ( i = 0; i <= dd1; i++ ) {
    SetForeground ( i, dd1, colour1, colour2 );
    t = (float)i/(float)dd1;
    mbs_multiBCHornerf ( n, 1, 2*(m+1), 0, (float*)pp, t, (float*)cp );
    mbs_ClipBC2f ( 4, pplanes, m, cp, xge_DrawBC2f );
  }

  for ( i = 0; i <= dd2; i++ ) {
    SetForeground ( i, dd2, colour1, colour2 );
    t = (float)i/(float)dd2;
    mbs_multiBCHornerf ( m, n+1, 2, 2*(m+1), (float*)pp, t, (float*)cp );
    mbs_ClipBC2f ( 4, pplanes, n, cp, xge_DrawBC2f );
  }

  pkv_SetScratchMemTop ( sp );
} /*DrawBezPatch2*/

void DrawBezPatch2a ( CameraRecf *PPos,
                      int n, int m, const point2f *p,
                      float u0, float u1, float v0, float v1,
                      int dd1, int dd2, int colour1, int colour2 )
{
  void    *sp;
  int     i, k;
  point2f *pp, *cp, *cq;
  float   t;

  sp = pkv_GetScratchMemTop ();
  k = max ( n, m ) + 1;
  cp = (point2f*)pkv_GetScratchMem ( 2*k*sizeof(point2f) );
  pp = (point2f*)pkv_GetScratchMem ( (n+1)*(m+1)*sizeof(point2f) );
  if ( !cp || !pp )
    exit ( 1 );
  cq = &cp[k];

  GetPClipRect ( PPos );

  for ( i = 0; i < (n+1)*(m+1); i++ )
    CameraProjectPoint2f ( PPos, &p[i], &pp[i] );

  for ( i = 0; i <= dd1; i++ ) {
    t = (float)i/(float)dd1;
    if ( t >= u0 && t <= u1 ) {
      SetForeground ( i, dd1, colour1, colour2 );
      mbs_multiBCHornerf ( n, 1, 2*(m+1), 0, (float*)pp, t, (float*)cp );
      mbs_DivideBC2f ( m, v0, cp, cq );
      mbs_DivideBC2f ( m, v1, cp, cq );
      mbs_ClipBC2f ( 4, pplanes, m, cq, xge_DrawBC2f );
    }
  }

  for ( i = 0; i <= dd2; i++ ) {
    t = (float)i/(float)dd2;
    if ( t >= v0 && t <= v1 ) {
      SetForeground ( i, dd2, colour1, colour2 );
      mbs_multiBCHornerf ( m, n+1, 2, 2*(m+1), (float*)pp, t, (float*)cp );
      mbs_DivideBC2f ( m, u0, cp, cq );
      mbs_DivideBC2f ( m, u1, cp, cq );
      mbs_ClipBC2f ( 4, pplanes, n, cq, xge_DrawBC2f );
    }
  }
  pkv_SetScratchMemTop ( sp );
} /*DrawBezPatch2a*/

