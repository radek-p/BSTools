
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#include "pkvaria.h"
#include "pknum.h"
#include "pkgeom.h"
#include "multibs.h"
#include "camera.h"
#include "psout.h"

#include "eg2holef.h"
#include "g2ps.h"

void DrawBSNet ( int k, point3f *Bpt, CameraRecf *CPos, int i, boolean last )
{
  int     ind[16];
  point3f q[16];
  int     j, l;

  ps_Set_Gray ( 0.0 );
  ps_Set_Line_Width ( 2.0 );
  gh_GetBspInd ( k, i, 0, ind );
  for ( j = 0; j < 16; j++ )
    CameraProjectPoint3f ( CPos, &Bpt[ind[j]], &q[j] );
  for ( j = 0; j <= 3; j++ )
    for ( l = 0; l < 3; l++ )
      ps_Draw_Line ( q[4*j+l].x, q[4*j+l].y, q[4*j+l+1].x, q[4*j+l+1].y );
  for ( j = 0; j < 3; j++ )
    for ( l = 0; l <= 3; l++ )
      ps_Draw_Line ( q[4*j+l].x, q[4*j+l].y, q[4*(j+1)+l].x, q[4*(j+1)+l].y );
  for ( j = 0; j <= 3; j++ )
    for ( l = 0; l <= 3; l++ )
      if ( j == 3 || l == 0 )
        ps_Mark_Circle ( q[4*j+l].x, q[4*j+l].y );
      else
        ps_Fill_Circle ( q[4*j+l].x, q[4*j+l].y, 12 );
  if ( last ) {
    i = (i+1) % k;
    gh_GetBspInd ( k, i, 0, ind );
    for ( j = 0; j < 16; j++ )
      CameraProjectPoint3f ( CPos, &Bpt[ind[j]], &q[j] );
    for ( j = 1; j <= 3; j++ )
      for ( l = 0; l < 3; l++ )
        ps_Draw_Line ( q[4*j+l].x, q[4*j+l].y, q[4*j+l+1].x, q[4*j+l+1].y );
    for ( j = 0; j < 3; j++ )
      for ( l = 0; l <= 3; l++ )
        ps_Draw_Line ( q[4*j+l].x, q[4*j+l].y, q[4*(j+1)+l].x, q[4*(j+1)+l].y );
    for ( j = 0; j <= 3; j++ )
      for ( l = 0; l <= 3; l++ )
        if ( j == 3 || l == 0 )
          ps_Mark_Circle ( q[4*j+l].x, q[4*j+l].y );
        else
          ps_Fill_Circle ( q[4*j+l].x, q[4*j+l].y, 12 );
  }
} /*DrawBSNet*/

static void DrawBezCurve3 ( CameraRecf *CPos, int n, point3f *cp )
{
#define DD 32
  void   *sp;
  int    i;
  point3f p, q;
  point2f *cc;

  sp = pkv_GetScratchMemTop ();
  cc = (point2f*)pkv_GetScratchMem ( (DD+1)*sizeof(point2f) );
  if ( !cc )
    exit ( 1 );

  for ( i = 0; i <= DD; i++ ) {
    mbs_BCHornerC3f ( n, cp, (float)i/(float)DD, &p );
    CameraProjectPoint3f ( CPos, &p, &q );
    memcpy ( &cc[i], &q, sizeof(point2f) );
  }
  ps_Draw_Polyline2f ( DD+1, cc );

  pkv_SetScratchMemTop ( sp );
#undef DD
} /*DrawBezCurve3*/

static void DrawBezPatch3 ( CameraRecf *CPos, int n, int m, point3f *cp )
{
#define D 8
  void    *sp;
  int     i;
  point3f *cc;

  sp = pkv_GetScratchMemTop ();
  cc = (point3f*)pkv_GetScratchMem ( (max(n,m)+1)*sizeof(point3f) );
  if ( !cc )
    exit ( 1 );

  for ( i = 0; i <= D; i++ ) {
    mbs_multiBCHornerf ( n, 1, 3*(m+1), 0, (float*)cp, (float)i/(float)D,
                         (float*)cc );
    if ( i == 0 || i == D )
      ps_Set_Line_Width ( 6.0 );
    else
      ps_Set_Line_Width ( 2.0 );
    DrawBezCurve3 ( CPos, m, cc );
  }
  for ( i = 0; i <= D; i++ ) {
    mbs_multiBCHornerf ( m, n+1, 3, 3*(m+1), (float*)cp, (float)i/(float)D,
                         (float*)cc );
    if ( i == 0 || i == D )
      ps_Set_Line_Width ( 6.0 );
    else
      ps_Set_Line_Width ( 2.0 );
    DrawBezCurve3 ( CPos, n, cc );
  }
  pkv_SetScratchMemTop ( sp );
#undef D
} /*DrawBezPatch3*/

void DrawBSPatch ( int k, float *knots, point3f *Bpt, CameraRecf *CPos, int i,
                   boolean first, boolean last )
{
  int     ind[16];
  point3f p[16], q[16];
  float   *ukn, *vkn;
  int     j, l;

  ps_Set_Gray ( 0.7 );
  if ( last ) {
    j = (i+1) % k;
    gh_GetBspInd ( k, j, 0, ind );
    for ( l = 0; l < 16; l++ )
      p[l] = Bpt[ind[l]];
    ukn = &knots[11*((j+k-1) % k)+3];
    vkn = &knots[11*j+0];
    mbs_BSPatchToBezf ( 3, 3, 7, ukn, 3, 7, vkn, 12, (float*)p,
                        NULL, NULL, NULL, NULL, NULL, NULL, 12, (float*)q );
    DrawBezPatch3 ( CPos, 3, 3, q );
  }
  if ( first )
    j = 0;
  else {
    j = 1;
    ps_Set_Gray ( 0.4 );
  }
  for ( ; j < 3; j++ ) {
    gh_GetBspInd ( k, i, j, ind );
    for ( l = 0; l < 16; l++ )
      p[l] = Bpt[ind[l]];
    ukn = &knots[11*((i+k-1) % k)+3];
    vkn = &knots[11*i+j];
    mbs_BSPatchToBezf ( 3, 3, 7, ukn, 3, 7, vkn, 12, (float*)p,
                        NULL, NULL, NULL, NULL, NULL, NULL, 12, (float*)q );
    DrawBezPatch3 ( CPos, 3, 3, q );
    ps_Set_Gray ( 0.4 );
  }
} /*DrawBSPatch*/

/* ////////////////////////////////////////////////////////////////////////// */
void DrawBSDomNet ( int k, point2f *Dompt, CameraRecf *PPos, int i, boolean last )
{
  int     ind[16];
  point2f q[16];
  int     j, l;

  ps_Set_Gray ( 0.0 );
  ps_Set_Line_Width ( 2.0 );
  gh_GetBspInd ( k, i, 0, ind );
  for ( j = 0; j < 16; j++ )
    CameraProjectPoint2f ( PPos, &Dompt[ind[j]], &q[j] );
  for ( j = 0; j <= 3; j++ )
    for ( l = 0; l < 3; l++ )
      ps_Draw_Line ( q[4*j+l].x, q[4*j+l].y, q[4*j+l+1].x, q[4*j+l+1].y );
  for ( j = 0; j < 3; j++ )
    for ( l = 0; l <= 3; l++ )
      ps_Draw_Line ( q[4*j+l].x, q[4*j+l].y, q[4*(j+1)+l].x, q[4*(j+1)+l].y );
  for ( j = 0; j <= 3; j++ )
    for ( l = 0; l <= 3; l++ )
      if ( j == 3 || l == 0 )
        ps_Mark_Circle ( q[4*j+l].x, q[4*j+l].y );
      else
        ps_Fill_Circle ( q[4*j+l].x, q[4*j+l].y, 12 );
  if ( last ) {
    i = (i+1) % k;
    gh_GetBspInd ( k, i, 0, ind );
    for ( j = 0; j < 16; j++ )
      CameraProjectPoint2f ( PPos, &Dompt[ind[j]], &q[j] );
    for ( j = 1; j <= 3; j++ )
      for ( l = 0; l < 3; l++ )
        ps_Draw_Line ( q[4*j+l].x, q[4*j+l].y, q[4*j+l+1].x, q[4*j+l+1].y );
    for ( j = 0; j < 3; j++ )
      for ( l = 0; l <= 3; l++ )
        ps_Draw_Line ( q[4*j+l].x, q[4*j+l].y, q[4*(j+1)+l].x, q[4*(j+1)+l].y );
    for ( j = 0; j <= 3; j++ )
      for ( l = 0; l <= 3; l++ )
        if ( j == 3 || l == 0 )
          ps_Mark_Circle ( q[4*j+l].x, q[4*j+l].y );
        else
          ps_Fill_Circle ( q[4*j+l].x, q[4*j+l].y, 12 );
  }
} /*DrawBSDomNet*/

void DrawBSDomNetNum ( int k, point2f *Dompt, CameraRecf *PPos )
{
  int     i;
  char    s[40];
  point2f q;

  for ( i = 0; i < k; i++ )
    DrawBSDomNet ( k, Dompt, PPos, i, false );
  ps_Write_Command ( "/Times-Roman findfont 50 scalefont setfont" );
  for ( i = 0; i <= 12*k; i++ ) {
    CameraProjectPoint2f ( PPos, &Dompt[i], &q );
    ps_Set_Gray ( 1.0 );
    ps_Fill_Circle ( q.x, q.y, 40.0 );
    ps_Set_Gray ( 0.0 );
    sprintf ( s, "%6.2f %6.2f moveto (%2d) show", q.x-25, q.y-15, i );
    ps_Write_Command ( s );
  }
} /*DrawBSDomNetNum*/

static void DrawBezCurve2 ( CameraRecf *PPos, int n, point2f *cp )
{
#define DD 32
  void   *sp;
  int    i;
  point2f p, *cc;

  sp = pkv_GetScratchMemTop ();
  cc = (point2f*)pkv_GetScratchMem ( (DD+1)*sizeof(point2f) );
  if ( !cc )
    exit ( 1 );

  for ( i = 0; i <= DD; i++ ) {
    mbs_BCHornerC2f ( n, cp, (float)i/(float)DD, &p );
    CameraProjectPoint2f ( PPos, &p, &cc[i] );
  }
  ps_Draw_Polyline2f ( DD+1, cc );

  pkv_SetScratchMemTop ( sp );
#undef DD
} /*DrawBezCurve2*/

static void DrawBezPatch2 ( CameraRecf *PPos, int n, int m, point2f *cp )
{
#define D 8
  void    *sp;
  int     i;
  point2f *cc;

  sp = pkv_GetScratchMemTop ();
  cc = (point2f*)pkv_GetScratchMem ( (max(n,m)+1)*sizeof(point2f) );
  if ( !cc )
    exit ( 1 );

  for ( i = 0; i <= D; i++ ) {
    mbs_multiBCHornerf ( n, 1, 2*(m+1), 0, (float*)cp, (float)i/(float)D,
                         (float*)cc );
    if ( i == 0 || i == D )
      ps_Set_Line_Width ( 6.0 );
    else
      ps_Set_Line_Width ( 2.0 );
    DrawBezCurve2 ( PPos, m, cc );
  }
  for ( i = 0; i <= D; i++ ) {
    mbs_multiBCHornerf ( m, n+1, 2, 2*(m+1), (float*)cp, (float)i/(float)D,
                         (float*)cc );
    if ( i == 0 || i == D )
      ps_Set_Line_Width ( 6.0 );
    else
      ps_Set_Line_Width ( 2.0 );
    DrawBezCurve2 ( PPos, n, cc );
  }
  pkv_SetScratchMemTop ( sp );
#undef D
} /*DrawBezPatch2*/

void DrawBSDomPatch ( int k, float *knots, point2f *Dompt, CameraRecf *PPos, int i,
                      boolean first, boolean last )
{
  int     ind[16];
  point2f p[16], q[16];
  float   *ukn, *vkn;
  int     j, l;

  ps_Set_Gray ( 0.7 );
  if ( last ) {
    j = (i+1) % k;
    gh_GetBspInd ( k, j, 0, ind );
    for ( l = 0; l < 16; l++ )
      p[l] = Dompt[ind[l]];
    ukn = &knots[11*((j+k-1) % k)+3];
    vkn = &knots[11*j+0];
    mbs_BSPatchToBezf ( 2, 3, 7, ukn, 3, 7, vkn, 8, (float*)p,
                        NULL, NULL, NULL, NULL, NULL, NULL, 8, (float*)q );
    DrawBezPatch2 ( PPos, 3, 3, q );
  }
  if ( first )
    j = 0;
  else {
    j = 1;
    ps_Set_Gray ( 0.4 );
  }
  for ( ; j < 3; j++ ) {
    gh_GetBspInd ( k, i, j, ind );
    for ( l = 0; l < 16; l++ )
      p[l] = Dompt[ind[l]];
    ukn = &knots[11*((i+k-1) % k)+3];
    vkn = &knots[11*i+j];
    mbs_BSPatchToBezf ( 2, 3, 7, ukn, 3, 7, vkn, 8, (float*)p,
                        NULL, NULL, NULL, NULL, NULL, NULL, 8, (float*)q );
    DrawBezPatch2 ( PPos, 3, 3, q );
    ps_Set_Gray ( 0.4 );
  }
} /*DrawBSDomPatch*/

/* ////////////////////////////////////////////////////////////////////////// */
void DrawBezCurve2a ( CameraRecf *PPos, int n, point2f *cp,
                      float t0, float t1 )
{
#define DD 32
  void   *sp;
  int    i;
  point2f p, *cc, *q;

  sp = pkv_GetScratchMemTop ();
  q = (point2f*)pkv_GetScratchMem ( (n+1)*sizeof(point2f) );
  cc = (point2f*)pkv_GetScratchMem ( (DD+1)*sizeof(point2f) );
  if ( !q || !cc )
    exit ( 1 );

  mbs_DivideBC2f ( n, t0, cp, q );
  mbs_DivideBC2f ( n, (float)((t1-t0)/(1.0-t0)), cp, q );
  for ( i = 0; i <= DD; i++ ) {
    mbs_BCHornerC2f ( n, q, (float)i/(float)DD, &p );
    CameraProjectPoint2f ( PPos, &p, &cc[i] );
  }
  ps_Draw_Polyline2f ( DD+1, cc );

  pkv_SetScratchMemTop ( sp );
#undef DD
} /*DrawBezCurve2a*/

void DrawBezPatch2a ( CameraRecf *PPos, int n, int m, const point2f *cp,   
                      float u0, float u1, float v0, float v1 )
{
#define D 8
  void    *sp;
  int     i;
  point2f *cc;
  float   t;

  sp = pkv_GetScratchMemTop ();
  cc = (point2f*)pkv_GetScratchMem ( (max(n,m)+1)*sizeof(point2f) );
  if ( !cc )
    exit ( 1 );

  for ( i = 0; i <= D; i++ ) {
    t = (float)i/(float)D;
    if ( u1 < 0.0 ) t = -t;
    mbs_multiBCHornerf ( n, 1, 2*(m+1), 0, (float*)cp, t, (float*)cc );
    if ( i == 0 || i == D )
      ps_Set_Line_Width ( 6.0 );
    else
      ps_Set_Line_Width ( 2.0 );
    if ( (t >= u0 && t <= u1) || (t >= u1 && t <= u0) )
      DrawBezCurve2a ( PPos, m, cc, v0, v1 );
  }
  for ( i = 0; i <= D; i++ ) {
    t = (float)i/(float)D;
    mbs_multiBCHornerf ( m, n+1, 2, 2*(m+1), (float*)cp, t, (float*)cc );
    if ( i == 0 || i == D )
      ps_Set_Line_Width ( 6.0 );
    else
      ps_Set_Line_Width ( 2.0 );
    if ( t >= v0 && t <= v1 )
      DrawBezCurve2a ( PPos, n, cc, u0, u1 );
  }
  pkv_SetScratchMemTop ( sp );
#undef D
} /*DrawBezPatch2a*/

/* ///////////////////////////////////////////////////////////////////////////// */
void DrawBezCurve3a ( CameraRecf *CPos, int n, point3f *cp,
                      float t0, float t1 )
{
#define DD 32
  void   *sp;
  int    i;
  point3f p, *q, s;
  point2f *cc;

  sp = pkv_GetScratchMemTop ();
  q = (point3f*)pkv_GetScratchMem ( (n+1)*sizeof(point3f) );
  cc = (point2f*)pkv_GetScratchMem ( (DD+1)*sizeof(point2f) );
  if ( !q || !cc )
    exit ( 1 );

  mbs_DivideBC3f ( n, t0, cp, q );
  mbs_DivideBC3f ( n, (float)((t1-t0)/(1.0-t0)), cp, q );
  for ( i = 0; i <= DD; i++ ) {
    mbs_BCHornerC3f ( n, q, (float)i/(float)DD, &p );
    CameraProjectPoint3f ( CPos, &p, &s );
    SetPoint2f ( &cc[i], s.x, s.y );
  }
  ps_Draw_Polyline2f ( DD+1, cc );

  pkv_SetScratchMemTop ( sp );
#undef DD
} /*DrawBezCurve3a*/

void DrawBezPatch3a ( CameraRecf *CPos, int n, int m, const point3f *cp,   
                      float u0, float u1, float v0, float v1 )
{
#define D 8
  void    *sp;
  int     i;
  point3f *cc;
  float   t;

  sp = pkv_GetScratchMemTop ();
  cc = (point3f*)pkv_GetScratchMem ( (max(n,m)+1)*sizeof(point3f) );
  if ( !cc )
    exit ( 1 );

  for ( i = 0; i <= D; i++ ) {
    t = (float)i/(float)D;
    if ( u1 < 0.0 ) t = -t;
    mbs_multiBCHornerf ( n, 1, 3*(m+1), 0, (float*)cp, t, (float*)cc );
    if ( i == 0 || i == D )
      ps_Set_Line_Width ( 6.0 );
    else
      ps_Set_Line_Width ( 2.0 );
    if ( (t >= u0 && t <= u1) || (t >= u1 && t <= u0) )
      DrawBezCurve3a ( CPos, m, cc, v0, v1 );
  }
  for ( i = 0; i <= D; i++ ) {
    t = (float)i/(float)D;
    mbs_multiBCHornerf ( m, n+1, 3, 3*(m+1), (float*)cp, t, (float*)cc );
    if ( i == 0 || i == D )
      ps_Set_Line_Width ( 6.0 );
    else
      ps_Set_Line_Width ( 2.0 );
    if ( t >= v0 && t <= v1 )
      DrawBezCurve3a ( CPos, n, cc, u0, u1 );
  }
  pkv_SetScratchMemTop ( sp );
#undef D
} /*DrawBezPatch3a*/

