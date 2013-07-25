
#include <math.h>  
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <sys/times.h>
#include <unistd.h>

#include "pkvaria.h"
#include "pknum.h"  
#include "pkgeom.h"   
#include "camera.h" 
#include "psout.h"  
#include "multibs.h"
#include "eg1holef.h"

#include "bslapf.h"
#include "drawitf.h"

CameraRecf   CPos, PPos;

void DrawBezCurve2f ( int n, const point2f *cp, float t0, float t1 )
{
#define DD 32
  void    *sp;
  int     i;
  point2f *cc, p;

  sp = pkv_GetScratchMemTop ();
  cc = (point2f*)pkv_GetScratchMem ( (DD+1)*sizeof(point2f) );
  if ( !cc )
    exit ( 1 );

  for ( i = 0; i <= DD; i++ ) {
    mbs_BCHornerC2f ( n, cp, (float)i/(float)DD, &p );
    CameraProjectPoint2f ( &PPos, &p, &cc[i] );
  }
  ps_Draw_Polyline2f ( DD+1, cc );

  pkv_SetScratchMemTop ( sp );
#undef DD
} /*DrawBezCurve2f*/

void DrawBezCurve3f ( int n, const point3f *cp, float t0, float t1 )
{
#define DD 32
  void    *sp;
  int     i, d;
  point3f p, q, *r, *s;
  point2f *cc;

  sp = pkv_GetScratchMemTop ();
  cc = (point2f*)pkv_GetScratchMem ( (DD+1)*sizeof(point2f) );
  r = (point3f*)pkv_GetScratchMem ( 2*(n+1)*sizeof(point3f) );
  if ( !cc || !r )
    exit ( 1 );
  s = &r[n+1];
  memcpy ( r, cp, (n+1)*sizeof(point3f) );
  mbs_DivideBC3f ( n, t0, r, s );
  mbs_DivideBC3f ( n, (t1-t0)/(1.0-t0), r, s );
  d = (int)((t1-t0)*DD+0.5);
  for ( i = 0; i <= d; i++ ) {
    mbs_BCHornerC3f ( n, s, (float)i/(float)d, &p );
    CameraProjectPoint3f ( &CPos, &p, &q );
    SetPoint2f ( &cc[i], q.x, q.y );
  }
  ps_Draw_Polyline2f ( d+1, cc );

  pkv_SetScratchMemTop ( sp );
#undef DD
} /*DrawBezCurve3f*/

void DrawBSCurve2f ( int n, int lkn, const float *knots, const point2f *cp )
{
#define DD 32
  void    *sp;
  int     i;
  float   t;
  point2f *cc, p;

  sp = pkv_GetScratchMemTop ();
  cc = (point2f*)pkv_GetScratchMem ( (DD+1)*sizeof(point2f) );
  if ( !cc )
    exit ( 1 );

  for ( i = 0; i <= DD; i++ ) {
    t = knots[n]+(float)i/(float)DD*(knots[lkn-n]-knots[n]);
    mbs_deBoorC2f ( n, lkn, knots, cp, t, &p );
    CameraProjectPoint2f ( &PPos, &p, &cc[i] );
  }
  ps_Draw_Polyline2f ( DD+1, cc );

  pkv_SetScratchMemTop ( sp );
#undef DD
} /*DrawBSCurve2f*/

void DrawBSCurve3f ( int n, int lkn, const float *knots, const point3f *cp )
{
#define DD 32
  void    *sp;
  int     i;
  float   t;
  point2f *cc;
  point3f p, q;

  sp = pkv_GetScratchMemTop ();
  cc = (point2f*)pkv_GetScratchMem ( (DD+1)*sizeof(point2f) );
  if ( !cc )
    exit ( 1 );

  for ( i = 0; i <= DD; i++ ) {
    t = knots[n]+(float)i/(float)DD*(knots[lkn-n]-knots[n]);
    mbs_deBoorC3f ( n, lkn, knots, cp, t, &p );
    CameraProjectPoint3f ( &PPos, &p, &q );
    SetPoint2f ( &cc[i], q.x, q.y );
  }
  ps_Draw_Polyline2f ( DD+1, cc );

  pkv_SetScratchMemTop ( sp );
#undef DD
} /*DrawBSCurve3f*/

void SetLW ( int i, int i0, int i1 )
{
  if ( i == i0 || i == i1 )
    ps_Set_Line_Width ( 4.0 );
  else
    ps_Set_Line_Width ( 2.0 );
} /*SetLW*/

void SetCLW ( int i, int i0, int i1 )
{
  if ( i == i0 || i == i1 ) {
    ps_Set_RGB ( 0.2, 0.0, 0.8 );
    ps_Set_Line_Width ( 1.5 );
  }
  else {
    ps_Set_Gray ( 0.5 );
    ps_Set_Line_Width ( 1.0 );
  }
} /*SetCLW*/

void DrawBezPatch2f ( int n, int m, const point2f *cp, int d0, int d1 )
{
  void    *sp;
  int     i;
  point2f *cc;

  sp = pkv_GetScratchMemTop ();
  cc = (point2f*)pkv_GetScratchMem ( (max(n,m)+1)*sizeof(point2f) );
  if ( !cc )
    exit ( 1 );

  for ( i = 0; i <= d0; i++ ) {
    mbs_multiBCHornerf ( n, 1, 2*(m+1), 0, (float*)cp, (float)i/(float)d0,
                         (float*)cc );
    SetLW ( i, 0, d0 );
    DrawBezCurve2f ( m, cc, 0.0, 1.0 );
  }
  for ( i = 0; i <= d1; i++ ) {
    mbs_multiBCHornerf ( m, n+1, 2, 2*(m+1), (float*)cp, (float)i/(float)d1,
                         (float*)cc );
    SetLW ( i, 0, d1 );
    DrawBezCurve2f ( n, cc, 0.0, 1.0 );
  }

  pkv_SetScratchMemTop ( sp );
} /*DrawBezPatch2f*/

void DrawBezPatch3f ( int n, int m, const point3f *cp, int d0, int d1 )
{
  void    *sp;
  int     i;
  point3f *cc;

  sp = pkv_GetScratchMemTop ();
  cc = (point3f*)pkv_GetScratchMem ( (max(n,m)+1)*sizeof(point3f) );
  if ( !cc )
    exit ( 1 );

  for ( i = 0; i <= d0; i++ ) {
    mbs_multiBCHornerf ( n, 1, 3*(m+1), 0, (float*)cp, (float)i/(float)d0,
                         (float*)cc );
    SetLW ( i, 0, d0 );
    DrawBezCurve3f ( m, cc, 0.0, 1.0 );
  }
  for ( i = 0; i <= d1; i++ ) {
    mbs_multiBCHornerf ( m, n+1, 3, 3*(m+1), (float*)cp, (float)i/(float)d1,
                         (float*)cc );
    SetLW ( i, 0, d1 );
    DrawBezCurve3f ( n, cc, 0.0, 1.0 );
  }

  pkv_SetScratchMemTop ( sp );
} /*DrawBezPatch3f*/

void DrawBezPatch3Rf ( int n, int m, const point3f *cp,
                       float u0, float u1, float v0, float v1,
                       int d0, int d1 )
{
  void    *sp;
  int     i;
  float   t;
  point3f *q;

  if ( u0 > u1 || v0 > v1 || u0 == 1.0 || v0 == 1.0 )
    return;
  sp = pkv_GetScratchMemTop ();
  i = max ( n, m ) + 1;
  q = (point3f*)pkv_GetScratchMem ( i*sizeof(point3f) );
  if ( !q )
    exit ( 1 );
  for ( i = 0; i <= d0; i++ ) {
    t = (float)i/(float)d0;
    if ( t >= u0 && t <= u1 ) {
      SetLW ( i, 0, d0 );
      mbs_multiBCHornerf ( n, 1, 3*(m+1), 0, (float*)cp, t, (float*)q );
      DrawBezCurve3f ( m, q, v0, v1 );
    }
  }
  for ( i = 0; i <= d1; i++ ) {
    t = (float)i/(float)d1;
    if ( t >= v0 && t <= v1 ) {
      SetLW ( i, 0, d1 );
      mbs_multiBCHornerf ( m, n+1, 3, 3*(m+1), (float*)cp, t, (float*)q );
      DrawBezCurve3f ( n, q, u0, u1 );
    }
  }
  pkv_SetScratchMemTop ( sp );
} /*DrawBezPatch3Rf*/

void DrawBezPatchNet2f ( int n, int m, const point2f *cp )
{
  void    *sp;
  int     i, j;
  point2f *cc;

  sp = pkv_GetScratchMemTop ();
  cc = (point2f*)pkv_GetScratchMem ( (max(n,m)+1)*sizeof(point2f) );
  ps_Set_Line_Width ( 2.0 );
  for ( i = 0; i <= n; i++ ) {
    for ( j = 0; j <= m; j++ )
      CameraProjectPoint2f ( &PPos, &cp[i*(m+1)+j], &cc[j] );
    ps_Draw_Polyline2f ( m+1, cc );
  }
  for ( j = 0; j <= m; j++ ) {
    for ( i = 0; i <= n; i++ )
      CameraProjectPoint2f ( &PPos, &cp[i*(m+1)+j], &cc[i] );
    ps_Draw_Polyline2f ( n+1, cc );
  }

  pkv_SetScratchMemTop ( sp );
} /*DrawBezPatchNet2f*/

void DrawBezPatchNet3f ( int n, int m, const point3f *cp )
{
  void    *sp;
  int     i, j;
  point3f q;
  point2f *r;

  sp = pkv_GetScratchMemTop ();
  r = pkv_GetScratchMem ( (n+1)*(m+1)*sizeof(point2f) );
  if ( !r )
    exit ( 1 );
  for ( i = 0; i < (n+1)*(m+1); i++ ) {
    CameraProjectPoint3f ( &CPos, &cp[i], &q );
    SetPoint2f ( &r[i], q.x, q.y );
  }
  ps_Set_Line_Width ( 1.0 );
  for ( i = 0; i <= n; i++ )
    ps_Draw_Polyline2f ( m+1, &r[i*(m+1)] );
  for ( j = 0; j <= m; j++ )
    for ( i = 0; i < n; i++ )
      ps_Draw_Line ( r[i*(m+1)+j].x, r[i*(m+1)+j].y,
                     r[(i+1)*(m+1)+j].x, r[(i+1)*(m+1)+j].y );
  for ( i = 0; i < (n+1)*(m+1); i++ )
    ps_Mark_Circle ( r[i].x, r[i].y );

  pkv_SetScratchMemTop ( sp );
} /*DrawBezPatchNet3f*/

void DrawBSPatch2f ( int n, int lknu, const float *knu,
                     int m, int lknv, const float *knv, const point2f *cp,
                     int d0, int d1 )
{
  void    *sp;
  int     ku, kv, pitch1, pitch2, pitch3;
  float   *b, *c;
  int     i, j, start;

  sp = pkv_GetScratchMemTop ();
  ku = mbs_NumKnotIntervalsf ( n, lknu, knu );
  kv = mbs_NumKnotIntervalsf ( m, lknv, knv );
  pitch1 = (lknv-m)*2;
  pitch2 = (m+1)*2*kv;
  pitch3 = (m+1)*2;
  b = pkv_GetScratchMemf ( pitch2*ku*(n+1) );
  c = pkv_GetScratchMemf ( (n+1)*pitch3 );
  if ( b && c ) {
    mbs_BSPatchToBezf ( 2, n, lknu, knu, m, lknv, knv, pitch1, (float*)cp,
                        &ku, NULL, NULL, &kv, NULL, NULL, pitch2, b );
    for ( i = 0; i < ku; i++ )
      for ( j = 0; j < kv; j++ ) {
        start = i*(n+1)*pitch2+j*pitch3;
        pkv_Selectf ( n+1, pitch3, pitch2, pitch3, &b[start], c );
        DrawBezPatch2f ( n, m, (point2f*)c, d0, d1 );
      }
  }

  pkv_SetScratchMemTop ( sp );
} /*DrawBSPatch2f*/

void DrawBSPatch3f ( int n, int lknu, const float *knu,
                     int m, int lknv, const float *knv, const point3f *cp,
                     int d0, int d1 )
{
  void    *sp;
  int     ku, kv, pitch1, pitch2, pitch3;
  float   *b, *c;
  int     i, j, start;

  sp = pkv_GetScratchMemTop ();
  ku = mbs_NumKnotIntervalsf ( n, lknu, knu );
  kv = mbs_NumKnotIntervalsf ( m, lknv, knv );
  pitch1 = (lknv-m)*3;
  pitch2 = (m+1)*3*kv;
  pitch3 = (m+1)*3;
  b = pkv_GetScratchMemf ( pitch2*ku*(n+1) );
  c = pkv_GetScratchMemf ( (n+1)*pitch3 );
  if ( b && c ) {
    mbs_BSPatchToBezf ( 3, n, lknu, knu, m, lknv, knv, pitch1, (float*)cp,
                        &ku, NULL, NULL, &kv, NULL, NULL, pitch2, b );
    for ( i = 0; i < ku; i++ )
      for ( j = 0; j < kv; j++ ) {
        start = i*(n+1)*pitch2+j*pitch3;
        pkv_Selectf ( n+1, pitch3, pitch2, pitch3, &b[start], c );
        DrawBezPatch3f ( n, m, (point3f*)c, d0, d1 );
      }
  }

  pkv_SetScratchMemTop ( sp );
} /*DrawBSPatch3f*/

void DrawBSPatch3Rf ( int n, int lknu, const float *knu,
                      int m, int lknv, const float *knv, const point3f *cp,
                      float u0, float u1, float v0, float v1,
                      int d0, int d1 )
{
  void    *sp;
  int     ku, kv, pitch1, pitch2, pitch3;
  float   *b, *c, *bknu, *bknv;
  float   ku0, ku1, kv0, kv1, uu0, uu1, vv0, vv1;
  int     i, j, start;

  if ( u0 > u1 || v0 > v1 )
    return;
  sp = pkv_GetScratchMemTop ();
  ku = mbs_NumKnotIntervalsf ( n, lknu, knu );
  kv = mbs_NumKnotIntervalsf ( m, lknv, knv );
  pitch1 = (lknv-m)*3;
  pitch2 = (m+1)*3*kv;
  pitch3 = (m+1)*3;
  b = pkv_GetScratchMemf ( pitch2*ku*(n+1) );
  c = pkv_GetScratchMemf ( (n+1)*pitch3 );
  bknu = pkv_GetScratchMemf ( (ku+1)*(n+1) );
  bknv = pkv_GetScratchMemf ( (kv+1)*(m+1) );
  if ( b && c && bknu && bknv ) {
    mbs_BSPatchToBezf ( 3, n, lknu, knu, m, lknv, knv, pitch1, (float*)cp,
                        &ku, &i, bknu, &kv, &j, bknv, pitch2, b );
    for ( i = 0; i < ku; i++ ) {
      ku0 = bknu[i*n+i+n];
      ku1 = bknu[i*n+i+n+1];
      if ( u0 < ku1 && u1 > ku0 ) {
        uu0 = max ( (u0-ku0)/(ku1-ku0), 0.0 );  uu0 = min ( uu0, 1.0 );
        uu1 = max ( (u1-ku0)/(ku1-ku0), 0.0 );  uu1 = min ( uu1, 1.0 );
        for ( j = 0; j < kv; j++ ) {
          kv0 = bknv[j*m+j+m];
          kv1 = bknv[j*m+j+m+1];
          if ( v0 < kv1 && v1 > kv0 ) {
            vv0 = max ( (v0-kv0)/(kv1-kv0), 0.0 );  vv0 = min ( vv0, 1.0 );
            vv1 = max ( (v1-kv0)/(kv1-kv0), 0.0 );  vv1 = min ( vv1, 1.0 );
            start = i*(n+1)*pitch2+j*pitch3;
            pkv_Selectf ( n+1, pitch3, pitch2, pitch3, &b[start], c );
            DrawBezPatch3Rf ( n, m, (point3f*)c,
                              uu0, uu1, vv0, vv1, d0, d1 );
          }
        }
      }
    }
  }

  pkv_SetScratchMemTop ( sp );
} /*DrawBSPatch3Rf*/

/* ////////////////////////////////////////////////////////////////////////// */
float FindBSPatchMaxLaplacianf ( int n, int lknu, const float *knu,
                                 int m, int lknv, const float *knv,
                                 const point3f *cp, int d0, int d1 )
{
  void    *sp;
  int     ku, kv, pitch1, pitch2, pitch3;
  float   *b, *c;
  int     i, j, k, l, start;
  float   lap, maxval;

  sp = pkv_GetScratchMemTop ();
  ku = mbs_NumKnotIntervalsf ( n, lknu, knu );
  kv = mbs_NumKnotIntervalsf ( m, lknv, knv );
  pitch1 = (lknv-m)*3;
  pitch2 = (m+1)*3*kv;
  pitch3 = (m+1)*3;
  b = pkv_GetScratchMemf ( pitch2*ku*(n+1) );
  c = pkv_GetScratchMemf ( (n+1)*pitch3 );
  maxval = 0.0;
  if ( b && c ) {
    mbs_BSPatchToBezf ( 3, n, lknu, knu, m, lknv, knv, pitch1, (float*)cp,
                        &ku, NULL, NULL, &kv, NULL, NULL, pitch2, b );
    for ( i = 0; i < ku; i++ )
      for ( j = 0; j < kv; j++ ) {
        start = i*(n+1)*pitch2+j*pitch3;
        pkv_Selectf ( n+1, pitch3, pitch2, pitch3, &b[start], c );
        for ( k = 0; k <= d0; k++ )
          for ( l = 0; l <= d1; l++ ) {
            lap = ComputeBezLaplacianf ( n, m, (point3f*)c,
                        (float)k/(float)d0, (float)l/(float)d1 );
            maxval = max ( maxval, fabs(lap) );
          }
      }
  }

  pkv_SetScratchMemTop ( sp );
  return maxval;
} /*FindBSPatchMaxLaplacianf*/

void DrawBSPatchLaplacianf ( int n, int lknu, const float *knu,
                             int m, int lknv, const float *knv,
                             const point3f *cp, int d0, int d1 )
{
#define DD 32
  void    *sp;
  int     ku, kv, pitch1, pitch2, pitch3;
  float   *b, *c, u, v;
  int     i, j, k, l, start;
  point3f p, q;
  point2f *cc;

  sp = pkv_GetScratchMemTop ();
  ku = mbs_NumKnotIntervalsf ( n, lknu, knu );
  kv = mbs_NumKnotIntervalsf ( m, lknv, knv );
  pitch1 = (lknv-m)*3;
  pitch2 = (m+1)*3*kv;
  pitch3 = (m+1)*3;
  b = pkv_GetScratchMemf ( pitch2*ku*(n+1) );
  c = pkv_GetScratchMemf ( (n+1)*pitch3 );
  cc = pkv_GetScratchMem ( (DD+1)*sizeof(point2f) );
  if ( b && c && cc ) {
    mbs_BSPatchToBezf ( 3, n, lknu, knu, m, lknv, knv, pitch1, (float*)cp,
                        &ku, NULL, NULL, &kv, NULL, NULL, pitch2, b );
    for ( i = 0; i < ku; i++ )
      for ( j = 0; j < kv; j++ ) {
        start = i*(n+1)*pitch2+j*pitch3;
        pkv_Selectf ( n+1, pitch3, pitch2, pitch3, &b[start], c );
        for ( k = 0; k <= d0; k++ ) {
          u = (float)k/(float)d0;
          for ( l = 0; l <= DD; l++ ) {
            v = (float)l/(float)DD;
            mbs_BCHornerP3f ( n, m, (point3f*)c, u, v, &p );
            p.z = ComputeBezLaplacianf ( n, m, (point3f*)c, u, v );
            CameraProjectPoint3f ( &CPos, &p, &q );
            SetPoint2f ( &cc[l], q.x, q.y );
          }
          SetCLW ( k, 0, d0 );
          ps_Draw_Polyline2f ( DD+1, cc );
        }
        for ( l = 0; l <= d1; l++ ) {
          v = (float)l/(float)d1;
          for ( k = 0; k <= DD; k++ ) {
            u = (float)k/(float)DD;
            mbs_BCHornerP3f ( n, m, (point3f*)c, u, v, &p );
            p.z = ComputeBezLaplacianf ( n, m, (point3f*)c, u, v );
            CameraProjectPoint3f ( &CPos, &p, &q );
            SetPoint2f ( &cc[k], q.x, q.y );
          }
          SetCLW ( l, 0, d1 );
          ps_Draw_Polyline2f ( DD+1, cc );
        }
      }
  }

  pkv_SetScratchMemTop ( sp );
#undef DD
} /*DrawBSPatchLaplacianf*/

/* ////////////////////////////////////////////////////////////////////////// */
float FindBezPatchMaxLaplacianf ( int n, int m, const point3f *cp, int d0, int d1 )
{
  int     k, l;
  float   lap, maxval;

  maxval = 0.0;
  for ( k = 0; k <= d0; k++ )
    for ( l = 0; l <= d1; l++ ) {
      lap = ComputeBezLaplacianf ( n, m, cp,
                      (float)k/(float)d0, (float)l/(float)d1 );
      maxval = max ( maxval, fabs(lap) );
    }
  return maxval;
} /*FindBezPatchMaxLaplacianf*/

void DrawBezPatchLaplacianf ( int n, int m, const point3f *cp, int d0, int d1 )
{
#define DD 32
  void    *sp;
  float   u, v;
  int     k, l;
  point3f p, q;
  point2f *cc;

  sp = pkv_GetScratchMemTop ();
  cc = pkv_GetScratchMem ( (DD+1)*sizeof(point2f) );
  if ( cc ) {
    for ( k = 0; k <= d0; k++ ) {
      u = (float)k/(float)d0;
      for ( l = 0; l <= DD; l++ ) {
        v = (float)l/(float)DD;
        mbs_BCHornerP3f ( n, m, cp, u, v, &p );
        p.z = ComputeBezLaplacianf ( n, m, cp, u, v );
        CameraProjectPoint3f ( &CPos, &p, &q );
        SetPoint2f ( &cc[l], q.x, q.y );
      }
      SetCLW ( k, 0, d0 );
      ps_Draw_Polyline2f ( DD+1, cc );
    }
    for ( l = 0; l <= d1; l++ ) {
      v = (float)l/(float)d1;
      for ( k = 0; k <= DD; k++ ) {
        u = (float)k/(float)DD;
        mbs_BCHornerP3f ( n, m, cp, u, v, &p );
        p.z = ComputeBezLaplacianf ( n, m, cp, u, v );
        CameraProjectPoint3f ( &CPos, &p, &q );
        SetPoint2f ( &cc[k], q.x, q.y );
      }
      SetCLW ( l, 0, d1 );
      ps_Draw_Polyline2f ( DD+1, cc );
    }
  }
  pkv_SetScratchMemTop ( sp );
#undef DD
} /*DrawBezPatchLaplacianf*/

