
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
#include "eg2holed.h"

#include "bslapd.h"
#include "drawitd.h"

CameraRecd   CPos, PPos;

void DrawBezCurve2d ( int n, const point2d *cp, double t0, double t1 )
{
#define DD 32
  void    *sp;
  int     i;
  point2d *cc, p;

  sp = pkv_GetScratchMemTop ();
  cc = (point2d*)pkv_GetScratchMem ( (DD+1)*sizeof(point2d) );
  if ( !cc )
    exit ( 1 );

  for ( i = 0; i <= DD; i++ ) {
    mbs_BCHornerC2d ( n, cp, (double)i/(double)DD, &p );
    CameraProjectPoint2d ( &PPos, &p, &cc[i] );
  }
  ps_Draw_Polyline2d ( DD+1, cc );

  pkv_SetScratchMemTop ( sp );
#undef DD
} /*DrawBezCurve2d*/

void DrawBezCurve3d ( int n, const point3d *cp, double t0, double t1 )
{
#define DD 32
  void    *sp;
  int     i, d;
  point3d p, q, *r, *s;
  point2d *cc;

  sp = pkv_GetScratchMemTop ();
  cc = (point2d*)pkv_GetScratchMem ( (DD+1)*sizeof(point2d) );
  r = (point3d*)pkv_GetScratchMem ( 2*(n+1)*sizeof(point3d) );
  if ( !cc || !r )
    exit ( 1 );
  s = &r[n+1];
  memcpy ( r, cp, (n+1)*sizeof(point3d) );
  mbs_DivideBC3d ( n, t0, r, s );
  mbs_DivideBC3d ( n, (t1-t0)/(1.0-t0), r, s );
  d = (int)((t1-t0)*DD+0.5);
  for ( i = 0; i <= d; i++ ) {
    mbs_BCHornerC3d ( n, s, (double)i/(double)d, &p );
    CameraProjectPoint3d ( &CPos, &p, &q );
    SetPoint2d ( &cc[i], q.x, q.y );
  }
  ps_Draw_Polyline2d ( d+1, cc );

  pkv_SetScratchMemTop ( sp );
#undef DD
} /*DrawBezCurve3d*/

void DrawBSCurve2d ( int n, int lkn, const double *knots, const point2d *cp )
{
#define DD 32
  void    *sp;
  int     i;
  double   t;
  point2d *cc, p;

  sp = pkv_GetScratchMemTop ();
  cc = (point2d*)pkv_GetScratchMem ( (DD+1)*sizeof(point2d) );
  if ( !cc )
    exit ( 1 );

  for ( i = 0; i <= DD; i++ ) {
    t = knots[n]+(double)i/(double)DD*(knots[lkn-n]-knots[n]);
    mbs_deBoorC2d ( n, lkn, knots, cp, t, &p );
    CameraProjectPoint2d ( &PPos, &p, &cc[i] );
  }
  ps_Draw_Polyline2d ( DD+1, cc );

  pkv_SetScratchMemTop ( sp );
#undef DD
} /*DrawBSCurve2d*/

void DrawBSCurve3d ( int n, int lkn, const double *knots, const point3d *cp )
{
#define DD 32
  void    *sp;
  int     i;
  double   t;
  point2d *cc;
  point3d p, q;

  sp = pkv_GetScratchMemTop ();
  cc = (point2d*)pkv_GetScratchMem ( (DD+1)*sizeof(point2d) );
  if ( !cc )
    exit ( 1 );

  for ( i = 0; i <= DD; i++ ) {
    t = knots[n]+(double)i/(double)DD*(knots[lkn-n]-knots[n]);
    mbs_deBoorC3d ( n, lkn, knots, cp, t, &p );
    CameraProjectPoint3d ( &PPos, &p, &q );
    SetPoint2d ( &cc[i], q.x, q.y );
  }
  ps_Draw_Polyline2d ( DD+1, cc );

  pkv_SetScratchMemTop ( sp );
#undef DD
} /*DrawBSCurve3d*/

void SetLW ( int i, int i0, int i1 )
{
  if ( i == i0 || i == i1 )
    ps_Set_Line_Width ( 1.5 );
  else
    ps_Set_Line_Width ( 1.0 );
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

void DrawBezPatch2d ( int n, int m, const point2d *cp, int d0, int d1 )
{
  void    *sp;
  int     i;
  point2d *cc;

  sp = pkv_GetScratchMemTop ();
  cc = (point2d*)pkv_GetScratchMem ( (max(n,m)+1)*sizeof(point2d) );
  if ( !cc )
    exit ( 1 );

  for ( i = 0; i <= d0; i++ ) {
    mbs_multiBCHornerd ( n, 1, 2*(m+1), 0, (double*)cp, (double)i/(double)d0,
                         (double*)cc );
    SetLW ( i, 0, d0 );
    DrawBezCurve2d ( m, cc, 0.0, 1.0 );
  }
  for ( i = 0; i <= d1; i++ ) {
    mbs_multiBCHornerd ( m, n+1, 2, 2*(m+1), (double*)cp, (double)i/(double)d1,
                         (double*)cc );
    SetLW ( i, 0, d1 );
    DrawBezCurve2d ( n, cc, 0.0, 1.0 );
  }

  pkv_SetScratchMemTop ( sp );
} /*DrawBezPatch2d*/

void DrawBezPatch3d ( int n, int m, const point3d *cp, int d0, int d1 )
{
  void    *sp;
  int     i;
  point3d *cc;

  sp = pkv_GetScratchMemTop ();
  cc = (point3d*)pkv_GetScratchMem ( (max(n,m)+1)*sizeof(point3d) );
  if ( !cc )
    exit ( 1 );

  for ( i = 0; i <= d0; i++ ) {
    mbs_multiBCHornerd ( n, 1, 3*(m+1), 0, (double*)cp, (double)i/(double)d0,
                         (double*)cc );
    SetLW ( i, 0, d0 );
    DrawBezCurve3d ( m, cc, 0.0, 1.0 );
  }
  for ( i = 0; i <= d1; i++ ) {
    mbs_multiBCHornerd ( m, n+1, 3, 3*(m+1), (double*)cp, (double)i/(double)d1,
                         (double*)cc );
    SetLW ( i, 0, d1 );
    DrawBezCurve3d ( n, cc, 0.0, 1.0 );
  }

  pkv_SetScratchMemTop ( sp );
} /*DrawBezPatch3d*/

void DrawBezPatch3Rd ( int n, int m, const point3d *cp,
                       double u0, double u1, double v0, double v1,
                       int d0, int d1 )
{
  void    *sp;
  int     i;
  double   t;
  point3d *q;

  if ( u0 > u1 || v0 > v1 || u0 == 1.0 || v0 == 1.0 )
    return;
  sp = pkv_GetScratchMemTop ();
  i = max ( n, m ) + 1;
  q = (point3d*)pkv_GetScratchMem ( i*sizeof(point3d) );
  if ( !q )
    exit ( 1 );
  for ( i = 0; i <= d0; i++ ) {
    t = (double)i/(double)d0;
    if ( t >= u0 && t <= u1 ) {
      SetLW ( i, 0, d0 );
      mbs_multiBCHornerd ( n, 1, 3*(m+1), 0, (double*)cp, t, (double*)q );
      DrawBezCurve3d ( m, q, v0, v1 );
    }
  }
  for ( i = 0; i <= d1; i++ ) {
    t = (double)i/(double)d1;
    if ( t >= v0 && t <= v1 ) {
      SetLW ( i, 0, d1 );
      mbs_multiBCHornerd ( m, n+1, 3, 3*(m+1), (double*)cp, t, (double*)q );
      DrawBezCurve3d ( n, q, u0, u1 );
    }
  }
  pkv_SetScratchMemTop ( sp );
} /*DrawBezPatch3Rd*/

void DrawBezPatchNet2d ( int n, int m, const point2d *cp )
{
  void    *sp;
  int     i, j;
  point2d *cc;

  sp = pkv_GetScratchMemTop ();
  cc = (point2d*)pkv_GetScratchMem ( (max(n,m)+1)*sizeof(point2d) );
  ps_Set_Line_Width ( 2.0 );
  for ( i = 0; i <= n; i++ ) {
    for ( j = 0; j <= m; j++ )
      CameraProjectPoint2d ( &PPos, &cp[i*(m+1)+j], &cc[j] );
    ps_Draw_Polyline2d ( m+1, cc );
  }
  for ( j = 0; j <= m; j++ ) {
    for ( i = 0; i <= n; i++ )
      CameraProjectPoint2d ( &PPos, &cp[i*(m+1)+j], &cc[i] );
    ps_Draw_Polyline2d ( n+1, cc );
  }

  pkv_SetScratchMemTop ( sp );
} /*DrawBezPatchNet2d*/

void DrawBezPatchNet3d ( int n, int m, const point3d *cp )
{
  void    *sp;
  int     i, j;
  point3d q;
  point2d *r;

  sp = pkv_GetScratchMemTop ();
  r = pkv_GetScratchMem ( (n+1)*(m+1)*sizeof(point2d) );
  if ( !r )
    exit ( 1 );
  for ( i = 0; i < (n+1)*(m+1); i++ ) {
    CameraProjectPoint3d ( &CPos, &cp[i], &q );
    SetPoint2d ( &r[i], q.x, q.y );
  }
  ps_Set_Line_Width ( 1.0 );
  for ( i = 0; i <= n; i++ )
    ps_Draw_Polyline2d ( m+1, &r[i*(m+1)] );
  for ( j = 0; j <= m; j++ )
    for ( i = 0; i < n; i++ )
      ps_Draw_Line ( r[i*(m+1)+j].x, r[i*(m+1)+j].y,
                     r[(i+1)*(m+1)+j].x, r[(i+1)*(m+1)+j].y );
  for ( i = 0; i < (n+1)*(m+1); i++ )
    ps_Mark_Circle ( r[i].x, r[i].y );

  pkv_SetScratchMemTop ( sp );
} /*DrawBezPatchNet3d*/

void DrawBSPatch2d ( int n, int lknu, const double *knu,
                     int m, int lknv, const double *knv, const point2d *cp,
                     int d0, int d1 )
{
  void    *sp;
  int     ku, kv, pitch1, pitch2, pitch3;
  double   *b, *c;
  int     i, j, start;

  sp = pkv_GetScratchMemTop ();
  ku = mbs_NumKnotIntervalsd ( n, lknu, knu );
  kv = mbs_NumKnotIntervalsd ( m, lknv, knv );
  pitch1 = (lknv-m)*2;
  pitch2 = (m+1)*2*kv;
  pitch3 = (m+1)*2;
  b = pkv_GetScratchMemd ( pitch2*ku*(n+1) );
  c = pkv_GetScratchMemd ( (n+1)*pitch3 );
  if ( b && c ) {
    mbs_BSPatchToBezd ( 2, n, lknu, knu, m, lknv, knv, pitch1, (double*)cp,
                        &ku, NULL, NULL, &kv, NULL, NULL, pitch2, b );
    for ( i = 0; i < ku; i++ )
      for ( j = 0; j < kv; j++ ) {
        start = i*(n+1)*pitch2+j*pitch3;
        pkv_Selectd ( n+1, pitch3, pitch2, pitch3, &b[start], c );
        DrawBezPatch2d ( n, m, (point2d*)c, d0, d1 );
      }
  }

  pkv_SetScratchMemTop ( sp );
} /*DrawBSPatch2d*/

void DrawBSPatch3d ( int n, int lknu, const double *knu,
                     int m, int lknv, const double *knv, const point3d *cp,
                     int d0, int d1 )
{
  void    *sp;
  int     ku, kv, pitch1, pitch2, pitch3;
  double   *b, *c;
  int     i, j, start;

  sp = pkv_GetScratchMemTop ();
  ku = mbs_NumKnotIntervalsd ( n, lknu, knu );
  kv = mbs_NumKnotIntervalsd ( m, lknv, knv );
  pitch1 = (lknv-m)*3;
  pitch2 = (m+1)*3*kv;
  pitch3 = (m+1)*3;
  b = pkv_GetScratchMemd ( pitch2*ku*(n+1) );
  c = pkv_GetScratchMemd ( (n+1)*pitch3 );
  if ( b && c ) {
    mbs_BSPatchToBezd ( 3, n, lknu, knu, m, lknv, knv, pitch1, (double*)cp,
                        &ku, NULL, NULL, &kv, NULL, NULL, pitch2, b );
    for ( i = 0; i < ku; i++ )
      for ( j = 0; j < kv; j++ ) {
        start = i*(n+1)*pitch2+j*pitch3;
        pkv_Selectd ( n+1, pitch3, pitch2, pitch3, &b[start], c );
        DrawBezPatch3d ( n, m, (point3d*)c, d0, d1 );
      }
  }

  pkv_SetScratchMemTop ( sp );
} /*DrawBSPatch3d*/

void DrawBSPatch3Rd ( int n, int lknu, const double *knu,
                      int m, int lknv, const double *knv, const point3d *cp,
                      double u0, double u1, double v0, double v1,
                      int d0, int d1 )
{
  void    *sp;
  int     ku, kv, pitch1, pitch2, pitch3;
  double   *b, *c, *bknu, *bknv;
  double   ku0, ku1, kv0, kv1, uu0, uu1, vv0, vv1;
  int     i, j, start;

  if ( u0 > u1 || v0 > v1 )
    return;
  sp = pkv_GetScratchMemTop ();
  ku = mbs_NumKnotIntervalsd ( n, lknu, knu );
  kv = mbs_NumKnotIntervalsd ( m, lknv, knv );
  pitch1 = (lknv-m)*3;
  pitch2 = (m+1)*3*kv;
  pitch3 = (m+1)*3;
  b = pkv_GetScratchMemd ( pitch2*ku*(n+1) );
  c = pkv_GetScratchMemd ( (n+1)*pitch3 );
  bknu = pkv_GetScratchMemd ( (ku+1)*(n+1) );
  bknv = pkv_GetScratchMemd ( (kv+1)*(m+1) );
  if ( b && c && bknu && bknv ) {
    mbs_BSPatchToBezd ( 3, n, lknu, knu, m, lknv, knv, pitch1, (double*)cp,
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
            pkv_Selectd ( n+1, pitch3, pitch2, pitch3, &b[start], c );
            DrawBezPatch3Rd ( n, m, (point3d*)c,
                              uu0, uu1, vv0, vv1, d0, d1 );
          }
        }
      }
    }
  }

  pkv_SetScratchMemTop ( sp );
} /*DrawBSPatch3Rd*/

/* ////////////////////////////////////////////////////////////////////////// */
double FindBSPatchMaxLaplaciand ( int n, int lknu, const double *knu,
                                 int m, int lknv, const double *knv,
                                 const point3d *cp, int d0, int d1 )
{
  void    *sp;
  int     ku, kv, pitch1, pitch2, pitch3;
  double   *b, *c;
  int     i, j, k, l, start;
  double   lap, maxval;

  sp = pkv_GetScratchMemTop ();
  ku = mbs_NumKnotIntervalsd ( n, lknu, knu );
  kv = mbs_NumKnotIntervalsd ( m, lknv, knv );
  pitch1 = (lknv-m)*3;
  pitch2 = (m+1)*3*kv;
  pitch3 = (m+1)*3;
  b = pkv_GetScratchMemd ( pitch2*ku*(n+1) );
  c = pkv_GetScratchMemd ( (n+1)*pitch3 );
  maxval = 0.0;
  if ( b && c ) {
    mbs_BSPatchToBezd ( 3, n, lknu, knu, m, lknv, knv, pitch1, (double*)cp,
                        &ku, NULL, NULL, &kv, NULL, NULL, pitch2, b );
    for ( i = 0; i < ku; i++ )
      for ( j = 0; j < kv; j++ ) {
        start = i*(n+1)*pitch2+j*pitch3;
        pkv_Selectd ( n+1, pitch3, pitch2, pitch3, &b[start], c );
        for ( k = 0; k <= d0; k++ )
          for ( l = 0; l <= d1; l++ ) {
            lap = ComputeBezLaplaciand ( n, m, (point3d*)c,
                        (double)k/(double)d0, (double)l/(double)d1 );
            maxval = max ( maxval, fabs(lap) );
          }
      }
  }

  pkv_SetScratchMemTop ( sp );
  return maxval;
} /*FindBSPatchMaxLaplaciand*/

void DrawBSPatchLaplaciand ( int n, int lknu, const double *knu,
                             int m, int lknv, const double *knv,
                             const point3d *cp, int d0, int d1 )
{
#define DD 32
  void    *sp;
  int     ku, kv, pitch1, pitch2, pitch3;
  double   *b, *c, u, v;
  int     i, j, k, l, start;
  point3d p, q;
  point2d *cc;

  sp = pkv_GetScratchMemTop ();
  ku = mbs_NumKnotIntervalsd ( n, lknu, knu );
  kv = mbs_NumKnotIntervalsd ( m, lknv, knv );
  pitch1 = (lknv-m)*3;
  pitch2 = (m+1)*3*kv;
  pitch3 = (m+1)*3;
  b = pkv_GetScratchMemd ( pitch2*ku*(n+1) );
  c = pkv_GetScratchMemd ( (n+1)*pitch3 );
  cc = pkv_GetScratchMem ( (DD+1)*sizeof(point2d) );
  if ( b && c && cc ) {
    mbs_BSPatchToBezd ( 3, n, lknu, knu, m, lknv, knv, pitch1, (double*)cp,
                        &ku, NULL, NULL, &kv, NULL, NULL, pitch2, b );
    for ( i = 0; i < ku; i++ )
      for ( j = 0; j < kv; j++ ) {
        start = i*(n+1)*pitch2+j*pitch3;
        pkv_Selectd ( n+1, pitch3, pitch2, pitch3, &b[start], c );
        for ( k = 0; k <= d0; k++ ) {
          u = (double)k/(double)d0;
          for ( l = 0; l <= DD; l++ ) {
            v = (double)l/(double)DD;
            mbs_BCHornerP3d ( n, m, (point3d*)c, u, v, &p );
            p.z = ComputeBezLaplaciand ( n, m, (point3d*)c, u, v );
            CameraProjectPoint3d ( &CPos, &p, &q );
            SetPoint2d ( &cc[l], q.x, q.y );
          }
          SetCLW ( k, 0, d0 );
          ps_Draw_Polyline2d ( DD+1, cc );
        }
        for ( l = 0; l <= d1; l++ ) {
          v = (double)l/(double)d1;
          for ( k = 0; k <= DD; k++ ) {
            u = (double)k/(double)DD;
            mbs_BCHornerP3d ( n, m, (point3d*)c, u, v, &p );
            p.z = ComputeBezLaplaciand ( n, m, (point3d*)c, u, v );
            CameraProjectPoint3d ( &CPos, &p, &q );
            SetPoint2d ( &cc[k], q.x, q.y );
          }
          SetCLW ( l, 0, d1 );
          ps_Draw_Polyline2d ( DD+1, cc );
        }
      }
  }

  pkv_SetScratchMemTop ( sp );
#undef DD
} /*DrawBSPatchLaplaciand*/

