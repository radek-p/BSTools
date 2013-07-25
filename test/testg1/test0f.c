
/* //////////////////////////////////////////////////// */
/* This file is a part of the BSTools procedure package */
/* written by Przemyslaw Kiciak.                        */
/* //////////////////////////////////////////////////// */

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
#include "datagenf.h"
#include "testgraphf.h"

/*
#define HOLE_K 5
#define FUNCA21
*/
#define NQUAD 16


#define HOLE_K 6
#define FUNCA25


char fn1[]  = "domsurrpf.ps";
char fn2[]  = "auxpatchesf.ps";
char fn3a[] = "jfunc1f.ps";
char fn4[]  = "dipatchesf.ps";
char fn5[]  = "partitionf.ps";
char fn6a[] = "basf0f.ps";
char fn6b[] = "basf1f.ps";
char fn6c[] = "basf2f.ps";
char fn6d[] = "basf3f.ps";
char fn6e[] = "basf4f.ps";
char fn6f[] = "basf5f.ps";
char fn6g[] = "basf6f.ps";
char fn7a[] = "basb0f.ps";
char fn7b[] = "basb1f.ps";
char fn7c[] = "basb2f.ps";
char fn7d[] = "basb3f.ps";
char fn7e[] = "basb4f.ps";
char fn7f[] = "basb5f.ps";
char fn7g[] = "basb6f.ps";
char fn7h[] = "basb7f.ps";
char fn10[] = "amatf.ps";

GHoleDomainf *domain;

CameraRecf   CPos, PPos;

/* ///////////////////////////////////////////////////////////////////////// */
void SetupPProj ( float w, float h, float x, float y, float diag )
{
  CameraInitFramef ( &PPos, true, true, w, h, x, y, 1.0, 4 );
  CameraInitPosf ( &PPos );
  PPos.vd.para.diag = diag;
  CameraSetMappingf ( &PPos );
} /*SetupPProj*/

void SetupCPos ( float w, float h, float x, float y )
{
  vector3f v;

  CameraInitFramef ( &CPos, false, true, w, h, x, y, 1.0, 4 );
  CameraInitPosf ( &CPos );
  SetVector3f ( &v, 0.0, 0.0, -46.0 );
  CameraMoveGf ( &CPos, &v );
  CameraRotXGf ( &CPos, -0.65*PI );
  CameraRotZGf ( &CPos, 0.05*PI );
  CameraZoomf ( &CPos, 11.0 );
} /*SetupCPos*/

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
  int     i;
  point3f p, q;
  point2f *cc;

  sp = pkv_GetScratchMemTop ();
  cc = (point2f*)pkv_GetScratchMem ( (DD+1)*sizeof(point2f) );
  if ( !cc )
    exit ( 1 );

  for ( i = 0; i <= DD; i++ ) {
    mbs_BCHornerC3f ( n, cp, (float)i/(float)DD, &p );
    CameraProjectPoint3f ( &CPos, &p, &q );
    SetPoint2f ( &cc[i], q.x, q.y );
  }
  ps_Draw_Polyline2f ( DD+1, cc );

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

void SetLW ( int i, int i0, int i1 )
{
    if ( i == i0 || i == i1 )
      ps_Set_Line_Width ( 4.0 );
    else
      ps_Set_Line_Width ( 2.0 );
} /*SetLW*/

void DrawBezPatch2f ( int n, int m, const point2f *cp )
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
    SetLW ( i, 0, D );
    DrawBezCurve2f ( m, cc, 0.0, 1.0 );
  }
  for ( i = 0; i <= D; i++ ) {
    mbs_multiBCHornerf ( m, n+1, 2, 2*(m+1), (float*)cp, (float)i/(float)D,
                         (float*)cc );
    SetLW ( i, 0, D );
    DrawBezCurve2f ( n, cc, 0.0, 1.0 );
  }

  pkv_SetScratchMemTop ( sp );
#undef D
} /*DrawBezPatch2f*/

void DrawBezPatch3f ( int n, int m, const point3f *cp )
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
    SetLW ( i, 0, D );
    DrawBezCurve3f ( m, cc, 0.0, 1.0 );
  }
  for ( i = 0; i <= D; i++ ) {
    mbs_multiBCHornerf ( m, n+1, 3, 3*(m+1), (float*)cp, (float)i/(float)D,
                         (float*)cc );
    SetLW ( i, 0, D );
    DrawBezCurve3f ( n, cc, 0.0, 1.0 );
  }

  pkv_SetScratchMemTop ( sp );
#undef D
} /*DrawBezPatch3f*/

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

void DrawPatchNet3f ( int n, int m, const point3f *cp )
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
} /*DrawPatchNet3f*/

void DrawBezPatchCPNum2f ( int n, int m, const point2f *cp )
{
  int     i;
  point2f p;
  char    s[30];

  ps_Write_Command ( "/Times-Roman findfont 40 scalefont setfont" );
  for ( i = 0; i < (m+1)*(n+1); i++ ) {
    CameraProjectPoint2f ( &PPos, &cp[i], &p );
    ps_MoveTo ( p.x, p.y );
    sprintf ( s, "(%d) show", i );
    ps_Write_Command ( s );
  }
} /*DrawBezPatchCPNum2f*/

/* ///////////////////////////////////////////////////////////////////////// */
void P1DrawBezPatch ( int n, int m, const point2f *cp )
{
  ps_Set_Gray ( 0.68 );
  DrawBezPatch2f ( n, m, cp );
  ps_Set_Gray ( 0.5 );
  DrawBezPatchNet2f ( n, m, cp );
  ps_Set_Gray ( 0.0 );
  DrawBezPatchCPNum2f ( n, m, cp );
} /*P1DrawBezPatch*/

void MakePicture1 ()
{
  ps_OpenFile ( fn1, 600 );
  ps_Write_Command ( "1 setlinecap" );
  SetupPProj ( 1800, 1500, 0.0, 0.0, 7.5 );
  gh_DrawDomSurrndPatchesf ( domain, P1DrawBezPatch );
  ps_CloseFile ();
  printf ( "%s\n", fn1 );
} /*MakePicture1*/

/* ///////////////////////////////////////////////////////////////////////// */
void MakePicture2 ()
{
  ps_OpenFile ( fn2, 600 );
  ps_Write_Command ( "1 setlinecap" );
  SetupPProj ( 1800, 1500, 0.0, 0.0, 7.5 );
  ps_Set_Gray ( 0.68 );
  gh_DrawDomSurrndPatchesf ( domain, DrawBezPatch2f );
  ps_Set_Gray ( 0.0 );
  g1h_DrawDomAuxPatchesf ( domain, DrawBezPatch2f );
  ps_CloseFile ();
  printf ( "%s\n", fn2 );
} /*MakePicture2*/

/* ///////////////////////////////////////////////////////////////////////// */
float gw, gh, gx, gy;

void DrawBPGraph ( int deg, const float *f )
{
#define D 32
  void    *sp;
  point2f *cc;
  float   x, y;
  int     i;

  sp = pkv_GetScratchMemTop ();
  cc = (point2f*)pkv_GetScratchMem ( (D+1)*sizeof(point2f) );
  if ( !cc )
    exit ( 1 );
  for ( i = 0; i <= D; i++ ) {
    x = (float)i/D;
    mbs_BCHornerC1f ( deg, f, x, &y );
    SetPoint2f ( &cc[i], gx+gw*x, gy+gh*y );
  }
  ps_Draw_Polyline2f ( D+1, cc );
  pkv_SetScratchMemTop ( sp );
#undef D
} /*DrawBPGraph*/

void DrawJFuncGraph ( int k, int l,
                      float w, float h, float x, float y, float gr )
{
  gw = w;  gh = h;  gx = x;  gy = y;
  ps_Set_Gray ( 0.0 );
  psl_SetLine ( x, y, x+w, y, 0.0, 1.0 );
  psl_Draw ( 0.0, 1.3, 2.0 );
  psl_Tick ( 0.0 );
  psl_Tick ( 1.0 );
  psl_Arrow ( 1.5, true );
  psl_SetLine ( x, y, x, y+h, 0.0, 1.0 );
  psl_Draw ( -1.2, 1.5, 2.0 );
  psl_Tick ( -1.0 );
  psl_Tick ( 0.0 );
  psl_Tick ( 1.0 );
  psl_Arrow ( 1.7, true );
  ps_Set_Line_Width ( 6.0 );
  ps_Set_RGB ( 1.0-gr, gr, 0.0 );
  g1h_DrawJFunctionf ( domain, k, l, DrawBPGraph );
} /*DrawJFuncGraph*/

void MakePicture3 ()
{
  int i;
  float gr;

  ps_OpenFile ( fn3a, 600 );
  ps_Write_Command ( "1 setlinecap" );
  for ( i = 0; i < HOLE_K; i++ ) {
    gr = (float)i/(float)(HOLE_K-1);
    DrawJFuncGraph ( i, 0, 200.0, 150.0, 50.0, 500.0*(HOLE_K-i), gr );    /*b01*/
    DrawJFuncGraph ( i, 1, 200.0, 150.0, 400.0, 500.0*(HOLE_K-i), gr );   /*c01*/
    DrawJFuncGraph ( i, 2, 200.0, 150.0, 750.0, 500.0*(HOLE_K-i), gr );   /*f01*/
    DrawJFuncGraph ( i, 3, 200.0, 150.0, 1100.0, 500.0*(HOLE_K-i), gr );  /*g01*/
    DrawJFuncGraph ( i, 4, 200.0, 150.0, 1450.0, 500.0*(HOLE_K-i), gr );  /*b11*/
    DrawJFuncGraph ( i, 5, 200.0, 150.0, 1800.0, 500.0*(HOLE_K-i), gr );  /*c11*/
    DrawJFuncGraph ( i, 6, 200.0, 150.0, 2150.0, 500.0*(HOLE_K-i), gr );  /*f11*/
    DrawJFuncGraph ( i, 7, 200.0, 150.0, 2500.0, 500.0*(HOLE_K-i), gr );  /*g11*/
  }
  ps_CloseFile ();
  printf ( "%s\n", fn3a );
} /*MakePicture3*/

/* ///////////////////////////////////////////////////////////////////////// */
void DumpPatch ( int n, int m, const point2f *cp )
{
  FILE *f;
  int i;

  if ( n > 3 ) {
    f = fopen ( "d0f.dat", "w+" );
    for ( i = 0; i < (n+1)*(m+1); i++ ) {
      fprintf ( f, "{%f,%f},", cp[i].x, cp[i].y );
      if ( (i+1) % (m+1) == 0 )
        fprintf ( f, "\n" );
    }
    fclose ( f );
    exit ( 0 );
  }
} /*DumpPatch*/

void MakePicture4 ()
{
/*
  g1h_DrawDiPatchesf ( domain, DumpPatch );
*/
  ps_OpenFile ( fn4, 600 );
  ps_Write_Command ( "1 setlinecap" );
  SetupPProj ( 1800, 1500, 0.0, 0.0, 7.5 );
  ps_Set_Gray ( 0.68 );
  gh_DrawDomSurrndPatchesf ( domain, DrawBezPatch2f );
  ps_Set_Gray ( 0.0 );
  g1h_DrawDiPatchesf ( domain, DrawBezPatch2f );
  ps_CloseFile ();
  printf ( "%s\n", fn4 );
} /*MakePicture4*/

/* ///////////////////////////////////////////////////////////////////////// */
void DrawPartition ( int hole_k, int hole_m,
                     float *partition, float *spart_alpha, float *spart_malpha,
                     float *spart_salpha, float *spart_knot, float alpha0,
                     boolean *spart_sgn, boolean *spart_both )
{
#define R 4.0
  int     i;
  point2f p, q, r, s;
  float   lgt;

  SetPoint2f ( &p, 0.0, 0.0 );
  CameraProjectPoint2f ( &PPos, &p, &q );
  ps_Set_Line_Width ( 2.0 );
  SetPoint2f ( &p, R*cos(alpha0), R*sin(alpha0) );
  CameraProjectPoint2f ( &PPos, &p, &r );
  ps_Set_Gray ( 0.68 );
  ps_Draw_Line ( q.x, q.y, r.x, r.y );
  ps_Set_Gray ( 0.0 );
  for ( i = 0; i < hole_k; i++ ) {
    SetPoint2f ( &p, R*cos(partition[i]), R*sin(partition[i]) );
    CameraProjectPoint2f ( &PPos, &p, &r );
    ps_Draw_Line ( q.x, q.y, r.x, r.y );
  }
  for ( i = 0; i < hole_k-hole_m; i++ )
    if ( spart_sgn[i] && !spart_both[i] ) {
      ps_Write_Command ( "[ 21 ] 0 setdash" );
      SetPoint2f ( &p, R*cos(spart_malpha[i]), R*sin(spart_malpha[i]) );
      CameraProjectPoint2f ( &PPos, &p, &r );
      ps_Draw_Line ( q.x, q.y, r.x, r.y );
      ps_Write_Command ( "[ ] 0 setdash" );
    }
  SetPoint2f ( &p, cos(alpha0), sin(alpha0) );
  CameraProjectPoint2f ( &PPos, &p, &r );
  SetPoint2f ( &p, cos(alpha0)+sin(alpha0), sin(alpha0)-cos(alpha0) );
  CameraProjectPoint2f ( &PPos, &p, &s );
  psl_SetLine ( r.x, r.y, s.x, s.y, 0.0, 1.0 );
  lgt = max( fabs(spart_knot[0]), fabs(spart_knot[hole_k-hole_m-1]) );
  psl_Draw ( -lgt-0.5, lgt+0.5, 2.0 );
  for ( i = 0; i < hole_k-hole_m; i++ )
    psl_Tick ( spart_knot[i] );

  ps_Mark_Circle ( q.x, q.y );
#undef R
} /*DrawPartition*/

void DrawBSpline2 ( int nkn, float alpha0, float *knots )
{
#define D 200
  void    *sp;
  float   *kn;
  point2f *d, *cc, c;
  int     lkn, i, j;
  float   t, sa, ca;

  if ( nkn < 4 )
    return;
  sp = pkv_GetScratchMemTop ();
  lkn = nkn+4;
  kn = pkv_GetScratchMemf ( lkn+1 );
  d = (point2f*)pkv_GetScratchMem ( (lkn-2)*sizeof(point2f) );
  cc = (point2f*)pkv_GetScratchMem ( (D+1)*sizeof(point2f) );
  if ( !kn || !d || !cc )
    exit ( 1 );
  kn[0] = kn[1] = kn[2] = knots[nkn-1];
  for ( i = nkn-2, j = 3;  i > 0;  i--, j++ )
    kn[j] = knots[i];
  kn[j] = kn[j+1] = kn[j+2] = knots[0];
  sa = sin ( alpha0 );  ca = cos ( alpha0 );
  ps_Set_Line_Width ( 4.0 );
  for ( j = 0; j < nkn-3; j++ ) {
    ps_Set_Gray ( 0.2+0.6*(float)j/(float)(lkn-2) );
    for ( i = 0; i < lkn-2; i++ ) {
      d[i].y = (kn[i+1]+kn[i+2])/2.0;
      d[i].x = 1.0;
    }
    d[2+j].x = 1.0+2.0;
    for ( i = 0; i < lkn-2; i++ ) {
      t      = ca*d[i].x - sa*d[i].y;
      d[i].y = sa*d[i].x + ca*d[i].y;
      d[i].x = t;
    }
    for ( i = 0; i <= D; i++ ) {
      t = kn[2] + (float)i/(float)D*(kn[lkn-2]-kn[2]);
      mbs_deBoorC2f ( 2, lkn, kn, d, t, &c );
      CameraProjectPoint2f ( &PPos, &c, &cc[i] );
    }
    ps_Draw_Polyline2f ( D+1, cc );
  }

  pkv_SetScratchMemTop ( sp );
#undef D
} /*DrawBSpline2*/

void MakePicture5 ()
{
  int     hole_k, hole_m;
  float   partition[8], part_delta[8], spart_alpha[8], spart_malpha[8],
          spart_salpha[8], spart_knot[8];
  float   alpha0;
  boolean spart_sgn[8], spart_both[8];

  ps_WriteBBox ( 1, 1, 310, 147 );
  ps_OpenFile ( fn5, 600 );
  ps_Write_Command ( "1 setlinecap" );

  g1h_ExtractPartitionf ( domain, &hole_k, &hole_m,
                          partition, part_delta, spart_alpha, spart_malpha,
                          spart_salpha, spart_knot, &alpha0,
                          spart_sgn, spart_both );
/*
  for ( i = 0; i < hole_k; i++ )
    printf ( "%6.3f %6.3f\n",
              partition[i], part_delta[i] );
  for ( i = 0; i < hole_k-hole_m; i++ )
    printf ( "%6.3f %6.3f %6.3f %6.3f %d %d\n",
              spart_alpha[i], spart_malpha[i], spart_salpha[i],
              spart_knot[i], spart_sgn[i], spart_both[i] );
*/
  if ( hole_k-hole_m >= 4 ) {
    SetupPProj ( 1200, 1200, 0.0, 0.0, 11.0 );
    DrawPartition ( hole_k, hole_m, partition, spart_alpha, spart_malpha,
                    spart_salpha, spart_knot, alpha0, spart_sgn, spart_both );
    DrawBSpline2 ( hole_k-hole_m, alpha0, spart_knot );
  }
  ps_CloseFile ();
  printf ( "%s\n", fn5 );
} /*MakePicture5*/

/* ///////////////////////////////////////////////////////////////////////// */
void DescribeGraph ( float s )
{
  char c[50];

  ps_Write_Command ( "/Times-Roman findfont 70 scalefont setfont" );
  ps_MoveTo ( CPos.xmin+10, CPos.ymin+CPos.height-80 );
  sprintf ( c, "(%7.3g) show", s );
  ps_Write_Command ( c );
} /*DescribeGraph*/

static boolean scaling;
static float graphscf;

void DrawBasisGraphBezPatch ( int n, int m, const point3f *cp )
{
  if ( n == 3 )
    ps_Set_Gray ( 0.68 );
  else
    ps_Set_Gray ( 0.25 );
  DrawBezPatch3f ( n, m, cp );
} /*DrawBasisGraphBezPatch*/

void DrawScBasisGraphBezPatch ( int n, int m, const point3f *cp )
{
  void    *sp;
  int     i;
  point3f *ccp;

  sp = pkv_GetScratchMemTop ();
  ccp = (point3f*)pkv_GetScratchMem ( (n+1)*(m+1)*sizeof(point3f) );
  memcpy ( ccp, cp, (n+1)*(m+1)*sizeof(point3f) );
  for ( i = 0; i < (n+1)*(m+1); i++ )
    ccp[i].z *= graphscf;
  if ( n == 3 )
    ps_Set_Gray ( 0.68 );
  else
    ps_Set_Gray ( 0.25 );
  DrawBezPatch3f ( n, m, ccp );
  pkv_SetScratchMemTop ( sp );
} /*DrawScBasisGraphBezPatch*/

void ScaleBasisGraphBezPatch ( int n, int m, const point3f *cp )
{
  void *sp;
  int  i;

  sp = pkv_GetScratchMemTop ();
  for ( i = 0; i < (n+1)*(m+1); i++ )
    graphscf = max ( graphscf, fabs(cp[i].z) );
  pkv_SetScratchMemTop ( sp );
} /*ScaleBasisGraphBezPatch*/

void DrawBasisAFuncGraph ( int fn, boolean scale )
{
  scaling = scale;
  ps_BeginDict ( 10 );
  ps_Write_Command ( "/setgr { setgray } bind def" );
  ps_Write_Command ( "/setgray { 0.5 mul 0.5 add setgr } def" );
  ps_Write_Command ( "0 setgray" );
  g1h_DrawBasCNetf ( domain, -1, DrawPatchNet3f );
  ps_EndDict ();
  ps_Write_Command ( "0 setgray" );
  if ( scaling ) {
    graphscf = 0.001;
    g1h_DrawBasAFunctionf ( domain, fn, ScaleBasisGraphBezPatch );
    graphscf = 0.5/graphscf;
    g1h_DrawBasAFunctionf ( domain, fn, DrawScBasisGraphBezPatch );
  }
  else
    g1h_DrawBasAFunctionf ( domain, fn, DrawBasisGraphBezPatch );
} /*DrawBasisAFuncGraph*/

void DrawBasisBFuncGraph ( int fn, boolean scale )
{
  scaling = scale;
  ps_BeginDict ( 10 );
  ps_Write_Command ( "/setgr { setgray } bind def" );
  ps_Write_Command ( "/setgray { 0.5 mul 0.5 add setgr } def" );
  ps_Write_Command ( "0 setgray" );
  g1h_DrawBasCNetf ( domain, fn, DrawPatchNet3f );
  ps_EndDict ();
  ps_Write_Command ( "0 setgray" );
  if ( scaling ) {
    graphscf = 0.001;
    g1h_DrawBasBFunctionf ( domain, fn, ScaleBasisGraphBezPatch );
    graphscf = 0.5/graphscf;
    g1h_DrawBasBFunctionf ( domain, fn, DrawScBasisGraphBezPatch );
  }
  else
    g1h_DrawBasBFunctionf ( domain, fn, DrawBasisGraphBezPatch );
} /*DrawBasisBFuncGraph*/

void DrawBasFuncGraph ( char ftype, int fn, float x, float y, boolean scale )
{
  char s[40];

  ps_GSave ();
  sprintf ( s, "%f %f translate", x, y );
  ps_Write_Command ( s );
  ps_Set_Gray ( 0.0 );
  ps_Set_Line_Width ( 1.0 );
  ps_Draw_Rect ( CPos.width, CPos.height, CPos.xmin, CPos.ymin );
  ps_Set_Clip_Rect ( CPos.width, CPos.height, CPos.xmin, CPos.ymin );
  switch ( ftype ) {
case 0:
    DrawBasisAFuncGraph ( fn, scale );
    break;
case 1:
    DrawBasisBFuncGraph ( fn, scale );
default:
    break;
  }
  ps_GRestore ();
} /*DrawBasFuncGraph*/

#define D   8
#define DD 48
void ScaleBasisGraphLaplacian ( int n, int m, const point3f *cp )
{
  void    *sp;
  int     i, j;
  point3f *lgp;

  if ( n > 3 ) {
    sp = pkv_GetScratchMemTop ();
    lgp = (point3f*)pkv_GetScratchMem ( (DD+1)*sizeof(point3f) );
    if ( !lgp )
      exit ( 1 );
    for ( i = 0; i <= D; i++ ) {
      TestTabLaplacianU ( n, m, cp, (float)i/(float)D, DD, lgp );
      for ( j = 0; j <= DD; j++ )
        graphscf = max ( graphscf, fabs(lgp[j].z) );
      TestTabLaplacianV ( n, m, cp, (float)i/(float)D, DD, lgp );
      for ( j = 0; j <= DD; j++ )
        graphscf = max ( graphscf, fabs(lgp[j].z) );
    }

    pkv_SetScratchMemTop ( sp );
  }
} /*ScaleBasisGraphLaplacian*/

void DrawBasisGraphLaplacian ( int n, int m, const point3f *cp )
{
  void    *sp;
  int     i, j;
  point3f *lgp, q;
  point2f *cc;

  sp = pkv_GetScratchMemTop ();
  lgp = (point3f*)pkv_GetScratchMem ( (DD+1)*sizeof(point3f) );
  cc = (point2f*)pkv_GetScratchMem ( (DD+1)*sizeof(point2f) );
  if ( !lgp || !cc )
    exit ( 0 );

  if ( m == 3 )
    ps_Set_Gray ( 0.68 );
  else
    ps_Set_Gray ( 0.25 );
  for ( i = 0; i <= D; i++ ) {
    TestTabLaplacianU ( n, m, cp, (float)i/(float)D, DD, lgp );
    for ( j = 0; j <= DD; j++ ) {
      lgp[j].z *= graphscf;
      CameraProjectPoint3f ( &CPos, &lgp[j], &q );
      cc[j].x = q.x;  cc[j].y = q.y;
    }
    SetLW ( i, 0, D );
    ps_Draw_Polyline2f ( DD+1, cc );
  }
  for ( i = 0; i <= D; i++ ) {
    TestTabLaplacianV ( n, m, cp, (float)i/(float)D, DD, lgp );
    for ( j = 0; j <= DD; j++ ) {
      lgp[j].z *= graphscf;
      CameraProjectPoint3f ( &CPos, &lgp[j], &q );
      cc[j].x = q.x;  cc[j].y = q.y;
    }
    SetLW ( i, 0, D );
    ps_Draw_Polyline2f ( DD+1, cc );
  }

  pkv_SetScratchMemTop ( sp );
} /*DrawBasisGraphLaplacian*/
#undef DD
#undef D

void DrawBasisALaplacian ( int fn )
{
  float s;

  graphscf = 0.001;
  g1h_DrawBasAFunctionf ( domain, fn, ScaleBasisGraphLaplacian );
  s = graphscf;
  graphscf = 0.25/graphscf;
  g1h_DrawBasAFunctionf ( domain, fn, DrawBasisGraphLaplacian );
  DescribeGraph ( s );
} /*DrawBasisALaplacian*/

void DrawBasisBLaplacian ( int fn )
{
  float s;

  graphscf = 0.001;
  g1h_DrawBasBFunctionf ( domain, fn, ScaleBasisGraphLaplacian );
  s = graphscf;
  graphscf = 0.25/graphscf;
  g1h_DrawBasBFunctionf ( domain, fn, DrawBasisGraphLaplacian );
  DescribeGraph ( s );
} /*DrawBasisBLaplacian*/

void DrawBasFuncLaplacian ( char ftype, int fn, float x, float y )
{
  char s[40];

  ps_GSave ();
  sprintf ( s, "%f %f translate", x, y );
  ps_Write_Command ( s );
  ps_Set_Gray ( 0.0 );
  ps_Set_Line_Width ( 1.0 );
  ps_Draw_Rect ( CPos.width, CPos.height, CPos.xmin, CPos.ymin );
  ps_Set_Clip_Rect ( CPos.width, CPos.height, CPos.xmin, CPos.ymin );
  switch ( ftype ) {
case 0:
    DrawBasisALaplacian ( fn );
    break;
case 1:
    DrawBasisBLaplacian ( fn );
default:
    break;
  }
  ps_GRestore ();
} /*DrawBasFuncLaplacian*/

char ssxy;

void MakePicture6x ( const char *fn, int ft, int ffn, int nf,
                     boolean sc0, boolean sc1, boolean sc2, boolean sc3 )
{
  nf -= ffn;
  if ( nf <= 0 )
    return;

  ps_OpenFile ( fn, 600 );
  ps_DenseScreen ();
  ps_Write_Command ( "1 setlinecap" );
  SetupCPos ( 1600.0, 1160.0, 0.0, 0.0 );

  DrawBasFuncGraph ( ft, ffn+0, 20.0, 3610.0, sc0 );
  DrawBasFuncLaplacian ( ft, ffn+0, 1660.0, 3610.0 );
  if ( nf > 1 ) {
    DrawBasFuncGraph ( ft, ffn+1, 20.0, 2410.0, sc1 );
    DrawBasFuncLaplacian ( ft, ffn+1, 1660.0, 2410.0 );
  }
  if ( nf > 2 ) {
    DrawBasFuncGraph ( ft, ffn+2, 20.0, 1210.0, sc2 );
    DrawBasFuncLaplacian ( ft, ffn+2, 1660.0, 1210.0);
  }
  if ( nf > 3 ) {
    DrawBasFuncGraph ( ft, ffn+3, 20.0,   10.0, sc3 );
    DrawBasFuncLaplacian ( ft, ffn+3, 1660.0,   10.0 );
  }

  ps_CloseFile ();
  printf ( "%s\n", fn );
} /*MakePicture6x*/

void MakePicture6 ()
{
  int nfa;

  nfa = g1h_V0SpaceDimf ( domain );
  MakePicture6x ( fn6a, 0,  0, nfa, false, false, false, true );
  MakePicture6x ( fn6b, 0,  4, nfa, true, true, true, true );
  MakePicture6x ( fn6c, 0,  8, nfa, true, true, true, true );
  MakePicture6x ( fn6d, 0, 12, nfa, true, true, true, true );
  MakePicture6x ( fn6e, 0, 16, nfa, true, true, true, true );
  MakePicture6x ( fn6f, 0, 20, nfa, true, true, true, true );
  MakePicture6x ( fn6g, 0, 24, nfa, true, true, true, true );
} /*MakePicture6*/

void MakePicture7 ()
{
  int nfb;

  nfb = 12*HOLE_K+1;
  MakePicture6x ( fn7a, 1,  0, nfb, false, false, false, false );
  MakePicture6x ( fn7b, 1,  4, nfb, false, false, false, false );
  MakePicture6x ( fn7c, 1,  8, nfb, false, false, false, false );
  MakePicture6x ( fn7d, 1, 12, nfb, false, false, false, false );
  MakePicture6x ( fn7e, 1, 16, nfb, false, false, false, false );
  MakePicture6x ( fn7f, 1, 20, nfb, false, false, false, false );
  MakePicture6x ( fn7g, 1, 24, nfb, false, false, false, false );
  MakePicture6x ( fn7h, 1, 28, nfb, false, false, false, false );
} /*MakePicture7*/

/* ///////////////////////////////////////////////////////////////////////// */
void PaintMatrix ( int n, const float *amat )
{
#define DMAX 40
  int   i, j, k;
  float amax, a;

  amax = fabs(amat[0]);
  for ( i = 1; i < n*(n+1)/2; i++ ) {
    a = fabs(amat[i]);
    if ( a > amax ) amax = a;
  }
  ps_Set_Gray ( 0.0 );
  ps_Set_Line_Width ( 2.0 );
  ps_Draw_Rect ( 2.0*DMAX*n+2.0, 2.0*DMAX*n+2.0, DMAX+1.0, DMAX+1.0 );
  for ( i = 0; i < n; i++ )
    for ( j = 0; j < n; j++ ) {
      k = pkn_SymMatIndex(i,j);
      a = fabs(amat[k])/amax;
      if ( a > 1.0e-10 ) {
        a = DMAX*pow(a,0.2);
        if ( amat[k] > 0 ) ps_Set_Gray ( 0.0 );
        else ps_Set_RGB ( 0.5, 1.0, 0.0 );
        ps_Fill_Circle ( 2.0*DMAX*(j+1), 2.0*DMAX*(n-i), a );
      }
    }
#undef DMAX
} /*PaintMatrix*/

void FindCondNumber ( int n, const float *a )
{
  void  *sp;
  float *b, *c, *x, *y, *z;
  float lambda_1, lambda_n, pr, yl, ca;
  int   i, j;

  sp = pkv_GetScratchMemTop ();
  /* power method of finding the greatest eigenvalue */

  j = -1;
  yl = 0.0;
  for ( i = 0; i < n*(n+1)/2; i++ )
    if ( fabs(a[i]) > yl ) {
      yl = max ( yl, fabs(a[i]) );
      j = i;
    }
  printf ( "a_max = %g\n", yl );

  b = pkv_GetScratchMemf ( n*(n+1)/2 );
  x = pkv_GetScratchMemf ( n );
  y = pkv_GetScratchMemf ( n );
  z = pkv_GetScratchMemf ( n );
  if ( !b || !x || !y || !z )
    exit ( 1 );
  c = &b[n*n];
  memcpy ( b, a, (n*(n+1)/2)*sizeof(float) );

  memset ( x, 0, n*sizeof(float) );
  x[0] = 1.0;
  for ( i = 0; i < 100; i++ ) {
    pkn_SymMatrixMultf ( n, b, 1, 1, x, 1, y );
    pr = pkn_ScalarProductf ( n, x, y );
    yl = pkn_SecondNormf ( n, y );
    ca = pr/yl;
    pkn_MultMatrixNumf ( 1, n, 0, y, 1.0/yl, 0, x );
    if ( ca > 0.9999999 )
      break;
  }
/*
  printf ( "ca = %f, yl = %g\n", ca, yl );
*/
  lambda_1 = yl;

  pkn_CholeskyDecompf ( n, b );
  memset ( x, 0, n*sizeof(float) );
  x[0] = 1.0;
  for ( i = 0; i < 100; i++ ) {
    pkn_LowerTrMatrixSolvef ( n, b, 1, 1, x, 1, z );
    pkn_UpperTrMatrixSolvef ( n, b, 1, 1, z, 1, y );
    pr = pkn_ScalarProductf ( n, x, y );
    yl = pkn_SecondNormf ( n, y );
    ca = pr/yl;
    pkn_MultMatrixNumf ( 1, n, 0, y, 1.0/yl, 0, x );
    if ( ca > 0.9999999 )
      break;
  }
/*
  printf ( "ca = %f, yl = %g\n", ca, yl );
*/
  lambda_n = 1.0/yl;
  printf ( "l_1 = %g, l_n = %g, cond A = %g\n",
           lambda_1, lambda_n, lambda_1/lambda_n );

  pkv_SetScratchMemTop ( sp );
} /*FindCondNumber*/

void DrawMatrices ( int nfa, int nfb, float *amat, float *bmat )
{
  PaintMatrix ( nfa, amat );
  FindCondNumber ( nfa, amat );
} /*DrawMatrices*/

void MakePicture10 ()
{
  ps_OpenFile ( fn10, 600 );
  g1h_DrawMatricesf ( domain, DrawMatrices );
  ps_CloseFile ();
  printf ( "%s\n", fn10 );
} /*MakePicture10*/

/* ///////////////////////////////////////////////////////////////////////// */
void MakePictures ()
{

  MakePicture1 ();
  MakePicture2 ();
  MakePicture3 ();
  MakePicture4 ();
  MakePicture5 ();
  MakePicture6 ();
  MakePicture7 ();
  MakePicture10 ();
} /*MakePictures*/

int main ()
{
  struct tms start, stop;
  float time;
  float param[2] = {0.0,0.0};

  pkv_InitScratchMem ( 2097152 );
  InitKnots ( HOLE_K );
  InitDomain ( HOLE_K, param );
  times ( &start );
  if ( (domain = gh_CreateDomainf ( HOLE_K, knots, Dompt )) ) {
    if ( g1h_ComputeBasisf ( domain ) ) {
      if ( g1h_ComputeFormMatrixf ( domain ) ) {
        times ( &stop );
        MakePictures ();
        time = (float)(stop.tms_utime-start.tms_utime)/
               (float)(sysconf(_SC_CLK_TCK));
        printf ( "time = %8.3f\n", time );
      }
    }
    gh_DestroyDomainf ( domain );
  }

  printf ( "Scratch memory used: %d bytes\n", pkv_MaxScratchTaken() );
  pkv_DestroyScratchMem ();
  exit ( 0 );
} /*main*/

