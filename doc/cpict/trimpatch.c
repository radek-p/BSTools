
/* //////////////////////////////////////////////////// */
/* This file is a part of the BSTools procedure package */
/* written by Przemyslaw Kiciak.                        */
/* //////////////////////////////////////////////////// */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "pkvaria.h"
#include "pkgeom.h"
#include "camera.h"
#include "pknum.h"
#include "multibs.h"
#include "psout.h"

char fn[] = "trimpatch.ps";

#define SCRATCHSIZE 65536

                             /* patch */
#define  n1   3
#define  NN1 10
float u1[NN1+1] = {-0.5, 0.0, 0.0, 0.0, 1.0, 2.0, 3.0, 4.0, 4.0, 4.0, 4.5};
#define  m1   3
#define  MM1 10
float v1[MM1+1] = {-0.5, 0.0, 0.0, 0.0, 1.0, 2.5, 2.5, 4.0, 4.0, 4.0, 4.5};
point3f cp1[NN1-n1][MM1-m1] =
  {{{-3.0,-3.0,-0.5},{-2.7,-3.0,-0.5},{-1.5,-3.0,-0.5},{0.0,-3.0,-0.5},{1.5,-3.0,-0.5},{2.7,-3.0,-0.5},{3.0,-3.0,-0.5}},
   {{-3.0,-2.7, 0.0},{-2.7,-2.7, 0.0},{-1.5,-2.7, 0.0},{0.0,-2.7, 0.0},{1.5,-2.7, 0.0},{2.7,-2.7, 0.0},{3.0,-2.7, 0.0}},
   {{-3.0,-1.5, 0.1},{-2.7,-1.5, 0.1},{-1.5,-1.5, 0.1},{0.0,-1.5, 0.1},{1.5,-1.5, 0.1},{2.7,-1.5, 0.1},{3.0,-1.5, 0.1}},
   {{-3.0, 0.0, 0.2},{-2.7, 0.0, 0.2},{-1.5, 0.0, 0.2},{0.0, 0.0, 0.2},{1.5, 0.0, 0.2},{2.7, 0.0, 0.2},{3.0, 0.0, 0.2}},
   {{-3.0,+1.5, 0.1},{-2.7,+1.5, 0.1},{-1.5,+1.5, 0.1},{0.0,+1.5, 0.1},{1.5,+1.5, 0.1},{2.7,+1.5, 0.1},{3.0,+1.5, 0.1}},
   {{-3.0,+2.7, 0.0},{-2.7,+2.7, 0.0},{-1.5,+2.7, 0.0},{0.0,+2.7, 0.0},{1.5,+2.7, 0.0},{2.7,+2.7, 0.0},{3.0,+2.7, 0.0}},
   {{-3.0,+3.0,-0.5},{-2.7,+3.0,-0.5},{-1.5,+3.0,-0.5},{0.0,+3.0,-0.5},{1.5,+3.0,-0.5},{2.7,+3.0,-0.5},{3.0,+3.0,-0.5}}};

                             /* patch boundary domain */
#define nt1a  3
#define NNt1a 10
float ut1a[NNt1a+1] =
  {-0.5, 0.0, 0.0, 0.0, 1.4, 2.8, 4.2, 5.6, 5.6, 5.6, 6.1};
point2f cpt1a[NNt1a-nt1a] =
  {{4.0,0.4},{3.5,0.4},{2.8,0.8},{2.6,2.0},{2.8,3.2},{3.5,3.6},{4.0,3.6}};

point2f cpd1a[3] =
  {{4.0,3.6},{4.0,4.0},{1.2,4.0}};

#define nt1b  3
#define NNt1b 7
point3f cpt1b[nt1b+1] =
  {{1.2,4.0,1.0},{1.5,3.0,1.0},{0.5,2.5,0.75},{0.0,2.5,1.0}};

point2f cpd1b[4] = {{0.0,2.5},{0.0,0.0},{4.0,0.0},{4.0,0.4}};

point2f cpd1c[5] = {{0.3,1.1},{1.6,1.6},{2.1,1.1},{1.6,0.6},{0.3,1.1}};

#define nt1c 3
point3f cpt1c[4] = {{2.0,3.3,1.0},{0.6,1.65,0.5},{0.6,1.25,0.5},{2.0,2.5,1.0}};

#define nt1d 3
point3f cpt1d[4] = {{2.0,2.5,1.0},{1.25,1.25,0.5},{1.25,1.5,0.5},{2.0,3.0,1.0}};

point2f cpd1e[2] = {{2.0,3.0},{2.0,3.3}};

polycurvef boundary1[8] =
  {{false,2,nt1a,NNt1a,&ut1a[0],(float*)&cpt1a[0]},   /* B-spline curve */
   {false,2,   1,    2/*2*/,    NULL,(float*)&cpd1a[0]},   /* polyline */
   {false,3,nt1b,   -1,    NULL,(float*)&cpt1b[0]},   /* Bezier curve */
   {true, 2,   1,    3,    NULL,(float*)&cpd1b[0]},   /* polyline */
   {true, 2,   1,    4,    NULL,(float*)&cpd1c[0]},   /* polyline */
   {false,3,nt1c,   -1,    NULL,(float*)&cpt1c[0]},   /* Bezier curve */
   {false,3,nt1d,   -1,    NULL,(float*)&cpt1d[0]},   /* Bezier curve */
   {true, 2,   1,    1,    NULL,(float*)&cpd1e[0]}};  /* line */

point2f frame[3] = {{100.0,100.0},{250.0,0.0},{0.0,250.0}};

CameraRecf CPos;

static void SetupCamera ( void )
{
  vector3f v;

  CameraInitFramef ( &CPos, false, true, 2000, 1500, 1200, 0, 1.0, 4 );
  CameraInitPosf ( &CPos );

  SetVector3f ( &v, 0.0, 0.0, -50.0 );
  CameraMoveGf ( &CPos, &v );
  CameraRotXCf ( &CPos, -0.65*PI );
  CameraRotZGf ( &CPos, -0.26*PI );
  CameraZoomf ( &CPos, 5.1 );
} /*SetupCamera*/

/* ///////////////////////////////////////////////////////////////////////// */
static void MapPoint ( const point2f *frame, const point2f *p, point2f *q )
{
  AddVector2Mf ( &frame[0], &frame[1], p->x, q );
  AddVector2Mf ( q, &frame[2], p->y, q );
} /*MapPoint*/

static void DrawLine1 ( point2f *p0, point2f *p1, int index )
{
  point2f q0, q1;

  if ( index & 1 ) {
    ps_Set_Line_Width ( (float)(2.0*abs(index)) );
    MapPoint ( frame, p0, &q0 );
    MapPoint ( frame, p1, &q1 );
    ps_Draw_Line ( q0.x, q0.y, q1.x, q1.y );
  }
} /*DrawLine1*/

static void DrawCurve1 ( int dim, int degree, const float *cp )
{
#define DENS 50
  int i, size;
  float t;
  point2f *c, p;

  ps_Set_Line_Width ( 6.0 );
  if ( degree == 1 ) {
    c = pkv_GetScratchMem ( size = 2*sizeof(point2f) );
    if ( c ) {
      if ( dim == 2 ) {
        MapPoint ( frame, (point2f*)cp, &c[0] );
        MapPoint ( frame, (point2f*)&cp[2], &c[1] );
      }
      else if ( dim == 3 ) {
        Point3to2f ( (point3f*)cp, &p );
        MapPoint ( frame, &p, &c[0] );
        Point3to2f ( (point3f*)&cp[3], &p );
        MapPoint ( frame, &p, &c[1] );
      }
      else goto out;
      ps_Draw_Line ( c[0].x, c[0].y, c[1].x, c[1].y );
      goto out;
    }
  }
  else if ( degree > 1 ) {
    c = pkv_GetScratchMem ( size = (DENS+1)*sizeof(point2f) );
    if ( c ) {
      if ( dim == 2 ) {
        for ( i = 0; i <= DENS; i++ ) {
          t = (float)i/(float)DENS;
          mbs_BCHornerC2f ( degree, cp, t, &p );
          MapPoint ( frame, &p, &c[i] );
        }
      }
      else if ( dim == 3 ) {
        for ( i = 0; i <= DENS; i++ ) {
          t = (float)i/(float)DENS;
          mbs_BCHornerC2Rf ( degree, (point3f*)cp, t, &p );
          MapPoint ( frame, &p, &c[i] );
        }
      }
      else goto out;
      ps_Draw_Polyline2f ( DENS+1, c );
out:
      pkv_FreeScratchMem ( size );
    }
  }
#undef DENS
} /*DrawCurve1*/

static void DrawLine2 ( point2f *p0, point2f *p1, int index )
{
#define LGT 0.05
  void     *sp;
  int      i, k;
  float    t, d;
  vector2f v;
  point2f  q, *c;
  point3f  p, r;

  if ( index & 1 ) {
    ps_Set_Line_Width ( (float)(2.0*abs(index)) );
    SubtractPoints2f ( p1, p0, &v );
    d = (float)sqrt ( DotProduct2f(&v,&v) );
    k = (int)(d/LGT+0.5);
    sp = pkv_GetScratchMemTop ();
    c = (point2f*)pkv_GetScratchMem ( (k+1)*sizeof(point2f) );
    for ( i = 0; i <= k; i++ ) {
      t = (float)i/(float)k;
      InterPoint2f ( p0, p1, t, &q );
      mbs_deBoorP3f ( n1, NN1, u1, m1, MM1, v1, 3*(MM1-m1), &cp1[0][0],
                      q.x, q.y, &p );
      CameraProjectPoint3f ( &CPos, &p, &r );
      c[i].x = r.x;  c[i].y = r.y;
    }
    ps_Draw_Polyline2f ( k+1, c );
    pkv_SetScratchMemTop ( sp );
  }
#undef LGT
} /*DrawLine2*/

static void DrawCurve2 ( int dim, int degree, const float *cp )
{
#define DENS 50
  int i, size;
  float t;
  point2f *c, q;
  point3f p, r;

  if ( degree >= 1 ) {
    ps_Set_Line_Width ( 6.0 );
    c = pkv_GetScratchMem ( size = (DENS+1)*sizeof(point2f) );
    if ( c ) {
      if ( dim == 2 ) {
        for ( i = 0; i <= DENS; i++ ) {
          t = (float)i/(float)DENS;
          mbs_BCHornerC2f ( degree, cp, t, &q );
          mbs_deBoorP3f ( n1, NN1, u1, m1, MM1, v1, 3*(MM1-m1), &cp1[0][0],
                          q.x, q.y, &p );
          CameraProjectPoint3f ( &CPos, &p, &r );
          c[i].x = r.x;  c[i].y = r.y;
        }
      }
      else if ( dim == 3 ) {
        for ( i = 0; i <= DENS; i++ ) {
          t = (float)i/(float)DENS;
          mbs_BCHornerC2Rf ( degree, (point3f*)cp, t, &q );
          mbs_deBoorP3f ( n1, NN1, u1, m1, MM1, v1, 3*(MM1-m1), &cp1[0][0],
                          q.x, q.y, &p );
          CameraProjectPoint3f ( &CPos, &p, &r );
          c[i].x = r.x;  c[i].y = r.y;
        }
      }
      else goto out;
      ps_Draw_Polyline2f ( DENS+1, c );
out:
      pkv_FreeScratchMem ( size );
    }
  }
#undef DENS
} /*DrawCurve2*/

static void test ( void )
{
  SetupCamera ();
  ps_WriteBBox ( 11, 12, 353, 140 );
  ps_OpenFile ( fn, 600 );

  ps_Write_Command ( "1 setlinecap" );
  mbs_DrawTrimBSPatchDomf ( n1, NN1, u1, m1, MM1, v1, 8, boundary1,
                            1, 2.0, 2.0, 1, 2.0, 2.0,
                            20, NULL, DrawLine1, DrawCurve1 );

  mbs_DrawTrimBSPatchDomf ( n1, NN1, u1, m1, MM1, v1, 8, boundary1,
                            4, 0.1, 2.0, 5, 0.3, 2.0,
                            20, NULL, DrawLine2, DrawCurve2 );

  ps_CloseFile ();
} /*test*/

int main ()
{
  pkv_InitScratchMem ( SCRATCHSIZE );
  test ();
  pkv_DestroyScratchMem ();
  exit ( 0 );
} /*main*/

