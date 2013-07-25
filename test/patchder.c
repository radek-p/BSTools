
/* //////////////////////////////////////////////////// */
/* This file is a part of the BSTools procedure package */
/* written by Przemyslaw Kiciak.                        */
/* //////////////////////////////////////////////////// */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "pkvaria.h"
#include "pkgeom.h"
#include "camera.h"
#include "pknum.h"
#include "multibs.h"
#include "psout.h"

char fn[] = "patchder.ps";

#define NN 3
#define MM 4

point3f cp[NN+1][MM+1] =
  {{{-1.0,-1.0,0.0},{-0.5,-1.0,0.0},{0.0,-1.0,0.0},{0.5,-1.0,0.0},{1.0,-1.0,0.0}},
   {{-1.0,-0.4,0.0},{-0.5,-0.4,-0.4},{0.4,-0.4,0.25},{0.5,-0.4,0.3},{1.0,-0.4,0.0}},
   {{-1.0, 0.4,0.0},{-0.5, 0.4,0.3},{0.0, 0.4,0.25},{0.5, 0.4,-0.4},{1.0, 0.4,0.0}},
   {{-1.0, 1.0,0.0},{-0.5, 1.0,0.0},{0.0, 1.0,0.0},{0.5, 1.0,0.0},{1.0, 1.0,0.0}}};

point3f cpd[(NN+2)*(MM+2)];

point3f cq0[3][3] = {{{-0.5,-0.5,1.0},{0.0,-0.5,0.0},{0.5,-0.5,0.0}},
                     {{-0.5, 0.0,1.0},{0.0, 0.0,0.0},{0.5, 0.0,0.0}},
                     {{-0.5, 0.5,1.0},{0.0, 0.5,0.0},{0.5, 0.5,0.0}}};
point3f cq1[2][3] = {{{-0.5,-0.5,1.0},{0.0,-0.5,0.0},{0.5,-0.5,0.0}},
                     {{-0.5, 0.5,1.0},{0.0, 0.5,0.0},{0.5, 0.5,0.0}}};
point3f cq2[3][2] = {{{-0.5,-0.5,1.0},{0.5,-0.5,1.0}},
                     {{-0.5, 0.0,0.0},{0.5, 0.0,0.0}},
                     {{-0.5, 0.5,0.0},{0.5, 0.5,0.0}}};
point3f cq3[2][2] = {{{-0.5,-0.5,1.0},{0.5,-0.5,0.0}},
                     {{-0.5, 0.5,0.0},{0.5, 0.5,0.0}}};

CameraRecf CPos;

static void SetupCamera ( void )
{
  vector3f v;

  CameraInitFramef ( &CPos, false, true, 1200, 900, 0, 0, 1.0, 4 );
  CameraInitPosf ( &CPos );

  SetVector3f ( &v, 0.0, 0.0, -20.0 );
  CameraMoveGf ( &CPos, &v );
  CameraRotXCf ( &CPos, -0.65*PI );
  CameraRotZGf ( &CPos, -0.14*PI );
  CameraZoomf ( &CPos, 5.0 );
} /*SetupCamera*/

static void DisplayPatch ( int n, int m, point3f *cp, int d1, int d2 )
{
  int     i, j;
  float   u, v;
  point2f *cc;
  point3f p, q;
  void    *sp;

  sp = pkv_GetScratchMemTop ();
  cc = pkv_GetScratchMem ( 4*max(d1,d2)*sizeof(point2f) );
  if ( !cc )
    exit ( 1 );

  ps_Set_Line_Width ( 6 );
  for ( i = 0; i <= d1; i++ ) {
    u = (float)i/(float)d1;
    for ( j = 0; j <= 4*d2; j++ ) {
      v = (float)j/(float)(4*d2);
      mbs_BCHornerP3f ( n, m, cp, u, v, &p );
      CameraProjectPoint3f ( &CPos, &p, &q );
      cc[j].x = q.x;  cc[j].y = q.y;
    }
    ps_Draw_Polyline2f ( 4*d2+1, cc );
    if ( i >= d1-1 )
      ps_Set_Line_Width ( 6 );
    else
      ps_Set_Line_Width ( 2 );
  }
  for ( i = 0; i <= d2; i++ ) {
    v = (float)i/(float)d2;
    for ( j = 0; j <= 4*d1; j++ ) {
      u = (float)j/(float)(4*d1);
      mbs_BCHornerP3f ( n, m, cp, u, v, &p );
      CameraProjectPoint3f ( &CPos, &p, &q );
      cc[j].x = q.x;  cc[j].y = q.y;
    }
    ps_Draw_Polyline2f ( 4*d1+1, cc );
    if ( i == d1-1 )
      ps_Set_Line_Width ( 6 );
    else
      ps_Set_Line_Width ( 2 );
  }

  pkv_SetScratchMemTop ( sp );
} /*DisplayPatch*/

static void DrawVect ( point3f *p, vector3f *v )
{
  point3f q, r, s;

  CameraProjectPoint3f ( &CPos, p, &r );
  AddVector3f ( p, v, &q );
  CameraProjectPoint3f ( &CPos, &q, &s );
  psl_SetLine ( r.x, r.y, s.x, s.y, 0.0, 1.0 );
  psl_Draw ( 0.0, 0.95, 4.0 );
  psl_Arrow ( 1.0, true );
} /*DrawVect*/

static void TestDer ( int n, int m, point3f *cp )
{
  point3f  p;
  vector3f du, dv, duu, duv, dvv;
  float    u, v;

  u = 0.3;  v = 0.4;
  ps_Set_RGB ( 1.0, 0.0, 0.0 );
  mbs_BCHornerDerP3f ( n, m, cp, u, v, &p, &du, &dv );
  DrawVect ( &p, &du );
  DrawVect ( &p, &dv );
  ps_Set_RGB ( 0.0, 0.7, 0.0 );
  mbs_BCHornerDer2P3f ( n, m, cp, u, v, &p, &du, &dv, &duu, &duv, &dvv );
  DrawVect ( &p, &du );
  DrawVect ( &p, &dv );
  ps_Set_RGB ( 0.0, 0.0, 1.0 );
  DrawVect ( &p, &duu );
  ps_Set_RGB ( 0.0, 0.0, 0.7 );
  DrawVect ( &p, &duv );
  ps_Set_RGB ( 0.0, 0.0, 0.5 );
  DrawVect ( &p, &dvv );
} /*TestDer*/

static void TestPatch ( int n, int m, point3f *cp )
{
  int nn, mm;

  ps_Write_Command ( "1 setlinecap" );
  ps_Set_Gray ( 0.68 );
  DisplayPatch ( n, m, &cp[0], 10, 10 );
  TestDer ( n, m, &cp[0] );
  ps_GSave ();
  ps_Write_Command ( "1200 0 translate" );
  mbs_BCDegElevP3f ( n, m, &cp[0], 1, 1, &nn, &mm, &cpd[0] );
  ps_Set_Gray ( 0.68 );
  DisplayPatch ( nn, mm, &cpd[0], 10, 10 );
  TestDer ( nn, mm, &cpd[0] );
  ps_GRestore ();
} /*TestPatch*/

static void OtherTest ()
{
  point3f  p;
  vector3f du, dv, duu, duv, dvv;
  char     format[] =
    "duu:(%5.2f,%5.2f,%5.2f), duv:(%5.2f,%5.2f,%5.2f), dvv:(%5.2f,%5.2f,%5.2f)\n";

  mbs_BCHornerDer2P3f ( 2, 2, cq0, 0.2, 0.3, &p, &du, &dv, &duu, &duv, &dvv );
  printf ( format, duu.x, duu.y, duu.z, duv.x, duv.y, duv.z, dvv.x, dvv.y, dvv.z );
  mbs_BCHornerDer2P3f ( 1, 2, cq1, 0.2, 0.3, &p, &du, &dv, &duu, &duv, &dvv );
  printf ( format, duu.x, duu.y, duu.z, duv.x, duv.y, duv.z, dvv.x, dvv.y, dvv.z );
  mbs_BCHornerDer2P3f ( 2, 1, cq2, 0.2, 0.3, &p, &du, &dv, &duu, &duv, &dvv );
  printf ( format, duu.x, duu.y, duu.z, duv.x, duv.y, duv.z, dvv.x, dvv.y, dvv.z );
  mbs_BCHornerDer2P3f ( 1, 1, cq3, 0.2, 0.3, &p, &du, &dv, &duu, &duv, &dvv );
  printf ( format, duu.x, duu.y, duu.z, duv.x, duv.y, duv.z, dvv.x, dvv.y, dvv.z );
} /*OtherTest*/

int main ()
{
  pkv_InitScratchMem ( 262144 );

  ps_WriteBBox ( 13, 22, 302, 88 );
  ps_OpenFile ( fn, 600 );
  SetupCamera ();
  TestPatch ( NN-1, MM, &cp[0][0] );
  ps_CloseFile ();

  OtherTest ();

  pkv_DestroyScratchMem ();
  exit ( 0 );
} /*main*/

