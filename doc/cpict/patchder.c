
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

CameraRecf CPos;

static void SetupCamera ( void )
{
  vector3f v;

  CameraInitFramef ( &CPos, false, true, 1600, 1200, 0, 0, 1.0, 4 );
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

static void TestDer ( int n, int m, point3f *cp, boolean second )
{
  point3f  p;
  vector3f du, dv, duu, duv, dvv;
  float    u, v;

  u = 0.3;  v = 0.4;
  ps_Set_Gray ( 0.4 );
  mbs_BCHornerDer2P3f ( n, m, cp, u, v, &p, &du, &dv, &duu, &duv, &dvv );
  if ( !second ) {
    DrawVect ( &p, &du );
    DrawVect ( &p, &dv );
  }
  if ( second ) {
    ps_Set_Gray ( 0.0 );
    DrawVect ( &p, &duu );
    DrawVect ( &p, &duv );
    DrawVect ( &p, &dvv );
  }
} /*TestDer*/

int main ()
{
  pkv_InitScratchMem ( 262144 );
  ps_WriteBBox ( 17, 29, 370, 120 );
  ps_OpenFile ( fn, 600 );
  ps_DenseScreen ();
  SetupCamera ();
  ps_Write_Command ( "1 setlinecap" );

  ps_Set_Gray ( 0.68 );
  DisplayPatch ( NN, MM, &cp[0][0], 10, 10 );
  TestDer ( NN, MM, &cp[0][0], false );

  ps_GSave ();
  ps_Write_Command ( "1600 0 translate" );
  ps_Set_Gray ( 0.68 );
  DisplayPatch ( NN, MM, &cp[0][0], 10, 10 );
  TestDer ( NN, MM, &cp[0][0], true );
  ps_GRestore ();

  ps_CloseFile ();
  pkv_DestroyScratchMem ();
  exit ( 0 );
} /*main*/

