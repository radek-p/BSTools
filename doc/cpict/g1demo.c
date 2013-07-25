
/* //////////////////////////////////////////////////// */
/* This file is a part of the BSTools procedure package */
/* written by Przemyslaw Kiciak.                        */
/* //////////////////////////////////////////////////// */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#include "pkvaria.h"
#include "pknum.h"
#include "pkgeom.h"
#include "camera.h"
#include "multibs.h"
#include "bookg1holef.h"
#include "datagen.h"
#include "psout.h"


char fn1[] = "g1patches1.ps";
char fn2[] = "g1patches2.ps";

point3f Bezpt[MAX_K][3][16];
point3f TBezpt[MAX_K][36];

CameraRecf CPos;

static point3f* GetSurrPatchCP ( int i, int j )
{
  return &Bezpt[i][j][0];
} /*GetSurrPatchCP*/

static void GetSurroundingPatch ( int i, int j, point3f *bcp )
{
  int     k;
  int     *ind;
  point3f *q;
  void    *sp;
  float   kn[8] = {0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 6.0};

  sp  = pkv_GetScratchMemTop ();
  ind = pkv_GetScratchMem ( 16*sizeof(int) );
  q   = pkv_GetScratchMem ( 16*sizeof(point3f) );
  if ( !ind || !q )
    exit ( 1 );

  GetBspInd ( i, j, ind );
  for ( k = 0; k < 16; k++ )
    q[k] = Bpt[ind[k]];
  mbs_BSPatchToBezf ( 3, 3, 7, kn, 3, 7, kn, 12, (float*)q,
                      NULL, NULL, NULL, NULL, NULL, NULL, 12, (float*)bcp );

  pkv_SetScratchMemTop ( sp );
} /*GetSurroundingPatch*/

static void SetupObject ( void )
{
  int i, j;

  InitHole ( hole_k = 5, 0.25, 0.0, 0.6 );
  for ( i = 0; i < hole_k; i++ )
    for ( j = 0; j < 3; j++ )
      GetSurroundingPatch ( i, j, &Bezpt[i][j][0] );
} /*SetupObject*/

static void SetupCamera ( void )
{
  vector3f v;

  CameraInitFramef ( &CPos, false, true, 1600, 1600, 0, 0, 1.0, 4 );
  CameraInitPosf ( &CPos );
  SetVector3f ( &v, 0.0, 0.0, -40.0 );
  CameraMoveGf ( &CPos, &v );
  CameraRotXGf ( &CPos, 0.15*PI );
  CameraZoomf ( &CPos, 4.0 );
} /*SetupCamera*/

static void DrawBezNet ( int n, int m, const point3f *cp, vector3f *tr )
{
  void    *sp;
  point2f *q;
  point3f p, r;
  int     i, j;

  sp = pkv_GetScratchMemTop ();
  q = pkv_GetScratchMem ( (n+1)*(m+1)*sizeof(point2f) );
  if ( !q )
    exit ( 1 );
  for ( i = 0; i < (n+1)*(m+1); i++ ) {
    AddVector3f ( &cp[i], tr, &r );
    CameraProjectPoint3f ( &CPos, &r, &p );
    SetPoint2f ( &q[i], p.x, p.y );
  }
  ps_Set_Line_Width ( 2.0 );
  for ( i = 0; i <= n; i++ )
    for ( j = 0; j < m; j++ )
      ps_Draw_Line ( q[(m+1)*i+j].x, q[(m+1)*i+j].y,
                     q[(m+1)*i+j+1].x, q[(m+1)*i+j+1].y );
  for ( i = 0; i < n; i++ )
    for ( j = 0; j <= m; j++ )
      ps_Draw_Line ( q[(m+1)*i+j].x, q[(m+1)*i+j].y,
                     q[(m+1)*(i+1)+j].x, q[(m+1)*(i+1)+j].y );
  for ( i = 0; i < (n+1)*(m+1); i++ )
    ps_Fill_Circle ( q[i].x, q[i].y, 6.0 );

  pkv_SetScratchMemTop ( sp );
} /*DrawBezNet*/

static void DrawBezPatch ( int n, int m, const point3f *cp )
{
#define D   6
#define DD 24
  void    *sp;
  int     i, j;
  float   u, v;
  point3f p, q;
  point2f *c;

  sp = pkv_GetScratchMemTop ();
  c = (point2f*)pkv_GetScratchMem ( (DD+1)*sizeof(point2f) );

  for ( i = 0; i <= D; i++ ) {
    u = (float)i/(float)D;
    for ( j = 0; j <= DD; j++ ) {
      v = (float)j/(float)DD;
      mbs_BCHornerP3f ( n, m, cp, u, v, &p );
      CameraProjectPoint3f ( &CPos, &p, &q );
      SetPoint2f ( &c[j], q.x, q.y );
    }
    if ( i == 0 || i == D )
      ps_Set_Line_Width ( 6.0 );
    else if ( i == 1 )
      ps_Set_Line_Width ( 2.0 );
    ps_Draw_Polyline2f ( DD+1, c );
  }
  for ( i = 0; i <= D; i++ ) {
    v = (float)i/(float)D;
    for ( j = 0; j <= DD; j++ ) {
      u = (float)j/(float)DD;
      mbs_BCHornerP3f ( n, m, cp, u, v, &p );
      CameraProjectPoint3f ( &CPos, &p, &q );
      SetPoint2f ( &c[j], q.x, q.y );
    }
    if ( i == 0 || i == D )
      ps_Set_Line_Width ( 6.0 );
    else if ( i == 1 )
      ps_Set_Line_Width ( 2.0 );
    ps_Draw_Polyline2f ( DD+1, c );
  }

  pkv_SetScratchMemTop ( sp );
#undef DD
#undef D
} /*DrawBezPatch*/

static void DrawPict1 ( void )
{
  int      i, j, k;
  point3f  p, q;
  char     s[30];
  vector3f tr;

  ps_WriteBBox ( 22, 3, 346, 152 );
  ps_OpenFile ( fn1, 600 );
  ps_DenseScreen ();
  ps_Write_Command ( "1 setlinecap" );
  ps_Write_Command ( "/Times-Roman findfont 50 scalefont setfont" );

  for ( i = 0; i < hole_k; i++ )
    for ( j = 0; j < 3; j++ ) {
      ps_Set_Gray ( 0.68 );
      DrawBezPatch ( 3, 3, &Bezpt[i][j][0] );
      ps_Set_Gray ( 0.0 );
      mbs_BCHornerP3f ( 3, 3, &Bezpt[i][j][0], 0.5, 0.5, &p );
      CameraProjectPoint3f ( &CPos, &p, &q );
      ps_MoveTo ( q.x-30, q.y-12 );
      sprintf ( s, "(%d,%d) show", i, j );
      ps_Write_Command ( s );
    }
  ps_GSave ();
  ps_BeginDict ( 10 );
  ps_Write_Command ( "1400 -80 translate" );
  ps_Write_Command ( "/setgr { setgray } bind def" );
  ps_Write_Command ( "/setgray { 1.0 add 0.5 mul setgr } def" );
  ps_Write_Command ( "0 setgray" );
  for ( i = 0; i < hole_k; i++ )
    for ( j = 1; j < 3; j++ ) {
      mbs_BCHornerP3f ( 3, 3, &Bezpt[i][j][0], 0.5, 0.5, &tr );
      MultVector3f ( 0.35, &tr, &tr );
      ps_Write_Command ( "0 setgray" );
      DrawBezNet ( 3, 3, &Bezpt[i][j][0], &tr );
      ps_Write_Command ( "0 setgr" );
      for ( k = 0; k < 8; k++ ) {
        AddVector3f ( &Bezpt[i][j][k], &tr, &p );
        CameraProjectPoint3f ( &CPos, &p, &q );
        if ( k <= 9 )
          ps_MoveTo ( q.x-10, q.y-12 );
        else
          ps_MoveTo ( q.x-25, q.y-12 );
        sprintf ( s, "(%d) show", k );
        ps_Write_Command ( s );
      }
    }
  ps_EndDict ();
  ps_GRestore ();

  printf ( "%s\n", fn1 );
  ps_CloseFile ();
} /*DrawPict1*/

static void DrawPict2 ( void )
{
  int i, j;

  ps_WriteBBox ( 22, 3, 172, 148 );
  ps_OpenFile ( fn2, 600 );

  ps_DenseScreen ();
  ps_Write_Command ( "1 setlinecap" );

  ps_Set_Gray ( 0.68 );
  for ( i = 0; i < hole_k; i++ )
    for ( j = 0; j < 3; j++ ) {
      DrawBezPatch ( 3, 3, &Bezpt[i][j][0] );
    }
  ps_Set_Gray ( 0.5 );
  for ( i = 0; i < hole_k; i++ )
    DrawBezPatch ( 5, 5, &TBezpt[i][0] );

  printf ( "%s\n", fn2 );
  ps_CloseFile ();
} /*DrawPict2*/

int main ( void )
{
  pkv_InitScratchMem ( 1048576 );

  SetupCamera ();
  SetupObject ();
  FillG1Holef ( hole_k, GetSurrPatchCP, 1.0, 1.0, &TBezpt[0][0] );
  DrawPict1 ();
  DrawPict2 ();

  pkv_DestroyScratchMem ();
  exit ( 0 );
} /*main*/

