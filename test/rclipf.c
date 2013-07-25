
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "pkvaria.h"
#include "pkgeom.h"
#include "camera.h"
#include "multibs.h"
#include "psout.h"

char fn[] = "rclipf.ps";

/* ///////////////////////////////////////////////////////////////////////// */
CameraRecf CPos;

point4f cpts[5] =
  {{7.0,0.0,0.0,1.0},{15.0,-20.0,0.0,0.5},{-4.0,1.0,0.0,1.0},
   {1.0,17.0,2.0,0.7},{-8.0,0.0,-2.0,1.0}};

static void SetupCamera ( void )
{
  vector3f v;

  CameraInitFramef ( &CPos, false, true, 1600, 1200, 200, 200, 1.0, 4 );
  CameraInitPosf ( &CPos );
  SetVector3f ( &v, 0.0, 0.0, -20.0 );
  CameraMoveGf ( &CPos, &v );
} /*SetupCamera*/

static void DrawCurve ( int degree, const point4f *cp )
{
#define DD 100
  void    *sp;
  point2f *cc;
  point3f p, q;
  int     i;

  sp = pkv_GetScratchMemTop ();
  cc = pkv_GetScratchMem ( (DD+1)*sizeof(point2f) );
  if ( cc ) {
    for ( i = 0; i <= DD; i++ ) {
      mbs_BCHornerC3Rf ( degree, cp, (float)i/(float)DD, &p );
      CameraProjectPoint3f ( &CPos, &p, &q );
      SetPoint2f ( &cc[i], q.x, q.y );
    }
    ps_Draw_Polyline2f ( DD+1, cc );
  }
  pkv_SetScratchMemTop ( sp );
#undef DD
} /*DrawCurve*/

static void MakePicture ( void )
{
  ps_OpenFile ( fn, 600 );
  ps_Set_Line_Width ( 1 );
  ps_Draw_Rect ( CPos.width, CPos.height, CPos.xmin, CPos.ymin );
  ps_Set_RGB ( 1.0, 0.0, 0.0 );
  ps_Set_Line_Width ( 2 );
  DrawCurve ( 4, cpts );

  ps_Set_RGB ( 0.0, 1.0, 0.0 );
  ps_Set_Line_Width ( 4 );
  mbs_ClipBC3Rf ( 4, CPos.cplane, 4, cpts, DrawCurve );
  
  ps_CloseFile ();
  printf ( "%s\n", fn );
} /*MakePicture*/

int main ()
{
  pkv_InitScratchMem ( 65536 );
  SetupCamera ();
  MakePicture ();
  printf ( "scratch memory used: %d\n", (int)pkv_MaxScratchTaken () );
  pkv_DestroyScratchMem ();
  exit ( 0 );
} /*main*/

