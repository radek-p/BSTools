
#include <pthread.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

#include "pkvaria.h"
#include "pkgeom.h"
#include "pknum.h"
#include "multibs.h"
#include "raybez.h"
#include "camera.h"
#include "psout.h"

char fn[] = "intersd.ps";

point4d f1[9] =
  {{0.0,0.0,0.7,1.0},{0.5,0.0,0.3,1.0},{1.0,0.0,-0.15,1.0},
   {0.0,0.5,0.5,1.0},{0.5,0.5,0.55,1.0},{1.0,0.5,-0.15,1.0},
   {0.0,1.0,0.05,1.0},{0.5,1.0,0.0,1.0},{1.0,1.0,-0.225,1.0}};
point4d f2[9] =
  {{0.0,0.0,-0.25,1.0},{0.5,0.0,-0.18,1.0},{1.0,0.0,-0.3,1.0},
   {0.0,0.5,0.025,1.0},{0.5,0.5,-0.24,1.0},{1.0,0.5,-0.26,1.0},
   {0.0,1.0,0.55,1.0},{0.5,1.0,0.66,1.0},{1.0,1.0,0.22,1.0}};

CameraRecd CPos;

double      x1, y1, x2, y2;

void SetCPos ( void )
{
  vector3d v;

  CameraInitFramed ( &CPos, false, true, 2780, 2100, 0, 0, 1.0, 4 );
  CameraInitPosd ( &CPos );
  SetVector3d ( &v, 0.0, 0.0, -12.0 );
  CameraMoveGd ( &CPos, &v );
  CameraRotZGd ( &CPos, 3*PI/5 );
  CameraRotXCd ( &CPos, -5*PI/8 );
  CameraZoomd ( &CPos, 4.5 );
} /*SetCPos*/

#define MAXARCS       10
#define MAXARCPOINTS 200
int      narc1;
int      nparc1[MAXARCS];
vector4d arc1[MAXARCPOINTS+1];

void ArcOut ( void *usrptr, rbiIntersArcd *arc, vector4d *ipt )
{
  int a;

  if ( narc1 < MAXARCS ) {
    if ( narc1 == 0 ) a = 0;
    else a = nparc1[narc1-1];
    if ( a+arc->npoints < MAXARCPOINTS ) {
      memcpy ( &arc1[a], ipt, arc->npoints*sizeof(vector4d) );
      nparc1[narc1] = a + arc->npoints;
    }
    narc1 ++;
  }
} /*ArcOut*/

#define density 48
/*#define density 24*/

void DrawGR ( double gr, point3d p[density+1][density+1] )
{
  point2d c[density+1];
  int     i, j;

  ps_Set_Gray ( gr );
  for ( i = 0; i <= density/4; i++ ) {
    if ( i == 0 || i == density/4 )
      ps_Set_Line_Width ( 6 );
    else
      ps_Set_Line_Width ( 2 );
    for ( j = 0; j <= density; j++ ) {
      c[j].x = p[4*i][j].x;
      c[j].y = p[4*i][j].y;
    }
    ps_Draw_Polyline2d ( density+1, c );
  }
  for ( j = 0; j <= density/4; j++ ) {
    if ( j == 0 || j == density/4 )
      ps_Set_Line_Width ( 6 );
    else
      ps_Set_Line_Width ( 2 );
    for ( i = 0; i <= density; i++ ) {
      c[i].x = p[i][4*j].x;
      c[i].y = p[i][4*j].y;
    }
    ps_Draw_Polyline2d ( density+1, c );
  }
} /*DrawGR*/

void DrawInsC ( double dx, double dy, int narc, int *nparc,
                vector4d *arc )
{
  int     i, j, k;
  point2d c[MAXARCPOINTS+1];
  point3d q;

  k = 0;
  for ( i = 0; i < narc; i++ ) {
    for ( j = k; j < nparc[i]; j++ ) {
      mbs_BCHornerP3Rd ( 2, 2, f1, arc[j].x, arc[j].y, &q );
      CameraProjectPoint3d ( &CPos, &q, &q );
      c[j-k].x = q.x+dx;
      c[j-k].y = q.y+dy;
    }
    ps_Draw_Polyline2d ( nparc[i]-k, c );
    k = nparc[i];
  }
} /*DrawInsC*/

void DrawGraphA ( double dx, double dy, void *f, double gr )
{
  int     i, j;
  point3d p[density+1][density+1];
  point3d q;

  for ( i = 0; i <= density; i++ )
    for ( j = 0; j <= density; j++ ) {
      mbs_BCHornerP3Rd ( 2, 2, f,
                         (double)i/(double)density, (double)j/(double)density, &q );
      CameraProjectPoint3d ( &CPos, &q, &p[i][j] );
      p[i][j].x = p[i][j].x+dx;
      p[i][j].y = p[i][j].y+dy;
    }

  DrawGR ( gr, p );

  ps_Set_Gray ( 0.0 );
  ps_Set_Line_Width ( 4 );
  DrawInsC ( dx, dy, narc1, nparc1, arc1 );
} /*DrawGraphA*/

int main ( void )
{
  setvbuf ( stdout, NULL, _IONBF, 0 );
  pkv_InitScratchMem ( 2*1048576 );
  SetCPos ();

  narc1 = 0;
  rbi_FindRBezIntersectiond ( 2, 2, f1, 2, 2, f2, 0.05, 20, ArcOut, NULL );

  ps_WriteBBox ( 0, 0, 235, 235  );
  ps_OpenFile ( fn, 600 );
  ps_Write_Command ( "1 setlinecap" );
  DrawGraphA ( -800.0, 0.0, f1, 0.72 );
  DrawGraphA ( -800.0, 0.0, f2, 0.5 );
  ps_CloseFile ();
  printf ( "%s\n", fn );

  pkv_DestroyScratchMem ();
  exit ( 0 );
} /*main*/

