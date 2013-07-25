
/* //////////////////////////////////////////////////// */
/* This file is a part of the BSTools procedure package */
/* written by Przemyslaw Kiciak.                        */
/* //////////////////////////////////////////////////// */

#include <stdlib.h>
#include <math.h>
#include <stdio.h>

#include "pkgeom.h"
#include "camera.h"
#include "psout.h"

char fn[] = "stereo.ps";

CameraRecf CPos;

static void SetupCamera ( void )
{
  CameraInitFramef ( &CPos, false, true, 1800, 1200, 0, 0, 1.0, 4 );
  SetPoint3f ( &CPos.position, 11.1,-4.99,-8.02 );
  CPos.psi = 1.83;  CPos.theta = 0.915;  CPos.phi = -1.987;
  CPos.vd.persp.f = 3.5;
  CameraSetMappingf ( &CPos );
} /*SetupCamera*/

static void DrawLine ( point3f *p1, point3f *p2, char *attr )
{
  point3f q1, q2;
  float l, t1, t2;

  CameraProjectPoint3f ( &CPos, p1, &q1 );
  CameraProjectPoint3f ( &CPos, p2, &q2 );
  l = (float)sqrt( (q2.x-q1.x)*(q2.x-q1.x) + (q2.y-q1.y)*(q2.y-q1.y) );
  psl_SetLine ( q1.x, q1.y, q2.x, q2.y, 0.0, l );
  if ( attr[0] == 'a' ) t1 = (float)(0.7*arrowl);  else t1 = 0.0;
  if ( attr[1] == 'a' ) t2 = (float)(l-0.7*arrowl);  else t2 = l;
  switch ( attr[2] ) {
case ' ':
    ps_Set_Gray ( 0.0 );
    psl_Draw ( t1, t2, 4.0 );
    break;
case 'd':
    ps_Set_Gray ( 0.5 );
    psl_Draw ( t1, t2, 1.4 );
    break;
case 'e':
    ps_Set_Gray ( 0.0 );
    ps_Write_Command ( "[ 21 ] 0 setdash" );
    psl_Draw ( t1, t2, 1.4 );
    ps_Write_Command ( "[] 0 setdash" );
    break;
default:
    break;
  }
  if ( attr[0] == 'a' ) psl_Arrow ( 0.0, false );
  else if ( attr[0] == 'b' )
    ps_Fill_Circle ( q1.x, q1.y, 10.0 );
  if ( attr[1] == 'a' ) psl_Arrow ( l, true );
  else if ( attr[1] == 'b' )
    ps_Fill_Circle ( q2.x, q2.y, 10.0 );
} /*DrawLine*/

static void dl ( float x1, float y1, float z1,
                 float x2, float y2, float z2, char *attr )
{
  point3f p1, p2;

  SetPoint3f ( &p1, x1, y1, z1 );
  SetPoint3f ( &p2, x2, y2, z2 );
  DrawLine ( &p1, &p2, attr );
} /*dl*/

int main ( void )
{
  ps_WriteBBox ( 3, 0, 268, 150 );
  ps_OpenFile ( fn, 600 );
  SetupCamera ();

/* osie XYW */
  dl ( -1.2,0.0,0.0, 1.2,0.0,0.0, " a " );
  dl ( -0.2,-0.8,0.0, -0.2,0.82,0.0, " a " );
  dl ( 0.2,-0.8,0.0, 0.2,0.82,0.0, " a " );
  dl ( -0.2,0.0,-0.2, -0.2,0.0,5.2, " a " );
  dl ( 0.2,0.0,-0.2, 0.2,0.0,5.2, " a " );

/* okno */
  dl ( -0.6,0.45,1.5, 0.6,0.45,1.5, "   " );
  dl ( 0.6,0.45,1.5, 0.6,-0.45,1.5, "   " );
  dl ( 0.6,-0.45,1.5, -0.6,-0.45,1.5, "   " );
  dl ( -0.6,-0.45,1.5, -0.6,0.45,1.5, "   " );

/* ostroslupy widzenia */
  dl ( 0.2,0,0, -1.4,0.9,3, "  d" );
  dl ( 0.2,0,0, 1.0,0.9,3, "  d" );
  dl ( 0.2,0,0, 1.0,-0.9,3, "  d" );
  dl ( 0.2,0,0, -1.4,-0.9,3, "b d" );

  dl ( -1.4,-0.9,3, -1.4,0.9,3, "  e" );
  dl ( -1.4,0.9,3, 1.0,0.9,3, "  e" );
  dl ( 1.0,0.9,3, 1.0,-0.9,3, "  e" );
  dl ( 1.0,-0.9,3, -1.4,-0.9,3, "  e" );

  dl ( -0.2,0,0, -1.0,0.9,3, "  d" );
  dl ( -0.2,0,0, 1.4,0.9,3, "  d" );
  dl ( -0.2,0,0, 1.4,-0.9,3, "  d" );
  dl ( -0.2,0,0, -1.0,-0.9,3, "b d" );

  dl ( -1.0,-0.9,3, -1.0,0.9,3, "  e" );
  dl ( -1.0,0.9,3, 1.4,0.9,3, "  e" );
  dl ( 1.4,0.9,3, 1.4,-0.9,3, "  e" );
  dl ( 1.4,-0.9,3, -1.0,-0.9,3, "  e" );

  ps_CloseFile ();
  printf ( fn );

  exit ( 0 );
} /*main*/

