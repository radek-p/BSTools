
/* //////////////////////////////////////////////////// */
/* This file is a part of the BSTools procedure package */
/* written by Przemyslaw Kiciak.                        */
/* //////////////////////////////////////////////////// */

#include <pthread.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <math.h>

#include "pkvaria.h"
#include "pknum.h"
#include "pkgeom.h"
#include "multibs.h"
#include "raybez.h"
#include "camera.h"
#include "psout.h"

/* ////////////////////////////////////////////////////////////////////////// */
#define WIDTH  640
#define HEIGHT 480
char fn[] = "bcoffsetf.ps";

CameraRecf CPos;
vector3f ld = { -0.5, -0.5, -1.0 };
BezCurveTreefp tree;

#define DEGREE 3
point3f cp[DEGREE+1] = {{-1.0},
    {-0.384426229508,0.0,0.677049180328},
    {0.338524590164,0.814754098361,0.0},
    {1.0}};

void SetupCamera ( void )
{
  CameraInitFramef ( &CPos, false, false, WIDTH, HEIGHT, 0, 0, 1.0, 4 );
  CameraInitPosf ( &CPos );
  SetPoint3f ( &CPos.position, -18.264203, 12.538737, 26.630790 );
  CPos.psi   = 1.373231;
  CPos.theta = 2.447704;
  CPos.phi   = 0.969175;
  CPos.vd.persp.f = 15.2;
  CameraSetMappingf ( &CPos );
  CameraRotYCf ( &CPos, PI );
} /*SetupCamera*/

static byte PixelColour ( ray3f *ray )
{
#define MAXLEVEL  26
#define MAXINTERS 10
  RayObjectIntersf inters[MAXINTERS];
  int ninters, i, j;
  float tmin;
  float gray, g0;

  rbez_FindRayBezcOffsetIntersf ( tree, ray, MAXLEVEL, MAXINTERS,
                                  &ninters, inters );
  if ( ninters ) {
    if ( ninters > 1 ) {
      j = 0;
      tmin = inters[0].t;
      for ( i = 0; i < ninters; i++ )
        if ( inters[i].t < tmin ) {
          tmin = inters[i].t;
          j = i;
        }
      if ( j )
        inters[0] = inters[j];
    }
    gray = 0.2;
    if ( DotProduct3f ( &inters[0].nv, &ray->v ) < 0.0 )
      g0 = DotProduct3f ( &ld, &inters[0].nv );
    else
      g0 = -DotProduct3f ( &ld, &inters[0].nv );
    if ( g0 > 0.0 )
      gray += 0.7*g0;
    if ( gray > 1.0 )
      gray = 1.0;
    return (byte)(255.0*gray);
  }
  else
    return 160;
#undef MAXINTERS
#undef MAXLEVEL
} /*PixelColour*/

void RayTrace ( void )
{
  int   x, y;
  ray3f ray;
  byte  buffer[WIDTH];

  CPos.upside = false;
  CameraSetMappingf ( &CPos );
  pkv_Tic ( NULL );
  ps_Init_BitmapP ( WIDTH, HEIGHT, 0, 0 );
  for ( y = 0; y < HEIGHT; y++ ) {
    printf ( "%4d\b\b\b\b", y );
    for ( x = 0; x < WIDTH; x++ ) {
      CameraRayOfPixelf ( &CPos, (float)x, (float)y, &ray );
      buffer[x] = PixelColour ( &ray );
    }
    ps_Out_LineP ( buffer );
  }
  printf ( "\ntime: %6.2f\n", pkv_Seconds ( pkv_Toc ( NULL ) ) );
} /*RayTrace*/

int main ( void )
{
  setvbuf ( stdout, NULL, _IONBF, 0 );
  pkv_InitScratchMem ( 655360 );
  SetupCamera ();
  NormalizeVector3f ( &ld );
  tree = rbez_NewBezCurveTreef ( 0, DEGREE, 0.0, 1.0, cp, 0.1 );
  ps_WriteBBox ( 0, 0, 154, 117 );
  ps_OpenFile ( fn, 300 );
  RayTrace ();
  ps_CloseFile ();

  pkv_DestroyScratchMem ();
  exit ( 0 );
} /*main*/

