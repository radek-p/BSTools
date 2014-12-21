
#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <sys/times.h>

#include "pkvaria.h"
#include "pknum.h"
#include "pkgeom.h"
#include "multibs.h"
#include "raybez.h"
#include "camera.h"
#include "psout.h"

char fn[] = "testraybsd.ps";

int n = 3, m = 3;
int lknu = 10, lknv = 10;
double knots_u[11] = {-0.5, 0.0, 0.0, 0.0, 1.0, 2.0, 3.0, 4.0, 4.0, 4.0, 4.5};
double knots_v[11] = {-0.5, 0.0, 0.0, 0.0, 1.0, 2.5, 2.5, 4.0, 4.0, 4.0, 4.5};
point3d cpoints[] =
  {{-3.0,-3.0,-0.5}, {-2.7,-3.0,-0.5}, {-1.5,-3.0,-0.5}, {0.0,-3.0,-0.5},
   {1.5,-3.0,-0.5}, {2.7,-3.0,-0.5}, {3.0,-3.0,-0.5},
   {-3.0,-2.7}, {-2.7,-2.7}, {-1.5,-2.7}, {0.0,-2.7}, {1.5,-2.7},
   {2.7,-2.7}, {3.0,-2.7},
   {-3.0,-1.5,0.1}, {-2.7,-1.5,0.1}, {-1.5,-1.5,0.1}, {0.0,-1.5,0.1},
   {1.5,-1.5,2.5}, {2.7,-1.5,0.1}, {3.0,-1.5,0.1},
   {-3.0,0.0,0.2}, {-2.7,0.0,0.2}, {-1.5,0.0,0.2}, {0.0,0.0,0.2},
   {1.5,0.0,0.2}, {2.7,0.0,0.2}, {3.0,0.0,0.2},
   {-3.0,1.5,0.1}, {-2.7,1.5,0.1}, {-1.5,1.5,0.1}, {0.0,1.5,0.1},
   {1.5,1.5,0.1}, {2.7,1.5,0.1}, {3.0,1.5,0.1},
   {-3.0,2.7}, {-2.7,2.7}, {-1.5,2.7}, {0.0,2.7}, {1.5,2.7}, {2.7,2.7},
   {3.0,2.7},
   {-3.0,3.0,-0.5}, {-2.7,3.0,-0.5}, {-1.5,3.0,-0.5}, {0.0,3.0,-0.5},
   {1.5,3.0,-0.5}, {2.7,3.0,-0.5}, {3.0,3.0,-0.5}};

#define WIDTH  640
#define HEIGHT 480

CameraRecd CPos;
vector3d ld;

BezPatchTreedp  tree;

static void SetCPos ( void )
{
  CameraInitFramed ( &CPos, false, false, WIDTH, HEIGHT, 0, 0, 1.0, 4 );
  CameraInitPosd ( &CPos );
  SetPoint3d ( &CPos.position, -15.31251039605, -72.93132382473, 49.231894764605 );
  CPos.psi   = 2.963570202665;
  CPos.theta = 2.1560152051;
  CPos.phi   = 2.934640744727;
  CPos.vd.persp.f = 10.0;
  CameraSetMappingd ( &CPos );
} /*SetCPos*/

static void SetLightDir ( void )
{
  SetVector3d ( &ld, 0.5, 0.5, -1.0 );
  NormalizeVector3d ( &ld );
} /*SetLightDir*/

static byte PixelColour ( ray3d *ray )
{
#define MAXLEVEL  16
#define MAXINTERS 10
  RayObjectIntersd inters[MAXINTERS];
  int ninters, i, j;
  double tmin;
  double gray, g0, l0;

  rbez_FindRayBezPatchIntersd ( tree, ray, MAXLEVEL, MAXINTERS,
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
    l0 = DotProduct3d ( &ray->v, &inters[0].nv );
    g0 = DotProduct3d ( &ld, &inters[0].nv );
    if ( g0*l0 > 0.0 ) {
      gray += 0.7*fabs ( g0 );
      if ( gray > 1.0 )
        gray = 1.0;
    }

    return (byte)(255.0*gray);
  }
  else
    return 160;

#undef MAXINTERS
#undef MAXLEVEL
} /*PixelColour*/

static void RayTrace ( void )
{
  int   x, y;
  ray3d ray;
  byte  buffer[WIDTH];

  CPos.upside = false;
  CameraSetMappingd ( &CPos );
  ps_Init_BitmapP ( WIDTH, HEIGHT, 0, 0 );
  for ( y = 0; y < HEIGHT; y++ ) {
    printf ( "%4d\b\b\b\b", y );
    for ( x = 0; x < WIDTH; x++ ) {
      CameraRayOfPixeld ( &CPos, x, y, &ray );
      buffer[x] = PixelColour ( &ray );
    }
    ps_Out_LineP ( buffer );
  }
} /*RayTrace*/

void DrawBezTreeVertexDivision ( BezPatchTreeVertexd *vertex,
                                 double w, double h, double x, double y )
{
  double w1, h1;

  if ( vertex->left ) {
    if ( vertex->ctlpoints )
      ps_Set_Line_Width ( 1.0 );
    else
      ps_Set_Line_Width ( 2.0 );
    if ( vertex->divdir == 0 ) {
      w1 = w*(vertex->left->u1-vertex->left->u0)/(vertex->u1-vertex->u0);
      ps_Draw_Line ( (float)(x+w1), (float)y, (float)(x+w1), (float)(y+h) );
      DrawBezTreeVertexDivision ( vertex->left, w1, h, x, y );
      DrawBezTreeVertexDivision ( vertex->right, w-w1, h, x+w1, y );
    }
    else {
      h1 = h*(vertex->left->v1-vertex->left->v0)/(vertex->v1-vertex->v0);
      ps_Draw_Line ( (float)x, (float)(y+h1), (float)(x+w), (float)(y+h1) );
      DrawBezTreeVertexDivision ( vertex->left, w, h1, x, y );
      DrawBezTreeVertexDivision ( vertex->right, w, h-h1, x, y+h1 );
    }
  }
} /*DrawBezTreeVertexDivision*/

void DrawBezTree ( BezPatchTreed *tree, double w, double x, double y )
{
  double h;

  h = w*(tree->root->v1-tree->root->v0)/(tree->root->u1-tree->root->u0);
  ps_Set_Gray ( 0.0 );
  ps_Set_Line_Width ( 2.0 );
  ps_Draw_Rect ( (float)w, (float)h, (float)x, (float)y );
  DrawBezTreeVertexDivision ( tree->root, w, h, x, y );
} /*DrawBezTree*/

int main ( void )
{
  struct tms start, stop;
  double time;

  setvbuf ( stdout, NULL, _IONBF, 0 );
  pkv_InitScratchMem ( 655360 );
  ps_WriteBBox ( -1, -1, 204, 88 );
  ps_OpenFile ( fn, 400 );
  SetCPos ();
  SetLightDir ();

  times ( &start );
  tree = rbez_NewBSPatchTreed ( 0, n, lknu, knots_u, m, lknv, knots_v,
                                lknv-m, cpoints );
  RayTrace ();
  DrawBezTree ( tree, 480, 650, 0 );

  rbez_DestroyBezPatchTreed ( tree );
  times ( &stop );
  time = (double)(stop.tms_utime-start.tms_utime)/(double)(sysconf(_SC_CLK_TCK));
  printf ( "time = %8.3f\n", time );
  ps_CloseFile ();

  pkv_DestroyScratchMem ();
  exit ( 0 );
} /*main*/

