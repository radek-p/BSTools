                                         
/* //////////////////////////////////////////////////// */
/* This file is a part of the BSTools procedure package */
/* written by Przemyslaw Kiciak.                        */
/* //////////////////////////////////////////////////// */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <sys/times.h>
#include <unistd.h>
#include <pthread.h>

#include "pkgeom.h"
#include "camera.h"
#include "multibs.h"
#include "raybez.h"
#include "psout.h"

char fn[] = "raybezp.ps";

#define n 5
#define m 5
point4f p[36] =
  {{1.23205    , 1.86602    , 0.00000    , 1.00000    },
   {1.08766    , 1.64734    ,-2.82807E-01, 8.82807E-01},
   {1.06547    , 1.45139    ,-5.24210E-01, 8.24210E-01},
   {1.16547    , 1.27819    ,-7.24210E-01, 8.24210E-01},
   {1.38766    , 1.12772    ,-8.82807E-01, 8.82807E-01},
   {1.73205    , 1.00000    ,-1.00000    , 1.00000    },
   {9.85640E-01, 1.72375    ,-2.45223E-15, 1.00000    },
   {8.69561E-01, 1.52082    ,-2.81499E-01, 8.82292E-01},
   {8.52375E-01, 1.32836    ,-5.24210E-01, 8.24210E-01},
   {9.32375E-01, 1.14361    ,-7.24210E-01, 8.24210E-01},
   {1.11013    , 9.67490E-01,-8.82807E-01, 8.82807E-01},
   {1.38564    , 8.00000E-01,-1.00000    , 1.00000    },
   {7.39230E-01, 1.57964    , 0.00000    , 1.00000    },
   {6.53471E-01, 1.39526    ,-2.90283E-01, 8.83009E-01},
   {6.39771E-01, 1.20972    ,-5.41115E-01, 8.22961E-01},
   {7.00206E-01, 1.01277    ,-7.30543E-01, 8.22911E-01},
   {8.34683E-01, 8.08461E-01,-8.83231E-01, 8.83231E-01},
   {1.03923    , 6.00000E-01,-1.00000    , 1.00000    },
   {4.92820E-01, 1.45862    , 0.00000    , 1.00000    },
   {4.36421E-01, 1.28919    ,-3.21810E-01, 8.83443E-01},
   {4.26850E-01, 1.10411    ,-5.79775E-01, 8.21558E-01},
   {4.67694E-01, 8.87191E-01,-7.47754E-01, 8.21536E-01},
   {5.58746E-01, 6.49149E-01,-8.83654E-01, 8.83654E-01},
   {6.92820E-01, 4.00000E-01,-1.00000    , 1.00000    },
   {2.46410E-01, 1.38564    , 0.00000    , 1.00000    },
   {2.17929E-01, 1.22225    ,-3.83472E-01, 8.82080E-01},
   {2.11223E-01, 1.02715    ,-6.54292E-01, 8.16225E-01},
   {2.31748E-01, 7.73389E-01,-7.80138E-01, 8.16885E-01},
   {2.78980E-01, 4.91655E-01,-8.83654E-01, 8.83654E-01},
   {3.46410E-01, 2.00000E-01,-1.00000    , 1.00000    },
   {0.00000    , 1.38564    , 0.00000    , 1.00000    },
   {0.00000    , 1.21577    ,-4.77404E-01, 8.77404E-01},
   {0.00000    , 9.63353E-01,-7.18808E-01, 8.18808E-01},
   {0.00000    , 6.57312E-01,-8.21509E-01, 8.21509E-01},
   {0.00000    , 3.26557E-01,-8.82807E-01, 8.82807E-01},
   {0.00000    , 0.00000    ,-1.00000    , 1.00000    }};

#define WIDTH  480
#define HEIGHT 360

CameraRecf CPos;
RBezPatchTreefp tree;
vector3f ld;

static void SetCPos ( void )
{
  vector3f v;

  CameraInitFramef ( &CPos, false, true, WIDTH, HEIGHT, 0, 0, 1.0, 4 );
  CameraInitPosf ( &CPos );
  SetPoint3f ( &CPos.position, 10.5, 16.5, -8.1 );
  CPos.psi   = (float)(4.5*PI/180.0);
  CPos.theta = (float)(67.83*PI/180.0);
  CPos.phi   = (float)(-34.0*PI/180.0);
  CPos.vd.persp.f = 7.5;
  CameraSetMappingf ( &CPos );
  CameraRotZCf ( &CPos, -3.0*PI/4.0 );
  CameraRotXCf ( &CPos, PI );
  CameraRotZCf ( &CPos, 0.5*PI );
  CameraRotYCf ( &CPos, 0.5*PI );
  SetVector3f ( &v, 2.0, 0.0, 0.0 );
  CameraMoveCf ( &CPos, &v );
} /*SetCPos*/

static void SetLightDir ( void )
{
  SetVector3f ( &ld, -0.5, +0.75, +0.5 );
  NormalizeVector3f ( &ld );
} /*SetLightDir*/

static void DrawTreeVertexDivision ( RBezPatchTreeVertexf *vertex,
                                     float w, float h, float x, float y )
{
  if ( vertex->left ) {
    if ( vertex->divdir ) {
      w *= 0.5;
      ps_Draw_Line ( x+w, y, x+w, y+h );
      DrawTreeVertexDivision ( vertex->left, w, h, x, y );
      DrawTreeVertexDivision ( vertex->right, w, h, x+w, y );
    }
    else {
      h *= 0.5;
      ps_Draw_Line ( x, y+h, x+w, y+h );
      DrawTreeVertexDivision ( vertex->left, w, h, x, y );
      DrawTreeVertexDivision ( vertex->right, w, h, x, y+h );
    }
  }
} /*DrawTreeVertexDivision*/

static void DrawRBezTree ( RBezPatchTreef *tree, float w, float h,
                           float x, float y )
{
  ps_Set_Gray ( 0.0 );
  ps_Set_Line_Width ( 2.0 );
  ps_Draw_Rect ( w, h, x, y );
  ps_Set_Line_Width ( 1.0 );
  DrawTreeVertexDivision ( tree->root, w, h, x, y );
} /*DrawRBezTree*/

static byte PixelColour ( ray3f *ray )
{
#define MAXLEVEL  16
#define MAXINTERS 10
  RayObjectIntersf inters[MAXINTERS];
  int ninters, i, j;
  float tmin;
  float gray, g0;

  rbez_FindRayRBezPatchIntersf ( tree, ray, MAXLEVEL, MAXINTERS,
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
    if ( ( g0 = (float)DotProduct3f ( &ld, &inters[0].nv ) ) > 0.0 )
      gray += (float)(0.7*g0);
    if ( gray > 1.0 )
      gray = 1.0;

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
  ray3f ray;
  byte  buffer[WIDTH];

  CPos.upside = false;
  CameraSetMappingf ( &CPos );
  ps_Init_BitmapP ( WIDTH, HEIGHT, 0, 0 );
  for ( y = 0; y < HEIGHT; y++ ) {
    printf ( "%4d\b\b\b\b", y );
    for ( x = 0; x < WIDTH; x++ ) {
      CameraRayOfPixelf ( &CPos, (float)x, (float)y, &ray );
      buffer[x] = PixelColour ( &ray );
    }
    ps_Out_LineP ( buffer );
  }
} /*RayTrace*/

int main ( void )
{
  struct tms start, stop;
  double time;

  setvbuf ( stdout, NULL, _IONBF, 0 );
  pkv_InitScratchMem ( 262144 );
  ps_WriteBBox ( 0, 0, 203, 88 );
  ps_OpenFile ( fn, 300 );
  SetCPos ();
  SetLightDir ();

  times ( &start );
  tree = rbez_NewRBezPatchTreef ( 0, n, m, 0.0, 1.0, 0.0, 1.0, p );

  RayTrace ();
  DrawRBezTree ( tree, 360.0, 360.0, 490.0, 0.0 );

  rbez_DestroyRBezPatchTreef ( tree );
  times ( &stop );
  time = (float)(stop.tms_utime-start.tms_utime)/(float)(sysconf(_SC_CLK_TCK));
  printf ( "time = %8.3f\n", time );

  ps_CloseFile ();
  pkv_DestroyScratchMem ();
  exit ( 0 );
} /*main*/

