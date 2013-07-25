
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

char fn[] = "raytestd.ps";

#define n 5
#define m 5
point4d p[36] =
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

#define WIDTH  320
#define HEIGHT 240

CameraRecd CPos;
RBezPatchTreedp tree;
vector3d ld;

static void SetCPos ( void )
{
  vector3d v;

  CameraInitFramed ( &CPos, false, false, WIDTH, HEIGHT, 0, 0, 1.0, 4 );
  CameraInitPosd ( &CPos );
  SetPoint3d ( &CPos.position, 10.5, 16.5, -8.1 );
  CPos.psi   = 4.5*PI/180.0;
  CPos.theta = 67.83*PI/180.0;
  CPos.phi   = -34.0*PI/180.0;
  CPos.vd.persp.f = 7.5;
  CameraSetMappingd ( &CPos );
  CameraRotZCd ( &CPos, -3*PI/4 );
  CameraRotXCd ( &CPos, PI );
  SetVector3d ( &v, 0.4, 0.0, 0.0 );
  CameraMoveCd ( &CPos, &v );
  CameraRotZCd ( &CPos, 0.5*PI );
} /*SetCPos*/

static void SetLightDir ( void )
{
  SetVector3d ( &ld, 0.5, 0.5, -1.0 );
  NormalizeVector3d ( &ld );
} /*SetLightDir*/

/*
static void DrawBezNet ( point4d *p )
{
  int    i;
  point3d q;
  point2d r[36];

  CPos.upside = true;
  CameraSetMappingd ( &CPos );
  for ( i = 0; i < 36; i++ ) {
    Point4to3d ( &p[i], &q );
    CameraProjectPoint3d ( &CPos, &q, &q );
    r[i].x = q.x;  r[i].y = q.y;
  }
  ps_Set_Line_Width ( 0.75 );
  for ( i = 0; i < 30; i++ )
    ps_Draw_Line ( r[i].x, r[i].y, r[i+6].x, r[i+6].y );
  for ( i = 0; i <= 5; i++ )
    ps_Draw_Polyline2d ( 6, &r[6*i] );
  for ( i = 0; i < 36; i++ )
    ps_Fill_Circle ( r[i].x, r[i].y, 2.0 );
}*/ /*DrawBezNet*/

static void DrawTreeVertexDivision ( RBezPatchTreeVertexd *vertex,
                                     double w, double h, double x, double y )
{
  if ( vertex->left ) {
    if ( vertex->divdir ) {
      w *= 0.5;
      ps_Draw_Line ( (float)(x+w), (float)y, (float)(x+w), (float)(y+h) );
      DrawTreeVertexDivision ( vertex->left, w, h, x, y );
      DrawTreeVertexDivision ( vertex->right, w, h, x+w, y );
    }
    else {
      h *= 0.5;
      ps_Draw_Line ( (float)x, (float)(y+h), (float)(x+w), (float)(y+h) );
      DrawTreeVertexDivision ( vertex->left, w, h, x, y );
      DrawTreeVertexDivision ( vertex->right, w, h, x, y+h );
    }
  }
} /*DrawTreeVertexDivision*/

static void DrawRBezTree ( RBezPatchTreed *tree, double w, double h, double x, double y )
{
  ps_Set_Gray ( 0.0 );
  ps_Set_Line_Width ( 1.0 );
  ps_Draw_Rect ( (float)w, (float)h, (float)x, (float)y );
  ps_Set_Line_Width ( 0.5 );
  DrawTreeVertexDivision ( tree->root, w, h, x, y );
} /*DrawRBezTree*/

static byte PixelColour ( ray3d *ray )
{
#define MAXLEVEL  16
#define MAXINTERS 10
  RayObjectIntersd inters[MAXINTERS];
  int ninters, i, j;
  double tmin;
  double gray, g0;

  rbez_FindRayRBezPatchIntersd ( tree, ray, MAXLEVEL, MAXINTERS,
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
    if ( ( g0 = DotProduct3d ( &ld, &inters[0].nv ) ) > 0.0 )
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
/*
static void TestTree ()
{
  RBezPatchTreeVertexdp a, b, c, d, e, f, g, h, i;

  a = tree->root;
  b = rbez_GetRBezLeftVertexd ( tree, a );
  c = rbez_GetRBezRightVertexd ( tree, a );

  d = rbez_GetRBezLeftVertexd ( tree, b );
  e = rbez_GetRBezRightVertexd ( tree, b );
  f = rbez_GetRBezLeftVertexd ( tree, c );
  g = rbez_GetRBezRightVertexd ( tree, c );

  h = rbez_GetRBezLeftVertexd ( tree, d );
  i = rbez_GetRBezRightVertexd ( tree, d );

  ps_Write_Command ( "0 0.5 0.5 setrgbcolor" );
  DrawBezNet ( h->ctlpoints );
  ps_Write_Command ( "1 0 0 setrgbcolor" );
  DrawBezNet ( i->ctlpoints );
  ps_Write_Command ( "0 0 1 setrgbcolor" );
  DrawBezNet ( e->ctlpoints );
  ps_Write_Command ( "0 0.7 0 setrgbcolor" );
  DrawBezNet ( f->ctlpoints );
  ps_Write_Command ( "0.5 0.5 0 setrgbcolor" );
  DrawBezNet ( g->ctlpoints );

  ps_Set_Gray ( 0.001 );
  DrawRBezTree ( tree, 240.0, 240.0, 322.0, 2.0 );
}*/ /*TestTree*/

int main ( void )
{
  struct tms start, stop;
  double time;

  setvbuf ( stdout, NULL, _IONBF, 0 );
  pkv_InitScratchMem ( 262144 );
  ps_OpenFile ( fn, 200 );
  SetCPos ();
  SetLightDir ();

  times ( &start );
  tree = rbez_NewRBezPatchTreed ( 0, n, m, 0.0, 1.0, 0.0, 1.0, p );

  RayTrace ();
  DrawRBezTree ( tree, 240.0, 240.0, 322.0, 2.0 );

  rbez_DestroyRBezPatchTreed ( tree );
  times ( &stop );
  time = (double)(stop.tms_utime-start.tms_utime)/(double)(sysconf(_SC_CLK_TCK));
  printf ( "time = %8.3f\n", time );

  ps_CloseFile ();
  pkv_DestroyScratchMem ();
  exit ( 0 );
} /*main*/

