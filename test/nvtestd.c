
/* //////////////////////////////////////////////////// */
/* This file is a part of the BSTools procedure package */
/* written by Przemyslaw Kiciak.                        */
/* //////////////////////////////////////////////////// */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "pkgeom.h"
#include "camera.h"
#include "multibs.h"
#include "psout.h"

char fn1[] = "nvtest1d.ps";
char fn2[] = "nvtest2d.ps";

point4d p4[36] =
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

point3d p3[36];

CameraRecd CPos;

static void SetupData ( void )
{
  int i;

  for ( i = 0; i < 36; i++ )
    Point4to3d ( &p4[i], &p3[i] );
} /*SetupData*/

static void SetCPos ( void )
{
  vector3d v;

  CameraInitFramed ( &CPos, false, true, 1200, 900, 0, 0, 1.0, 4 );
  CameraInitPosd ( &CPos );
  SetVector3d ( &v, -0.6, -0.5, -20.0 );
  CameraMoveCd ( &CPos, &v );
  CameraRotZCd ( &CPos, PI );
  CameraRotYCd ( &CPos, -PI/9.0 );
  CameraRotXCd ( &CPos, -PI/9.0 );
  CameraRotZCd ( &CPos, 0.2*PI );
  CameraZoomd ( &CPos, 5.0 );
} /*SetCPos*/

static void DrawZero ( int n, int m, vector3d *p )
{
  point3d q, r, s;

  SetVector3d ( &q, 0.0, 0.0, 0.0 );
  CameraProjectPoint3d ( &CPos, &q, &r );

  mbs_BCHornerP3d ( n, m, p, 0.0, 0.0, &q );
  CameraProjectPoint3d ( &CPos, &q, &s );
  psl_SetLine ( (float)r.x, (float)r.y, (float)s.x, (float)s.y, 0.0, 1.0 );
  psl_Draw ( 0.0, 0.95, 3.0 );
  psl_Arrow ( 1.0, true );

  mbs_BCHornerP3d ( n, m, p, 1.0, 0.0, &q );
  CameraProjectPoint3d ( &CPos, &q, &s );
  psl_SetLine ( (float)r.x, (float)r.y, (float)s.x, (float)s.y, 0.0, 1.0 );
  psl_Draw ( 0.0, 0.95, 3.0 );
  psl_Arrow ( 1.0, true );

  mbs_BCHornerP3d ( n, m, p, 0.0, 1.0, &q );
  CameraProjectPoint3d ( &CPos, &q, &s );
  psl_SetLine ( (float)r.x, (float)r.y, (float)s.x, (float)s.y, 0.0, 1.0 );
  psl_Draw ( 0.0, 0.95, 3.0 );
  psl_Arrow ( 1.0, true );
  
  mbs_BCHornerP3d ( n, m, p, 1.0, 1.0, &q );
  CameraProjectPoint3d ( &CPos, &q, &s );
  psl_SetLine ( (float)r.x, (float)r.y, (float)s.x, (float)s.y, 0.0, 1.0 );
  psl_Draw ( 0.0, 0.95, 3.0 );
  psl_Arrow ( 1.0, true );

  ps_Fill_Circle ( (float)r.x, (float)r.y, 10.0 );
} /*DrawZero*/

static void DrawBezNet3 ( int n, int m, point3d *p )
{
  int     i, j, ncp;
  int    size_r;
  point3d q;
  point2d *r;

  ps_Set_Line_Width ( 2.0 );
  ncp = (n+1)*(m+1);
  r = (point2d*)pkv_GetScratchMem ( size_r = ncp*sizeof(point2d) );
  for ( i = 0; i < ncp; i++ ) {
    CameraProjectPoint3d ( &CPos, &p[i], &q );
    r[i].x = q.x;  r[i].y = q.y;
  }
  for ( i = 0; i < ncp; i += m+1 ) {
    ps_Draw_Polyline2d ( m+1, &r[i] );
    if ( i < ncp-m-1 )
      for ( j = 0; j <= m; j++ )
        ps_Draw_Line ( (float)r[i+j].x, (float)r[i+j].y, (float)r[i+j+m+1].x, (float)r[i+j+m+1].y );
  }
  for ( i = 0; i < ncp; i++ )
    ps_Mark_Circle ( (float)r[i].x, (float)r[i].y );
  pkv_FreeScratchMem ( size_r );
} /*DrawBezNet3*/

static void DrawBezNet4 ( int n, int m, point4d *p )
{
  int     i, j, ncp;
  int    size_r;
  point3d q;
  point2d *r;

  ps_Set_Line_Width ( 2.0 );
  ncp = (n+1)*(m+1);
  r = (point2d*)pkv_GetScratchMem ( size_r = ncp*sizeof(point2d) );
  for ( i = 0; i < ncp; i++ ) {
    Point4to3d ( &p[i], &q );
    CameraProjectPoint3d ( &CPos, &q, &q );
    r[i].x = q.x;  r[i].y = q.y;
  }
  for ( i = 0; i < ncp; i += m+1 ) {
    ps_Draw_Polyline2d ( m+1, &r[i] );
    if ( i < ncp-m-1 )
      for ( j = 0; j <= m; j++ )
        ps_Draw_Line ( (float)r[i+j].x, (float)r[i+j].y, (float)r[i+j+m+1].x, (float)r[i+j+m+1].y );
  }
  for ( i = 0; i < ncp; i++ )
    ps_Mark_Circle ( (float)r[i].x, (float)r[i].y );
  pkv_FreeScratchMem ( size_r );
} /*DrawBezNet4*/

int main ()
{
#define nn 4

  int      ndegu, ndegv, size_ncp;
  vector3d *ncp;

  pkv_InitScratchMem ( 262144 );
  SetupData ();
  SetCPos ();

  ps_OpenFile ( fn1, 600 );
  mbs_BezP3RNormalDeg ( nn, 5, &ndegu, &ndegv );
  printf ( "%d, %d\n", ndegu, ndegv );
  ncp = (vector3d*)pkv_GetScratchMem (
                      size_ncp = (ndegu+1)*(ndegv+1)*sizeof(vector3d) );
  mbs_BezP3RNormald ( nn, 5, p4, &ndegu, &ndegv, ncp );
  DrawBezNet4 ( nn, 5, p4 );
  ps_Write_Command ( "900 -100 translate" );
  DrawZero ( ndegu, ndegv, ncp );
  DrawBezNet3 ( ndegu, ndegv, ncp );

  pkv_FreeScratchMem ( size_ncp );
  ps_CloseFile ();

  ps_OpenFile ( fn2, 600 );
  mbs_BezP3NormalDeg ( nn, 5, &ndegu, &ndegv );
  printf ( "%d, %d\n", ndegu, ndegv );
  ncp = (vector3d*)pkv_GetScratchMem (
                      size_ncp = (ndegu+1)*(ndegv+1)*sizeof(vector3d) );
  mbs_BezP3Normald ( nn, 5, p3, &ndegu, &ndegv, ncp );
  DrawBezNet3 ( nn, 5, p3 );
  ps_Write_Command ( "900 700 translate" );
  DrawBezNet3 ( ndegu, ndegv, ncp );
  DrawZero ( ndegu, ndegv, ncp );

  pkv_FreeScratchMem ( size_ncp );
  ps_CloseFile ();

  pkv_DestroyScratchMem ();
  exit ( 0 );

#undef nn
} /*main*/

