
/* //////////////////////////////////////////////////// */
/* This file is a part of the BSTools procedure package */
/* written by Przemyslaw Kiciak.                        */
/* //////////////////////////////////////////////////// */

#include <stdarg.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "pkvaria.h"
#include "pknum.h"
#include "pkgeom.h"
#include "camera.h"
#include "multibs.h"
#include "psout.h"

char fn[] = "bscoonsd.ps";

#define n  3
#define N  8
#define m  3
#define M  9
double ukn[N+1] = {0.0,0.5,0.5,0.5,1.0,2.0,2.0,2.0,2.5};
double vkn[M+1] = {-0.5,0.0,0.0,0.2,1.0,2.0,3.0,3.0,3.0,3.5};
point3d cp[(N-n)*(M-m)] =
  {{0.0,0.0,0.2},{1.0,0.0,0.2},{2.0,0.0,0.0},{3.0,0.0,0.0},{4.0,0.0,0.0},{5.0,0.0,0.0},
   {0.0,1.0,0.2},{1.0,0.2,0.4},{2.0,1.0,0.2},{3.0,1.0,0.0},{4.0,1.0,0.0},{5.0,1.0,0.0},
   {0.0,2.0,0.2},{1.0,2.0,0.4},{2.0,2.0,0.2},{3.0,2.0,0.0},{4.0,2.0,0.0},{5.0,2.0,0.0},
   {0.0,3.0,0.0},{1.0,3.0,1.0},{2.0,3.0,0.0},{3.0,3.0,0.0},{4.8,3.8,1.0},{5.0,3.0,0.0},
   {0.0,4.0,0.0},{1.0,4.0,0.0},{2.0,4.0,0.0},{3.0,4.0,0.0},{4.0,4.0,0.0},{5.0,4.0,0.0}};

CameraRecd CPos;

point3d c00[N-n], c01[N-n], c02[N-n], c10[N-n], c11[N-n], c12[N-n];
point3d d00[M-m], d01[M-m], d02[M-m], d10[M-m], d11[M-m], d12[M-m];

int bscdegu, bscdegv, bsclknu, bsclknv;
double bscknu[40], bscknv[40];
point3d bsccp[200];


static void PrepareData ( void )
{
  mbs_multideBoorDer2d ( m, M, vkn, N-n, 3, 3*(M-m), (double*)cp, vkn[m],
                         (double*)c00, (double*)c01, (double*)c02 );
  mbs_multideBoorDer2d ( m, M, vkn, N-n, 3, 3*(M-m), (double*)cp, vkn[M-m],
                         (double*)c10, (double*)c11, (double*)c12 );
  mbs_multideBoorDer2d ( n, N, ukn, 1, 3*(M-m), 0, (double*)cp, ukn[n],
                         (double*)d00, (double*)d01, (double*)d02 );
  mbs_multideBoorDer2d ( n, N, ukn, 1, 3*(M-m), 0, (double*)cp, ukn[N-n],
                         (double*)d10, (double*)d11, (double*)d12 );
} /*PrepareData*/

static void SetupCamera ( void )
{
  point3d  p;
  vector3d v;

  CameraInitFramed ( &CPos, false, true, 1200, 900, 400, 0, 1.0, 4 );
  CameraInitPosd ( &CPos );

  SetPoint3d ( &p, 2.5, 2.0, 0.0 );
  CameraSetRotCentred ( &CPos, &p, true, true );
  SetVector3d ( &v, 2.5, 2.0, -25.0 );
  CameraMoveGd ( &CPos, &v );
  CameraRotXCd ( &CPos, -0.65*PI );
  CameraRotZGd ( &CPos, -0.14*PI );
  CameraZoomd ( &CPos, 5.0 );
} /*SetupCamera*/

static void DrawBSPatch ( int udeg, int lknu, double *knu,
                          int vdeg, int lknv, double *knv,
                          point3d *cp,
                          int du, int dv, int dd, double wd1, double wd2 )
{
  double   *cpa, *cpb;
  int     i, j, k, l, ku, kv, NNa, MMa;
  int     pitch;
  double   t;
  point3d p;
  point2d *c;
  void    *st;
 
  st = pkv_GetScratchMemTop ();

  ku = mbs_NumKnotIntervalsd ( udeg, lknu, knu );
  kv = mbs_NumKnotIntervalsd ( vdeg, lknv, knv );
  pitch = 3*kv*(vdeg+1);
  cpa = pkv_GetScratchMemd ( 3*ku*kv*(udeg+1)*(vdeg+1) );
  mbs_BSPatchToBezd ( 3, udeg, lknu, knu, vdeg, lknv, knv, 3*(lknv-vdeg),
                      (double*)cp,
                      &ku, &NNa, NULL, &kv, &MMa, NULL, pitch, cpa );

  cpb = pkv_GetScratchMemd ( pitch*sizeof(double) );
            c = (point2d*)pkv_GetScratchMem ( (max(du,dv)*dd+1)*sizeof(point2d) );
        /* draw lines of constant parameter u */
  for ( k = 0; k < ku; k++ ) {
    ps_Set_Line_Width ( (float)wd1 );
    for ( i = 0; i < du; i++ ) {
      t = (double)i/(double)du;
      mbs_multiBCHornerd ( udeg, 1, pitch, 0, &cpa[k*(udeg+1)*pitch], t, cpb );
      for ( l = 0; l < kv; l++ ) {
        for ( j = 0; j <= dv*dd; j++ ) {
          t = (double)j/(double)(dv*dd);
          mbs_BCHornerC3d ( vdeg, (point3d*)&cpb[l*(vdeg+1)*3], t, &p );
          CameraProjectPoint3d ( &CPos, &p, (point3d*)&c[j] );
        }
        ps_Draw_Polyline2d ( dv*dd+1, c );
      }
      ps_Set_Line_Width ( (float)wd2 );
    }
  }
  ps_Set_Line_Width ( (float)wd1 );
  mbs_multiBCHornerd ( udeg, 1, pitch, 0, &cpa[(ku-1)*(udeg+1)*pitch], 1.0, cpb );
  for ( l = 0; l < kv; l++ ) {
    for ( j = 0; j <= dv*dd; j++ ) {
      t = (double)j/(double)(dv*dd);
      mbs_BCHornerC3d ( vdeg, (point3d*)&cpb[l*(vdeg+1)*3], t, &p );
      CameraProjectPoint3d ( &CPos, &p, (point3d*)&c[j] );
    }
    ps_Draw_Polyline2d ( dv*dd+1, c );
  }
                             /* draw lines of constant parameter v */
  for ( l = 0; l < kv; l++ ) {
    ps_Set_Line_Width ( (float)wd1 );
    for ( j = 0; j < dv; j++ ) {
      t = (double)j/(double)dv;
      mbs_multiBCHornerd ( vdeg, NNa-udeg, 3, pitch, &cpa[l*(vdeg+1)*3], t, cpb );
      for ( k = 0; k < ku; k++ ) {
        for ( i = 0; i <= du*dd; i++ ) {
          t = (double)i/(double)(du*dd);
          mbs_BCHornerC3d ( udeg, (point3d*)&cpb[k*(udeg+1)*3], t, &p );
          CameraProjectPoint3d ( &CPos, &p, (point3d*)&c[i] );
        }
        ps_Draw_Polyline2d ( du*dd+1, c );
      }
      ps_Set_Line_Width ( (float)wd2 );
    }
  }
  ps_Set_Line_Width ( (float)wd1 );
  mbs_multiBCHornerd ( vdeg, NNa-udeg, 3, pitch, &cpa[(kv-1)*(vdeg+1)*3], 1.0, cpb );
  for ( k = 0; k < ku; k++ ) {
    for ( i = 0; i <= du*dd; i++ ) {
      t = (double)i/(double)(du*dd);  
      mbs_BCHornerC3d ( udeg, (point3d*)&cpb[k*(udeg+1)*3], t, &p );
      CameraProjectPoint3d ( &CPos, &p, (point3d*)&c[i] );
    }
    ps_Draw_Polyline2d ( du*dd+1, c );
  }
        /* deallocate buffers */
  pkv_SetScratchMemTop ( st );
} /*DrawBSPatch*/

static void DrawPicture ( void )
{
  ps_OpenFile ( fn, 600 );
  ps_Write_Command ( "1 setlinecap" );
  SetupCamera ();

  ps_GSave ();
  ps_Write_Command ( "0 900 translate" );
  DrawBSPatch ( n, N, ukn, m, M, vkn, cp, 6, 6, 4, 6.0, 2.0 );

  mbs_BSC1CoonsToBSd ( 3,
      n, N, ukn, (double*)c00, n, N, ukn, (double*)c01,
      n, N, ukn, (double*)c10, n, N, ukn, (double*)c11,
      m, M, vkn, (double*)d00, m, M, vkn, (double*)d01,
      m, M, vkn, (double*)d10, m, M, vkn, (double*)d11,
      &bscdegu, &bsclknu, bscknu, &bscdegv, &bsclknv, bscknv, (double*)bsccp );

  ps_Set_RGB ( 0.0, 1.0, 1.0 );
  DrawBSPatch ( bscdegu, bsclknu, bscknu, bscdegv, bsclknv, bscknv, bsccp,
                6, 6, 4, 6.0, 2.0 );

  ps_GRestore ();

  ps_GSave ();
  DrawBSPatch ( n, N, ukn, m, M, vkn, cp, 6, 6, 4, 6.0, 2.0 );

  mbs_BSC2CoonsToBSd ( 3,
      n, N, ukn, (double*)c00, n, N, ukn, (double*)c01, n, N, ukn, (double*)c02,
      n, N, ukn, (double*)c10, n, N, ukn, (double*)c11, n, N, ukn, (double*)c12,
      m, M, vkn, (double*)d00, m, M, vkn, (double*)d01, m, M, vkn, (double*)d02,
      m, M, vkn, (double*)d10, m, M, vkn, (double*)d11, m, M, vkn, (double*)d12,
      &bscdegu, &bsclknu, bscknu, &bscdegv, &bsclknv, bscknv, (double*)bsccp );

  ps_Set_RGB ( 0.0, 1.0, 1.0 );
  DrawBSPatch ( bscdegu, bsclknu, bscknu, bscdegv, bsclknv, bscknv, bsccp,
                6, 6, 4, 6.0, 2.0 );
  ps_GRestore ();

  ps_CloseFile ();
  printf ( "%s\n", fn );
} /*DrawPicture*/

int main ( void )
{
  pkv_InitScratchMem ( 262144 );
  PrepareData ();
  DrawPicture ();
  pkv_DestroyScratchMem ();
  exit ( 0 );
} /*main*/

