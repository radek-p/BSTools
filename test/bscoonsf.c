
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

char fn[] = "bscoonsf.ps";

#define n  3
#define N  8
#define m  3
#define M  9
float ukn[N+1] = {0.0,0.5,0.5,0.5,1.0,2.0,2.0,2.0,2.5};
float vkn[M+1] = {-0.5,0.0,0.0,0.2,1.0,2.0,3.0,3.0,3.0,3.5};
point3f cp[(N-n)*(M-m)] =
  {{0.0,0.0,0.2},{1.0,0.0,0.2},{2.0,0.0,0.0},{3.0,0.0,0.0},{4.0,0.0,0.0},{5.0,0.0,0.0},
   {0.0,1.0,0.2},{1.0,0.2,0.4},{2.0,1.0,0.2},{3.0,1.0,0.0},{4.0,1.0,0.0},{5.0,1.0,0.0},
   {0.0,2.0,0.2},{1.0,2.0,0.4},{2.0,2.0,0.2},{3.0,2.0,0.0},{4.0,2.0,0.0},{5.0,2.0,0.0},
   {0.0,3.0,0.0},{1.0,3.0,1.0},{2.0,3.0,0.0},{3.0,3.0,0.0},{4.8,3.8,1.0},{5.0,3.0,0.0},
   {0.0,4.0,0.0},{1.0,4.0,0.0},{2.0,4.0,0.0},{3.0,4.0,0.0},{4.0,4.0,0.0},{5.0,4.0,0.0}};

CameraRecf CPos;

point3f c00[N-n], c01[N-n], c02[N-n], c10[N-n], c11[N-n], c12[N-n];
point3f d00[M-m], d01[M-m], d02[M-m], d10[M-m], d11[M-m], d12[M-m];

int bscdegu, bscdegv, bsclknu, bsclknv;
float bscknu[40], bscknv[40];
point3f bsccp[200];


static void PrepareData ( void )
{
  mbs_multideBoorDer2f ( m, M, vkn, N-n, 3, 3*(M-m), (float*)cp, vkn[m],
                         (float*)c00, (float*)c01, (float*)c02 );
  mbs_multideBoorDer2f ( m, M, vkn, N-n, 3, 3*(M-m), (float*)cp, vkn[M-m],
                         (float*)c10, (float*)c11, (float*)c12 );
  mbs_multideBoorDer2f ( n, N, ukn, 1, 3*(M-m), 0, (float*)cp, ukn[n],
                         (float*)d00, (float*)d01, (float*)d02 );
  mbs_multideBoorDer2f ( n, N, ukn, 1, 3*(M-m), 0, (float*)cp, ukn[N-n],
                         (float*)d10, (float*)d11, (float*)d12 );
} /*PrepareData*/

static void SetupCamera ( void )
{
  point3f  p;
  vector3f v;

  CameraInitFramef ( &CPos, false, true, 1200, 900, 400, 0, 1.0, 4 );
  CameraInitPosf ( &CPos );

  SetPoint3f ( &p, 2.5, 2.0, 0.0 );
  CameraSetRotCentref ( &CPos, &p, true, true );
  SetVector3f ( &v, 2.5, 2.0, -25.0 );
  CameraMoveGf ( &CPos, &v );
  CameraRotXCf ( &CPos, -0.65*PI );
  CameraRotZGf ( &CPos, -0.14*PI );
  CameraZoomf ( &CPos, 5.0 );
} /*SetupCamera*/

static void DrawBSPatch ( int udeg, int lknu, float *knu,
                          int vdeg, int lknv, float *knv,
                          point3f *cp,
                          int du, int dv, int dd, float wd1, float wd2 )
{
  float   *cpa, *cpb;
  int     i, j, k, l, ku, kv, NNa, MMa;
  int     pitch;
  float   t;
  point3f p;
  point2f *c;
  void    *st;
 
  st = pkv_GetScratchMemTop ();

  ku = mbs_NumKnotIntervalsf ( udeg, lknu, knu );
  kv = mbs_NumKnotIntervalsf ( vdeg, lknv, knv );
  pitch = 3*kv*(vdeg+1);
  cpa = pkv_GetScratchMemf ( 3*ku*kv*(udeg+1)*(vdeg+1) );
  mbs_BSPatchToBezf ( 3, udeg, lknu, knu, vdeg, lknv, knv, 3*(lknv-vdeg),
                      (float*)cp,
                      &ku, &NNa, NULL, &kv, &MMa, NULL, pitch, cpa );

  cpb = pkv_GetScratchMemf ( pitch*sizeof(float) );
            c = (point2f*)pkv_GetScratchMem ( (max(du,dv)*dd+1)*sizeof(point2f) );
        /* draw lines of constant parameter u */
  for ( k = 0; k < ku; k++ ) {
    ps_Set_Line_Width ( wd1 );
    for ( i = 0; i < du; i++ ) {
      t = (float)i/(float)du;
      mbs_multiBCHornerf ( udeg, 1, pitch, 0, &cpa[k*(udeg+1)*pitch], t, cpb );
      for ( l = 0; l < kv; l++ ) {
        for ( j = 0; j <= dv*dd; j++ ) {
          t = (float)j/(float)(dv*dd);
          mbs_BCHornerC3f ( vdeg, (point3f*)&cpb[l*(vdeg+1)*3], t, &p );
          CameraProjectPoint3f ( &CPos, &p, (point3f*)&c[j] );
        }
        ps_Draw_Polyline2f ( dv*dd+1, c );
      }
      ps_Set_Line_Width ( wd2 );
    }
  }
  ps_Set_Line_Width ( wd1 );
  mbs_multiBCHornerf ( udeg, 1, pitch, 0, &cpa[(ku-1)*(udeg+1)*pitch], 1.0, cpb );
  for ( l = 0; l < kv; l++ ) {
    for ( j = 0; j <= dv*dd; j++ ) {
      t = (float)j/(float)(dv*dd);
      mbs_BCHornerC3f ( vdeg, (point3f*)&cpb[l*(vdeg+1)*3], t, &p );
      CameraProjectPoint3f ( &CPos, &p, (point3f*)&c[j] );
    }
    ps_Draw_Polyline2f ( dv*dd+1, c );
  }
                             /* draw lines of constant parameter v */
  for ( l = 0; l < kv; l++ ) {
    ps_Set_Line_Width ( wd1 );
    for ( j = 0; j < dv; j++ ) {
      t = (float)j/(float)dv;
      mbs_multiBCHornerf ( vdeg, NNa-udeg, 3, pitch, &cpa[l*(vdeg+1)*3], t, cpb );
      for ( k = 0; k < ku; k++ ) {
        for ( i = 0; i <= du*dd; i++ ) {
          t = (float)i/(float)(du*dd);
          mbs_BCHornerC3f ( udeg, (point3f*)&cpb[k*(udeg+1)*3], t, &p );
          CameraProjectPoint3f ( &CPos, &p, (point3f*)&c[i] );
        }
        ps_Draw_Polyline2f ( du*dd+1, c );
      }
      ps_Set_Line_Width ( wd2 );
    }
  }
  ps_Set_Line_Width ( wd1 );
  mbs_multiBCHornerf ( vdeg, NNa-udeg, 3, pitch, &cpa[(kv-1)*(vdeg+1)*3], 1.0, cpb );
  for ( k = 0; k < ku; k++ ) {
    for ( i = 0; i <= du*dd; i++ ) {
      t = (float)i/(float)(du*dd);  
      mbs_BCHornerC3f ( udeg, (point3f*)&cpb[k*(udeg+1)*3], t, &p );
      CameraProjectPoint3f ( &CPos, &p, (point3f*)&c[i] );
    }
    ps_Draw_Polyline2f ( du*dd+1, c );
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

  mbs_BSC1CoonsToBSf ( 3,
      n, N, ukn, (float*)c00, n, N, ukn, (float*)c01,
      n, N, ukn, (float*)c10, n, N, ukn, (float*)c11,
      m, M, vkn, (float*)d00, m, M, vkn, (float*)d01,
      m, M, vkn, (float*)d10, m, M, vkn, (float*)d11,
      &bscdegu, &bsclknu, bscknu, &bscdegv, &bsclknv, bscknv, (float*)bsccp );

  ps_Set_RGB ( 0.0, 1.0, 1.0 );
  DrawBSPatch ( bscdegu, bsclknu, bscknu, bscdegv, bsclknv, bscknv, bsccp,
                6, 6, 4, 6.0, 2.0 );

  ps_GRestore ();

  ps_GSave ();
  DrawBSPatch ( n, N, ukn, m, M, vkn, cp, 6, 6, 4, 6.0, 2.0 );

  mbs_BSC2CoonsToBSf ( 3,
      n, N, ukn, (float*)c00, n, N, ukn, (float*)c01, n, N, ukn, (float*)c02,
      n, N, ukn, (float*)c10, n, N, ukn, (float*)c11, n, N, ukn, (float*)c12,
      m, M, vkn, (float*)d00, m, M, vkn, (float*)d01, m, M, vkn, (float*)d02,
      m, M, vkn, (float*)d10, m, M, vkn, (float*)d11, m, M, vkn, (float*)d12,
      &bscdegu, &bsclknu, bscknu, &bscdegv, &bsclknv, bscknv, (float*)bsccp );

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

