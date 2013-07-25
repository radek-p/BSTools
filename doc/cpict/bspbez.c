
/* //////////////////////////////////////////////////// */
/* This file is a part of the BSTools procedure package */
/* written by Przemyslaw Kiciak.                        */
/* //////////////////////////////////////////////////// */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#include "pkvaria.h"
#include "pknum.h"
#include "pkgeom.h"
#include "camera.h"
#include "multibs.h"
#include "psout.h"

#define fn "bspbez.ps"
#define SCRATCHMEMSIZE 262144

#define  n1   3
#define  NN1 10
int n = n1;
int NN = NN1;
float u[NN1+1] = {-0.5, -0.4, 0.0, 0.0, 1.0, 2.0, 3.0, 4.0, 4.0, 4.4, 4.5};
#define  m1   3
#define  MM1 10
int m = m1;
int MM = MM1;
float v[MM1+1] = {-0.5, -0.4, 0.0, 0.0, 1.0, 2.5, 2.5, 4.0, 4.2, 4.3, 4.5};
point3f cp[NN1-n1][MM1-m1] =
  {{{-3.0,-3.0,-0.5},{-2.2,-3.0,-0.5},{-1.3,-3.0,-0.5},{0.0,-3.0,-0.5},{1.5,-3.0,-0.5},{2.3,-3.0,-0.5},{3.0,-3.0,-0.5}},
   {{-3.0,-2.7, 0.0},{-2.2,-2.7, 0.0},{-1.3,-2.7, 0.0},{0.0,-2.7, 0.0},{1.5,-2.7, 0.0},{2.3,-2.7, 0.0},{3.0,-2.7, 0.0}},
   {{-3.0,-1.5, 0.1},{-2.2,-1.5, 0.1},{-1.3,-1.5, 0.1},{0.0,-1.5, 0.1},{1.5,-1.5, 0.1},{2.3,-1.5, 0.1},{3.0,-1.5, 0.1}},
   {{-3.0, 0.0, 0.2},{-2.2, 0.0, 0.2},{-1.3, 0.0, 0.2},{0.0, 0.0, 0.2},{1.5, 0.0, 0.2},{2.3, 0.0, 0.2},{3.0, 0.0, 0.2}},
   {{-3.0,+1.5, 0.1},{-2.2,+1.5, 0.1},{-1.3,+1.5, 0.1},{0.0,+1.5, 0.1},{1.5,+1.5, 0.1},{2.3,+1.5, 0.1},{3.0,+1.5, 0.1}},
   {{-3.0,+2.7, 0.0},{-2.2,+2.7, 0.0},{-1.3,+2.7, 0.0},{0.0,+2.7, 0.0},{1.5,+2.7, 0.0},{2.3,+2.7, 0.0},{3.0,+2.7, 0.0}},
   {{-3.0,+3.0,-0.5},{-2.2,+3.0,-0.5},{-1.3,+3.0,-0.5},{0.0,+3.0,-0.5},{1.5,+3.0,-0.5},{2.3,+3.0,-0.5},{3.0,+3.0,-0.5}}};

point2f frame[3] = {{300,-300},{400,-300},{300,-200}};

CameraRecf CPos;


static void SetupCamera ( void )
{
  vector3f v;

  CameraInitFramef ( &CPos, false, true, 1200, 900, 400, 0, 1.0, 4 );
  CameraInitPosf ( &CPos );

  SetVector3f ( &v, 0.0, 0.0, -50.0 );
  CameraMoveGf ( &CPos, &v );
  CameraRotXCf ( &CPos, -0.65*PI );
  CameraRotZGf ( &CPos, -0.14*PI );
  CameraZoomf ( &CPos, 5.0 );
} /*SetupCamera*/

static void DisplayBSCNet ( int n, int NN, int m, int MM,
                            const float *cp, float wd )
{
  int         i, j;
  point3f     p, q;
  point2f     *buf; 
  const float *acp;
  int         size_buf;

  ps_Set_Line_Width ( wd );
  buf = (point2f*)pkv_GetScratchMem ( size_buf =
                    max( NN-n, MM-m )*sizeof(point2f) + sizeof(float) );
  for ( i = 0; i < NN-n; i++ ) {
    acp = &cp[3*i*(MM-m)];
    for ( j = 0; j < MM-m; j++ ) {
      memcpy ( &p, &acp[3*j], 3*sizeof(float) );
      CameraProjectPoint3f ( &CPos, &p, (point3f*)&buf[j] );
    }
    ps_Draw_Polyline2f ( MM-m, buf );
  }
  for ( j = 0; j < MM-m; j++ ) {
    acp = &cp[3*j];
    for ( i = 0; i < NN-n; i++ ) {
      memcpy ( &p, &acp[3*i*(MM-m)], 3*sizeof(float) );
      CameraProjectPoint3f ( &CPos, &p, (point3f*)&buf[i] );
    }
    ps_Draw_Polyline2f ( NN-n, buf );
  }
  for ( i = 0; i < NN-n; i++ ) {
    acp = &cp[3*i*(MM-m)];
    for ( j = 0; j < MM-m; j++ ) {
      memcpy ( &p, &acp[3*j], 3*sizeof(float) );
      CameraProjectPoint3f ( &CPos, &p, &q );
      ps_Mark_Circle ( q.x, q.y );
    }
  }
  
  pkv_FreeScratchMem ( size_buf );
} /*DisplayBSCNet*/

static void DisplayBSPatch ( int n, int NN, const float *u,
                             int m, int MM, const float *v,
                             const float *cp,
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

  ku = mbs_NumKnotIntervalsf ( n, NN, u );
  kv = mbs_NumKnotIntervalsf ( m, MM, v );
  pitch = 3*kv*(m+1);
  cpa = pkv_GetScratchMemf ( 3*ku*kv*(n+1)*(m+1) );
  mbs_BSPatchToBezf ( 3, n, NN, u, m, MM, v, 3*(MM-m), cp,
                      &ku, &NNa, NULL, &kv, &MMa, NULL, pitch, cpa );

  cpb = pkv_GetScratchMemf ( pitch*sizeof(float) );
  c = (point2f*)pkv_GetScratchMem ( (max(du,dv)*dd+1)*sizeof(point2f) );
                             /* draw lines of constant parameter u */
  for ( k = 0; k < ku; k++ ) {
    ps_Set_Line_Width ( wd1 );
    for ( i = 0; i < du; i++ ) {
      t = (float)i/(float)du;
      mbs_multiBCHornerf ( n, 1, pitch, 0, &cpa[k*(n+1)*pitch], t, cpb );
      for ( l = 0; l < kv; l++ ) {
        for ( j = 0; j <= dv*dd; j++ ) {
          t = (float)j/(float)(dv*dd);
          mbs_BCHornerC3f ( m, (point3f*)&cpb[l*(m+1)*3], t, &p );
          CameraProjectPoint3f ( &CPos, &p, (point3f*)&c[j] );
        }
        ps_Draw_Polyline2f ( dv*dd+1, c );
      }
      ps_Set_Line_Width ( wd2 );  
    }
  }
  ps_Set_Line_Width ( wd1 );
  mbs_multiBCHornerf ( n, 1, pitch, 0, &cpa[(ku-1)*(n+1)*pitch], 1.0, cpb );
  for ( l = 0; l < kv; l++ ) {
    for ( j = 0; j <= dv*dd; j++ ) {
      t = (float)j/(float)(dv*dd);
      mbs_BCHornerC3f ( m, (point3f*)&cpb[l*(m+1)*3], t, &p );
      CameraProjectPoint3f ( &CPos, &p, (point3f*)&c[j] );
    }
    ps_Draw_Polyline2f ( dv*dd+1, c );
  }
                             /* draw lines of constant parameter v */
  for ( l = 0; l < kv; l++ ) {
    ps_Set_Line_Width ( wd1 );
    for ( j = 0; j < dv; j++ ) {
      t = (float)j/(float)dv;
      mbs_multiBCHornerf ( m, NNa-n, 3, pitch, &cpa[l*(m+1)*3], t, cpb );
      for ( k = 0; k < ku; k++ ) {
        for ( i = 0; i <= du*dd; i++ ) {
          t = (float)i/(float)(du*dd);
          mbs_BCHornerC3f ( n, (point3f*)&cpb[k*(n+1)*3], t, &p );
          CameraProjectPoint3f ( &CPos, &p, (point3f*)&c[i] );
        }
        ps_Draw_Polyline2f ( du*dd+1, c );
      }
      ps_Set_Line_Width ( wd2 );
    }
  }
  ps_Set_Line_Width ( wd1 );
  mbs_multiBCHornerf ( m, NNa-n, 3, pitch, &cpa[(kv-1)*(m+1)*3], 1.0, cpb );
  for ( k = 0; k < ku; k++ ) {
    for ( i = 0; i <= du*dd; i++ ) {
      t = (float)i/(float)(du*dd);
      mbs_BCHornerC3f ( n, (point3f*)&cpb[k*(n+1)*3], t, &p );
      CameraProjectPoint3f ( &CPos, &p, (point3f*)&c[i] );
    }
    ps_Draw_Polyline2f ( du*dd+1, c );
  }
                             /* deallocate buffers */
  pkv_SetScratchMemTop ( st );
} /*DisplayBSPatch*/

static void DrawDDots ( float x, float dx, float y, float dy, int n )
{
  int i;
  float r;

  r = (float)sqrt(dx*dx+dy*dy);
  dx *= (float)(20.0/r);
  dy *= (float)(20.0/r);

  x += dy;
  y -= dx;
  for ( i = 0; i < n; i++, x += dy, y -= dx )
    ps_Fill_Circle ( x, y, 6.0 );
} /*DrawDDots*/

static void DisplayBSDomain ( int n, int NN, const float *u,
                              int m, int MM, const float *v,
                              point2f *frame )
{
  vector2f fv[2], fr[3];
  point2f  corner[5], c;
  point2f  aux[2];
  int      i, k;

                             /* compute the corners of the domain image */
  SubtractPoints2f ( &frame[1], &frame[0], &fv[0] );
  SubtractPoints2f ( &frame[2], &frame[0], &fv[1] );
  AddVector2Mf ( &frame[0], &fv[0], u[n], &corner[0] );
  AddVector2Mf ( &corner[0], &fv[1], v[m], &corner[0] );
  AddVector2Mf ( &frame[0], &fv[0], u[NN-n], &corner[1] );
  AddVector2Mf ( &corner[1], &fv[1], v[m], &corner[1] );
  AddVector2Mf ( &frame[0], &fv[0], u[n], &corner[3] );
  AddVector2Mf ( &corner[3], &fv[1], v[MM-m], &corner[3] );
  AddVector2Mf ( &frame[0], &fv[0], u[NN-n], &corner[2] );
  AddVector2Mf ( &corner[2], &fv[1], v[MM-m], &corner[2] );
  corner[4] = corner[0];

  AddVector2Mf ( &frame[0], &fv[0], u[0], &fr[0] );
  AddVector2Mf ( &fr[0], &fv[1], v[0], &fr[0] );
  AddVector2Mf ( &frame[0], &fv[0], u[NN], &fr[1] );
  AddVector2Mf ( &fr[1], &fv[1], v[0], &fr[1] );
  AddVector2Mf ( &frame[0], &fv[0], u[0], &fr[2] );
  AddVector2Mf ( &fr[2], &fv[1], v[MM], &fr[2] );
  ps_Set_Line_Width ( 2.0 );
  ps_Draw_Line ( fr[0].x, fr[0].y, fr[1].x, fr[1].y );
  ps_Draw_Line ( fr[0].x, fr[0].y, fr[2].x, fr[2].y );

                             /* draw the domain outline */
  ps_Set_Line_Width ( 3.0 );
  ps_Draw_Polyline2f ( 5, corner );
                             /* draw polynomial junction lines */
  ps_Set_Line_Width ( 1.0 );
  for ( i = n+1; i < NN-n; i++ ) {
    if ( u[i] > u[i-1] ) {
      AddVector2Mf ( &frame[0], &fv[0], u[i], &aux[0] );  aux[1] = aux[0];
      AddVector2Mf ( &aux[0], &fv[1], v[m], &aux[0] );
      AddVector2Mf ( &aux[1], &fv[1], v[MM-m], &aux[1] );
      ps_Draw_Line ( aux[0].x, aux[0].y, aux[1].x, aux[1].y );
    }
  }
  for ( i = m+1; i < MM-m; i++ ) {
    if ( v[i] > v[i-1] ) {
      AddVector2Mf ( &frame[0], &fv[1], v[i], &aux[0] );  aux[1] = aux[0];
      AddVector2Mf ( &aux[0], &fv[0], u[n], &aux[0] );
      AddVector2Mf ( &aux[1], &fv[0], u[NN-n], &aux[1] );
      ps_Draw_Line ( aux[0].x, aux[0].y, aux[1].x, aux[1].y );
    }
  }
                             /* draw knots */
  psl_SetLine ( frame[0].x, fr[0].y, frame[0].x+fv[0].x, fr[0].y, 0.0, 1.0 );
  for ( i = 1, k = 1; i <= NN; i++ ) {
    if ( u[i] > u[i-1] ) {
      if ( i < NN )
        psl_Tick ( u[i] );
      if ( i > 1 && k > 0 ) {
        AddVector2Mf ( &frame[0], &fv[0], u[i-1], &c );
        AddVector2Mf ( &c, &fv[1], v[m], &c );
        DrawDDots ( c.x, fv[0].x, fr[0].y, fv[0].y, k );
      }
      k = 1;
    }
    else
      k++;
  }
  psl_SetLine ( fr[0].x, frame[0].y, fr[0].x, frame[0].y+fv[1].y, 0.0, 1.0 );
  for ( i = 1, k = 1; i <= MM; i++ ) {
    if ( v[i] > v[i-1] ) {
      if ( i < MM )
        psl_Tick ( v[i] );
      if ( i > 1 && k > 0 ) {
        AddVector2Mf ( &frame[0], &fv[1], v[i-1], &c );
        AddVector2Mf ( &c, &fv[0], u[n], &c );
        DrawDDots ( fr[0].x, -fv[1].x, c.y, -fv[1].y, k );
      }
      k = 1;
    }
    else
      k++;
  }
} /*DisplayBSDomain*/

int main ()
{
  int   NNa, MMa, ku, kv, pitch;
  float *ua, *va, *cpa;

  pkv_InitScratchMem ( SCRATCHMEMSIZE );
  ps_WriteBBox ( 24, 10, 371, 145 );
  ps_OpenFile ( fn, 600 );
  ps_DenseScreen ();
  ps_Write_Command ( "1 setlinecap" );
  SetupCamera ();

  NNa = mbs_LastknotMaxInsf ( n, NN, u, &ku );
  MMa = mbs_LastknotMaxInsf ( m, MM, v, &kv );
  pitch = 3*kv*(m+1);
  ua = (float*)pkv_GetScratchMem ( (NNa+1)*sizeof(float) );
  va = (float*)pkv_GetScratchMem ( (MMa+1)*sizeof(float) );
  cpa = pkv_GetScratchMemf ( 3*ku*kv*(n+1)*(m+1) );
  mbs_BSPatchToBezf ( 3, n, NN, u, m, MM, v, 3*(MM-m), (float*)cp,
                      &ku, &NNa, ua, &kv, &MMa, va, pitch, cpa );

  ps_GSave ();
  ps_Write_Command ( "0 500 translate" );
  ps_Set_Gray ( 0.68 );
  DisplayBSPatch ( n, NN, u, m, MM, v, (float*)cp, 6, 6, 4, 6.0, 2.0 );
  ps_Set_Gray ( 0.0 );
  DisplayBSCNet ( n, NN, m, MM, (float*)cp, 4.0 );
  DisplayBSDomain ( n, NN, u, m, MM, v, frame );
  ps_GRestore ();

  ps_GSave ();
  ps_Write_Command ( "1500 500 translate" );
  ps_Set_Gray ( 0.68 );
  DisplayBSPatch ( n, NN, u, m, MM, v, (float*)cp, 6, 6, 4, 6.0, 2.0 );
  ps_Set_Gray ( 0.0 );
  DisplayBSCNet ( n, NNa, m, MMa, (float*)cpa, 4.0 );
  DisplayBSDomain ( n, NNa, ua, m, MMa, va, frame );
  ps_GRestore ();

  ps_CloseFile ();
  printf ( "%s\n", fn );
  pkv_DestroyScratchMem ();
  exit ( 0 );
} /*main*/

