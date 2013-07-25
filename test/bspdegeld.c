
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

#define fn "bspdegeld.ps"
#define SCRATCHMEMSIZE 262144

#define  n1   3
#define  NN1 10
int n = n1;
int NN = NN1;
double u[NN1+1] = {-0.5, -0.4, 0.0, 0.0, 1.0, 2.0, 3.0, 4.0, 4.0, 4.4, 4.5};
#define  m1   3
#define  MM1 10
int m = m1;
int MM = MM1;
double v[MM1+1] = {-0.5, -0.4, 0.0, 0.0, 1.0, 2.5, 2.5, 4.0, 4.2, 4.3, 4.5};
point3d cp[NN1-n1][MM1-m1] =
  {{{-3.0,-3.0,-0.5},{-2.7,-3.0,-0.5},{-1.5,-3.0,-0.5},{0.0,-3.0,-0.5},{1.5,-3.0,-0.5},{2.7,-3.0,-0.5},{3.0,-3.0,-0.5}},
   {{-3.0,-2.7, 0.0},{-2.7,-2.7, 0.0},{-1.5,-2.7, 0.0},{0.0,-2.7, 0.0},{1.5,-2.7, 0.0},{2.7,-2.7, 0.0},{3.0,-2.7, 0.0}},
   {{-3.0,-1.5, 0.1},{-2.7,-1.5, 0.1},{-1.5,-1.5, 0.1},{0.0,-1.5, 0.1},{1.5,-1.5, 0.1},{2.7,-1.5, 0.1},{3.0,-1.5, 0.1}},
   {{-3.0, 0.0, 0.2},{-2.7, 0.0, 0.2},{-1.5, 0.0, 0.2},{0.0, 0.0, 0.2},{1.5, 0.0, 0.2},{2.7, 0.0, 0.2},{3.0, 0.0, 0.2}},
   {{-3.0,+1.5, 0.1},{-2.7,+1.5, 0.1},{-1.5,+1.5, 0.1},{0.0,+1.5, 0.1},{1.5,+1.5, 0.1},{2.7,+1.5, 0.1},{3.0,+1.5, 0.1}},
   {{-3.0,+2.7, 0.0},{-2.7,+2.7, 0.0},{-1.5,+2.7, 0.0},{0.0,+2.7, 0.0},{1.5,+2.7, 0.0},{2.7,+2.7, 0.0},{3.0,+2.7, 0.0}},
   {{-3.0,+3.0,-0.5},{-2.7,+3.0,-0.5},{-1.5,+3.0,-0.5},{0.0,+3.0,-0.5},{1.5,+3.0,-0.5},{2.7,+3.0,-0.5},{3.0,+3.0,-0.5}}};

point2d frame[3] = {{300,-300},{400,-300},{300,-200}};

CameraRecd CPos;


static void SetupCamera ( void )
{
  vector3d v;

  CameraInitFramed ( &CPos, false, true, 1200, 900, 400, 0, 1.0, 4 );
  CameraInitPosd ( &CPos );

  SetVector3d ( &v, 0.0, 0.0, -50.0 );
  CameraMoveGd ( &CPos, &v );
  CameraRotXCd ( &CPos, -0.65*PI );
  CameraRotZGd ( &CPos, -0.14*PI );
  CameraZoomd ( &CPos, 5.0 );
} /*SetupCamera*/

static void DisplayBSCNet ( int n, int NN, int m, int MM,
                            const double *cp, double wd )
{
  int         i, j;
  point3d     p, q;
  point2d     *buf; 
  const double *acp;
  int         size_buf;

  ps_Set_Line_Width ( (float)wd );
  buf = (point2d*)pkv_GetScratchMem ( size_buf =
                    max( NN-n, MM-m )*sizeof(point2d) + sizeof(double) );
  for ( i = 0; i < NN-n; i++ ) {
    acp = &cp[3*i*(MM-m)];
    for ( j = 0; j < MM-m; j++ ) {
      memcpy ( &p, &acp[3*j], 3*sizeof(double) );
      CameraProjectPoint3d ( &CPos, &p, (point3d*)&buf[j] );
    }
    ps_Draw_Polyline2d ( MM-m, buf );
  }
  for ( j = 0; j < MM-m; j++ ) {
    acp = &cp[3*j];
    for ( i = 0; i < NN-n; i++ ) {
      memcpy ( &p, &acp[3*i*(MM-m)], 3*sizeof(double) );
      CameraProjectPoint3d ( &CPos, &p, (point3d*)&buf[i] );
    }
    ps_Draw_Polyline2d ( NN-n, buf );
  }
  for ( i = 0; i < NN-n; i++ ) {
    acp = &cp[3*i*(MM-m)];
    for ( j = 0; j < MM-m; j++ ) {
      memcpy ( &p, &acp[3*j], 3*sizeof(double) );
      CameraProjectPoint3d ( &CPos, &p, &q );
      ps_Mark_Circle ( (float)q.x, (float)q.y );
    }
  }
  
  pkv_FreeScratchMem ( size_buf );
} /*DisplayBSCNet*/

static void DisplayBSPatch ( int n, int NN, const double *u,
                             int m, int MM, const double *v,
                             const double *cp,
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

  ku = mbs_NumKnotIntervalsd ( n, NN, u );
  kv = mbs_NumKnotIntervalsd ( m, MM, v );
  pitch = 3*kv*(m+1);
  cpa = pkv_GetScratchMemd ( 3*ku*kv*(n+1)*(m+1) );
  mbs_BSPatchToBezd ( 3, n, NN, u, m, MM, v, 3*(MM-m), cp,
                      &ku, &NNa, NULL, &kv, &MMa, NULL, pitch, cpa );

  cpb = pkv_GetScratchMemd ( pitch*sizeof(double) );
  c = (point2d*)pkv_GetScratchMem ( (max(du,dv)*dd+1)*sizeof(point2d) );
                             /* draw lines of constant parameter u */
  for ( k = 0; k < ku; k++ ) {
    ps_Set_Line_Width ( (float)wd1 );
    for ( i = 0; i < du; i++ ) {
      t = (double)i/(double)du;
      mbs_multiBCHornerd ( n, 1, pitch, 0, &cpa[k*(n+1)*pitch], t, cpb );
      for ( l = 0; l < kv; l++ ) {
        for ( j = 0; j <= dv*dd; j++ ) {
          t = (double)j/(double)(dv*dd);
          mbs_BCHornerC3d ( m, (point3d*)&cpb[l*(m+1)*3], t, &p );
          CameraProjectPoint3d ( &CPos, &p, (point3d*)&c[j] );
        }
        ps_Draw_Polyline2d ( dv*dd+1, c );
      }
      ps_Set_Line_Width ( (float)wd2 );  
    }
  }
  ps_Set_Line_Width ( (float)wd1 );
  mbs_multiBCHornerd ( n, 1, pitch, 0, &cpa[(ku-1)*(n+1)*pitch], 1.0, cpb );
  for ( l = 0; l < kv; l++ ) {
    for ( j = 0; j <= dv*dd; j++ ) {
      t = (double)j/(double)(dv*dd);
      mbs_BCHornerC3d ( m, (point3d*)&cpb[l*(m+1)*3], t, &p );
      CameraProjectPoint3d ( &CPos, &p, (point3d*)&c[j] );
    }
    ps_Draw_Polyline2d ( dv*dd+1, c );
  }
                             /* draw lines of constant parameter v */
  for ( l = 0; l < kv; l++ ) {
    ps_Set_Line_Width ( (float)wd1 );
    for ( j = 0; j < dv; j++ ) {
      t = (double)j/(double)dv;
      mbs_multiBCHornerd ( m, NNa-n, 3, pitch, &cpa[l*(m+1)*3], t, cpb );
      for ( k = 0; k < ku; k++ ) {
        for ( i = 0; i <= du*dd; i++ ) {
          t = (double)i/(double)(du*dd);
          mbs_BCHornerC3d ( n, (point3d*)&cpb[k*(n+1)*3], t, &p );
          CameraProjectPoint3d ( &CPos, &p, (point3d*)&c[i] );
        }
        ps_Draw_Polyline2d ( du*dd+1, c );
      }
      ps_Set_Line_Width ( (float)wd2 );
    }
  }
  ps_Set_Line_Width ( (float)wd1 );
  mbs_multiBCHornerd ( m, NNa-n, 3, pitch, &cpa[(kv-1)*(m+1)*3], 1.0, cpb );
  for ( k = 0; k < ku; k++ ) {
    for ( i = 0; i <= du*dd; i++ ) {
      t = (double)i/(double)(du*dd);
      mbs_BCHornerC3d ( n, (point3d*)&cpb[k*(n+1)*3], t, &p );
      CameraProjectPoint3d ( &CPos, &p, (point3d*)&c[i] );
    }
    ps_Draw_Polyline2d ( du*dd+1, c );
  }
                             /* deallocate buffers */
  pkv_SetScratchMemTop ( st );
} /*DisplayBSPatch*/

static void DrawDDots ( double x, double dx, double y, double dy, int n )
{
  int i;
  double r;
  r = sqrt(dx*dx+dy*dy);
  dx *= 20.0/r;
  dy *= 20.0/r;

  x += dy;
  y -= dx;
  for ( i = 0; i < n; i++, x += dy, y -= dx )
    ps_Fill_Circle ( (float)x, (float)y, 6.0 );
} /*DrawDDots*/

static void DisplayBSDomain ( int n, int NN, const double *u,
                              int m, int MM, const double *v,
                              point2d *frame )
{
  vector2d fv[2], fr[3];
  point2d  corner[5], c;
  point2d  aux[2];
  int      i, k;

                             /* compute the corners of the domain image */
  SubtractPoints2d ( &frame[1], &frame[0], &fv[0] );
  SubtractPoints2d ( &frame[2], &frame[0], &fv[1] );
  AddVector2Md ( &frame[0], &fv[0], u[n], &corner[0] );
  AddVector2Md ( &corner[0], &fv[1], v[m], &corner[0] );
  AddVector2Md ( &frame[0], &fv[0], u[NN-n], &corner[1] );
  AddVector2Md ( &corner[1], &fv[1], v[m], &corner[1] );
  AddVector2Md ( &frame[0], &fv[0], u[n], &corner[3] );
  AddVector2Md ( &corner[3], &fv[1], v[MM-m], &corner[3] );
  AddVector2Md ( &frame[0], &fv[0], u[NN-n], &corner[2] );
  AddVector2Md ( &corner[2], &fv[1], v[MM-m], &corner[2] );
  corner[4] = corner[0];

  AddVector2Md ( &frame[0], &fv[0], u[0], &fr[0] );
  AddVector2Md ( &fr[0], &fv[1], v[0], &fr[0] );
  AddVector2Md ( &frame[0], &fv[0], u[NN], &fr[1] );
  AddVector2Md ( &fr[1], &fv[1], v[0], &fr[1] );
  AddVector2Md ( &frame[0], &fv[0], u[0], &fr[2] );
  AddVector2Md ( &fr[2], &fv[1], v[MM], &fr[2] );
  ps_Set_Line_Width ( 2.0 );
  ps_Draw_Line ( (float)fr[0].x, (float)fr[0].y, (float)fr[1].x, (float)fr[1].y );
  ps_Draw_Line ( (float)fr[0].x, (float)fr[0].y, (float)fr[2].x, (float)fr[2].y );

                             /* draw the domain outline */
  ps_Set_Line_Width ( 3.0 );
  ps_Draw_Polyline2d ( 5, corner );
                             /* draw polynomial junction lines */
  ps_Set_Line_Width ( 1.0 );
  for ( i = n+1; i < NN-n; i++ ) {
    if ( u[i] > u[i-1] ) {
      AddVector2Md ( &frame[0], &fv[0], u[i], &aux[0] );  aux[1] = aux[0];
      AddVector2Md ( &aux[0], &fv[1], v[m], &aux[0] );
      AddVector2Md ( &aux[1], &fv[1], v[MM-m], &aux[1] );
      ps_Draw_Line ( (float)aux[0].x, (float)aux[0].y, (float)aux[1].x, (float)aux[1].y );
    }
  }
  for ( i = m+1; i < MM-m; i++ ) {
    if ( v[i] > v[i-1] ) {
      AddVector2Md ( &frame[0], &fv[1], v[i], &aux[0] );  aux[1] = aux[0];
      AddVector2Md ( &aux[0], &fv[0], u[n], &aux[0] );
      AddVector2Md ( &aux[1], &fv[0], u[NN-n], &aux[1] );
      ps_Draw_Line ( (float)aux[0].x, (float)aux[0].y, (float)aux[1].x, (float)aux[1].y );
    }
  }
                             /* draw knots */
  psl_SetLine ( (float)frame[0].x, (float)fr[0].y,
                (float)(frame[0].x+fv[0].x), (float)fr[0].y, 0.0, 1.0 );
  for ( i = 1, k = 1; i <= NN; i++ ) {
    if ( u[i] > u[i-1] ) {
      if ( i < NN )
        psl_Tick ( (float)u[i] );
      if ( i > 1 && k > 0 ) {
        AddVector2Md ( &frame[0], &fv[0], u[i-1], &c );
        AddVector2Md ( &c, &fv[1], v[m], &c );
        DrawDDots ( c.x, fv[0].x, fr[0].y, fv[0].y, k );
      }
      k = 1;
    }
    else
      k++;
  }
  psl_SetLine ( (float)fr[0].x, (float)frame[0].y, (float)fr[0].x,
                (float)(frame[0].y+fv[1].y), 0.0, 1.0 );
  for ( i = 1, k = 1; i <= MM; i++ ) {
    if ( v[i] > v[i-1] ) {
      if ( i < MM )
        psl_Tick ( (float)v[i] );
      if ( i > 1 && k > 0 ) {
        AddVector2Md ( &frame[0], &fv[1], v[i-1], &c );
        AddVector2Md ( &c, &fv[0], u[n], &c );
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
  int   na, NNa, ma, MMa, ku, kv;
  double *ua, *va, *cpa;

  pkv_InitScratchMem ( SCRATCHMEMSIZE );
  ps_WriteBBox ( 24, 9, 374, 288 );
  ps_OpenFile ( fn, 600 );
  ps_DenseScreen ();
  SetupCamera ();

  ku = mbs_NumKnotIntervalsd ( n, NN, u );
  kv = mbs_NumKnotIntervalsd ( m, MM, v );
  ua = pkv_GetScratchMem ( (NN+2+ku)*sizeof(double) );
  va = pkv_GetScratchMem ( (MM+2+kv)*sizeof(double) );
  cpa = pkv_GetScratchMem ( (NN-n+ku)*(MM-m+kv)*sizeof(point3d) );

  ps_GSave ();
  ps_Write_Command ( "0 1700 translate" );
  ps_Set_Gray ( 0.68 );
  DisplayBSPatch ( n, NN, u, m, MM, v, (double*)cp, 6, 6, 4, 6.0, 2.0 );
  ps_Set_Gray ( 0.0 );
  DisplayBSCNet ( n, NN, m, MM, (double*)cp, 4.0 );
  DisplayBSDomain ( n, NN, u, m, MM, v, frame );
  ps_GRestore ();

  ps_GSave ();
  ps_Write_Command ( "1500 1700 translate" );
  mbs_multiBSDegElevd ( 1, 3*(MM-m), n, NN, u, 0, (double*)cp, 1,
                        &na, &NNa, ua, 0, cpa, false );
  ps_Set_Gray ( 0.68 );
  DisplayBSPatch ( na, NNa, ua, m, MM, v, (double*)cpa, 6, 6, 4, 6.0, 2.0 );
  ps_Set_Gray ( 0.0 );
  DisplayBSCNet ( na, NNa, m, MM, (double*)cpa, 4.0 );
  DisplayBSDomain ( na, NNa, ua, m, MM, v, frame );
  ps_GRestore ();

  ps_GSave ();
  ps_Write_Command ( "0 500 translate" );
  mbs_multiBSDegElevd ( NN-n, 3, m, MM, v, 3*(MM-m), (double*)cp, 1,
                        &ma, &MMa, va, 3*(MM-m+kv), cpa, false );
  ps_Set_Gray ( 0.68 );
  DisplayBSPatch ( n, NN, u, ma, MMa, va, (double*)cpa, 6, 6, 4, 6.0, 2.0 );
  ps_Set_Gray ( 0.0 );
  DisplayBSCNet ( n, NN, ma, MMa, (double*)cpa, 4.0 );
  DisplayBSDomain ( n, NN, u, ma, MMa, va, frame );
  ps_GRestore ();

  ps_GSave ();
  ps_Write_Command ( "1500 500 translate" );
  mbs_multiBSDegElevd ( NN-n, 3, m, MM, v, 3*(MM-m), (double*)cp, 1,
                        &ma, &MMa, va, 3*(MM-m+kv), cpa, false );
  mbs_multiBSDegElevd ( 1, 3*(MMa-ma), n, NN, u, 0, NULL, 1,
                        &na, &NNa, ua, 0, cpa, false );
  ps_Set_Gray ( 0.68 );
  DisplayBSPatch ( na, NNa, ua, ma, MMa, va, (double*)cpa, 6, 6, 4, 6.0, 2.0 );
  ps_Set_Gray ( 0.0 );
  DisplayBSCNet ( na, NNa, ma, MMa, (double*)cpa, 4.0 );
  DisplayBSDomain ( na, NNa, ua, ma, MMa, va, frame );
  ps_GRestore ();

  ps_CloseFile ();
  pkv_DestroyScratchMem ();
  exit ( 0 );
} /*main*/

