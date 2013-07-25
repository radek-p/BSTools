
/* //////////////////////////////////////////////////// */
/* This file is a part of the BSTools procedure package */
/* written by Przemyslaw Kiciak.                        */
/* //////////////////////////////////////////////////// */

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <sys/times.h>
#include <unistd.h>

#include "pkvaria.h"
#include "pknum.h"
#include "pkgeom.h"
#include "camera.h"
#include "psout.h"
#include "multibs.h"  

#include "eg1holed.h"
#include "datagend.h"
#include "testgraphd.h"

#define HOLE_K 6

char fn1[]  = "g1hnlconstrd.ps";

CameraRecd   CPos;

GHoleDomaind *domain;

int nfinalp;
point3d FinalCPoints[2*HOLE_K][121];

double iparams[3] = {0.5,0.0,0.5};
point3d cnt[2];

#define DBDIM 4

/* ///////////////////////////////////////////////////////////////////////// */
void SetupCamera ( double w, double h, double x, double y )
{
  CameraInitFramed ( &CPos, false, true, w, h, x, y, 1.0, 4 );
  CameraInitPosd ( &CPos );
  SetPoint3d ( &CPos.position, 20.402914, -7.444557, -26.986362 );
  CPos.psi = -0.173305;  CPos.theta = 0.677664;  CPos.phi = -1.920662;
  CPos.vd.persp.f = 4.0;
  CameraSetMappingd ( &CPos );
} /*SetupCamera*/

void DrawBezCurve3d ( int n, const point3d *cp, double t0, double t1 )
{
#define DD 32
  void    *sp;
  int     i;
  point3d p, q;
  point2d *cc;

  sp = pkv_GetScratchMemTop ();
  cc = (point2d*)pkv_GetScratchMem ( (DD+1)*sizeof(point2d) );
  if ( !cc )
    exit ( 1 );

  for ( i = 0; i <= DD; i++ ) {
    mbs_BCHornerC3d ( n, cp, (double)i/(double)DD, &p );
    CameraProjectPoint3d ( &CPos, &p, &q );
    SetPoint2d ( &cc[i], q.x, q.y );
  }
  ps_Draw_Polyline2d ( DD+1, cc );

  pkv_SetScratchMemTop ( sp );
#undef DD
} /*DrawBezCurve3d*/

void SetLW ( int i, int i0, int i1 )
{
    if ( i == i0 || i == i1 )
      ps_Set_Line_Width ( 4.0 );
    else
      ps_Set_Line_Width ( 2.0 );
} /*SetLW*/

void DrawBezPatch3d ( int n, int m, const point3d *cp )
{
#define D 8
  void    *sp;
  int     i;
  point3d *cc;

  sp = pkv_GetScratchMemTop ();
  cc = (point3d*)pkv_GetScratchMem ( (max(n,m)+1)*sizeof(point3d) );
  if ( !cc )
    exit ( 1 );

  for ( i = 0; i <= D; i++ ) {
    mbs_multiBCHornerd ( n, 1, 3*(m+1), 0, (double*)cp, (double)i/(double)D,
                         (double*)cc );
    SetLW ( i, 0, D );
    DrawBezCurve3d ( m, cc, 0.0, 1.0 );
  }
  for ( i = 0; i <= D; i++ ) {
    mbs_multiBCHornerd ( m, n+1, 3, 3*(m+1), (double*)cp, (double)i/(double)D,
                         (double*)cc );
    SetLW ( i, 0, D );
    DrawBezCurve3d ( n, cc, 0.0, 1.0 );
  }

  pkv_SetScratchMemTop ( sp );
#undef D
} /*DrawBezPatch3d*/

void DrawPatchNet3d ( int n, int m, const point3d *cp )
{
  void    *sp;
  int     i, j;
  point3d q;
  point2d *r;

  sp = pkv_GetScratchMemTop ();
  r = pkv_GetScratchMem ( (n+1)*(m+1)*sizeof(point2d) );
  if ( !r )
    exit ( 1 );
  for ( i = 0; i < (n+1)*(m+1); i++ ) {
    CameraProjectPoint3d ( &CPos, &cp[i], &q );
    SetPoint2d ( &r[i], q.x, q.y );
  }
  ps_Set_Line_Width ( 1.0 );
  for ( i = 0; i <= n; i++ )
    ps_Draw_Polyline2d ( m+1, &r[i*(m+1)] );
  for ( j = 0; j <= m; j++ )
    for ( i = 0; i < n; i++ )
      ps_Draw_Line ( r[i*(m+1)+j].x, r[i*(m+1)+j].y,
                     r[(i+1)*(m+1)+j].x, r[(i+1)*(m+1)+j].y );
  for ( i = 0; i < (n+1)*(m+1); i++ )
    ps_Mark_Circle ( r[i].x, r[i].y );

  pkv_SetScratchMemTop ( sp );
} /*DrawPatchNet3d*/

void ShowPoint ( point3d *p )
{
  point3d q;

  CameraProjectPoint3d ( &CPos, p, &q );
  ps_Mark_Circle ( q.x, q.y );
} /*ShowPoint*/

/* ///////////////////////////////////////////////////////////////////////// */
void GetHoleSurrndPatch ( int i, int j, point3d *dombcp )
{
  int     k;
  int     *ind;
  point3d *q;
  void    *sp;
  double   *ukn, *vkn;

  sp  = pkv_GetScratchMemTop ();
  ind = pkv_GetScratchMem ( 16*sizeof(int) );
  q   = pkv_GetScratchMem ( 16*sizeof(point3d) );
  if ( !ind || !q )
    exit ( 1 );

  gh_GetBspInd ( HOLE_K, i, j, ind );
  for ( k = 0; k < 16; k++ )
    q[k] = Bpt[ind[k]];

  ukn = &knots[11*((i+HOLE_K-1) % HOLE_K)+3];
  vkn = &knots[11*i+j];
  mbs_BSPatchToBezd ( 3, 3, 7, ukn, 3, 7, vkn, 12, (double*)q,
                      NULL, NULL, NULL, NULL, NULL, NULL, 12, (double*)dombcp );

  pkv_SetScratchMemTop ( sp );
} /*GetHoleSurrndPatch*/

void DrawBezPatches ()
{
  int     i, j;
  point3d cp[16];

  ps_Set_Gray ( 0.72 );
  for ( i = 0; i < HOLE_K; i++ )
    for ( j = 0; j < 3; j++ ) {
      GetHoleSurrndPatch ( i, j, cp );
      DrawBezPatch3d ( 3, 3, cp );
    }
} /*DrawBezPatches*/

void DrawFinalPatches ( int p0, int p1 )
{
  int i;

  ps_Set_Gray ( 0.25 );
  for ( i = p0; i < p1; i++ )
    DrawBezPatch3d ( G1H_FINALDEG, G1H_FINALDEG, FinalCPoints[i] );
} /*DrawFinalPatches*/

void MakePicture1 ()
{
  ps_OpenFile ( fn1, 600 );
  ps_Write_Command ( "1 setlinecap" );
  SetupCamera ( 1800, 1460, 0.0, 0.0 );

  DrawBezPatches ();
  DrawFinalPatches ( 0, HOLE_K );
  ShowPoint ( &cnt[0] );

  ps_GSave ();
  ps_Write_Command ( "1800 0 translate" );
  DrawBezPatches ();
  DrawFinalPatches ( HOLE_K, 2*HOLE_K );
  ShowPoint ( &cnt[0] );
  ps_GRestore ();

  ps_CloseFile ();
  printf ( "%s\n", fn1 );
} /*MakePicture1*/

void MakePictures ()
{
  MakePicture1 ();
} /*MakePictures*/

/* ///////////////////////////////////////////////////////////////////////// */
void OutPatch ( int n, int m, const point3d *cp, void *usrptr )
{
  if ( n != G1H_FINALDEG || m != G1H_FINALDEG || nfinalp >= 2*HOLE_K )
    exit ( 1 );
  memcpy ( FinalCPoints[nfinalp], cp, 121*sizeof(point3d) );
  nfinalp ++;
} /*OutPatch*/

int main ()
{
  struct tms start, stop;
  double  time;
  double  param[2] = {0.0,0.0};
  double  *bfder, *cmat;
  int    nfunc_a, nfunc_c;
  double  fval, efval, *acoeff, *bcoeff;

  pkv_InitScratchMem ( 2097152 );
  acoeff = pkv_GetScratchMemd ( 90 );
  bcoeff = pkv_GetScratchMemd ( 474 );
  if ( !acoeff || !bcoeff )
    goto failure;
  InitKnots ( HOLE_K );
  InitDomain ( HOLE_K, param );
  times ( &start );
  if ( (domain = gh_CreateDomaind ( HOLE_K, knots, Dompt )) ) {
    if ( g1h_ComputeBasisd ( domain ) ) {
      nfinalp = 0;
      InitHole ( HOLE_K, iparams );

                 /* setup the constraints */
      nfunc_a = g1h_V0SpaceDimd ( domain );
      nfunc_c = HOLE_K*DBDIM;
      bfder = pkv_GetScratchMemd ( nfunc_a*HOLE_K*3 );
      if ( !bfder )
        goto failure;
      if ( !g1h_GetBPDerivativesd ( domain, 0, bfder ) )
        goto failure;
      cmat = pkv_GetScratchMemd ( 2*(nfunc_a+nfunc_c) );
      if ( !cmat )
        goto failure;
      pkv_Selectd ( nfunc_a, 1, 3, 1, bfder, cmat );
      pkv_Selectd ( nfunc_a, 1, 3, 1, &bfder[1], &cmat[nfunc_a] );
      if ( !g1h_SetConstraintMatrixd ( domain, 2, cmat ) )
        goto failure;
      SetPoint3d ( &cnt[0], 0.15, 0.05, -0.6 );
      SetPoint3d ( &cnt[1], 1.0, 0.0, 0.5 );
      if ( !g1h_NLFillHoleConstrd ( domain, Bpt,
                                    2, cnt, acoeff, NULL, OutPatch ) )
        goto failure;

      memset ( cmat, 0, 2*(nfunc_a+nfunc_c)*sizeof(double) );
      pkv_Selectd ( nfunc_a, 1, 3, 1, bfder, &cmat[nfunc_c] );
      pkv_Selectd ( nfunc_a, 1, 3, 1, &bfder[1], &cmat[nfunc_a+nfunc_c+nfunc_c] );
      if ( !g1h_SetExtConstraintMatrixd ( domain, 2, cmat ) )
        goto failure;
      if ( !g1h_NLExtFillHoleConstrd ( domain, Bpt,
                                       2, cnt, bcoeff, NULL, OutPatch ) )
        goto failure;

      times ( &stop );
      fval = g1h_FunctionalValued ( domain, 3, (double*)Bpt, acoeff );
      efval = g1h_ExtFunctionalValued ( domain, 3, (double*)Bpt, bcoeff );
      printf ( "Functional values = %f, %f\n", fval, efval );
      MakePictures ();
      time = (double)(stop.tms_utime-start.tms_utime)/
             (double)(sysconf(_SC_CLK_TCK));
      printf ( "time = %8.3f\n", time );
    }
    gh_DestroyDomaind ( domain );
  }

  printf ( "Scratch memory used: %d bytes\n", pkv_MaxScratchTaken() );
  pkv_DestroyScratchMem ();
  exit ( 0 );

failure:
  exit ( 1 );
} /*main*/

