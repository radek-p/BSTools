
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

#include "eg2holed.h"
#include "datagend.h"
#include "testgraphd.h"
#include "readdatd.h"

#define HOLE_K 8

char infn[] = "d8a.dat";
char fn1[]  = "g2hconstr1d.ps";
char fn2[]  = "g2hconstr2d.ps";

int hole_k;

CameraRecd   CPos;

GHoleDomaind *domain;

int nfinalp;
point3d FinalCPoints[2*HOLE_K][121];

double iparams[3] = {0.5,0.0,0.5};

#define DBDIM 16

/* ///////////////////////////////////////////////////////////////////////// */
void SetupCamera ( double w, double h, double x, double y )
{
  CameraInitFramed ( &CPos, false, true, w, h, x, y, 1.0, 4 );
  CameraInitPosd ( &CPos );
  SetPoint3d ( &CPos.position, 20.402914, -7.444557, -26.986362 );
  CPos.psi = -0.173305;  CPos.theta = 0.677664;  CPos.phi = -1.920662;
  CPos.vd.persp.f = 3.0;
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

  ps_Set_RGB ( 1.0, 0.0, 0.0 );
  GetHoleSurrndPatch ( 0, 1, cp );
  DrawBezPatch3d ( 3, 3, cp );
  ShowPoint ( &cp[3] );
  ps_Set_RGB ( 0.0, 1.0, 0.0 );
  GetHoleSurrndPatch ( 4, 1, cp );
  DrawBezPatch3d ( 3, 3, cp );
  ShowPoint ( &cp[3] );
} /*DrawBezPatches*/

void DrawFinalPatches ( int p0, int p1 )
{
  int i;

  ps_Set_Gray ( 0.25 );
  for ( i = p0; i < p1; i++ )
    DrawBezPatch3d ( G2H_FINALDEG, G2H_FINALDEG, FinalCPoints[i] );
} /*DrawFinalPatches*/

void MakePictures ( char *fn )
{
  ps_OpenFile ( fn, 600 );
  ps_Write_Command ( "1 setlinecap" );
  SetupCamera ( 1800, 1460, 0.0, 0.0 );

  DrawBezPatches ();
  DrawFinalPatches ( 0, HOLE_K );

  ps_GSave ();
  ps_Write_Command ( "1800 0 translate" );
  DrawBezPatches ();
  DrawFinalPatches ( HOLE_K, 2*HOLE_K );
  ps_GRestore ();

  ps_CloseFile ();
  printf ( "%s\n", fn );
} /*MakePictures*/

/* ///////////////////////////////////////////////////////////////////////// */
void OutPatch ( int n, int m, const double *cp, void *usrptr )
{
  if ( n != G2H_FINALDEG || m != G2H_FINALDEG || nfinalp >= 2*HOLE_K )
    exit ( 1 );
  memcpy ( FinalCPoints[nfinalp], cp, 121*sizeof(point3d) );
  nfinalp ++;
} /*OutPatch*/

void SetupConstraints ( point3d *centp, vector3d *nvect )
{
  void *sp;
  point3d *cp, q, r, s;

  sp = pkv_GetScratchMemTop ();
  cp = (point3d*)pkv_GetScratchMem ( 16*sizeof(point3d) );
  if ( !cp )
    exit ( 1 );
  GetHoleSurrndPatch ( 0, 1, cp );
  q = cp[3];
  GetHoleSurrndPatch ( 4, 1, cp );
  r = cp[3];
  GetHoleSurrndPatch ( 2, 1, cp );
  s = cp[3];
  MidPoint3d ( &q, &r, centp );
  SubtractPoints3d ( &q, &r, &r );
  SubtractPoints3d ( &q, &s, &s );
  CrossProduct3d ( &s, &r, nvect );
  NormalizeVector3d ( nvect );
  printf ( "centp = %f, %f, %f\n", centp->x, centp->y, centp->z );
  printf ( "nvect = %f, %f, %f\n", nvect->x, nvect->y, nvect->z );
  pkv_SetScratchMemTop ( sp );
} /*SetupConstraints*/

void MakeTest1 ()
{
#define NCONSTR 11
  struct   tms start, stop;
  double   time;
  double   *bfder, *cmat;
  vector3d *dmat;
  int      nfunc_a, nfunc_c, i, j, k;
  point3d  centp;
  vector3d nvect;

  ReadDatFile ( infn, &hole_k, knots, Dompt, Bpt, NULL );
  SetupConstraints ( &centp, &nvect );

  times ( &start );
  if ( (domain = gh_CreateDomaind ( HOLE_K, knots, Dompt )) ) {
    if ( g2h_ComputeBasisd ( domain ) ) {
      nfinalp = 0;
                 /* setup the constraints */
      nfunc_a = g2h_V0SpaceDimd ( domain );
      nfunc_c = HOLE_K*DBDIM;
      bfder = pkv_GetScratchMemd ( nfunc_a*5 );
      cmat = pkv_GetScratchMemd ( NCONSTR*(nfunc_a+nfunc_c) );
      dmat = pkv_GetScratchMem ( NCONSTR*sizeof(vector3d) );
      if ( !bfder || !cmat || !dmat )
        goto failure;

      memset ( cmat, 0, NCONSTR*nfunc_a*sizeof(double) );
      if ( !g2h_GetBPDerivativesd ( domain, 0, bfder ) )
        goto failure;
      pkv_Selectd ( nfunc_a, 1, 5, 1, bfder, cmat );
      for ( i = 0, k = 1;  i < hole_k;  i += 2 ) {
        if ( !g2h_GetBPDerivativesd ( domain, i, bfder ) )
          goto failure;
        if ( i == 0 || i == 2 )
          for ( j = 2;  j <= 4;  j++, k++ )
            pkv_Selectd ( nfunc_a, 1, 5, 1, &bfder[j], &cmat[k*nfunc_a] );
        else
          for ( j = 3;  j <= 4;  j++, k++ )
            pkv_Selectd ( nfunc_a, 1, 5, 1, &bfder[j], &cmat[k*nfunc_a] );
      }
      if ( !g2h_SetConstraintMatrixd ( domain, k, cmat ) )
        goto failure;
      memset ( dmat, 0, NCONSTR*sizeof(vector3d) );
      dmat[0] = centp;
      if ( !g2h_FillHoleConstrd ( domain, 3, (double*)Bpt,
                                  k, &dmat[0].x, NULL, NULL, OutPatch ) )
        goto failure;


      memset ( cmat, 0, NCONSTR*(nfunc_a+nfunc_c)*sizeof(double) );
      pkv_Selectd ( nfunc_a, 1, 5, 1, bfder, &cmat[nfunc_c] );
      for ( i = 0, k = 1;  i < hole_k;  i += 2 ) {
        if ( !g2h_GetBPDerivativesd ( domain, i, bfder ) )
          goto failure;
        if ( i == 0 || i == 2 )
          for ( j = 2;  j <= 4;  j++, k++ )
            pkv_Selectd ( nfunc_a, 1, 5, 1, &bfder[j],
                          &cmat[k*(nfunc_c+nfunc_a)+nfunc_c] );
        else
          for ( j = 3;  j <= 4;  j++, k++ )
            pkv_Selectd ( nfunc_a, 1, 5, 1, &bfder[j],
                          &cmat[k*(nfunc_c+nfunc_a)+nfunc_c] );
      }
      if ( !g2h_SetExtConstraintMatrixd ( domain, k, cmat ) )
        goto failure;
      memset ( dmat, 0, NCONSTR*sizeof(vector3d) );
      dmat[0] = centp;
      if ( !g2h_ExtFillHoleConstrd ( domain, 3, (double*)Bpt,
                                     k, &dmat[0].x, NULL, NULL, OutPatch ) )
        goto failure;

      times ( &stop );
      MakePictures ( fn1 );
      time = (double)(stop.tms_utime-start.tms_utime)/
             (double)(sysconf(_SC_CLK_TCK));
      printf ( "time = %8.3f\n", time );

    }
    gh_DestroyDomaind ( domain );
  }
  return;

failure:
  exit ( 1 );
} /*MakeTest1*/

void MakeTest2 ()
{
#define NCONSTR 11
  struct   tms start, stop;
  double   time;
  double   *bfder, *cmat, *dmat;
  int      nfunc_a, nfunc_c, i, j, k;
  point3d  centp;
  vector3d nvect;

  ReadDatFile ( infn, &hole_k, knots, Dompt, Bpt, NULL );
  SetupConstraints ( &centp, &nvect );

  times ( &start );
  if ( (domain = gh_CreateDomaind ( HOLE_K, knots, Dompt )) ) {
    if ( g2h_ComputeBasisd ( domain ) ) {
      nfinalp = 0;
                 /* setup the constraints */
      nfunc_a = g2h_V0SpaceDimd ( domain );
      nfunc_c = HOLE_K*DBDIM;
      bfder = pkv_GetScratchMemd ( nfunc_a*5 );
      cmat = pkv_GetScratchMemd ( 3*NCONSTR*(nfunc_a+nfunc_c) );
      dmat = pkv_GetScratchMemd ( NCONSTR );
      if ( !bfder || !cmat || !dmat )
        goto failure;

      memset ( cmat, 0, 3*NCONSTR*nfunc_a*sizeof(double) );
      if ( !g2h_GetBPDerivativesd ( domain, 0, bfder ) )
        goto failure;
      pkv_Selectd ( nfunc_a, 1, 5, 1, bfder, &cmat[2*nfunc_a] );
      for ( i = 0, k = 1;  i < hole_k;  i += 2 ) {
        if ( !g2h_GetBPDerivativesd ( domain, i, bfder ) )
          goto failure;
        if ( i == 0 || i == 2 )
          for ( j = 2;  j <= 4;  j++, k++ )
            pkv_Selectd ( nfunc_a, 1, 5, 1, &bfder[j], &cmat[(3*k+2)*nfunc_a] );
        else
          for ( j = 3;  j <= 4;  j++, k++ )
            pkv_Selectd ( nfunc_a, 1, 5, 1, &bfder[j], &cmat[(3*k+2)*nfunc_a] );
      }

      if ( !g2h_SetAltConstraintMatrixd ( domain, 3, k, cmat ) )
        goto failure;
      memset ( dmat, 0, NCONSTR*sizeof(double) );
      dmat[0] = DotProduct3d ( &centp, &nvect );
      if ( !g2h_FillHoleAltConstrd ( domain, 3, (double*)Bpt,
                                     k, dmat, NULL, NULL, OutPatch ) )
        goto failure;


      memset ( cmat, 0, 3*NCONSTR*(nfunc_a+nfunc_c)*sizeof(double) );
      pkv_Selectd ( nfunc_a, 1, 5, 1, bfder, &cmat[2*(nfunc_c+nfunc_a)+nfunc_c] );
      for ( i = 0, k = 1;  i < hole_k;  i += 2 ) {
        if ( !g2h_GetBPDerivativesd ( domain, i, bfder ) )
          goto failure;
        if ( i == 0 || i == 2 )
          for ( j = 2;  j <= 4;  j++, k++ )
            pkv_Selectd ( nfunc_a, 1, 5, 1, &bfder[j],
                          &cmat[(3*k+2)*(nfunc_c+nfunc_a)+nfunc_c] );
        else
          for ( j = 3;  j <= 4;  j++, k++ )
            pkv_Selectd ( nfunc_a, 1, 5, 1, &bfder[j],
                          &cmat[(3*k+2)*(nfunc_c+nfunc_a)+nfunc_c] );
      }
      if ( !g2h_SetExtAltConstraintMatrixd ( domain, 3, k, cmat ) )
        goto failure;
      memset ( dmat, 0, NCONSTR*sizeof(double) );
      dmat[0] = DotProduct3d ( &centp, &nvect );
      if ( !g2h_ExtFillHoleAltConstrd ( domain, 3, (double*)Bpt,
                                        k, dmat, NULL, NULL, OutPatch ) )
        goto failure;

      times ( &stop );
/*
      fval = g2h_FunctionalValued ( domain, 3, (double*)Bpt, 0, NULL );
      cfval = g2h_FunctionalValued ( domain, 3, (double*)Bpt, 2, &cnt[0].x );
      efval = g2h_ExtFunctionalValued ( domain, 3, (double*)Bpt, 0, NULL );  
      ecfval = g2h_ExtFunctionalValued ( domain, 3, (double*)Bpt, 2, &cnt[0].x );
      printf ( "Functional values = %f, %f, %f, %f\n",
               fval, cfval, efval, ecfval );
*/
      MakePictures ( fn2 );
      time = (double)(stop.tms_utime-start.tms_utime)/
             (double)(sysconf(_SC_CLK_TCK));
      printf ( "time = %8.3f\n", time );

    }
    gh_DestroyDomaind ( domain );
  }
  return;

failure:
  exit ( 1 );
} /*MakeTest2*/

int main ()
{
  pkv_InitScratchMem ( 3145728 );
  MakeTest1 ();
  MakeTest2 ();
  printf ( "Scratch memory used: %d bytes\n", pkv_MaxScratchTaken() );
  pkv_DestroyScratchMem ();
  exit ( 0 );
} /*main*/

