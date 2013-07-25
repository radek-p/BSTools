
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
#include "drawitd.h"

#define HOLE_K 6

char fn1[]  = "g1hq2altcd.ps";

CameraRecd   CPos;

GHoleDomaind *domain;

int nfinalp;
point3d FinalCPoints[2*HOLE_K][121];

double iparams[3] = {0.5,0.0,0.5};

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

void DrawBezPatches ( void )
{
  int     i, j;
  point3d cp[16];

  ps_Set_Gray ( 0.72 );
  for ( i = 0; i < HOLE_K; i++ )
    for ( j = 0; j < 3; j++ ) {
      GetHoleSurrndPatch ( i, j, cp );
      DrawBezPatch3d ( 3, 3, cp, 8, 8 );
    }
} /*DrawBezPatches*/

void DrawFinalPatches ( int p0, int p1 )
{
  int i;

  ps_Set_Gray ( 0.25 );
  for ( i = p0; i < p1; i++ )
    DrawBezPatch3d ( G1H_FINALDEG, G1H_FINALDEG, FinalCPoints[i], 8, 8 );
} /*DrawFinalPatches*/

void MakePicture1 ( void )
{
  ps_OpenFile ( fn1, 600 );
  ps_Write_Command ( "1 setlinecap" );
  SetupCamera ( 1800, 1460, 0.0, 0.0 );

  DrawBezPatches ();
  DrawFinalPatches ( 0, HOLE_K );
  printf ( "central point: x = %f, y = %f, z = %f\n", FinalCPoints[0][0].x,
           FinalCPoints[0][0].y, FinalCPoints[0][0].z );

  ps_GSave ();
  ps_Write_Command ( "1800 0 translate" );
  DrawBezPatches ();
  DrawFinalPatches ( HOLE_K, 2*HOLE_K );
  printf ( "central point: x = %f, y = %f, z = %f\n", FinalCPoints[HOLE_K][0].x,
           FinalCPoints[HOLE_K][0].y, FinalCPoints[HOLE_K][0].z );
  ps_GRestore ();

  ps_CloseFile ();
  printf ( "%s\n", fn1 );
} /*MakePicture1*/

void MakePictures ( void )
{
  MakePicture1 ();
} /*MakePictures*/

/* ///////////////////////////////////////////////////////////////////////// */
void OutPatch ( int n, int m, const double *cp, void *usrptr )
{
  if ( n != G1H_FINALDEG || m != G1H_FINALDEG || nfinalp >= 2*HOLE_K )
    exit ( 1 );
  memcpy ( FinalCPoints[nfinalp], cp, 121*sizeof(point3d) );
  nfinalp ++;
} /*OutPatch*/

int main ( void )
{
  struct tms start, stop;
  double  time;
  double  param[2] = {0.0,0.0};
  double  *bfder, *cmat;
  int    nfunc_a, nfunc_c;
  double  fval, efval, cnstr[2], *acoeff, *bcoeff;

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
      cmat = pkv_GetScratchMemd ( 6*(nfunc_a+nfunc_c) );
      if ( !cmat )
        goto failure;
      memset ( cmat, 0, 6*(nfunc_a+nfunc_c)*sizeof(double) );
      pkv_Selectd ( nfunc_a, 1, 3, 1, bfder, &cmat[2*nfunc_a] );
      pkv_Selectd ( nfunc_a, 1, 3, 1, &bfder[1], &cmat[(2+3)*nfunc_a] );
      if ( !g1h_SetAltConstraintMatrixd ( domain, 3, 2, cmat ) )
        goto failure;
      cnstr[0] = -0.55;
      cnstr[1] = -0.0;
      if ( !g1h_Q2FillHoleAltConstrd ( domain, 3, (double*)Bpt,
                                       2, cnstr, acoeff, NULL, OutPatch ) )
        goto failure;
      memset ( cmat, 0, 6*(nfunc_a+nfunc_c)*sizeof(double) );
      pkv_Selectd ( nfunc_a, 1, 3, 1, bfder,
                    &cmat[2*(nfunc_c+nfunc_a)+nfunc_c] );
      pkv_Selectd ( nfunc_a, 1, 3, 1, &bfder[1],
                    &cmat[(2+3)*(nfunc_a+nfunc_c)+nfunc_c] );
      if ( !g1h_SetExtAltConstraintMatrixd ( domain, 3, 2, cmat ) )
        goto failure;
      if ( !g1h_Q2ExtFillHoleAltConstrd ( domain, 3, (double*)Bpt,
                                          2, cnstr, bcoeff, NULL, OutPatch ) )
        goto failure;

      times ( &stop );
      fval = g1h_FunctionalValued ( domain, 3, (double*)Bpt, acoeff );
      efval = g1h_ExtFunctionalValued ( domain, 3, (double*)Bpt, bcoeff );
      printf ( "Functional values = %f, %f\n",
               fval, efval );
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

