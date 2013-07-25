
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

#define HOLE_K 6

char fn1[] = "g2holeextd.ps";
char fn2[] = "g2extmatrixd.ps";

CameraRecd   CPos;

GHoleDomaind *domain;

int nfinalp;
point3d FinalCPoints[HOLE_K][121];

double iparams[3] = {0.5,0.0,0.5};

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

/* ///////////////////////////////////////////////////////////////////////// */
void GetHoleSurrndPatch ( int i, int j, point3d *dombcp )
{
  int     k;
  int     *ind;
  point3d *q;
  void    *sp;
  double  *ukn, *vkn;

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

void DrawFinalPatches ()
{
  int i;

  ps_Set_Gray ( 0.25 );
  for ( i = 0; i < nfinalp; i++ )
    DrawBezPatch3d ( G2H_FINALDEG, G2H_FINALDEG, FinalCPoints[i] );
} /*DrawFinalPatches*/

void MakePicture1 ()
{
  ps_OpenFile ( fn1, 600 );
  ps_Write_Command ( "1 setlinecap" );
  SetupCamera ( 1800, 1460, 0.0, 0.0 );

  DrawBezPatches ();
  DrawFinalPatches ();

  ps_CloseFile ();
  printf ( "%s\n", fn1 );
} /*MakePicture1*/

/* ///////////////////////////////////////////////////////////////////////// */
double amax;

double CRadius ( double a, double scf )
{
  if ( a < 0.0 )
    ps_Set_Gray ( 0.68 );
  else
    ps_Set_Gray ( 0.0 );
  return 0.66*scf*pow( fabs(a)/amax, 0.2 );
} /*CRadius*/

void DrawSymBlock ( int size, double scf, double posx, double posy, double *coeff )
{
  int    i, j, k;
  double r;

  for ( i = k = 0;  i < size; i++ ) {
    for ( j = 0;  j < i;  j++, k++ ) {
      r = CRadius ( coeff[k], scf );
      if ( r > 0.05 ) {
        ps_Fill_Circle ( posx + i*scf + 0.5*scf, posy - j*scf - 0.5*scf, r );
        ps_Fill_Circle ( posx + j*scf + 0.5*scf, posy - i*scf - 0.5*scf, r );
      }
    }
    r = CRadius ( coeff[k], scf );
    if ( r > 0.05 )
      ps_Fill_Circle ( posx + i*scf + 0.5*scf, posy - i*scf - 0.5*scf, r );
    k++;
  }
} /*DrawSymBlock*/

void DrawFullBlock ( int nrows, int ncols, double scf,
                     double posx, double posy, double *coeff )
{
  int    i, j, k;
  double r;

  for ( i = k = 0;  i < nrows;  i++ )
    for ( j = 0;  j < ncols;  j++, k++ ) {
      r = CRadius ( coeff[k], scf );
      if ( r > 0.05 )
        ps_Fill_Circle ( posx + j*scf + 0.5*scf, posy - i*scf - 0.5*scf, r );
    }
} /*DrawFullBlock*/

void DrawFullBlockT ( int nrows, int ncols, double scf,
                      double posx, double posy, double *coeff )
{
  int    i, j, k;
  double r;

  for ( i = k = 0;  i < nrows;  i++ )
    for ( j = 0;  j < ncols;  j++, k++ ) {
      r = CRadius ( coeff[k], scf );
      if ( r > 0.05 )
        ps_Fill_Circle ( posx + i*scf + 0.5*scf, posy - j*scf - 0.5*scf, r );
    }
} /*DrawFullBlockT*/

void DrawMatrix ( int k, int r, int s, double *Aii, double *Bi )
{
  int    i, n;
  double scf;
  double *Aki, *Akk;

  Aki = &Aii[pkn_Block1FindBlockPos( k, r, s, k, 0 )];
  Akk = &Aii[pkn_Block1FindBlockPos( k, r, s, k, k )];
  n = s+k*r;
  printf ( "n = %d\n", n );
  scf = 25.0;

  amax = fabs(Aii[0]);
  for ( i = 1; i < k*r*(r+1)/2; i++ )
    amax = max ( amax, fabs(Aii[i]) );
  for ( i = 0; i < s*(s+1)/2; i++ )
    amax = max ( amax, fabs(Akk[i]) );
  for ( i = 0; i < k*r*s; i++ )
    amax = max ( amax, fabs(Aki[i]) );

  for ( i = 0; i < k; i++ ) {
    DrawSymBlock ( r, scf, i*r*scf, (n-i*r)*scf, &Aii[r*(r+1)/2] );
    DrawFullBlock ( s, r, scf, i*r*scf, s*scf, &Aki[i*r*s] );
    DrawFullBlockT ( s, r, scf, k*r*scf, (n-i*r)*scf, &Aki[i*r*s] );
  }
  DrawSymBlock ( s, scf, k*r*scf, s*scf, Akk );
} /*DrawMatrix*/

void ExtCondNumber ( int k, int r, int s, double *Aii, double *Bi )
{
#define ITER 50
  void   *sp;
  int    n, i;
  double *Lii, *x, *y;
  double ca, pr, yl, lmin, lmax;

  sp = pkv_GetScratchMemTop ();
  n = k*r+s;
  Lii = pkv_GetScratchMemd ( pkn_Block1ArraySize ( k, r, s ) );
  x = pkv_GetScratchMemd ( n );
  y = pkv_GetScratchMemd ( n );
  if ( Lii && x && y ) {
    memcpy ( Lii, Aii, pkn_Block1ArraySize ( k, r, s )*sizeof(double) );
    if ( !pkn_Block1CholeskyDecompMd ( k, r, s, Lii ) )
      exit ( 1 );

    memset ( x, 0, n*sizeof(double) );
    x[0] = 1.0;
    for ( i = 0; i < ITER; i++ ) {
      pkn_Block1SymMatrixMultd ( k, r, s, Aii, 1, 1, x, 1, y );
      pr = pkn_ScalarProductd ( n, x, y );
      yl = pkn_SecondNormd ( n, y );
      ca = pr/yl;pkn_MultMatrixNumd ( 1, n, 0, y, 1.0/yl, 0, x );
      if  (ca > 0.999999 )
        break;
    }
    lmax = yl;

    memset ( x, 0, n*sizeof(double) );
    x[0] = 1.0;
    for ( i = 0; i < ITER; i++ ) {
      memcpy ( y, x, n*sizeof(double) );
      pkn_Block1LowerTrMSolved ( k, r, s, Lii, 1, 1, y );
      pkn_Block1UpperTrMSolved ( k, r, s, Lii, 1, 1, y );
      pr = pkn_ScalarProductd ( n, x, y );
      yl = pkn_SecondNormd ( n, y );
      ca = pr/yl;pkn_MultMatrixNumd ( 1, n, 0, y, 1.0/yl, 0, x );
      if  (ca > 0.999999 )
        break;
    }
    lmin = 1.0/yl;
    printf ( "lambda_max = %f,  lambda_min = %f,  cond_2(A) = %f\n",
             lmax, lmin, lmax*yl );
  }

  pkv_SetScratchMemTop ( sp );
#undef ITER
} /*ExtCondNumber*/

void MakePicture2 ()
{
  ps_OpenFile ( fn2, 600 );
  g2h_DrawExtMatricesd ( domain, DrawMatrix );
  g2h_DrawExtMatricesd ( domain, ExtCondNumber );
  ps_CloseFile ();
  printf ( "%s\n", fn2 );
} /*MakePicture2*/

/* ///////////////////////////////////////////////////////////////////////// */
void MakePictures ()
{
  MakePicture1 ();
  MakePicture2 ();
} /*MakePictures*/

/* ///////////////////////////////////////////////////////////////////////// */
void OutPatch ( int n, int m, const double *cp, void *usrptr )
{
  if ( n != G2H_FINALDEG || m != G2H_FINALDEG || nfinalp >= HOLE_K )
    exit ( 1 );
  memcpy ( FinalCPoints[nfinalp], cp, 121*sizeof(point3d) );
  nfinalp ++;
} /*OutPatch*/

int main ()
{
  struct tms start, stop;
  double time;
  double param[2] = {0.0,0.0};
  double funcval, *acoeff;

  pkv_InitScratchMem ( 4194304 );
  acoeff = pkv_GetScratchMemd ( 474 );
  if ( !acoeff )
    goto failure;
  InitKnots ( HOLE_K );
  InitDomain ( HOLE_K, param );
  times ( &start );
  if ( (domain = gh_CreateDomaind ( HOLE_K, knots, Dompt )) ) {
    if ( g2h_ComputeBasisd ( domain ) ) {
      if ( g2h_ComputeExtFormMatrixd ( domain ) ) {
        nfinalp = 0;
        InitHole ( HOLE_K, iparams );
        if ( !g2h_ExtFillHoled ( domain, 3, (double*)Bpt, acoeff, NULL, OutPatch ) )
          goto failure;
        times ( &stop );
        funcval = g2h_ExtFunctionalValued ( domain, 3, (double*)Bpt, acoeff );       
        printf ( "Functional value = %f\n", funcval );
        MakePictures ();
        time = (double)(stop.tms_utime-start.tms_utime)/
               (double)(sysconf(_SC_CLK_TCK));
        printf ( "time = %8.3f\n", time );
      }
    }
    gh_DestroyDomaind ( domain );
  }

  printf ( "Scratch memory used: %d bytes\n", pkv_MaxScratchTaken () );
  pkv_DestroyScratchMem ();
  exit ( 0 );

failure:
  exit ( 1 );
} /*main*/

