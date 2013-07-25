
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

char fn1[]  = "g1hnonlind.ps";

CameraRecd   CPos;

GHoleDomaind *domain;

int nfinalp;
point3d FinalCPoints[2*HOLE_K][(G1H_FINALDEG+1)*(G1H_FINALDEG+1)];

double iparams[3] = {0.5,0.0,0.5};

vector3d nlnvect;

/* ///////////////////////////////////////////////////////////////////////// */
void SetupCamera ( double w, double h, double x, double y )
{
  CameraInitFramed ( &CPos, false, true, w, h, x, y, 1.0, 4 );
  CameraInitPosd ( &CPos );
  SetPoint3d ( &CPos.position, 20.402914, -7.444557, -26.986362 );
  CPos.psi = -0.173305;  CPos.theta = 0.677664;  CPos.phi = -1.920662;
  CPos.vd.persp.f = 5.0;
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

void DrawFinalPatches ( int first )
{
  int i;

  ps_Set_Gray ( 0.25 );
  for ( i = first; i < first+HOLE_K; i++ )
    DrawBezPatch3d ( G1H_FINALDEG, G1H_FINALDEG, FinalCPoints[i] );
} /*DrawFinalPatches*/

void DrawNLNV ()
{
  point3d p0, p1;

  SetPoint3d ( &p1, 0.0, 0.0, 0.0 );
  CameraProjectPoint3d ( &CPos, &p1, &p0 );
  CameraProjectPoint3d ( &CPos, &nlnvect, &p1 );
  ps_Set_Gray ( 0.0 );
  psl_SetLine ( p0.x, p0.y, p1.x, p1.y, 0.0, 1.0 );
  psl_Draw ( 0.0, 0.9, 6.0 );
  psl_Arrow ( 1.0, true );
} /*DrawNLNV*/

void MakePicture1 ()
{
  ps_OpenFile ( fn1, 600 );
  ps_Write_Command ( "1 setlinecap" );

  SetupCamera ( 2400, 1800, 0.0, 0.0 );
  DrawBezPatches ();
  DrawFinalPatches ( 0 );

  SetupCamera ( 2400, 1800, 2400.0, 0.0 );
  DrawBezPatches ();
  DrawFinalPatches ( HOLE_K );

  ps_CloseFile ();
  printf ( "%s\n", fn1 );
} /*MakePicture1*/

void MakePictures ()
{
  MakePicture1 ();
} /*MakePictures*/

/* ///////////////////////////////////////////////////////////////////////// */
double amax;

double CRadius ( double a, double scf )
{
  if ( a < 0.0 )
    ps_Set_RGB ( 1.0, 0.0, 0.0 );
  else
    ps_Set_Gray ( 0.0 );
  return 0.66*scf*pow( fabs(a)/amax, 0.2 );
} /*CRadius*/

void DrawSymBlock ( int size, double scf, double posx, double posy, double *coeff )
{
  int   i, j, k;
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
  int   i, j, k;
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
  int   i, j, k;
  double r;

  for ( i = k = 0;  i < nrows;  i++ )
    for ( j = 0;  j < ncols;  j++, k++ ) {
      r = CRadius ( coeff[k], scf );
      if ( r > 0.05 )
        ps_Fill_Circle ( posx + i*scf + 0.5*scf, posy - j*scf - 0.5*scf, r );
    }
} /*DrawFullBlockT*/

void DrawMatrix ( int k, int r, int s, double *Aii, double *Aki, double *Akk )
{
  int    i, n;
  double scf;

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

void ExtCondNumber ( int k, int r, int s, double *Aii )
{
#define ITER 50
  void  *sp;
  int   n, i;
  double *Lii, *x, *y;
  double ca, pr, yl, lmin, lmax;

  sp = pkv_GetScratchMemTop ();
  n = k*r+s;
  Lii = pkv_GetScratchMemd ( pkn_Block1ArraySize ( k, r, s ));
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
      ca = pr/yl;
      pkn_MultMatrixNumd ( 1, n, 0, y, 1.0/yl, 0, x );
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

void DrawExtMatrix ( int k, int r, int s, double *Aii, double *Aki, double *Akk )
{
  ps_OpenFile ( "nlextmatrix.ps", 600 );
  DrawMatrix ( k, r, s, Aii, Aki, Akk );
  ps_CloseFile ();
  exit ( 0 );
} /*DrawExtMatrix*/

/* ///////////////////////////////////////////////////////////////////////// */
void OutPatch ( int n, int m, const point3d *cp, void *usrptr )
{
  if ( n != G1H_FINALDEG || m != G1H_FINALDEG || nfinalp >= 2*HOLE_K )
    exit ( 1 );
  memcpy ( FinalCPoints[nfinalp], cp,
          (G1H_FINALDEG+1)*(G1H_FINALDEG+1)*sizeof(point3d) );
  nfinalp ++;
} /*OutPatch*/

/* ///////////////////////////////////////////////////////////////////////// */
void OutputNLNV ( vector3d *nv )
{
  nlnvect = *nv;
} /*OutputNLNV*/

void SetupCameraG ( double w, double h, double x, double y )
{
  vector3d v;

  CameraInitFramed ( &CPos, false, true, w, h, x, y, 1.0, 4 );
  CameraInitPosd ( &CPos );
  SetVector3d ( &v, 0.0, 0.0, -46.0 );
  CameraMoveGd ( &CPos, &v );
  CameraRotXGd ( &CPos, -0.65*PI );
  CameraRotZGd ( &CPos, 0.05*PI );
  CameraZoomd ( &CPos, 10.0 );
} /*SetupCameraG*/

void DrawFGraph ( int k, int n, point2d *xy, double *z )
{
/*#define SF 0.03125*/
/*#define SF 0.0625*/
#define SF 0.025
  int i, j, l;
  point3d p, q, r;

  ps_OpenFile ( "grnl.ps", 600 );
  SetupCameraG ( 1800, 1350, 100, 100 );
  SetPoint3d ( &p, 0.0, 0.0, 0.0 );
  CameraProjectPoint3d ( &CPos, &p, &r );
  SetPoint3d ( &p, 1.0, 0.0, 0.0 );
  CameraProjectPoint3d ( &CPos, &p, &q );
  psl_SetLine ( r.x, r.y, q.x, q.y, 0.0, 1.0  );
  ps_Set_RGB ( 1.0, 0.0, 0.0 );
  psl_Draw ( 0.0, 0.95, 4.0 );
  psl_Arrow ( 1.0, true );
  ps_Set_RGB ( 0.0, 1.0, 0.0 );
  SetPoint3d ( &p, 0.0, 1.0, 0.0 );
  CameraProjectPoint3d ( &CPos, &p, &q );
  psl_SetLine ( r.x, r.y, q.x, q.y, 0.0, 1.0  );
  psl_Draw ( 0.0, 0.95, 4.0 );
  psl_Arrow ( 1.0, true );

/*
  ps_Set_Gray ( 0.68 );
  for ( i = 0; i < n; i++ ) {
    SetPoint3d ( &p, xy[i].x, xy[i].y, 0.0 );
    CameraProjectPoint3d ( &CPos, &p, &q );
    ps_Fill_Circle ( q.x, q.y, 3.0 );
  }
*/
  ps_Set_Line_Width ( 2.0 );
  for ( i = 0; i < k; i++ ) {
    ps_Set_Gray ( 0.3+0.1*i );
    for ( j = 0; j < n; j++ ) {
      SetPoint3d ( &p, xy[i*n*n+j*n].x, xy[i*n*n+j*n].y, z[i*n*n+j*n]*SF );
      CameraProjectPoint3d ( &CPos, &p, &q );
      for ( l = 1; l < n; l++ ) {
        r = q;
        SetPoint3d ( &p, xy[i*n*n+j*n+l].x, xy[i*n*n+j*n+l].y, z[i*n*n+j*n+l]*SF );
        CameraProjectPoint3d ( &CPos, &p, &q );
        ps_Draw_Line ( r.x, r.y, q.x, q.y );
      }
    }
    for ( j = 0; j < n; j++ ) {
      SetPoint3d ( &p, xy[i*n*n+j].x, xy[i*n*n+j].y, z[i*n*n+j]*SF );
      CameraProjectPoint3d ( &CPos, &p, &q );
      for ( l = 1; l < n; l++ ) {
        r = q;
        SetPoint3d ( &p, xy[i*n*n+l*n+j].x, xy[i*n*n+l*n+j].y, z[i*n*n+l*n+j]*SF );
        CameraProjectPoint3d ( &CPos, &p, &q );
        ps_Draw_Line ( r.x, r.y, q.x, q.y );
      }
    }
  }
  ps_CloseFile ();
  printf ( "grnl.ps\n" );
} /*DrawFGraph*/

/* ///////////////////////////////////////////////////////////////////////// */
int main ()
{
  struct tms start, stop;
  double time;
  double funcval[3];
  double param[2] = {0.0,0.0};
  vector3d *acoeff;

  pkv_InitScratchMem ( 2097152 );
  acoeff = (vector3d*)pkv_GetScratchMem ( 158*sizeof(vector3d) );
  if ( !acoeff )
    goto failure; 
  InitKnots ( HOLE_K );
  InitDomain ( HOLE_K, param );
  times ( &start );
  if ( (domain = gh_CreateDomaind ( HOLE_K, knots, Dompt )) ) {
    if ( g1h_ComputeBasisd ( domain ) ) {
      InitHole ( HOLE_K, iparams );
      nfinalp = 0;
      if ( !g1h_NLFillHoled ( domain, Bpt, (double*)acoeff, NULL, OutPatch ) )
        goto failure;
      if ( !g1h_NLFunctionalValued ( domain, Bpt, acoeff, &funcval[1] ) )
        goto failure;
      printf ( "E_a = %f, E_b = %f\n", funcval[1], funcval[2] );
      if ( !g1h_NLExtFillHoled ( domain, Bpt, (double*)acoeff, NULL, OutPatch ) )
        goto failure;
      if ( !g1h_NLExtFunctionalValued ( domain, Bpt, acoeff, &funcval[1] ) )
        goto failure;
      printf ( "E_a = %f, E_b = %f\n", funcval[1], funcval[2] );
      times ( &stop );
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

