
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

#include "eg1holef.h"
#include "datagenf.h"
#include "testgraphf.h"
#include "drawitf.h"

#define HOLE_K 6

char fn1[]  = "g1hq2altcf.ps";

CameraRecf   CPos;

GHoleDomainf *domain;

int nfinalp;
point3f FinalCPoints[2*HOLE_K][121];

float iparams[3] = {0.5,0.0,0.5};

#define DBDIM 4

/* ///////////////////////////////////////////////////////////////////////// */
void SetupCamera ( float w, float h, float x, float y )
{
  CameraInitFramef ( &CPos, false, true, w, h, x, y, 1.0, 4 );
  CameraInitPosf ( &CPos );
  SetPoint3f ( &CPos.position, 20.402914, -7.444557, -26.986362 );
  CPos.psi = -0.173305;  CPos.theta = 0.677664;  CPos.phi = -1.920662;
  CPos.vd.persp.f = 4.0;
  CameraSetMappingf ( &CPos );
} /*SetupCamera*/

void DrawPatchNet3f ( int n, int m, const point3f *cp )
{
  void    *sp;
  int     i, j;
  point3f q;
  point2f *r;

  sp = pkv_GetScratchMemTop ();
  r = pkv_GetScratchMem ( (n+1)*(m+1)*sizeof(point2f) );
  if ( !r )
    exit ( 1 );
  for ( i = 0; i < (n+1)*(m+1); i++ ) {
    CameraProjectPoint3f ( &CPos, &cp[i], &q );
    SetPoint2f ( &r[i], q.x, q.y );
  }
  ps_Set_Line_Width ( 1.0 );
  for ( i = 0; i <= n; i++ )
    ps_Draw_Polyline2f ( m+1, &r[i*(m+1)] );
  for ( j = 0; j <= m; j++ )
    for ( i = 0; i < n; i++ )
      ps_Draw_Line ( r[i*(m+1)+j].x, r[i*(m+1)+j].y,
                     r[(i+1)*(m+1)+j].x, r[(i+1)*(m+1)+j].y );
  for ( i = 0; i < (n+1)*(m+1); i++ )
    ps_Mark_Circle ( r[i].x, r[i].y );

  pkv_SetScratchMemTop ( sp );
} /*DrawPatchNet3f*/

void ShowPoint ( point3f *p )
{
  point3f q;

  CameraProjectPoint3f ( &CPos, p, &q );
  ps_Mark_Circle ( q.x, q.y );
} /*ShowPoint*/

/* ///////////////////////////////////////////////////////////////////////// */
void GetHoleSurrndPatch ( int i, int j, point3f *dombcp )
{
  int     k;
  int     *ind;
  point3f *q;
  void    *sp;
  float   *ukn, *vkn;

  sp  = pkv_GetScratchMemTop ();
  ind = pkv_GetScratchMem ( 16*sizeof(int) );
  q   = pkv_GetScratchMem ( 16*sizeof(point3f) );
  if ( !ind || !q )
    exit ( 1 );

  gh_GetBspInd ( HOLE_K, i, j, ind );
  for ( k = 0; k < 16; k++ )
    q[k] = Bpt[ind[k]];

  ukn = &knots[11*((i+HOLE_K-1) % HOLE_K)+3];
  vkn = &knots[11*i+j];
  mbs_BSPatchToBezf ( 3, 3, 7, ukn, 3, 7, vkn, 12, (float*)q,
                      NULL, NULL, NULL, NULL, NULL, NULL, 12, (float*)dombcp );

  pkv_SetScratchMemTop ( sp );
} /*GetHoleSurrndPatch*/

void DrawBezPatches ( void )
{
  int     i, j;
  point3f cp[16];

  ps_Set_Gray ( 0.72 );
  for ( i = 0; i < HOLE_K; i++ )
    for ( j = 0; j < 3; j++ ) {
      GetHoleSurrndPatch ( i, j, cp );
      DrawBezPatch3f ( 3, 3, cp, 8, 8 );
    }
} /*DrawBezPatches*/

void DrawFinalPatches ( int p0, int p1 )
{
  int i;

  ps_Set_Gray ( 0.25 );
  for ( i = p0; i < p1; i++ )
    DrawBezPatch3f ( G1H_FINALDEG, G1H_FINALDEG, FinalCPoints[i], 8, 8 );
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
void OutPatch ( int n, int m, const float *cp, void *usrptr )
{
  if ( n != G1H_FINALDEG || m != G1H_FINALDEG || nfinalp >= 2*HOLE_K )
    exit ( 1 );
  memcpy ( FinalCPoints[nfinalp], cp, 121*sizeof(point3f) );
  nfinalp ++;
} /*OutPatch*/

int main ( void )
{
  struct tms start, stop;
  float  time;
  float  param[2] = {0.0,0.0};
  float  *bfder, *cmat;
  int    nfunc_a, nfunc_c;
  float  fval, efval, cnstr[2], *acoeff, *bcoeff;

  pkv_InitScratchMem ( 2097152 );
  acoeff = pkv_GetScratchMemf ( 90 );
  bcoeff = pkv_GetScratchMemf ( 474 );
  if ( !acoeff || !bcoeff )
    goto failure;
  InitKnots ( HOLE_K );
  InitDomain ( HOLE_K, param );
  times ( &start );
  if ( (domain = gh_CreateDomainf ( HOLE_K, knots, Dompt )) ) {
    if ( g1h_ComputeBasisf ( domain ) ) {
      nfinalp = 0;
      InitHole ( HOLE_K, iparams );

                 /* setup the constraints */
      nfunc_a = g1h_V0SpaceDimf ( domain );
      nfunc_c = HOLE_K*DBDIM;
      bfder = pkv_GetScratchMemf ( nfunc_a*HOLE_K*3 );
      if ( !bfder )
        goto failure;
      if ( !g1h_GetBPDerivativesf ( domain, 0, bfder ) )
        goto failure;
      cmat = pkv_GetScratchMemf ( 6*(nfunc_a+nfunc_c) );
      if ( !cmat )
        goto failure;
      memset ( cmat, 0, 6*(nfunc_a+nfunc_c)*sizeof(float) );
      pkv_Selectf ( nfunc_a, 1, 3, 1, bfder, &cmat[2*nfunc_a] );
      pkv_Selectf ( nfunc_a, 1, 3, 1, &bfder[1], &cmat[(2+3)*nfunc_a] );
      if ( !g1h_SetAltConstraintMatrixf ( domain, 3, 2, cmat ) )
        goto failure;
      cnstr[0] = -0.55;
      cnstr[1] = -0.0;
      if ( !g1h_Q2FillHoleAltConstrf ( domain, 3, (float*)Bpt,
                                       2, cnstr, acoeff, NULL, OutPatch ) )
        goto failure;
      memset ( cmat, 0, 6*(nfunc_a+nfunc_c)*sizeof(float) );
      pkv_Selectf ( nfunc_a, 1, 3, 1, bfder,
                    &cmat[2*(nfunc_c+nfunc_a)+nfunc_c] );
      pkv_Selectf ( nfunc_a, 1, 3, 1, &bfder[1],
                    &cmat[(2+3)*(nfunc_a+nfunc_c)+nfunc_c] );
      if ( !g1h_SetExtAltConstraintMatrixf ( domain, 3, 2, cmat ) )
        goto failure;
      if ( !g1h_Q2ExtFillHoleAltConstrf ( domain, 3, (float*)Bpt,
                                          2, cnstr, bcoeff, NULL, OutPatch ) )
        goto failure;

      times ( &stop );
      fval = g1h_FunctionalValuef ( domain, 3, (float*)Bpt, acoeff );
      efval = g1h_ExtFunctionalValuef ( domain, 3, (float*)Bpt, bcoeff );
      printf ( "Functional values = %f, %f\n",
               fval, efval );
      MakePictures ();
      time = (float)(stop.tms_utime-start.tms_utime)/
             (float)(sysconf(_SC_CLK_TCK));
      printf ( "time = %8.3f\n", time );
    }
    gh_DestroyDomainf ( domain );
  }

  printf ( "Scratch memory used: %d bytes\n", pkv_MaxScratchTaken() );
  pkv_DestroyScratchMem ();
  exit ( 0 );

failure:
  exit ( 1 );
} /*main*/

