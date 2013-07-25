
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

#include "eg2holef.h"
#include "datagenf.h"

/* number of sides of the hole - must be 3, 5, 6 or 8 */
#define HOLE_K 5

char fn1[]  = "testg2f.ps";

CameraRecf   CPos;
GHoleDomainf *domain;
int nfinalp;
point3f FinalCPoints[HOLE_K][121];


static void OutPatch ( int n, int m, const float *cp, void *usrptr )
{
  if ( n != G2H_FINALDEG || m != G2H_FINALDEG || nfinalp >= HOLE_K )
    exit ( 1 );
  memcpy ( FinalCPoints[nfinalp], cp, 121*sizeof(point3f) );
  nfinalp ++;
} /*OutPatch*/

/* ///////////////////////////////////////////////////////////////////////// */
static void SetupCamera ( short w, short h, short x, short y )
{
  CameraInitFramef ( &CPos, false, true, w, h, x, y, 1.0, 4 );
  CameraInitPosf ( &CPos );
  SetPoint3f ( &CPos.position, 20.402914, -7.444557, -26.986362 );
  CPos.psi = -0.173305;  CPos.theta = 0.677664;  CPos.phi = -1.920662;
  CPos.vd.persp.f = 4.0;
  CameraSetMappingf ( &CPos );
} /*SetupCamera*/

static void DrawBezCurve3f ( int n, const point3f *cp )
{
#define DD 32
  void    *sp;
  int     i;
  point3f p, q;
  point2f *cc;

  sp = pkv_GetScratchMemTop ();
  cc = (point2f*)pkv_GetScratchMem ( (DD+1)*sizeof(point2f) );
  if ( !cc )
    exit ( 1 );

  for ( i = 0; i <= DD; i++ ) {
    mbs_BCHornerC3f ( n, cp, (float)i/(float)DD, &p );
    CameraProjectPoint3f ( &CPos, &p, &q );
    SetPoint2f ( &cc[i], q.x, q.y );
  }
  ps_Draw_Polyline2f ( DD+1, cc );

  pkv_SetScratchMemTop ( sp );
#undef DD
} /*DrawBezCurve3f*/

static void SetLW ( int i, int i0, int i1 )
{
    if ( i == i0 || i == i1 )
      ps_Set_Line_Width ( 4.0 );
    else
      ps_Set_Line_Width ( 2.0 );
} /*SetLW*/

static void DrawBezPatch3f ( int n, int m, const point3f *cp )
{
#define D 8
  void    *sp;
  int     i;
  point3f *cc;

  sp = pkv_GetScratchMemTop ();
  cc = (point3f*)pkv_GetScratchMem ( (max(n,m)+1)*sizeof(point3f) );
  if ( !cc )
    exit ( 1 );

  for ( i = 0; i <= D; i++ ) {
    mbs_multiBCHornerf ( n, 1, 3*(m+1), 0, (float*)cp, (float)i/(float)D,
                         (float*)cc );
    SetLW ( i, 0, D );
    DrawBezCurve3f ( m, cc );
  }
  for ( i = 0; i <= D; i++ ) {
    mbs_multiBCHornerf ( m, n+1, 3, 3*(m+1), (float*)cp, (float)i/(float)D,
                         (float*)cc );
    SetLW ( i, 0, D );
    DrawBezCurve3f ( n, cc );
  }

  pkv_SetScratchMemTop ( sp );
#undef D
} /*DrawBezPatch3f*/

/* ///////////////////////////////////////////////////////////////////////// */
static void GetHoleSurrndPatch ( int i, int j, point3f *dombcp )
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

static void DrawBezPatches ( void )
{
  int     i, j;
  point3f cp[16];

  ps_Set_Gray ( 0.72 );
  for ( i = 0; i < HOLE_K; i++ )
    for ( j = 0; j < 3; j++ ) {
      GetHoleSurrndPatch ( i, j, cp );
      DrawBezPatch3f ( 3, 3, cp );
    }
} /*DrawBezPatches*/

static void DrawFinalPatches ( void )
{
  int i;

  ps_Set_Gray ( 0.25 );
  for ( i = 0; i < nfinalp; i++ )
    DrawBezPatch3f ( G2H_FINALDEG, G2H_FINALDEG, FinalCPoints[i] );
} /*DrawFinalPatches*/

static void MakePicture ( void )
{
  ps_OpenFile ( fn1, 600 );
  ps_Write_Command ( "1 setlinecap" );
  SetupCamera ( 1800, 1460, 0.0, 0.0 );

  DrawBezPatches ();
  DrawFinalPatches ();

  ps_CloseFile ();
  printf ( "%s\n", fn1 );
} /*MakePicture*/

int main ( void )
{
  float dparams[2] = {0.0,0.0};
  float sparams[3] = {0.5,0.0,0.5};

    /* Initialize the scratch memory pool. 2MB should be enough even for */
    /* double precision */
  pkv_InitScratchMem ( 2097152 );

    /* Prepare data. Knots are equidistant; this may be changed, but */
    /* the knot sequences must satisfy the condition (12) */
  InitKnots ( HOLE_K );
    /* Prepare the domain control points. The parameters should be */
    /* between 0 and 1 */
  InitDomain ( HOLE_K, dparams );
    /* Prepare the surface control points. They may be arbitrarily changed */
  InitHole ( HOLE_K, sparams );

    /* Now the proper construction. The first part is creating an object */
    /* which will represent the basis and contain all data necessary to */
    /* fill the hole. */
  if ( (domain = gh_CreateDomainf ( HOLE_K, knots, Dompt )) ) {
    if ( g2h_ComputeBasisf ( domain ) ) {
      nfinalp = 0;
      if ( !g2h_FillHolef ( domain, 3, (float*)Bpt, NULL, NULL, OutPatch ) )
        goto failure;
      MakePicture ();
    }
    gh_DestroyDomainf ( domain );
  }

  printf ( "Scratch memory used: %d bytes\n", (int)pkv_MaxScratchTaken() );
  pkv_DestroyScratchMem ();
  exit ( 0 );

failure:
  exit ( 1 );
} /*main*/

