
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
#include "drawitd.h"
#include "writdatd.h"

#define _DRAW_BASIS
#define _DRAW_LAPLACIAN

#define WIDTH  920.0
#define HEIGHT 690.0

#define QUAD_FACTOR 10  /* 10 */

#define NK 1
#define M1 1
#define M2 1

#define HOLE_K 5
#define FUNCA21

char fn[] = "eg2hnlspld.ps";
char datfn1[] = "eg2hnlspld.dat";

GHoleDomaind *domain;

int nfpatches = 0;
int fpdeg, fplkn;
double *fpknots;
point3d *fpcp[HOLE_K];

int nfbpatches = 0;
int bfpdeg;
point3d *bfpcp[2*HOLE_K];

/* ///////////////////////////////////////////////////////////////////////// */
void SetupCPos ( double w, double h, double x, double y )
{
  vector3d v;

  CameraInitFramed ( &CPos, false, true, w, h, x, y, 1.0, 4 );
  CameraInitPosd ( &CPos );

  SetVector3d ( &v, 0.0, 0.0, -60.0 );
  CameraMoveGd ( &CPos, &v );
  CameraRotXGd ( &CPos, 0.25*PI );
  CameraRotZGd ( &CPos, 0.2*PI );
  CameraZoomd ( &CPos, 12.0 );
/*
  SetVector3d ( &v, 0.0, 0.0, -46.0 );
  CameraMoveGd ( &CPos, &v );
  CameraRotXGd ( &CPos, -0.65*PI );
  CameraRotZGd ( &CPos, 0.05*PI ); 
  CameraZoomd ( &CPos, 11.0 );
*/
} /*SetupCPos*/

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

  ps_Set_Gray ( 0.5 );
  for ( i = 0; i < HOLE_K; i++ )
    for ( j = 0; j < 3; j++ ) {
      GetHoleSurrndPatch ( i, j, cp );
      DrawBezPatch3d ( 3, 3, cp, 6, 6 );
    }
} /*DrawBezPatches*/

void DrawFinalBezPatches ( void )
{
  int i;

  ps_Set_RGB ( 1.0, 0.0, 0.0 );
  for ( i = 0; i < nfpatches; i++ )
    DrawBezPatch3d ( fpdeg, fpdeg, bfpcp[i], 2*(NK+1), 2*(NK+1) );
} /*DrawFinalBezPatches*/

void DrawFinalPatches ( void )
{
  int i;

  ps_Set_Gray ( 0.25 );
  for ( i = 0; i < nfpatches; i++ )
    DrawBSPatch3d ( fpdeg, fplkn, fpknots, fpdeg, fplkn, fpknots,
                    fpcp[i], 4, 4 );
} /*DrawFinalPatches*/

void MakePicture1 ( void )
{
  ps_OpenFile ( fn, 600 );
  ps_Write_Command ( "1 setlinecap" );
  SetupCPos ( 2*WIDTH, 2*HEIGHT, 0.0, 0.0 );
  ps_Set_Line_Width ( 1.0 );
  ps_Draw_Rect ( CPos.width, CPos.height, 0, 0 );
  ps_Set_Clip_Rect ( CPos.width, CPos.height, 0, 0 );
  DrawBezPatches ();
  DrawFinalBezPatches ();
  DrawFinalPatches ();

  ps_CloseFile ();
  printf ( "%s\n", fn );
} /*MakePicture1*/

/* ///////////////////////////////////////////////////////////////////////// */
void MakePictures ( void )
{
  MakePicture1 ();
} /*MakePictures*/

/* ///////////////////////////////////////////////////////////////////////// */
void WriteData ( void )
{
  int i;
  vector3d ld = {0.2,-0.3,1.0};

  OpenDataFile ( datfn1 );
  SetupCPos ( WIDTH, HEIGHT, 0.0, 0.0 );
  WriteCamera ( &CPos );
  NormalizeVector3d ( &ld );
  WriteLightDir ( &ld );
  WriteKnots ( HOLE_K, (void*)knots );
  WriteDomCP ( 12*HOLE_K+1, Dompt );
  WriteSurfCP ( 12*HOLE_K+1, Bpt );
  WriteSurfPatches ( HOLE_K, (void*)knots, Bpt );
  for ( i = 0; i < nfpatches; i++ )
    WriteBSPatch (  fpdeg, fplkn, fpknots, fpdeg, fplkn, fpknots, fpcp[i] );
  CloseDataFile ();
  printf ( "%s\n", datfn1 );
/*
  OpenDataFile ( datfn2 );
  SetupCPos ( WIDTH, HEIGHT, 0.0, 0.0 );
  WriteCamera ( &CPos );
  WriteLightDir ( &ld );
  WriteKnots ( HOLE_K, (void*)knots );
  WriteDomCP ( 12*HOLE_K+1, Dompt );
  WriteSurfCP ( 12*HOLE_K+1, Bpt );
  WriteSurfPatches ( HOLE_K, (void*)knots, Bpt );
  for ( i = 0; i < HOLE_K; i++ )
    WriteBezPatch (  fpdeg, fpdeg, &bfpcp[HOLE_K+i]->x );
  CloseDataFile ();
  printf ( "%s\n", datfn2 );

  OpenDataFile ( datfn3 );
  SetupCPos ( WIDTH, HEIGHT, 0.0, 0.0 );
  WriteCamera ( &CPos );
  WriteLightDir ( &ld );
  WriteKnots ( HOLE_K, (void*)knots );
  WriteDomCP ( 12*HOLE_K+1, Dompt );
  WriteSurfCP ( 12*HOLE_K+1, Bpt );
  WriteSurfPatches ( HOLE_K, (void*)knots, Bpt );
  for ( i = 0; i < HOLE_K; i++ )
    WriteBezPatch (  fpdeg, fpdeg, &bfpcp[i]->x );
  CloseDataFile ();
  printf ( "%s\n", datfn3 );
*/
} /*WriteData*/

/* ///////////////////////////////////////////////////////////////////////// */
void OutputSplPatch ( int degu, int lknu, const double *knotsu,
                      int degv, int lknv, const double *knotsv,
                      const point3d *cp, void *usrptr )
{
  int s;

  if ( nfpatches == 0 ) {
    fpknots = malloc ( (lknu+1)*sizeof(double) );
    if ( !fpknots )
      exit ( 1 );
    memcpy ( fpknots, knotsu, (lknu+1)*sizeof(double) );
    fpdeg = degu;
    fplkn = lknu;
    printf ( "%d %d %d %d\n", degu, lknu, degv, lknv );
  }
  s = (lknu-degu)*(lknu-degu);
  fpcp[nfpatches] = malloc ( s*sizeof(point3d) );
  if ( !fpcp[nfpatches] )
    exit ( 1 );
  memcpy ( fpcp[nfpatches], cp, s*sizeof(point3d) );
  nfpatches ++;
} /*OutputSplPatch*/

void OutputBezPatch ( int degu, int degv, const double *cp, void *usrptr )
{
  bfpdeg = degu;
  bfpcp[nfbpatches] = malloc ( (degu+1)*(degu+1)*sizeof(point3d) );
  if ( !bfpcp[nfbpatches] )
    exit ( 1 );
  memcpy ( bfpcp[nfbpatches], cp, (degu+1)*(degu+1)*sizeof(point3d) );
  nfbpatches ++;
} /*OutputBezPatch*/

int main ( void )
{
  struct tms start, stop;
  double time;
  double param[2] = { 0.0, 0.0 };
  double iparams[3] = /*{ 0.0, 0.0, 0.0 }*/ { 0.7, 0.1, 1.0 };
  double *acoeff, fval;
  int nfunc_a, nfunc_b, nfunc_c, nfunc_d;

  setvbuf ( stdout, NULL, _IONBF, 0 );
  pkv_InitScratchMem ( 134217728 ); /* 128MB */
  acoeff = pkv_GetScratchMemd ( 262144 );
  if ( !acoeff )
    exit ( 1 );
  InitKnots ( HOLE_K );
  ModifyKnots ( HOLE_K );
  InitDomain ( HOLE_K, param );
  InitHole ( HOLE_K, iparams );
  times ( &start );
  if ( (domain = gh_CreateDomaind ( HOLE_K, knots, Dompt )) ) {
    if ( g2h_ComputeSplBasisd ( domain, NK, M1, M2 ) ) {
      g2h_DrawSplBasFuncNumd ( domain, &nfunc_a, &nfunc_b, &nfunc_c, &nfunc_d );
      printf ( "a=%d, b=%d, c=%d, d=%d\n", nfunc_a, nfunc_b, nfunc_c, nfunc_d );
      printf ( "    A    B    C    D\n" );
      printf ( "%5d%5d%5d%5d%5d\n", 0, nfunc_a, nfunc_a+nfunc_b,
               nfunc_a+nfunc_b+nfunc_c, nfunc_a+nfunc_b+nfunc_c+nfunc_d );
      printf ( "%5d     %5d%5d\n", nfunc_c+nfunc_d, 0, nfunc_c );
      nfbpatches = 0;
      if ( !g2h_FillHoled ( domain, 3, (double*)Bpt, acoeff, NULL, OutputBezPatch ) )
        exit ( 1 );
      fval = g2h_FunctionalValued ( domain, 3, (double*)Bpt, acoeff );
      printf ( "fval = %f\n", fval );
      if ( !g2h_ExtFillHoled ( domain, 3, (double*)Bpt, acoeff, NULL, OutputBezPatch ) )
        exit ( 1 );
      fval = g2h_ExtFunctionalValued ( domain, 3, (double*)Bpt, acoeff );
      printf ( "fval = %f\n", fval );
      nfpatches = 0;
      if ( g2h_NLSplFillHoled ( domain, Bpt, NULL, NULL, OutputSplPatch ) ) {
        times ( &stop );
        time = (double)(stop.tms_utime-start.tms_utime)/
               (double)(sysconf(_SC_CLK_TCK));
        printf ( "time = %8.3f\n", time );
        MakePictures ();
        WriteData ();
      }
    }
    gh_DestroyDomaind ( domain );
  }

  printf ( "Scratch memory used: %d bytes\a\n", pkv_MaxScratchTaken() );
  pkv_DestroyScratchMem ();
  exit ( 0 );
} /*main*/

