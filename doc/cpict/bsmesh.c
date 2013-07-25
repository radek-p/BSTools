
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
#include "bsmesh.h"
#include "bsfile.h"
#include "psout.h"

char ifn[] = "m1.bs";
char fn[] = "bsmesh.ps";

#define MAXNV     500
#define MAXNHE   2000
#define MAXNFAC   500

int         nv1, nhe1, nfac1;
int         *mvhei1, *mfhei1;
point2d     *mvp1;
BSMvertex   *mv1;
BSMhalfedge *mhe1;
BSMfacet    *mfac1;
boolean     m1ok = false;

CameraRecd CPos;

void InputMesh ( void *userdata,
                 const char *name, int degree,
                 int nv, const BSMvertex *mv, const int *mvhei,
                 const point4d *vc, int nhe, const BSMhalfedge *mhe,
                 int nfac, const BSMfacet *mfac, const int *mfhei,
                 int spdimen, boolean rational, byte *mkv )
{
  int i;

  nv1   = nv;
  nhe1  = nhe;
  nfac1 = nfac;
  mv1    = malloc ( nv*sizeof(BSMvertex) );
  mvhei1 = malloc ( nhe*sizeof(int) );
  mhe1   = malloc ( nhe*sizeof(BSMhalfedge) );
  mfac1  = malloc ( nfac*sizeof(BSMfacet) );
  mfhei1 = malloc ( nhe*sizeof(int) );
  mvp1   = malloc ( nv*sizeof(point2d) );
  if ( !mv1 || !mvhei1 || !mhe1 || !mfac1 || !mfhei1 || !mvp1 )
    exit ( 1 );
  memcpy ( mv1, mv, nv*sizeof(BSMvertex) );
  memcpy ( mvhei1, mvhei, nhe*sizeof(int) );
  memcpy ( mhe1, mhe, nhe*sizeof(BSMhalfedge) );
  memcpy ( mfac1, mfac, nfac*sizeof(BSMfacet) );
  memcpy ( mfhei1, mfhei, nhe*sizeof(int) );
  for ( i = 0; i < nv; i++ )
    Point4to2d ( &vc[i], &mvp1[i] );
  m1ok = true;
} /*InputMesh*/

void SetupCamera ( void )
{
  vector3d v;

  CameraInitFramed ( &CPos, true, false, 2000, 2000, 0, 0, 1.0, 4.0 );
  CameraInitPosd ( &CPos );
  CPos.vd.para.wdt = 1.6;
  CPos.vd.para.dim_case = 1;
  SetVector3d ( &v, 0.0, 0.25, 0.0 );
  CameraMoveGd ( &CPos, &v );
} /*SetupCamera*/

void DrawMesh ( int nv, BSMvertex *mv, int *mvhei, point2d *mcp,
                int nhe, BSMhalfedge *mhe,
                int nfac, BSMfacet *mfac, int *mfhei )
{
  void    *sp;
  int     i, j, fhe, d;
  point2d *mcq, fc, ev;
  char    s[60];

  sp = pkv_GetScratchMemTop ();
  mcq = pkv_GetScratchMem ( nv*sizeof(point3d) );
  if ( !mcq )
    exit ( 1 );
  for ( i = 0; i < nv; i++ )
    CameraProjectPoint2d ( &CPos, &mcp[i], &mcq[i] );
  ps_Set_Line_Width ( 1.0 );
  for ( i = 0; i < nhe; i++ )
    if ( mhe[i].otherhalf < 0 || mhe[i].otherhalf > i )
      ps_Draw_Line ( mcq[mhe[i].v0].x, mcq[mhe[i].v0].y,
                     mcq[mhe[i].v1].x, mcq[mhe[i].v1].y );
  ps_Write_Command ( "/Times-Roman findfont 100 scalefont setfont" );
  for ( i = 0; i < nv; i++ ) {
    ps_Set_Gray ( 1.0 );
    ps_Fill_Circle ( mcq[i].x, mcq[i].y, 80 );
    ps_Set_Gray ( 0.0 );
    ps_Draw_Circle ( mcq[i].x, mcq[i].y, 80 );
    sprintf ( s, "%f %f moveto (%d) show", mcq[i].x-30, mcq[i].y-40, i );
    ps_Write_Command ( s );
  }
  for ( i = 0; i < nfac; i++ ) {
    d = mfac[i].degree;
    fhe = mfac[i].firsthalfedge;
    fc = mcq[mhe[mfhei[fhe]].v0];
    for ( j = 1; j < d; j++ )
      AddVector2d ( &fc, &mcq[mhe[mfhei[fhe+j]].v0], &fc );
    MultVector2d ( 1.0/(double)d, &fc, &fc );
    ps_Set_Gray ( 0.5 );
    ps_Draw_Circle ( fc.x, fc.y, 80 );
    sprintf ( s, "%f %f moveto (%d) show", fc.x-30, fc.y-40, i );
    ps_Write_Command ( s );
  }
  ps_Set_Gray ( 0.0 );
  for ( i = 0; i < nhe; i++ ) {
    SubtractPoints2d ( &mcq[mhe[i].v1], &mcq[mhe[i].v0], &ev );
    NormalizeVector2d ( &ev );
    MidPoint2d ( &mcq[mhe[i].v1], &mcq[mhe[i].v0], &fc );
    fc.x -= 20.0*ev.y;
    fc.y += 20.0*ev.x;
    psl_SetLine ( fc.x-280.0*ev.x, fc.y-280.0*ev.y,
                  fc.x+280.0*ev.x, fc.y+280.0*ev.y, 0.0, 1.0 );
    psl_Draw ( 0.0, 0.95, 6.0 );
    psl_Arrow ( 1.0, true );
    fc.x -= 50.0*ev.y;
    fc.y += 50.0*ev.x;
    sprintf ( s, "%f %f moveto (%d) show", fc.x-30, fc.y-40, i );
    ps_Write_Command ( s );
  }

  pkv_SetScratchMemTop ( sp );
} /*DrawMesh*/

int main ( void )
{
  bsf_UserReaders rr;

  pkv_InitScratchMem ( 655360 );
  bsf_ClearReaders ( &rr );
  bsf_BSM4ReadFuncd ( &rr, InputMesh, 10, MAXNV, MAXNHE, MAXNFAC );
  bsf_ReadBSFiled ( ifn, &rr );
  if ( !m1ok )
    exit ( 1 );
  SetupCamera ();
  ps_WriteBBox ( -4, -6, 218, 180 );
  ps_OpenFile ( fn, 600 );
  DrawMesh ( nv1, mv1, mvhei1, mvp1, nhe1, mhe1, nfac1, mfac1, mfhei1 );
  ps_CloseFile ();
  printf ( "%s\n", fn );
  pkv_DestroyScratchMem ();
  exit ( 0 );
} /*main*/

