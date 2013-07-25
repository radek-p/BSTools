
/* //////////////////////////////////////////////////// */
/* This file is a part of the BSTools procedure package */
/* written by Przemyslaw Kiciak                         */
/* //////////////////////////////////////////////////// */

#include <stdio.h>    
#include <stdlib.h>   
#include <math.h>
#include <malloc.h>  
#include <string.h>     

#include <X11/Xlib.h>
#include <X11/Xutil.h>

#include "pkvaria.h"
#include "pknum.h"  
#include "pkgeom.h"   
#include "camerad.h"
#include "multibs.h"
#include "eg1holed.h"
#include "eg2holed.h"
#include "xgedit.h"   

#include "render.h"
#include "edpolicz.h"
#include "edwidgets.h"
#include "drawbezd.h"
#include "splhole.h"
#include "datagend.h"
  

/* ////////////////////////////////////////////////////////////////////////// */
void FindBoundingBox ( Box3d *bbox )
{
  int    i;
  double d;

  bbox->x0 = bbox->x1 = hole_cp[0].x;
  bbox->y0 = bbox->y1 = hole_cp[0].y;
  bbox->z0 = bbox->z1 = hole_cp[0].z;
  for ( i = 1; i < nctrlp; i++ ) {
    if ( hole_cp[i].x > bbox->x1 )      bbox->x1 = hole_cp[i].x;
    else if ( hole_cp[i].x < bbox->x0 ) bbox->x0 = hole_cp[i].x;
    if ( hole_cp[i].y > bbox->y1 )      bbox->y1 = hole_cp[i].y;
    else if ( hole_cp[i].y < bbox->y0 ) bbox->y0 = hole_cp[i].y;
    if ( hole_cp[i].z > bbox->z1 )      bbox->z1 = hole_cp[i].z;
    else if ( hole_cp[i].z < bbox->z0 ) bbox->z0 = hole_cp[i].z;
  }
        /* extend it */
  d = 0.025*(bbox->x1-bbox->x0);  bbox->x0 -= d;  bbox->x1 +=d;
  d = 0.025*(bbox->y1-bbox->y0);  bbox->y0 -= d;  bbox->y1 +=d;
  d = 0.025*(bbox->z1-bbox->z0);  bbox->z0 -= d;  bbox->z1 +=d;
} /*FindBoundingBox*/

void FindDomainBoundingBox ( Box2d *bbox )
{
  int    i;
  double d;

  bbox->x0 = bbox->x1 = domain_cp[0].x;
  bbox->y0 = bbox->y1 = domain_cp[0].y;
  for ( i = 1; i < nctrlp; i++ ) {
    if ( domain_cp[i].x > bbox->x1 )      bbox->x1 = domain_cp[i].x;
    else if ( domain_cp[i].x < bbox->x0 ) bbox->x0 = domain_cp[i].x;
    if ( domain_cp[i].y > bbox->y1 )      bbox->y1 = domain_cp[i].y;
    else if ( domain_cp[i].y < bbox->y0 ) bbox->y0 = domain_cp[i].y;
  }
        /* extend it */
  d = 0.025*(bbox->x1-bbox->x0);  bbox->x0 -= d;  bbox->x1 +=d;
  d = 0.025*(bbox->y1-bbox->y0);  bbox->y0 -= d;  bbox->y1 +=d;
} /*FindDomainBoundingBox*/

void ProjectSurfaceNet ( void )
{
  int     i, j;
  point3d p;

  for ( i = 0; i < 4; i++ )
    for ( j = 0; j < nctrlp; j++ ) {
      CameraProjectPoint3d ( &swind.CPos[i], &hole_cp[j], &p );
      SetPoint2d ( &rpoints[i][j], p.x, p.y );
    }
} /*ProjectSurfaceNet*/

void ProjectDomainNet ( void )
{
  int i;

  for ( i = 0; i < nctrlp; i++ )
    CameraProjectPoint2d ( &domwind.CPos, &domain_cp[i], &rdpoints[i] );
} /*ProjectDomainNet*/

/* ////////////////////////////////////////////////////////////////////////// */
void MyDrawLined ( point2d *p0, point2d *p1 )
{
  xgeDrawLine ( (int)p0->x, (int)p0->y, (int)p1->x, (int)p1->y );
} /*MyDrawLined*/

void MyMarkPointd ( point2d *p )
{
  xgeFillRectangle ( 3, 3, (int)p->x-1, (int)p->y-1 );
} /*MyMarkPointd*/

void MyWriteNumber ( short x, short y, int n )
{
  char s[10];

  sprintf ( s, "%d", n );
  xgeDrawString ( s, x, y );
} /*MyWriteNumber*/

/* ////////////////////////////////////////////////////////////////////////// */
void DrawGHControlNet ( int id )
{
  void *sp;
  int  i, j, k;
  int  *ind;

  sp = pkv_GetScratchMemTop ();
  id &= 0x03;
  ind = pkv_GetScratchMem ( 16*sizeof(int) );
  if ( ind ) {
    xgeSetForeground ( xgec_Green1 );
    for ( i = 0; i < hole_k; i++ ) {
      gh_GetBspInd ( hole_k, i, 0, ind );
      for ( j = 0; j < 4; j++ )
        for ( k = 0; k < 3; k++ )
          MyDrawLined ( &rpoints[id][ind[4*j+k]], &rpoints[id][ind[4*j+k+1]] );
      for ( j = 0; j < 3; j++ )
        for ( k = 0; k < 3; k++ )
          MyDrawLined ( &rpoints[id][ind[4*j+k]], &rpoints[id][ind[4*(j+1)+k]] );
    }
  }
  xgeSetForeground ( xgec_Yellow );
  for ( i = 0; i < nctrlp; i++ )
    if ( !mkhcp[i] )
      MyMarkPointd ( &rpoints[id][i] );
  xgeSetForeground ( xgec_OrangeRed );
  for ( i = 0; i < nctrlp; i++ )
    if ( mkhcp[i] )
      MyMarkPointd ( &rpoints[id][i] );
  pkv_SetScratchMemTop ( sp );
} /*DrawGHControlNet*/

boolean GetSurfBezPatch ( int i, int j, point3d *bcp )
{
#define GetKnotSequence(a,b) &knots[11*(a)+b];
  void    *sp;
  int     *ind, k;
  point3d *cp;
  double  *ukn, *vkn;

  sp = pkv_GetScratchMemTop ();
  ind = pkv_GetScratchMem ( 16*sizeof(int) );
  cp = pkv_GetScratchMem ( 16*sizeof(point3d) );
  if ( ind && cp ) {
    ukn = GetKnotSequence ( (i+hole_k-1)%hole_k, 3 );
    vkn = GetKnotSequence ( i, j );
    gh_GetBspInd ( hole_k, i, j, ind );
    for ( k = 0; k < 16; k++ )
      cp[k] = hole_cp[ind[k]];
    mbs_BSPatchToBezd ( 3, 3, 7, ukn, 3, 7, vkn, 12, (double*)cp,
                        NULL, NULL, NULL, NULL, NULL, NULL, 12, (double*)bcp );
    pkv_SetScratchMemTop ( sp );
    return true;
  }
  else {
    pkv_SetScratchMemTop ( sp );
    return false;
  }
#undef GetKnotSequence
} /*GetSurfBezPatch*/

void DrawGHSurfPatches ( int id )
{
  void    *sp;
  int     i, j;
  point3d *bcp;

  sp = pkv_GetScratchMemTop ();
  id &= 0x03;
  if ( (bcp = pkv_GetScratchMem ( 16*sizeof(point3d) )) ) {
    for ( i = 0; i < hole_k; i++ )
      for ( j = 0; j < 3; j++ ) {
        if ( GetSurfBezPatch ( i, j, bcp ) )
          DrawBezPatch3d ( &swind.CPos[id], 3, 3, bcp, 6, 6,
                           xgec_White, xgec_Grey3  );
      }
  }
  pkv_SetScratchMemTop ( sp );
} /*DrawGHSurfPatches*/

void DrawGHSurfFillingPatches ( int id, int sn )
{
  int i, size;

  id &= 0x03;
  switch ( sn ) {
case 1:
    if ( final_cp1 ) {
      if ( options1.spline ) {
        size = final_lkn1-final_deg1;  size *= size;
        for ( i = 0; i < final_np1; i++ )
          DrawBSPatch3d ( &swind.CPos[id], final_deg1, final_lkn1, final_knots1,
                          final_deg1, final_lkn1, final_knots1,
                          &final_cp1[i*size], 6, 6, xgec_White, xgec_Grey2 );
      }
      else {
        size = final_deg1+1;  size *= size;
        for ( i = 0; i < final_np1; i++ )
          DrawBezPatch3d ( &swind.CPos[id], final_deg1, final_deg1,
                           &final_cp1[i*size], 6, 6, xgec_White, xgec_Grey2 );
      }
    }
    break;
case 2:
    if ( final_cp2 ) {
      if ( options2.spline ) {
        size = final_lkn2-final_deg2;  size *= size;
        for ( i = 0; i < final_np2; i++ )
          DrawBSPatch3d ( &swind.CPos[id], final_deg2, final_lkn2, final_knots2,
                          final_deg2, final_lkn2, final_knots2,
                          &final_cp2[i*size], 6, 6, xgec_Yellow, xgec_Gold );
      }
      else {
        size = final_deg2+1;  size *= size;
        for ( i = 0; i < final_np2; i++ )
          DrawBezPatch3d ( &swind.CPos[id], final_deg2, final_deg2,
                           &final_cp2[i*size], 6, 6, xgec_Yellow, xgec_Gold );
      }
    }
    break;
default:
    break;
  }
} /*DrawGHSurfFillingPatches*/

void DrawGHSurfNumbers ( int id )
{
  void    *sp;
  int     i;
  point3d *bcp;

  sp = pkv_GetScratchMemTop ();
  id &= 0x03;
  if ( (bcp = pkv_GetScratchMem ( 32*sizeof(point3d) )) ) {
    xgeSetForeground ( xgec_Cyan );
    for ( i = 0; i < hole_k; i++ ) {
      if ( GetSurfBezPatch ( i, 2, bcp ) )
        if ( CameraClipPoint3d ( &swind.CPos[id], &bcp[0], &bcp[1] ) )
          MyWriteNumber ( (short)bcp[1].x, (short)bcp[1].y, i );
    }
  }
  pkv_SetScratchMemTop ( sp );
#undef GetKnotSequence
} /*DrawGHSurfNumbers*/

/* ////////////////////////////////////////////////////////////////////////// */
void DrawGHDomainControlNet ( void )
{
  void *sp;
  int  i, j, k;
  int  *ind;

  sp = pkv_GetScratchMemTop ();
  ind = pkv_GetScratchMem ( 16*sizeof(int) );
  if ( ind ) {
    xgeSetForeground ( xgec_Green1 );
    for ( i = 0; i < hole_k; i++ ) {
      gh_GetBspInd ( hole_k, i, 0, ind );
      for ( j = 0; j < 4; j++ )
        for ( k = 0; k < 3; k++ )
          MyDrawLined ( &rdpoints[ind[4*j+k]], &rdpoints[ind[4*j+k+1]] );
      for ( j = 0; j < 3; j++ )
        for ( k = 0; k < 3; k++ )
          MyDrawLined ( &rdpoints[ind[4*j+k]], &rdpoints[ind[4*(j+1)+k]] );
    }
  }
  xgeSetForeground ( xgec_Yellow );
  for ( i = 0; i < nctrlp; i++ )
    if ( !mkdcp[i] )
      MyMarkPointd ( &rdpoints[i] );
  xgeSetForeground ( xgec_OrangeRed );
  for ( i = 0; i < nctrlp; i++ )
    if ( mkdcp[i] )
      MyMarkPointd ( &rdpoints[i] );
  pkv_SetScratchMemTop ( sp );
} /*DrawDomainControlNet*/

boolean GetDomBezPatch ( int i, int j, point2d *bcp )
{
#define GetKnotSequence(a,b) &knots[11*(a)+b];
  void    *sp;
  int     *ind, k;
  point2d *cp;
  double  *ukn, *vkn;

  sp = pkv_GetScratchMemTop ();
  ind = pkv_GetScratchMem ( 16*sizeof(int) );
  cp = pkv_GetScratchMem ( 16*sizeof(point2d) );
  if ( ind && cp ) {
    ukn = GetKnotSequence ( (i+hole_k-1)%hole_k, 3 );
    vkn = GetKnotSequence ( i, j );
    gh_GetBspInd ( hole_k, i, j, ind );
    for ( k = 0; k < 16; k++ )
      cp[k] = domain_cp[ind[k]];
    mbs_BSPatchToBezd ( 2, 3, 7, ukn, 3, 7, vkn, 8, (double*)cp,
                        NULL, NULL, NULL, NULL, NULL, NULL, 8, (double*)bcp );
    pkv_SetScratchMemTop ( sp );
    return true;
  }
  else {
    pkv_SetScratchMemTop ( sp );
    return false;
  }
#undef GetKnotSequence
} /*GetDomBezPatch*/

void DrawDomSurrndPatch ( int n, int m, const point2d *cp )
{
  DrawBezPatch2d ( &domwind.CPos, n, m, cp, 6, 6, xgec_White, xgec_Grey3 );
} /*DrawDomSurrndPatch*/

void DrawDomainSurrPatches ( void )
{
  if ( domain1 )
    gh_DrawDomSurrndPatchesd ( domain1, DrawDomSurrndPatch );
} /*DrawDomainSurrPatches*/

void DrawDomainPatches ( int dn )
{
  int i, size;

  switch ( dn ) {
case 1:
    if ( domain1 && domain_bcp1 && domain_np1 == hole_k ) {
      size = (domain_deg1+1)*(domain_deg1+1);
      for ( i = 0; i < hole_k; i++ )
        DrawBezPatch2d ( &domwind.CPos, domain_deg1, domain_deg1,
                         &domain_bcp1[i*size], 6, 6, xgec_White, xgec_Grey2 );
    }
    break;
case 2:
    if ( domain2 && domain_bcp2 && domain_np2 == hole_k ) {
      size = (domain_deg2+1)*(domain_deg2+1);
      for ( i = 0; i < hole_k; i++ )
        DrawBezPatch2d ( &domwind.CPos, domain_deg2, domain_deg2,
                         &domain_bcp2[i*size], 6, 6, xgec_Yellow, xgec_Gold );
    }
    break;
default:
    break;
  }
} /*DrawDomainPatches*/

void DrawDomainNumbers ( void )
{
  void    *sp;
  int     i;
  point2d *bcp;

  sp = pkv_GetScratchMemTop ();
  if ( (bcp = pkv_GetScratchMem ( 32*sizeof(point3d) )) ) {
    xgeSetForeground ( xgec_Cyan );
    for ( i = 0; i < hole_k; i++ ) {
      if ( GetDomBezPatch ( i, 2, bcp ) ) {
        CameraProjectPoint2d ( &domwind.CPos, &bcp[0], &bcp[1] );
        MyWriteNumber ( (short)bcp[1].x, (short)bcp[1].y, i );
      }
    }
  }
  pkv_SetScratchMemTop ( sp );
#undef GetKnotSequence
} /*DrawDomainNumbers*/

/* ////////////////////////////////////////////////////////////////////////// */
boolean FindNearestPoint ( int npts, point2d *pts,
                           short x, short y, short mindist,
                           int *current_point )
{
  int   i, j;
  short d, e;

  e = (short)(mindist+1);
  j = -1;
  for ( i = 0; i < npts; i++ ) {
    d = (short)(fabs(pts[i].x-(double)x)+fabs(pts[i].y-(double)y));
    if ( d < e ) {
      e = d;
      j = i;
    }
  }
  if ( j >= 0 ) {
    *current_point = j;
    return true;
  }
  else
    return false;
} /*FindNearestPoint*/

/* ////////////////////////////////////////////////////////////////////////// */
boolean FindNearestSurfCPoint ( int id, short x, short y, short mindist )
{
  id &= 0x03;
  if ( view_surf_cp )
    return FindNearestPoint ( nctrlp, rpoints[id], x, y, mindist,
                              &swind.current_point );
  else
    return false;
} /*FindNearestSurfCPoint*/

void SetSurfCPoint ( int id, short x, short y )
{
  int     i;
  point3d p;

  id &= 0x03;
  swind_picture = false;
  i = swind.current_point;
  CameraProjectPoint3d ( &swind.CPos[id], &hole_cp[i], &p );
  p.x = (double)x;
  p.y = (double)y;
  CameraUnProjectPoint3d ( &swind.CPos[id], &p, &hole_cp[i] );
  for ( id = 0; id < 4; id++ ) {
    CameraProjectPoint3d ( &swind.CPos[id], &hole_cp[i], &p );  
    SetPoint2d ( &rpoints[id][i], p.x, p.y );
  }
  InitConstraintFrame ( 1 );
  InitConstraintFrame ( 2 );
} /*SetSurfCPoint*/

void SelectCPoints ( Box2s *sel_rect, int npts, point2d *rpts,
                     boolean *mkpts, boolean value )
{
  int    i, j;
  double d, e;

  j = -1;
  if ( sel_rect->x1-sel_rect->x0 < 2 && sel_rect->y1-sel_rect->y0 < 2 ) {
    d = 4.0;
    for ( i = 0; i < npts; i++ ) {
      e = fabs(sel_rect->x0-rpts[i].x)+fabs(sel_rect->y0-rpts[i].y);
      if ( e < d ) {
        d = e;
        j = i;
      }
    }
    if ( d <= 3.0 )
      mkpts[j] = value;
  }
  else {
    for ( i = 0; i < npts; i++ )
      if ( rpts[i].x >= sel_rect->x0 && rpts[i].x <= sel_rect->x1 &&
           rpts[i].y >= sel_rect->y0 && rpts[i].y <= sel_rect->y1 )
        mkpts[i] = value;
  }
} /*SelectCPoints*/

void SelectSurfCPoints ( int id )
{
  id &= 0x03;
  if ( view_surf_cp )
    SelectCPoints ( &swind.selection_rect, nctrlp, rpoints[id], mkhcp, true );
} /*SelectSurfCPoints*/

void UnselectSurfCPoints ( int id )
{
  id &= 0x03;
  if ( view_surf_cp )
    SelectCPoints ( &swind.selection_rect, nctrlp, rpoints[id], mkhcp, false );
} /*UnselectSurfCPoints*/

void SaveSurfCPoints ( void )
{
  memcpy ( saved_cp, hole_cp, nctrlp*sizeof(point3d) );
} /*SaveSurfCPoints*/

void TransformSurfCPoints ( void )
{
  int i;

  swind_picture = false;
  for ( i = 0; i < nctrlp; i++ )
    if ( mkhcp[i] )
      TransPoint3d ( &swind.gwtrans, &saved_cp[i], &hole_cp[i] );
  ProjectSurfaceNet ();
  InitConstraintFrame ( 1 );
  InitConstraintFrame ( 2 );
} /*TransformSurfCPoints*/

/* ////////////////////////////////////////////////////////////////////////// */
boolean FindNearestDomainCPoint ( short x, short y, short mindist )
{
  if ( view_dom_cp )
    return FindNearestPoint ( nctrlp, rdpoints, x, y, mindist,
                              &domwind.current_point );
  else
    return false;
} /*FindNearestDomainCPoint*/

void SetDomainCPoint ( short x, short y )
{
  int i;

  swind_picture = false;
  i = domwind.current_point;
  SetPoint2d ( &rdpoints[i], x, y );
  CameraUnProjectPoint2d ( &domwind.CPos, &rdpoints[i], &domain_cp[i] );
  UpdateDomains ();
} /*SetDomainCPoint*/

void SelectDomainCPoints ( void )
{
  if ( view_dom_cp )
    SelectCPoints ( &domwind.selection_rect, nctrlp, rdpoints, mkdcp, true );
} /*SelectDomainCPoints*/

void UnselectDomainCPoints ( void )
{
  if ( view_dom_cp )
    SelectCPoints ( &domwind.selection_rect, nctrlp, rdpoints, mkdcp, false );
} /*UnselectDomainCPoints*/

void SaveDomainCPoints ( void )
{
  memcpy ( saved_cp, domain_cp, nctrlp*sizeof(point2d) );
} /*SaveDomainCPoints*/

void TransformDomainCPoints ( void )
{
  double *scp;
  int    i;

  swind_picture = false;
  scp = &saved_cp[0].x;
  for ( i = 0; i < nctrlp; i++ )
    if ( mkdcp[i] )
      TransPoint2d ( &domwind.gwtrans, (point2d*)&scp[2*i], &domain_cp[i] );
  UpdateDomains ();
  ProjectDomainNet ();
} /*TransformDomainCPoints*/

