
/* //////////////////////////////////////////////////// */
/* This file is a part of the BSTools procedure package */
/* written by Przemyslaw Kiciak                         */
/* //////////////////////////////////////////////////// */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
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


void SetBCColour ( int d, int i, xgecolour_int c0, xgecolour_int c1 )
{
  if ( i == 0 || i == d )
    xgeSetForeground ( c0 );
  else if ( i == 1 )
    xgeSetForeground ( c1 );
} /*SetBCColour*/

void Get2dClipLines ( CameraRecd *CPos, vector3d *cliplines )
{
  SetVector3d ( &cliplines[0], 1.0, 0.0,  (double)CPos->xmin );
  SetVector3d ( &cliplines[1], 0.0, 1.0,  (double)CPos->ymin );
  SetVector3d ( &cliplines[2], -1.0, 0.0, (double)(CPos->xmin+CPos->width) );
  SetVector3d ( &cliplines[3], 0.0, -1.0, (double)(CPos->ymin+CPos->height) );
} /*Get2dClipLines*/

void DrawBezPatch2d ( CameraRecd *CPos,
                      int n, int m, const point2d *cp, 
                      int d0, int d1, xgecolour_int c0, xgecolour_int c1 )
{
  void     *sp;
  int      i, j;
  point2d  *p, *q;
  double   t;
  vector3d cliplines[4];

  sp = pkv_GetScratchMemTop ();
  i = max(n,m)+1;
  if ( (p = pkv_GetScratchMem ( 2*i*sizeof(point2d) )) ) {
    q = &p[i];
    Get2dClipLines ( CPos, cliplines );
    for ( i = 0; i <= d0; i++ ) {
      SetBCColour ( d0, i, c0, c1 );
      t = (double)i/(double)d0;
      mbs_multiBCHornerd ( n, 1, 2*(m+1), 0, (double*)cp, t, (double*)p );
      for ( j = 0; j <= m; j++ )
        CameraProjectPoint2d ( CPos, &p[j], &q[j] );
      mbs_ClipBC2d ( 4, cliplines, m, q, xge_DrawBC2d );
    }
    for ( i = 0; i <= d1; i++ ) {
      SetBCColour ( d1, i, c0, c1 );
      t = (double)i/(double)d1;
      mbs_multiBCHornerd ( m, n+1, 2, 2*(m+1), (double*)cp, t, (double*)p );
      for ( j = 0; j <= n; j++ )
        CameraProjectPoint2d ( CPos, &p[j], &q[j] );
      mbs_ClipBC2d ( 4, cliplines, n, q, xge_DrawBC2d );
    }
  }
  pkv_SetScratchMemTop ( sp );
} /*DrawBezPatch2d*/

void ProjectBezCP3d ( CameraRecd *CPos, point3d *p, point3d *q )
{
  TransPoint3d ( &CPos->CTr, p, q );
  if ( CPos->parallel )
    q->z = 1.0;
  else {
    q->x += CPos->vd.persp.xi0*q->z;
    q->y += CPos->vd.persp.eta0*q->z;
  }
} /*ProjectBezCP3d*/

void DrawBezPatch3d ( CameraRecd *CPos,
                      int n, int m, const point3d *cp,
                      int d0, int d1, xgecolour_int c0, xgecolour_int c1 )
{
  void    *sp;
  int     i, j;
  point3d *p, *q;
  double   t;
  vector3d cliplines[4];

  sp = pkv_GetScratchMemTop ();
  i = max(n,m)+1;
  if ( (p = pkv_GetScratchMem ( 2*i*sizeof(point3d) )) ) {
    q = &p[i];
    Get2dClipLines ( CPos, cliplines );
    for ( i = 0; i <= d0; i++ ) {
      SetBCColour ( d0, i, c0, c1 );
      t = (double)i/(double)d0;
      mbs_multiBCHornerd ( n, 1, 3*(m+1), 0, (double*)cp, t, (double*)p );
      for ( j = 0; j <= m; j++ )
        ProjectBezCP3d ( CPos, &p[j], &q[j] );
      mbs_ClipBC2Rd ( 4, cliplines, m, q, xge_DrawBC2Rd );
    }
    for ( i = 0; i <= d1; i++ ) {
      SetBCColour ( d1, i, c0, c1 );
      t = (double)i/(double)d1;
      mbs_multiBCHornerd ( m, n+1, 3, 3*(m+1), (double*)cp, t, (double*)p );
      for ( j = 0; j <= n; j++ )
        ProjectBezCP3d ( CPos, &p[j], &q[j] );
      mbs_ClipBC2Rd ( 4, cliplines, n, q, xge_DrawBC2Rd );
    }
    
  }
  pkv_SetScratchMemTop ( sp );
} /*DrawBezPatch3d*/

void DrawBSPatch3d ( CameraRecd *CPos,
                     int n, int lknu, const double *knu,  
                     int m, int lknv, const double *knv, const point3d *cp,  
                     int d0, int d1, xgecolour_int c0, xgecolour_int c1 )
{
  void   *sp;
  int    ku, kv, pitch, md, i, j;
  double *b, *c;

  sp = pkv_GetScratchMemTop ();
  ku = mbs_NumKnotIntervalsd ( n, lknu, knu );
  kv = mbs_NumKnotIntervalsd ( m, lknv, knv );
  md = 3*(m+1);
  pitch = md*kv;
  b = pkv_GetScratchMemd ( pitch*(n+1)*ku );
  c = pkv_GetScratchMemd ( md*(n+1) );
  if ( b && c ) {
    mbs_BSPatchToBezd ( 3, n, lknu, knu, m, lknv, knv, 3*(lknv-m), (double*)cp,
                        &ku, NULL, NULL, &kv, NULL, NULL, pitch, b );
    for ( i = 0; i < ku; i++ )
      for ( j = 0; j < kv; j++ ) {
        pkv_Selectd ( n+1, md, pitch, md, &b[(n+1)*i*pitch+md*j], c );
        DrawBezPatch3d ( CPos, n, m, (point3d*)c, d0, d1, c0, c1 );
      }
  }
  pkv_SetScratchMemTop ( sp );
} /*DrawBSPatch3d*/

