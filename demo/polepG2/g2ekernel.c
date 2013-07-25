
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2005, 2007                            */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <malloc.h>
#include <pthread.h>

#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <X11/cursorfont.h>

#include "pkvaria.h"
#include "pknum.h"
#include "pkgeom.h"
#include "camera.h"
#include "multibs.h"
#include "raybez.h"
#include "eg2holef.h"

#include "oldxgedit.h"
#include "datagenf.h"
#include "edg2hole.h"
#include "g2ekernel.h"
#include "drawbez.h"
#include "render.h"


CameraRecf PPos[3], DomPPos;
CameraRecf   CPos;

int hole_k;
int hole_np;
int nconstr;

vector2f firstpartition[GH_MAX_K];

float surfparams[4] = {0.0, 0.0, 0.0, 0.0};
float domparams[3] = {0.0, 0.0, 0.0};

float    knots[GH_MAX_K][11];
point2f  domcp[12*GH_MAX_K+1];
point3f  surfcp[12*GH_MAX_K+1];
point3f  constrcp[4*GH_MAX_K+1];
vector3f constrder[4*GH_MAX_K+1];

vector3f ConstraintNormal = { 0.0, 0.0, 1.0 };

boolean swCoonsPatches = true;
boolean swBezierPatches = false;
boolean swConstraintsOn = false;
boolean swNormalConstr = false;
boolean swConstraint[4*GH_MAX_K+1] =
  { false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false };
boolean swZConstraint[4*GH_MAX_K+1] =
  { false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false };
boolean swPointMode = false;
boolean swZeroDer   = false;

GHoleDomainf *domain = NULL;

point2f  central_point;
vector2f domcder[GH_MAX_K];

boolean CentralPointValid     = false;
boolean FormMatrixValid       = false;
boolean FinalSurfValid        = false;
boolean ExtFormMatrixValid    = false;
boolean ConstraintsValid      = false;
boolean ConstrMatrixValid     = false;
boolean ExtConstrMatrixValid  = false;
boolean NConstrMatrixValid    = false;
boolean ExtNConstrMatrixValid = false;
boolean NLFinalSurfValid      = false;

point2f knoti[GH_MAX_K][11];
point2f domcpi[12*GH_MAX_K+1];
point2f surfcpi[4][12*GH_MAX_K+1];
point2f constrcpi[4][4*GH_MAX_K+1];
point2f constrnvi[4];

point2f central_pointi;
point2f domcderi[GH_MAX_K];

static int colour1, colour2;
static int dens1, dens2;

Box3f RefBBox;

point3f FinalCP[2*GH_MAX_K][(G2H_FINALDEG+1)*(G2H_FINALDEG+1)];
int finaldegu, finaldegv;
int nfinal;

/* queues for switching the constraints for derivatives of order 1 and 2 */
char constrql1 = 0;
static char ConstrQueue1[2] = { -1, -1 };
char constrql2 = 0;
static char ConstrQueue2[3] = { -1, -1, -1 };

/* ////////////////////////////////////////////////////////////////////////// */
static int G2OptionProc ( GHoleDomainf *domain, int query, int qn,
                          int *ndata, int **idata, float **fdata )
{
  switch ( query ) {
case G2HQUERY_CENTRAL_POINT:
    if ( CentralPointValid && swDisplayCentralPoint ) {
      *ndata = 2;
      *idata = NULL;
      *fdata = &central_point.x;
      return G2H_CENTRAL_POINT_GIVEN;
    }
    else
/*      return G2H_DEFAULT;*/
      return G2H_CENTRAL_POINT_ALT;

case G2HQUERY_CENTRAL_DERIVATIVES1:
    if ( CentralPointValid && swUseDerivatives1 ) {
      *ndata = 2*hole_k;
      *idata = NULL;
      *fdata = &domcder[0].x;
      return G2H_CENTRAL_DERIVATIVES1_GIVEN;
    }
    else
      return G2H_DEFAULT;

case G2HQUERY_DOMAIN_CURVES:
    if ( swAltDomCurves )
      return G2H_DOMAIN_CURVES_DEG4;
    else
      return G2H_DEFAULT;

case G2HQUERY_BASIS:
    if ( swRestrictBasis )
      return G2H_USE_RESTRICTED_BASIS;
    else
      return G2H_DEFAULT;

default:
    return G2H_DEFAULT;
  }
} /*G2OptionProc*/

/* ////////////////////////////////////////////////////////////////////////// */
void RecreateDomain ()
{
  if ( domain )
    gh_DestroyDomainf ( domain );
  domain = gh_CreateDomainf ( hole_k, &knots[0][0], domcp );

  if ( !domain ) {
    printf ( "error: cannot create domain" );
    exit ( 1 );
  }
  g2h_SetOptionProcf ( domain, G2OptionProc );
  if ( g2h_ComputeBasisf ( domain ) )
    if ( !CentralPointValid || !swUseDerivatives1 ) {
      g2h_ExtractCentralPointf ( domain, &central_point, domcder );
      CentralPointValid = true;
    }
  FormMatrixValid = ExtFormMatrixValid = FinalSurfValid =
    NLFinalSurfValid = ConstrMatrixValid = ExtConstrMatrixValid = false;
  swDisplayNLFinalPatches = PictureIsOn = false;
} /*RecreateDomain*/

boolean UpdateFormMatrix ()
{
  if ( !FormMatrixValid )
    FormMatrixValid = g2h_ComputeFormMatrixf ( domain );
  return FormMatrixValid;
} /*UpdateFormMatrix*/

boolean UpdateExtFormMatrix ()
{
  if ( !ExtFormMatrixValid )
    ExtFormMatrixValid = g2h_ComputeExtFormMatrixf ( domain );
  return ExtFormMatrixValid;
} /*UpdateExtFormMatrix*/

void OutFinalPatch ( int n, int m, const float *cp, void *usrptr )
{
  if ( nfinal < 2*GH_MAX_K && n == G2H_FINALDEG && m == G2H_FINALDEG ) {
    memcpy ( FinalCP[nfinal], cp, (n+1)*(m+1)*sizeof(point3f) );
    finaldegu = n;
    finaldegv = m;
    nfinal++;
  }
} /*OutFinalPatch*/

boolean UpdateFinalSurf ()
{
  void     *sp;
  vector3f *rhs;
  float    *rhsf;

  if ( swCoonsPatches ) {
    if ( FormMatrixValid ) {
      nfinal = 0;
      if ( !swConstraintsOn || !nconstr || !ConstraintsValid ) {
        FinalSurfValid = g2h_FillHolef ( domain, 3, (float*)surfcp,
                                         NULL, NULL, OutFinalPatch );
        if ( !ConstraintsValid )
          GetConstraintsHandle ();
      }
      else {
        sp = pkv_GetScratchMemTop ();
        if ( !swNormalConstr ) {
          if ( !ConstrMatrixValid )
            if ( !SetupConstraintsMatrix () )
              return false;
          rhs = pkv_GetScratchMem ( nconstr*sizeof(point3f) );
          if ( !rhs )
            return false;
          SetupConstraintsRHS ( rhs );
          FinalSurfValid = g2h_FillHoleConstrf ( domain, 3, (float*)surfcp,
                                   nconstr, (float*)rhs, NULL, NULL,
                                   OutFinalPatch );
        }
        else {
          if ( !NConstrMatrixValid )
            if ( !SetupNConstraintsMatrix () )
              return false;
          rhsf = pkv_GetScratchMemf ( nconstr );
          if ( !rhsf )
            return false;
          SetupNConstraintsRHS ( rhsf );
          FinalSurfValid = g2h_FillHoleAltConstrf ( domain, 3, (float*)surfcp,
                                   nconstr, rhsf, NULL, NULL, OutFinalPatch );
        }
        pkv_SetScratchMemTop ( sp );
      }
      return FinalSurfValid;
    }
    else
      return false;
  }
  else {
    if ( ExtFormMatrixValid ) {
      nfinal = 0;
      if ( !swConstraintsOn || !nconstr || !ConstraintsValid ) {
        FinalSurfValid = g2h_ExtFillHolef ( domain, 3, (float*)surfcp,
                                            NULL, NULL, OutFinalPatch );
        if ( !ConstraintsValid )
          GetConstraintsHandle ();
      }
      else {
        sp = pkv_GetScratchMemTop ();
        if ( !swNormalConstr ) {
          if ( !ExtConstrMatrixValid )
            if ( !SetupExtConstraintsMatrix () )
              return false;
          rhs = pkv_GetScratchMem ( nconstr*sizeof(point3f) );
          if ( !rhs )
            return false;
          SetupConstraintsRHS ( rhs );
          FinalSurfValid = g2h_ExtFillHoleConstrf ( domain, 3, (float*)surfcp,
                                    nconstr, (float*)rhs, NULL, NULL,
                                    OutFinalPatch );
        }
        else {
          if ( !ExtNConstrMatrixValid )
            if ( !SetupExtNConstraintsMatrix () )
              return false;
          rhsf = pkv_GetScratchMemf ( nconstr );
          if ( !rhsf )
            return false;
          SetupNConstraintsRHS ( rhsf );
          FinalSurfValid = g2h_ExtFillHoleAltConstrf ( domain, 3, (float*)surfcp,
                                  nconstr, rhsf, NULL, NULL,
                                  OutFinalPatch );
        }
        pkv_SetScratchMemTop ( sp );
      }
      return FinalSurfValid;
    }
    else
      return false;
  }
} /*UpdateFinalSurf*/

boolean UpdateNLFinalSurf ( void )
{
  typedef void (*outfunc3) (int n, int m, const point3f *cp, void *usrptr );
  void     *sp;
  vector3f *rhs;
  float    *rhsf;

  sp = pkv_GetScratchMemTop ();
  nfinal = hole_k;
  if ( swCoonsPatches ) {
    if ( !swConstraintsOn || !nconstr )
      NLFinalSurfValid = g2h_NLFillHolef ( domain, surfcp, NULL,
                                           NULL, (outfunc3)OutFinalPatch );
    else {
      if ( !swNormalConstr ) {
        if ( !ConstrMatrixValid )
          if ( !SetupConstraintsMatrix () )
            return false;
        rhs = pkv_GetScratchMem ( nconstr*sizeof(point3f) );
        if ( !rhs )
          return false;
        SetupConstraintsRHS ( rhs );
        NLFinalSurfValid =
          g2h_NLFillHoleConstrf ( domain, surfcp, nconstr, rhs,
                                  NULL, NULL, (outfunc3)OutFinalPatch );
      }
      else {
        if ( !NConstrMatrixValid )
          if ( !SetupNConstraintsMatrix () )
            return false;
        rhsf = pkv_GetScratchMemf ( nconstr );
        if ( !rhsf )
          return false;
        SetupNConstraintsRHS ( rhsf );
        NLFinalSurfValid = g2h_NLFillHoleAltConstrf ( domain, surfcp, nconstr,
                                      rhsf, NULL, NULL,
                                      (outfunc3)OutFinalPatch  ); 
      }
    }
  }
  else {
    if ( !swConstraintsOn || !nconstr )
      NLFinalSurfValid = g2h_NLExtFillHolef ( domain, surfcp, NULL,
                                              NULL, (outfunc3)OutFinalPatch );
    else {
      if ( !swNormalConstr ) {
        if ( !ExtConstrMatrixValid )
          if ( !SetupExtConstraintsMatrix () )
            return false;
        rhs = pkv_GetScratchMem ( nconstr*sizeof(point3f) );
        if ( !rhs )
          return false;
        SetupConstraintsRHS ( rhs );
        NLFinalSurfValid =
          g2h_NLExtFillHoleConstrf ( domain, surfcp, nconstr, rhs,
                                     NULL, NULL, (outfunc3)OutFinalPatch );
      }
      else {
        if ( !ExtNConstrMatrixValid )
          if ( !SetupExtNConstraintsMatrix () )
            return false;
          rhsf = pkv_GetScratchMemf ( nconstr );
          if ( !rhsf )
            return false;
          SetupNConstraintsRHS ( rhsf );
          NLFinalSurfValid = g2h_NLExtFillHoleAltConstrf ( domain, surfcp,
                                 nconstr, rhsf, NULL, NULL,
                                 (outfunc3)OutFinalPatch );
      }
    }
  }
  pkv_SetScratchMemTop ( sp );
  return NLFinalSurfValid;
} /*UpdateNLFinalSurf*/

/* ////////////////////////////////////////////////////////////////////////// */
void SetupRefBox ( float x0, float x1, float y0, float y1, float z0, float z1 )
{
  RefBBox.x0 = x0;  RefBBox.x1 = x1;
  RefBBox.y0 = y0;  RefBBox.y1 = y1;
  RefBBox.z0 = z0;  RefBBox.z1 = z1;
} /*SetupRefBox*/

void ProjectDomainCP ()
{
  int i;

  for ( i = 0; i < hole_np; i++ )
    CameraProjectPoint2f ( &DomPPos, &domcp[i], &domcpi[i] );
} /*ProjectDomainCP*/

void ProjectDomainCentralPoint ()
{
  int     i;
  point2f q;

  CameraProjectPoint2f ( &DomPPos, &central_point, &central_pointi );
  for ( i = 0; i < hole_k; i++ ) {
    AddVector2f ( &central_point, &domcder[i], &q );
    CameraProjectPoint2f ( &DomPPos, &q, &domcderi[i] );
  }
} /*ProjectDomainCentralPoint*/

void ProjectScene ( int id )
{
  int     i;
  point3f q;

  if ( id < 3 )
    for ( i = 0; i < hole_np; i++ ) {
      CameraProjectPoint3f ( &PPos[id], &surfcp[i], &q );
      memcpy ( &surfcpi[id][i], &q, sizeof(point2f) );
    }
  else
    for ( i = 0; i < hole_np; i++ ) {
      CameraProjectPoint3f ( &CPos, &surfcp[i], &q );
      memcpy ( &surfcpi[3][i], &q, sizeof(point2f) );
    }
} /*ProjectScene*/

void ProjectSurfCP ()
{
  int id;

  for ( id = 0; id < 4; id++ )
    ProjectScene ( id );
} /*ProjectSurfCP*/

void ResetCPos ()
{
  float    bx, by, bz, r;
  vector3f v;

  CameraInitPosf ( &CPos );
  SetPoint3f ( &CPos.g_centre, 0.5*(RefBBox.x0+RefBBox.x1),
           0.5*(RefBBox.y0+RefBBox.y1), 0.5*(RefBBox.z0+RefBBox.z1) );
  CPos.c_centre = CPos.g_centre;
  bx = RefBBox.x0-RefBBox.x1;
  by = RefBBox.y0-RefBBox.y1;
  bz = RefBBox.z0-RefBBox.z1;
  r = sqrt ( bx*bx + by*by + bz*bz );
  SetVector3f ( &v, 0.0, 0.0, -10.0*r );
  CameraMoveGf ( &CPos, &v );
  CameraZoomf ( &CPos, 5.0 );
} /*ResetCPos*/

void SetupPerspProj ( ed_rect *edr, boolean reset )
{
  CameraInitFramef ( &CPos, false, false,
                     edr->w, edr->h, edr->x, edr->y, 1.0, 4 );
  if ( reset )
    ResetCPos ();
  else
    CameraSetMappingf ( &CPos );
} /*SetupCPos*/

void SetupParProj ( ed_rect *edr )
{
  float dx, dy, dz, bx, by, bz, s;
  float aspect = 1.0;
  int   i;
  point3f centre;

  dx = min ( edr[0].w, edr[2].w );
  dy = min ( edr[1].w, edr[2].h/aspect );
  dz = min ( edr[0].h, edr[1].h )/aspect;
  bx = RefBBox.x1-RefBBox.x0;
  by = RefBBox.y1-RefBBox.y0;
  bz = RefBBox.z1-RefBBox.z0;
  s = min ( (dx/bx), (dy/by) );
  s = min ( s, (dz/bz) );
  SetPoint3f ( &centre, 0.5*(RefBBox.x0+RefBBox.x1),
          0.5*(RefBBox.y0+RefBBox.y1), 0.5*(RefBBox.z0+RefBBox.z1) );
  PPos[0].vd.para.wdt = edr[0].w/s;
  PPos[1].vd.para.wdt = edr[1].w/s;
  PPos[2].vd.para.wdt = edr[2].w/s;
  for ( i = 0; i < 3; i++ ) {
    CameraInitFramef ( &PPos[i], true, false,
                         edr[i].w, edr[i].h,
                         edr[i].x, edr[i].y, aspect, 4 );
    PPos[i].position = centre;
    PPos[i].vd.para.dim_case = 1;
    CameraSetMappingf ( &PPos[i] );
    CameraTurnGf ( &PPos[i], 0.0, PI, 0.0 );
  }
  CameraTurnGf ( &PPos[0], 0.0, -0.5*PI, 0.0 );
  CameraTurnGf ( &PPos[1], 0.0, -0.5*PI, -0.5*PI );
} /*SetupParProj*/

void SetupDomPPos ( int w, int h, int x, int y, float d )
{
  CameraInitFramef ( &DomPPos, true, false, w, h, x, y, 1.0, 4 );
  DomPPos.vd.para.diag = d;
  CameraSetMappingf ( &DomPPos );
  CameraTurnGf ( &DomPPos, 0.0, PI, 0.0 );
} /*SetupDomPPos*/

void InitProjections ()
{
  int i;

  /* The window dimensions here are irrelevant, as this */
  /* procedure is called in order to initialise the viewer position. */
  /* During window initialization and after each size change */
  /* the proper dimensions will be set/ */

  CameraInitFramef ( &CPos, false, false, 100.0, 100.0, 0.0, 0.0, 1.0, 4 );
  CameraInitPosf ( &CPos );

  for ( i = 0; i < 4; i++ ) {
    CameraInitFramef ( &PPos[i], true, false, 100.0, 100.0, 0.0, 0.0, 1.0, 4 );
    CameraInitPosf ( &PPos[i] );
  }
  CameraInitFramef ( &DomPPos, true, false, 100.0, 100.0, 0.0, 0.0, 1.0, 4 );
  CameraInitPosf ( &DomPPos );
  DomPPos.vd.para.diag = 12.0;
  CameraSetMappingf ( &DomPPos );
} /*InitProjections*/

/* ////////////////////////////////////////////////////////////////////////// */
void MyDrawLine ( point2f *p1, point2f *p2 )
{
  XDrawLine ( thedisplay, thepixmap, thegc,
              (int)p1->x, (int)p1->y, (int)p2->x, (int)p2->y );
} /*MyDrawLine*/

void MyMarkRect ( point2f *p )
{
  int x, y;

  x = (int)p->x;  y = (int)p->y;
  XFillRectangle ( thedisplay, thepixmap, thegc, x-1, y-1, 3, 3 );
} /*MyMarkRect*/

void RedrawDomainCP ( boolean bright )
{
  void *sp;
  int  i, j, k;
  int  *ind;

  sp = pkv_GetScratchMemTop ();
  ind = pkv_GetScratchMem ( 16*sizeof(int) );
  if ( !ind )
    exit ( 1 );
  if ( bright )
    XSetForeground ( thedisplay, thegc, c_green );
  else
    XSetForeground ( thedisplay, thegc, c_blue );
  for ( i = 0; i < hole_k; i++ ) {
    gh_GetBspInd ( hole_k, i, 0, ind );
    for ( j = 0; j < 4; j++ )
      for ( k = 0; k < 3; k++ )
        MyDrawLine ( &domcpi[ind[4*j+k]], &domcpi[ind[4*j+k+1]] );
    for ( j = 0; j < 3; j++ )
      for ( k = 0; k < 3; k++ )
        MyDrawLine ( &domcpi[ind[4*j+k]], &domcpi[ind[4*(j+1)+k]] );
  }
  if ( bright )
    XSetForeground ( thedisplay, thegc, c_yellow );
  for ( i = 0; i < hole_np; i++ )
    MyMarkRect ( &domcpi[i] );
  pkv_SetScratchMemTop ( sp );
} /*RedrawDomainCP*/

static void DrawDomBezPatch ( int n, int m, const point2f *cp )
{
  DrawBezPatch2 ( &DomPPos, n, m, cp, dens1, dens2, colour1, colour2 );
} /*DrawDomBezPatch*/

static void DrawDomAuxPatch ( int n, int m, const point2f *cp )
{
  DrawBezPatch2a ( &DomPPos, n, m, cp, 0.0, 0.33, 0.0, 1.0,
                   dens1, dens2, colour1, colour2 );
} /*DrawDomAuxPatch*/

void RedrawDomainSurrPatches ()
{
  colour1 = c_lt_grey;
  colour2 = c_dk_grey;
  dens1 = dens2 = 8;
  gh_DrawDomSurrndPatchesf ( domain, DrawDomBezPatch );
} /*RedrawDomainSurrPatches*/

void RedrawDomainCurves ()
{
  colour1 = c_white;
  colour2 = c_dk_grey;
  dens1 = dens2 = 1;
  g2h_DrawDiPatchesf ( domain, DrawDomBezPatch );
} /*RedrawDomainCurves*/

void RedrawDomainAuxPatches ()
{
  colour1 = c_yellow;
  colour2 = c_red;
  dens1 = dens2 = 8;
  g2h_DrawDomAuxPatchesf ( domain, DrawDomAuxPatch );
} /*RedrawDomainAuxPatches*/

void RedrawDomainPatches ()
{
  colour1 = c_white;
  colour2 = c_lt_grey;
  dens1 = dens2 = 8;
  g2h_DrawDiPatchesf ( domain, DrawDomBezPatch );
} /*RedrawDomainPatches*/

void RedrawDomainPartition ()
{
} /*RedrawDomainPartition*/

static int pointnum;

static void DrawDomNumbers ( int n, int m, const point2f *cp )
{
  point2f p;
  char    s[8];

  CameraProjectPoint2f ( &DomPPos, &cp[m], &p );
  sprintf ( s, "%d", pointnum );
  pointnum++;
  XDrawString ( thedisplay, thepixmap, thegc,
                (int)p.x-2, (int)p.y+5, s, strlen(s) );
} /*DrawDomNumbers*/

void RedrawDomainNumbers ()
{
  XSetForeground ( thedisplay, thegc, c_cyan );
  pointnum = 0;
  g2h_DrawDomAuxPatchesf ( domain, DrawDomNumbers );
} /*RedrawDomainNumbers*/

/* ////////////////////////////////////////////////////////////////////////// */
int FindNearestDomCP ( int x, int y )
{
#define TOL 10.0
  int i, j;
  float dmin, d;

  dmin = TOL+1.0;
  j = -1;
  for ( i = 0; i < hole_np; i++ ) {
    d = fabs(domcpi[i].x-x)+fabs(domcpi[i].y-y);
    if ( d <= dmin ) {
      dmin = d;
      j = i;
    }
  }
  if ( dmin > TOL )
    j = -1;
  return j;
#undef TOL
} /*FindNearestDomCP*/

void SetDomCP ( int np, int x, int y )
{
  if ( np >= 0 && np < hole_np ) {
    SetPoint2f ( &domcpi[np], (float)x, (float)y );
    CameraUnProjectPoint2f ( &DomPPos, &domcpi[np], &domcp[np] );
    RecreateDomain ();
  }
} /*SetDomCP*/

/* ////////////////////////////////////////////////////////////////////////// */
void RedrawHoleKnots ( int w, int h, int x, int y )
{
  int     i, j;
  int     uy;
  point2f p;

  XSetForeground ( thedisplay, thegc, c_lt_green );
  for ( i = 0; i < hole_k; i++ ) {
    uy = y + 4 +i*6;
    XDrawLine ( thedisplay, thepixmap, thegc, x+10, uy, x+w-10, uy );
  }
  XSetForeground ( thedisplay, thegc, c_white );
  for ( i = 0; i < hole_k; i++ ) {
    p.y = (float)(y + 4 + i*6);
    for ( j = 1; j < 10; j++ ) {
      p.x = (float)(x+10) + knots[i][j]*(float)(w-20);
      MyMarkRect ( &p );
    }
  }
} /*RedrawHoleKnots*/

boolean FindNearestKnot ( int w, int h, int xc, int yc,
                          int x, int y, int *ink, int *jnk )
{
#define TOL 7
  int xu, yu, d, dmin;
  int i, ii = 0, j;

  dmin = TOL+1;
  for ( i = 0; i < hole_k; i++ ) {
    yu = yc + 4 + i*6;
    d = abs ( yu-y );
    if ( d < dmin ) {
      *ink = ii = i;  dmin = d;
    }
  }
  if ( dmin <= TOL ) {
    dmin = TOL+1;    
    for ( j = 1; j < 10; j++ ) {
      xu = xc+10 + (int)(knots[ii][j]*(float)(w-20));
      d = abs ( xu-x );
      if ( d < dmin ) {
        *jnk = j;  dmin = d;
      }
    }
  }
  return ( dmin <= TOL );
#undef TOL
} /*FindNearestKnot*/

boolean SetKnot ( int w, int h, int xc, int yc, int i, int j, int x )
{
  int   ii, jj, k;
  float u, v;

  if ( i >= 0 && i < hole_k && j > 0 && j < 10 ) {
    u = (float)(x-xc-10)/(float)(w-20);
    v = 1.0-u;
    jj = 10-j;
    switch ( j ) {
  case 2:
  case 3:
      k = (i+hole_k-2) % hole_k;
      goto proceed1;

  case 7:
  case 8:
      k = (i+2) % hole_k;
proceed1:
      if ( u > knots[i][j-1] && u < knots[i][j+1] &&
           v > knots[k][jj-1] && v < knots[k][jj+1] ) {
        knots[i][j] = u;
        knots[k][jj] = v;
        break;
      }
      else
        return false;

default:
      if ( hole_k == 8 ) {         /* in general - divisible by 4 */
        if ( u > knots[i][j-1] && u < knots[i][j+1] ) {
          knots[i][j] = u;
          knots[(i+2) % hole_k][10-j] = v;
          knots[(i+4) % hole_k][j] = u;
          knots[(i+6) % hole_k][10-j] = v;
          break;
        }
        else
          return false;
      }
      else if ( j != 5 ) {
        if ( hole_k & 0x01 ) {  /* odd */
          for ( ii = 0; ii < hole_k; ii++ )
            if ( u <= knots[ii][j-1] || u >= knots[ii][j+1] ||
                 v <= knots[ii][jj-1] || v >= knots[ii][jj+1] )
              return false;
          for ( ii = 0; ii < hole_k; ii++ ) {
            knots[ii][j] = u;
            knots[ii][jj] = v;
          }
        }
        else {                       /* even, indivisible by 4 */
          for ( ii = i & 0x01; ii < hole_k; ii += 2 )
            if ( u <= knots[ii][j-1] || u >= knots[ii][j+1] ||
                 v <= knots[ii][jj-1] || v >= knots[ii][jj+1] )
              return false;
          for ( ii = i & 0x01; ii < hole_k; ii += 2 ) {
            knots[ii][j] = u;
            if ( jj != 5 )
              knots[ii][jj] = v;
          }
        }
      }
    }
  } /* switch */

  RecreateDomain ();
  return true;
} /*SetKnot*/

/* ////////////////////////////////////////////////////////////////////////// */
void InitFirstPartition ( int k )
{
  int   i;
  float a;

  for ( i = 0; i < hole_k; i++ ) {
    a = 2.0*PI*(float)i/(float)hole_k;
    SetVector2f ( &firstpartition[i], 3.0*cos(a), 3.0*sin(a) );
  }
} /*InitFirstPartition*/

void RedrawFirstPartition ( boolean bright )
{
  void    *sp;
  int     i;
  point2f p0, *p;

  sp = pkv_GetScratchMemTop ();
  p = (point2f*)pkv_GetScratchMem ( hole_k*sizeof(point2f) );
  if ( !p )
    exit ( 1 );
  SetPoint2f ( p, 0.0, 0.0 );
  CameraProjectPoint2f ( &DomPPos, p, &p0 );
  for ( i = 0; i < hole_k; i++ )
    CameraProjectPoint2f ( &DomPPos, &firstpartition[i], &p[i] );

  if ( bright )
    XSetForeground ( thedisplay, thegc, c_lt_green );
  else
    XSetForeground ( thedisplay, thegc, c_blue );
  for ( i = 0; i < hole_k; i++ )
    MyDrawLine ( &p0, &p[i] );
  MyMarkRect ( &p0 );
  if ( bright )
    XSetForeground ( thedisplay, thegc, c_lt_yellow );
  for ( i = 0; i < hole_k; i++ )
    MyMarkRect ( &p[i] );

  pkv_SetScratchMemTop ( sp );
} /*RedrawFirstPartition*/

int FindNearestFirstPartitionVector ( int x, int y )
{
#define TOL 10.0
  int     i, j;
  point2f p;
  float   d, e;

  d = TOL+1.0;
  j = -1;
  for ( i = 0; i < hole_k; i++ ) {
    CameraProjectPoint2f ( &DomPPos, &firstpartition[i], &p );
    e = fabs(p.x-(float)x)+fabs(p.y-(float)y);
    if ( e < d ) {
      d = e;
      j = i;
    }
  }
  if ( d > TOL )
    j = -1;
  return j;
#undef TOL
} /*FindNearestFirstPartitionVector*/

boolean SetFirstPartitionVector ( int nfp, int x, int y )
{
  int     i;
  point2f p, q;
  float   d, a0, a1, a2;

  SetPoint2f ( &p, x, y );
  CameraUnProjectPoint2f ( &DomPPos, &p, &q );
  d = sqrt ( DotProduct2f ( &q, &q ) );
  if ( d == 0.0 )
    return false;
  if ( d < 1.0 )
    NormalizeVector2f ( &q );
  else if ( d > 5.0 ) {
    NormalizeVector2f ( &q );
    MultVector2f ( 5.0, &q, &q );
  }
  a1 = atan2 ( q.y, q.x );
  i = (nfp+1) % hole_k;
  a2 = atan2 ( firstpartition[i].y, firstpartition[i].x );
  i = (nfp+hole_k-1) % hole_k;
  a0 = atan2 ( firstpartition[i].y, firstpartition[i].x );
  a0 = TrimAnglef ( a1-a0 );
  a1 = TrimAnglef ( a2-a1 );
  swDisplayFinalPatches   = FinalSurfValid     = false;
  swDisplayNLFinalPatches = NLFinalSurfValid   = false;
  PictureIsOn = false;
  if ( a0 > 0.0 && a0 < PI && a1 > 0.0 && a1 < PI ) {
    firstpartition[nfp] = q;
    return true;
  }
  else
    return false;
} /*SetFirstPartitionVector*/

/* ////////////////////////////////////////////////////////////////////////// */
int FindNearestCPointPart ( int x, int y )
{
#define TOL 10.0
  int i, j;
  float dmin, d;

  j = -1;
  d = dmin = fabs(central_pointi.x-x)+fabs(central_pointi.y-y);
  if ( d <= TOL ) j = 0;
  for ( i = 0; i < hole_k; i++ ) {
    d = fabs(domcderi[i].x-x)+fabs(domcderi[i].y-y);
    if ( d < dmin ) {
      dmin = d;
      j = i+1;
    }
  }
  return j;
#undef TOL
} /*FindNearestCPointPart*/

boolean SetCPointPart ( int nfp, int x, int y )
{
  point2f q, r;

/* !!! there is no validation of the partition of the full angle !!! */
  SetPoint2f ( &q, (float)x, (float)y );
  CameraUnProjectPoint2f ( &DomPPos, &q, &r );
  if ( nfp > 0 ) {
    nfp--;
    domcderi[nfp] = q;
    SubtractPoints2f ( &r, &central_point, &domcder[nfp] );
  }
  else {
    central_point = r;
  }
  swDisplayFinalPatches   = FinalSurfValid     = false;
  swDisplayNLFinalPatches = NLFinalSurfValid   = false;
  PictureIsOn = false;
  return true;
} /*SetCPointPart*/

void RedrawDomainCentralPoint ()
{
  int i;

  if ( !CentralPointValid )
    return;
  ProjectDomainCentralPoint ();
  XSetForeground ( thedisplay, thegc, c_cyan );
  for ( i = 0; i < hole_k; i++ )
    MyDrawLine ( &central_pointi, &domcderi[i] );
  XSetForeground ( thedisplay, thegc, c_lt_red );
  MyMarkRect ( &central_pointi );
  if ( !swUseDerivatives1 )
    XSetForeground ( thedisplay, thegc, c_cyan );
  for ( i = 0; i < hole_k; i++ )
    MyMarkRect ( &domcderi[i] );
} /*RedrawCentralPoint*/

/* ////////////////////////////////////////////////////////////////////////// */
void RedrawSurfaceCP ( int id, boolean bright )
{
  void *sp;
  int  i, j, k;
  int  *ind;

  sp = pkv_GetScratchMemTop ();
  ind = pkv_GetScratchMem ( 16*sizeof(int) );
  if ( !ind )
    return;
  if ( bright )
    XSetForeground ( thedisplay, thegc, c_green );
  else
    XSetForeground ( thedisplay, thegc, c_blue );
  for ( i = 0; i < hole_k; i++ ) {
    gh_GetBspInd ( hole_k, i, 0, ind );
    for ( j = 0; j < 4; j++ )
      for ( k = 0; k < 3; k++ )
        MyDrawLine ( &surfcpi[id][ind[4*j+k]], &surfcpi[id][ind[4*j+k+1]] );
    for ( j = 0; j < 3; j++ )
      for ( k = 0; k < 3; k++ )
        MyDrawLine ( &surfcpi[id][ind[4*j+k]], &surfcpi[id][ind[4*(j+1)+k]] );
  }
  if ( bright )
    XSetForeground ( thedisplay, thegc, c_yellow );
  for ( i = 0; i < hole_np; i++ )
    MyMarkRect ( &surfcpi[id][i] );
  pkv_SetScratchMemTop ( sp );
} /*RedrawSurfaceCP*/

static float *GetKnotSequencef ( int i )
{
  if ( i < 0 ) i += hole_k;
  else if ( i >= hole_k ) i -= hole_k;

  return knots[i];
} /*GetKnotSequencef*/

void RedrawSurface ( int id )
{
  void    *sp;
  int     i, j, k;
  int     *ind;
  point3f *cp, *cq;
  float   *ukn, *vkn;

  sp = pkv_GetScratchMemTop ();
  ind = pkv_GetScratchMem ( 16*sizeof(int) );
  cp = pkv_GetScratchMem ( 32*sizeof(point3f) );
  if ( ind && cp ) {
    cq = &cp[16];
    for ( i = 0; i < hole_k; i++ )
      for ( j = 0; j < 3; j++ ) {
        ukn = GetKnotSequencef ( i-1 );  ukn += 3;
        vkn = GetKnotSequencef ( i );    vkn += j;
        gh_GetBspInd ( hole_k, i, j, ind );
        for ( k = 0; k < 16; k++ )
          cp[k] = surfcp[ind[k]];
        mbs_BSPatchToBezf ( 3, 3, 7, ukn, 3, 7, vkn, 12, (float*)cp,
                      NULL, NULL, NULL, NULL, NULL, NULL, 12, (float*)cq );
        if ( id < 3 )
          DrawBezPatch3b ( &PPos[id], 3, 3, cq, 8, 8, c_lt_grey, c_dk_grey );
        else
          DrawBezPatch3a ( &CPos, 3, 3, cq, 8, 8, c_lt_grey, c_dk_grey );
      }
  }
  pkv_SetScratchMemTop ( sp );
} /*RedrawSurface*/

void RedrawFinalPatches ( int id )
{
  int i;

  if ( !FinalSurfValid )
    FinalSurfValid = UpdateFinalSurf ();
  if ( FinalSurfValid ) {
    if ( id < 3 ) {
      for ( i = 0; i < hole_k; i++ )
        DrawBezPatch3b ( &PPos[id], finaldegu, finaldegv, FinalCP[i],
                         8, 8, c_white, c_lt_grey );
    }
    else {
      for ( i = 0; i < hole_k; i++ )
        DrawBezPatch3a ( &CPos, finaldegu, finaldegv, FinalCP[i],
                         8, 8, c_white, c_lt_grey );
    }
  }
} /*RedrawFinalPatches*/

void RedrawNLFinalPatches ( int id )
{
  int i;

  if ( NLFinalSurfValid ) {
    if ( id < 3 ) {
      for ( i = hole_k; i < 2*hole_k; i++ )
        DrawBezPatch3b ( &PPos[id], finaldegu, finaldegv, FinalCP[i],
                         8, 8, c_lt_yellow, c_dk_yellow );
    }
    else {
      for ( i = hole_k; i < 2*hole_k; i++ )
        DrawBezPatch3a ( &CPos, finaldegu, finaldegv, FinalCP[i],
                         8, 8, c_lt_yellow, c_dk_yellow );
    }
  }
} /*RedrawNLFinalPatches*/

void RedrawSurfNumbers ( int id )
{
  void    *sp;
  int     i, k, *ind;
  point3f *cp, *cq;
  float   *ukn, *vkn;
  char    s[8];

  sp = pkv_GetScratchMemTop ();
  ind = pkv_GetScratchMem ( 16*sizeof(int) );
  cp = pkv_GetScratchMem ( 32*sizeof(point3f) );
  if ( ind && cp ) {
    XSetForeground ( thedisplay, thegc, c_cyan );
    cq = &cp[16];
    for ( i = 0; i < hole_k; i++ ) {
      ukn = GetKnotSequencef ( i-1 ) + 3;
      vkn = GetKnotSequencef ( i ) + 2;
      gh_GetBspInd ( hole_k, i, 2, ind );
      for ( k = 0; k < 16; k++ )
        cp[k] = surfcp[ind[k]];
      mbs_BSPatchToBezf ( 3, 3, 7, ukn, 3, 7, vkn, 12, (float*)cp,
                          NULL, NULL, NULL, NULL, NULL, NULL, 12, (float*)cq );
      if ( id < 3 )
        CameraProjectPoint3f ( &PPos[id], &cq[0], &cp[0] );
      else
        CameraProjectPoint3f ( &CPos, &cq[0], &cp[0] );
      sprintf ( s, "%d", i );
      XDrawString ( thedisplay, thepixmap, thegc,
                    (int)cp[0].x-2, (int)cp[0].y+5, s, strlen(s) );
    }
  }
  pkv_SetScratchMemTop ( sp );
} /*RedrawSurfNumbers*/

void SetCPoint ( int id, int np, int x, int y )
{
  point3f r, s;

  if ( eswEdRendering )
    SetRendPoint ( id, np, x, y );
  else {
    if ( swConstraintsOn ) {
      if ( np == 4*hole_k+1 ) {
        AddVector3f ( &constrcp[0], &ConstraintNormal, &s );
        CameraProjectPoint3f ( &PPos[id], &s, &r );
        r.x = (float)x;  r.y = (float)y;
        CameraUnProjectPoint3f ( &PPos[id], &r, &s );
        SubtractPoints3f ( &s, &constrcp[0], &ConstraintNormal );
        NConstrMatrixValid = ExtNConstrMatrixValid = false;
      }
      else {
        CameraProjectPoint3f ( &PPos[id], &constrcp[np], &r );
        r.x = (float)x;  r.y = (float)y;
        CameraUnProjectPoint3f ( &PPos[id], &r, &s );
        UpdateDerConstraintsHandle ( np, &s );
      }
      for ( id = 0; id < 4; id++ )
        ProjectConstraintsHandle ( id );
    }
    else {
      CameraProjectPoint3f ( &PPos[id], &surfcp[np], &r );
      r.x = (float)x;  r.y = (float)y;
      CameraUnProjectPoint3f ( &PPos[id], &r, &surfcp[np] );
      ProjectSurfCP ();
    }
    FinalSurfValid = NLFinalSurfValid = PictureIsOn = false;
    swDisplayNLFinalPatches = false;
  }
} /*SetCPoint*/

int FindNearestPoint ( int id, int x, int y )
{
#define TOL 10.0
  int   i, j;
  float d, e;

  if ( eswEdRendering )
    return FindNearestRendPoint ( id, x, y );

  d = TOL+1.0;
  j = -1;
  if ( swConstraintsOn ) {
    if ( swNormalConstr ) {
      d = fabs((float)x-constrnvi[id].x) + fabs((float)y-constrnvi[id].y);
      j = 4*hole_k+1;
    }
    for ( i = 0; i <= 4*hole_k; i++ ) {
      if ( swConstraint[i] ) {
        e = fabs((float)x-constrcpi[id][i].x) + fabs((float)y-constrcpi[id][i].y);
        if ( e < d ) {
          j = i;
          d = e;
        }
      }
    }
  }
  else {
    for ( i = 0; i < hole_np; i++ ) {
      e = fabs((float)x-surfcpi[id][i].x) + fabs((float)y-surfcpi[id][i].y);
      if ( e < d ) {
        j = i;
        d = e;
      }
    }
  }
  if ( d > TOL )
    j = -1;
  return j;
#undef TOL
} /*FindNearestPoint*/

/* ////////////////////////////////////////////////////////////////////////// */
void InitSurface ( int k )
{
  surfparams[0] = surfparams[1] = surfparams[2] = 0.0;
  domparams[0] = domparams[1] = domparams[2] = 0.0;

  InitKnots ( k, &knots[0][0] );
  hole_np = InitHole ( k, surfparams, surfcp );
  InitFirstPartition ( k );
  InitDomain ( k, firstpartition, domparams, domcp );
  CentralPointValid = false;
  RecreateDomain ();
} /*InitSurface*/

void FlattenSurface ()
{
  int i;

  for ( i = 0; i < hole_np; i++ )
    surfcp[i].z = 0.0;
  swDisplayFinalPatches   = FinalSurfValid     = false;
  swDisplayNLFinalPatches = NLFinalSurfValid   = false;
  PictureIsOn = false;
} /*FlattenSurface*/

void FitSurfaceToDomain ()
{
  pkv_Selectf ( hole_np, 2, 2, 3, (float*)domcp, (float*)surfcp );
  swDisplayFinalPatches   = FinalSurfValid     = false;
  swDisplayNLFinalPatches = NLFinalSurfValid   = false;
  PictureIsOn = false;  
} /*FitSurfaceToDomain*/

/* ////////////////////////////////////////////////////////////////////////// */
void FitLSQPlanef ( int n, const point3f *p, point3f *p0, vector3f *nv )
{
#define MAXITER 30
#define TOL  1.0e-7
  void     *sp;
  vector3f v, w;
  float    *a, *aa, /*l,*/ m;
  int      i;

  sp = pkv_GetScratchMemTop ();
  a = pkv_GetScratchMemf ( 15 );
  if ( !a )
    exit ( 1 );
  aa = &a[9];

  memcpy ( p0, p, sizeof(point3f) );
  for ( i = 1; i < n; i++ )
    AddVector3f ( p0, &p[i], p0 );
  MultVector3f ( 1.0/(float)n, p0, p0 );

  memset ( a, 0, 9*sizeof(float) );
  for ( i = 0; i < n; i++ ) {
    SubtractPoints3f ( &p[i], p0, &v );
    a[0] += v.x*v.x;
    a[3] += v.x*v.y;  a[4] += v.y*v.y;
    a[6] += v.x*v.z;  a[7] += v.y*v.z;  a[8] += v.z*v.z;
  }
  a[1] = a[3];  a[2] = a[6];  a[5] = a[7];
  a[0] += TOL;  a[4] += TOL;  a[8] += TOL;

  pkn_QRDecomposeMatrixf ( 3, 3, a, aa );
  SetVector3f ( nv, 0.3, -0.4, 0.5 );
  NormalizeVector3f ( nv );
  for ( i = 0; i < MAXITER; i++ ) {
    w = v = *nv;
    pkn_multiReflectVectorf ( 3, 3, a, aa, 1, 1, (float*)&w );
    pkn_multiMultInvUTVectorf ( 3, a, 1, 1, (float*)&w, 1, (float*)nv );
/*
    l = DotProduct3f ( &v, nv );
*/
    NormalizeVector3f ( nv );
    m = fabs(nv->x-v.x) + fabs(nv->y-v.y) + fabs(nv->z-v.z);
    if ( m < TOL )
      break;
/*
printf ( "l = %f,  m = %e,  nv = (%f,%f,%f)\n", l, m, nv->x, nv->y, nv->z );
*/
  }

  pkv_SetScratchMemTop ( sp );
#undef TOL
#undef MAXITER
} /*FitLSQPlanef*/

void ProjectPointsOntoPlanef ( const point3f *p0, const vector3f *nv,
                               int n, const point3f *p, point2f *q )
{
  int      i;
  vector3f w, v;
  float    g, h;

  w = *nv;
  NormalizeVector3f ( &w );
  if ( w.z > 0.0 )
    w.z += 1.0;
  else
    w.z -= 1.0;
  g = -2.0/DotProduct3f ( &w, &w );

  for ( i = 0; i < n; i++ ) {
    SubtractPoints3f ( &p[i], p0, &v );
    h = DotProduct3f ( &w, &v );
    AddVector3Mf ( &v, &w, g*h, &v );
    SetPoint2f ( &q[i], v.x, v.y );
/*
printf ( "q = (%f,%f)\n", q[i].x, q[i].y );
*/
  }
} /*ProjectPointsOntoPlanef*/

void FitDomainToSurface ()
{
  void *sp;
  point3f  p0;
  vector3f nv;
  unsigned char *bfcpn;
  point3f  *bfcp;
  int      i, n;

  sp = pkv_GetScratchMemTop ();

  n = 6*hole_k+1;
  bfcpn = pkv_GetScratchMem ( n*sizeof(unsigned char) );
  bfcp = pkv_GetScratchMem ( n*sizeof(point3f) );
  if ( !bfcpn || !bfcp )
    exit ( 1 );
  gh_DrawBFcpn ( hole_k, bfcpn );
  for ( i = 0; i < n; i++ )
    bfcp[i] = surfcp[bfcpn[i]];
  /* instead of fitting the plane to the surface control points, */
  /* and old option commented out below, the plane normal vector */
  /* is computed with the method used for the nonlinear optimization */
/*
  FitLSQPlanef ( n, bfcp, &p0, &nv );
*/
  SetPoint3f ( &p0, 0.0, 0.0, 0.0 );
  g2h_ComputeNLNormalf ( domain, surfcp, &nv );

  ProjectPointsOntoPlanef ( &p0, &nv, hole_np, surfcp, domcp );
  CentralPointValid = false;
  RecreateDomain ();
  ProjectDomainCP ();

  pkv_SetScratchMemTop ( sp );
} /*FitDomainToSurface*/

/* ////////////////////////////////////////////////////////////////////////// */
void ProjectConstraintsHandle ( int id )
{
  int     i;
  point3f q, r;

  if ( id < 3 ) {
    for ( i = 0; i <= 4*hole_k; i++ ) {
      CameraProjectPoint3f ( &PPos[id], &constrcp[i], &q );
      memcpy ( &constrcpi[id][i], &q, sizeof(point2f) );
    }
    AddVector3f ( &constrcp[0], &ConstraintNormal, &q );
    CameraProjectPoint3f ( &PPos[id], &q, &r );
    memcpy ( &constrnvi[id], &r, sizeof(point2f) );
  }
  else {
    for ( i = 0; i <= 4*hole_k; i++ ) {
      CameraProjectPoint3f ( &CPos, &constrcp[i], &q );
      memcpy ( &constrcpi[3][i], &q, sizeof(point2f) );
    }
    AddVector3f ( &constrcp[0], &ConstraintNormal, &q );
    CameraProjectPoint3f ( &CPos, &q, &r );
    memcpy ( &constrnvi[3], &r, sizeof(point2f) );
  }
} /*ProjectConstraintsHandle*/

void GetConstraintsHandle ()
{
#define NFCP ((G2H_FINALDEG+1)*(G2H_FINALDEG+1))
  void    *sp;
  point3f *cp, *fcp, v1, v2;
  int     i, j, k, l, m, n;

  sp = pkv_GetScratchMemTop ();
  if ( FinalSurfValid ) {
    cp = (point3f*)pkv_GetScratchMem ( 5*sizeof(point3f) );
    if ( !cp )
      goto failure;

    if ( FinalSurfValid )
      fcp = &FinalCP[0][0];
    else
      fcp = &FinalCP[hole_k][0];
    constrcp[0] = constrder[0] = fcp[0];
    for ( i = 0; i < hole_k; i++ ) {
      pkv_Selectc ( 5, sizeof(point3f), (G2H_FINALDEG+1)*sizeof(point3f),
                    sizeof(point3f), (char*)&fcp[i*NFCP], (char*)cp );

            /* rescale derivatives */
      for ( j = 1; j <= 4; j++ )
        for ( k = 4; k >= j; k-- )
          SubtractPoints3f ( &cp[k], &cp[k-1], &cp[k] );
      
      for ( j = 1, k = l = G2H_FINALDEG, m = n = 4;
            j <= 4;
            j++, l *= --k, n *= --m ) {
        MultVector3f ( (double)l, &cp[j], &constrder[4*i+j] );
        MultVector3f ( (double)l/(double)n, &cp[j], &cp[j] );
      }
      for ( j = 4; j >= 1; j-- )
        for ( k = j; k <= 4; k++ )
          AddVector3f ( &cp[k], &cp[k-1], &cp[k] );

      memcpy ( &constrcp[4*i+1], &cp[1], 4*sizeof(point3f) );
    }
    SubtractPoints3f ( &fcp[1], &fcp[0], &v1 );
    SubtractPoints3f ( &fcp[NFCP+1], &fcp[0], &v2 );
    CrossProduct3f ( &v1, &v2, &ConstraintNormal );
    NormalizeVector3f ( &ConstraintNormal );
    ConstraintsValid = true;
  }
  else {
failure:
    ConstraintsValid = false;
  }
  pkv_SetScratchMemTop ( sp );
#undef NFCP
} /*GetConstraintHandle*/

static int ConstrHLK ( int i )
{
  int j, k;

  for ( j = 1, k = 0;  j <= 4;  j++ )
    if ( swConstraint[4*i+j] )
       k = j;
  return k;
} /*ConstrHLK*/

static void SetConstrHColour ( int i )
{
  if ( !swConstraint[i] )
    XSetForeground ( thedisplay, thegc, c_cyan );
  else if ( !swZConstraint[i] )
    XSetForeground ( thedisplay, thegc, c_lt_red );
  else
    XSetForeground ( thedisplay, thegc, c_dk_red );
} /*SetConstrHColour*/

void RedrawConstraintsHandle ( int id, boolean bright )
{
  int i, j, k;

  if ( ConstraintsValid ) {
    if ( bright )
      XSetForeground ( thedisplay, thegc, c_cyan );
    else
      XSetForeground ( thedisplay, thegc, c_blue );
    ProjectConstraintsHandle ( id );
    for ( i = 0; i < hole_k; i++ ) {
      if ( (k = ConstrHLK ( i )) ) {
        MyDrawLine ( &constrcpi[id][0], &constrcpi[id][4*i+1] );
        for ( j = 1; j < k; j++ )
          MyDrawLine ( &constrcpi[id][4*i+j], &constrcpi[id][4*i+j+1] );
      }
    }
    SetConstrHColour ( 0 );
    MyMarkRect ( &constrcpi[id][0] );
    for ( i = 0; i < hole_k; i++ ) {
      if ( (k = ConstrHLK ( i )) )
        for ( j = 1; j <= k; j++ ) {
          SetConstrHColour ( 4*i+j );
          MyMarkRect ( &constrcpi[id][4*i+j] );
        }
    }
    if ( swNormalConstr ) {
      if ( bright )
        XSetForeground ( thedisplay, thegc, c_lt_green );
      MyDrawLine ( &constrcpi[id][0], &constrnvi[id] );
      if ( bright )
        XSetForeground ( thedisplay, thegc, c_lt_red );
      MyMarkRect ( &constrnvi[id] );
    }
  }
} /*RedrawConstraintHandle*/

int GetConstraintsNumber ()
{
  int i, n;

  n = 0;
  for ( i = 0; i <= 4*hole_k; i++ )
    if ( swConstraint[i] )
      n++;
  return n;
} /*GetConstraintsNumber*/

boolean SetupConstraintsMatrix ()
{
  void  *sp;
  float *constr, *drval;
  int   dimV0;
  int   i, j, s, ncn;

  if ( !(nconstr = GetConstraintsNumber ()) )
    return true;

  sp = pkv_GetScratchMemTop ();
  if ( !(dimV0 = g2h_V0SpaceDimf ( domain )) )
    goto failure;

  constr = pkv_GetScratchMemf ( nconstr*dimV0 );
  drval  = pkv_GetScratchMemf ( 5*dimV0 );
  if ( !constr || !drval )
    goto failure;

  ncn = 0;
  s   = 1;
  for ( i = 0; i < hole_k; i++ ) {
    if ( !g2h_GetBPDerivativesf ( domain, i, drval ) )
      goto failure;
    if ( i == 0 && swConstraint[0] ) {
      pkv_Selectf ( dimV0, 1, 5, 1, drval, constr );
      ncn = 1;
    }
    for ( j = 1;  j <= 4;  j++, s++ )
      if ( swConstraint[s] ) {
        pkv_Selectf ( dimV0, 1, 5, 1, &drval[j], &constr[ncn*dimV0] );
        ncn++;
      }
  }
  pkv_FreeScratchMemf ( 5*dimV0 );
  if ( !g2h_SetConstraintMatrixf ( domain, nconstr, constr ) )
    goto failure;

  pkv_SetScratchMemTop ( sp );
  ConstrMatrixValid = true;
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  ConstrMatrixValid = false;
  return false;
} /*SetupConstraintsMatrix*/

boolean SetupExtConstraintsMatrix ()
{
  void  *sp;
  float *constr, *drval;
  int   dimV0, dimExtV0, dd;
  int   i, j, s, ncn;

  if ( !(nconstr = GetConstraintsNumber ()) )
    return true;

  sp = pkv_GetScratchMemTop ();
  if ( !(dimV0 = g2h_V0SpaceDimf ( domain )) )
    goto failure;
  dimExtV0 = g2h_ExtV0SpaceDimf ( domain );
  dd = dimExtV0-dimV0;

  constr = pkv_GetScratchMemf ( nconstr*dimExtV0 );
  drval  = pkv_GetScratchMemf ( 5*dimV0 );
  if ( !constr || !drval )
    goto failure;

  memset ( constr, 0, nconstr*dimExtV0*sizeof(float) );
  ncn = 0;
  s   = 1;
  for ( i = 0; i < hole_k; i++ ) {
    if ( !g2h_GetBPDerivativesf ( domain, i, drval ) )
      goto failure;
    if ( i == 0 && swConstraint[0] ) {
      pkv_Selectf ( dimV0, 1, 5, 1, drval, &constr[dd] );
      ncn = 1;
    }
    for ( j = 1;  j <= 4;  j++, s++ )
      if ( swConstraint[s] ) {
        pkv_Selectf ( dimV0, 1, 5, 1, &drval[j], &constr[ncn*dimExtV0+dd] );
        ncn++;
      }
  }
  pkv_FreeScratchMemf ( 5*dimV0 );
  if ( !g2h_SetExtConstraintMatrixf ( domain, nconstr, constr ) )
    goto failure;

  pkv_SetScratchMemTop ( sp );
  ExtConstrMatrixValid = true;
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  ExtConstrMatrixValid = false;
  return false;
} /*SetupExtConstraintsMatrix*/

void SetupConstraintsRHS ( point3f *rhs )
{
  vector3f cder[5];
  int      i, j, k, l, s, ncn;

  cder[0] = constrcp[0];
  if ( swConstraint[0] ) {
    rhs[0] = cder[0];
    ncn = 1;
  }
  else
    ncn = 0;
  s = 1;
  for ( i = 0; i < hole_k; i++ ) {
        /* compute the derivatives */
    memcpy ( &cder[1], &constrcp[4*i+1], 4*sizeof(point3f) );
    for ( j = 1; j <= 4; j++ )
      for ( k = 4; k >= j; k-- )
        SubtractPoints3f ( &cder[k], &cder[k-1], &cder[k] );
    for ( j = 1, k = l = 4;  j <= 4;  j++, k *= --l )
      MultVector3f ( (float)k, &cder[j], &cder[j] );
        /* now select the constraints */
    for ( j = 1;  j <= 4;  j++, s++ )
      if ( swConstraint[s] ) {
        if ( j >= 2 && swZConstraint[s] )
          SetVector3f ( &rhs[ncn], 0.0, 0.0, 0.0 );
        else
          rhs[ncn] = cder[j];
        ncn++;
      }
  }
} /*SetupConstraintsRHS*/

boolean SetupNConstraintsMatrix ()
{
  void  *sp;
  float *constr, *drval;
  int   dimV0;
  int   i, j, s, ncn;

  if ( !(nconstr = GetConstraintsNumber ()) )
    return true;

  sp = pkv_GetScratchMemTop ();
  if ( !(dimV0 = g2h_V0SpaceDimf ( domain )) )
    goto failure;

  constr = pkv_GetScratchMemf ( 3*nconstr*dimV0 );
  drval  = pkv_GetScratchMemf ( 5*dimV0 );
  if ( !constr || !drval )
    goto failure;

  ncn = 0;
  s = 1;
  for ( i = 0; i < hole_k; i++ ) {
    if ( !g2h_GetBPDerivativesf ( domain, i, drval ) )
      goto failure;
    if ( i == 0 && swConstraint[0] ) {
      pkv_Selectf ( dimV0, 1, 5, 1, drval, constr );
      ncn = 1;
    }
    for ( j = 1;  j <= 4;  j++, s++ )
      if ( swConstraint[s] ) {
        pkv_Selectf ( dimV0, 1, 5, 1, &drval[j], &constr[3*ncn*dimV0] );
        ncn++;
      }
  }
  pkv_FreeScratchMemf ( 5*dimV0 );
  pkv_Movef ( nconstr, dimV0, 3*dimV0, dimV0, constr );
  pkv_Movef ( nconstr, dimV0, 3*dimV0, 2*dimV0, constr );
  pkn_MultMatrixNumf ( nconstr, dimV0, 3*dimV0, constr,
                       ConstraintNormal.x, 3*dimV0, constr );
  pkn_MultMatrixNumf ( nconstr, dimV0, 3*dimV0, &constr[dimV0],
                       ConstraintNormal.y, 3*dimV0, &constr[dimV0] );
  pkn_MultMatrixNumf ( nconstr, dimV0, 3*dimV0, &constr[2*dimV0],
                       ConstraintNormal.z, 3*dimV0, &constr[2*dimV0] );
  if ( !g2h_SetAltConstraintMatrixf ( domain, 3, nconstr, constr ) )
    goto failure;

  pkv_SetScratchMemTop ( sp );
  NConstrMatrixValid = true;
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  NConstrMatrixValid = false;
  return false;
} /*SetupNConstraintsMatrix*/

boolean SetupExtNConstraintsMatrix ()
{
  void  *sp;
  float *constr, *drval;
  int   dimV0, dimExtV0, dd;
  int   i, j, s, ncn;

  if ( !(nconstr = GetConstraintsNumber ()) )
    return true;

  sp = pkv_GetScratchMemTop ();
  if ( !(dimV0 = g2h_V0SpaceDimf ( domain )) )
    goto failure;
  dimExtV0 = g2h_ExtV0SpaceDimf ( domain );
  dd = dimExtV0-dimV0;

  constr = pkv_GetScratchMemf ( 3*nconstr*dimExtV0 );
  drval  = pkv_GetScratchMemf ( 5*dimV0 );
  if ( !constr || !drval )
    goto failure;

  memset ( constr, 0, 3*nconstr*dimExtV0*sizeof(float) );
  ncn = 0;
  s   = 1;
  for ( i = 0; i < hole_k; i++ ) {
    if ( !g2h_GetBPDerivativesf ( domain, i, drval ) )
      goto failure;
    if ( i == 0 && swConstraint[0] ) {
      pkv_Selectf ( dimV0, 1, 5, 1, drval, &constr[dd] );
      ncn = 1;
    }
    for ( j = 1;  j <= 4;  j++, s++ )
      if ( swConstraint[s] ) {
        pkv_Selectf ( dimV0, 1, 5, 1, &drval[j], &constr[3*ncn*dimExtV0+dd] );
        ncn++;
      }
  }
  pkv_FreeScratchMemf ( 5*dimV0 );
  pkv_Movef ( nconstr, dimV0, 3*dimExtV0, dimExtV0, &constr[dd] );
  pkv_Movef ( nconstr, dimV0, 3*dimExtV0, 2*dimExtV0, &constr[dd] );
  pkn_MultMatrixNumf ( nconstr, dimV0, 3*dimExtV0, &constr[dd],
                       ConstraintNormal.x, 3*dimExtV0, &constr[dd] );
  pkn_MultMatrixNumf ( nconstr, dimV0, 3*dimExtV0, &constr[dimExtV0+dd],
                       ConstraintNormal.y, 3*dimExtV0, &constr[dimExtV0+dd] );
  pkn_MultMatrixNumf ( nconstr, dimV0, 3*dimExtV0, &constr[2*dimExtV0+dd],
                       ConstraintNormal.z, 3*dimExtV0, &constr[2*dimExtV0+dd] );
  if ( !g2h_SetExtAltConstraintMatrixf ( domain, 3, nconstr, constr ) )
    goto failure;

  pkv_SetScratchMemTop ( sp );
  ExtNConstrMatrixValid = true;
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  ExtNConstrMatrixValid = false;
  return false;
} /*SetupExtNConstraintsMatrix*/

void SetupNConstraintsRHS ( float *rhs )
{
  vector3f cder[5];
  int      i, j, k, l, s, ncn;

  cder[0] = constrcp[0];
  if ( swConstraint[0] ) {
    rhs[0] = DotProduct3f ( &cder[0], &ConstraintNormal );
    ncn = 1;
  }
  else
    ncn = 0;
  s = 1;
  for ( i = 0; i < hole_k; i++ ) {
        /* compute the derivatives */
    memcpy ( &cder[1], &constrcp[4*i+1], 4*sizeof(point3f) );
    for ( j = 1; j <= 4; j++ )
      for ( k = 4; k >= j; k-- )
        SubtractPoints3f ( &cder[k], &cder[k-1], &cder[k] );
    for ( j = 1, k = l = 4;  j <= 4;  j++, k *= --l )
      MultVector3f ( (float)k, &cder[j], &cder[j] );
        /* now select the constrianits */
    for ( j = 1;  j <= 4;  j++, s++ )
      if ( swConstraint[s] ) {
        if ( swZConstraint[s] )
          rhs[ncn] = 0.0;
        else
          rhs[ncn] = DotProduct3f ( &cder[j], &ConstraintNormal );
        ncn++;
      }
  }
} /*SetupNConstraintsRHS*/

void SetDefaultConstraints ()
{
  GetConstraintsHandle ();
  NConstrMatrixValid = ExtNConstrMatrixValid = false;
} /*SetDefaultConstraints*/

void TurnConstraints ()
{
  if ( swConstraintsOn ) {
    if ( !ConstraintsValid )
      GetConstraintsHandle ();
    NConstrMatrixValid = ExtNConstrMatrixValid = false;
  }
  FinalSurfValid = NLFinalSurfValid = false;
} /*TurnConstraints*/

static void InsertSWQ ( char *queue, int ql, char *qn, int nc, int nd )
{
  int i, j, cn, cnn;

      /* remove from the queue all switches turned off */
  cnn = 4*nc+nd;
  for ( i = j = 0; i < *qn; i++ ) {
    cn = 4*queue[i]+nd;
    if ( swConstraint[cn] && cn != cnn )
      queue[j++] = queue[i];
  }
  if ( j >= ql ) {
       /* turn off the first switch in the queue */
    cn = 4*queue[0]+nd;
    swConstraint[cn] = false;
    memmove ( queue, &queue[1], ql-1 );
    j = ql-1;
  }
  queue[j] = nc;
  *qn = j+1;
} /*InsertSWQ*/

void UpdateDerConstraintsHandle ( int np, point3f *hp )
{
  int      i, j, k, l, nc, nd;
  vector3f v, d[5];

  if ( !swPointMode ) {
    SubtractPoints3f ( hp, &constrcp[np], &v );
    if ( np == 0 ) {
      constrder[0] = d[0] = *hp;
      for ( i = 0; i < hole_k; i++ ) {
        memcpy ( &d[1], &constrder[4*i+1], 4*sizeof(vector3f) );
        for ( j = 1, k = l = 4;  j <= 4;  j++, l *= --k )
          MultVector3f ( 1.0/(double)l, &d[j], &d[j] );
        for ( j = 4; j >= 1; j-- )
          for ( k = j; k <= 4; k++ )
            AddVector3f ( &d[k], &d[k-1], &d[k] );
        memcpy ( &constrcp[4*i+1], &d[1], 4*sizeof(point3f) );
      }
    }
    else {
      nc = (np-1) / 4;
      nd = (np-1) % 4 + 1;
      d[0] = constrcp[0];
      memcpy ( &d[1], &constrcp[4*nc+1], 4*sizeof(vector3f) );
      d[nd] = *hp;
      for ( j = 1; j <= 4; j++ )
        for ( k = 4; k >= j; k-- )
          SubtractPoints3f ( &d[k], &d[k-1], &d[k] );
      for ( j = 1, k = l = 4;  j <= 4;  j++, l *= --k )
        MultVector3f ( (double)l, &d[j], &d[j] );
      if ( nd < 4 )
        memcpy ( &d[nd+1], &constrcp[4*nc+nd+1], (4-nd)*sizeof(point3f) );
      memcpy ( &constrcp[4*nc+1], &d[1], 4*sizeof(point3f) );
      for ( j = 1, k = l = 4;  j <= 4;  j++, l *= --k )
        MultVector3f ( 1.0/(double)l, &d[j], &d[j] );
      for ( j = 4; j >= 1; j-- )
        for ( k = j; k <= 4; k++ )
          AddVector3f ( &d[k], &d[k-1], &d[k] );
      memcpy ( &constrcp[4*nc+1], &d[1], 4*sizeof(vector3f) );
    }
  }
  constrcp[np] = *hp;
} /*UpdateDerConstraintsHandle*/

void SwitchAConstraint ( int cno )
{
  int nc, nd;

  if ( swZeroDer ) {
    if ( cno == 0 )
      swZConstraint[0] = false;
    if ( swConstraintsOn )
      FinalSurfValid = NLFinalSurfValid = false;
  }
  else {
    if ( swConstraint[cno] && cno > 0 ) {
      nc = (cno-1) / 4;
      nd = (cno-1) % 4 + 1;
      switch ( nd ) {
    case 1:
        InsertSWQ ( ConstrQueue1, 2, &constrql1, nc, 1 );
        break;
    case 2:
        InsertSWQ ( ConstrQueue2, 3, &constrql2, nc, 2 );
        break;
   default:
        break;
      }
    }
    ConstrMatrixValid = ExtConstrMatrixValid =
    NConstrMatrixValid = ExtNConstrMatrixValid = false;
    if ( swConstraintsOn ) {
      if ( swCoonsPatches ) {
        if ( !swNormalConstr )
          SetupConstraintsMatrix ();
        else
          SetupNConstraintsMatrix ();
      }
      else {
        if ( !swNormalConstr )
          SetupExtConstraintsMatrix ();
        else
          SetupExtNConstraintsMatrix ();
      }
      FinalSurfValid = NLFinalSurfValid = false;
    }
  }
} /*SwitchAConstraint*/

