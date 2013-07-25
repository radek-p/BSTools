
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
static int nccurves;
static GHoptions *constropt;

static void OutConstrCurve ( int n, const double *cp )
{
  memcpy ( &constropt->constrcp[nccurves*(n+1)], cp, (n+1)*sizeof(point3d) );
  nccurves ++;
} /*OutConstrCurve*/

static void OutSplConstrCurve ( int n, int lkn, const double *kn,
                                const double *cp )
{
  memcpy ( &constropt->constrcp[nccurves*(lkn-n)], cp,
           (lkn-n)*sizeof(point3d) );
  nccurves ++;
} /*OutSplConstrCurve*/

static void IgnorePatch ( int n, int m, const double *cp, void *usrptr )
{
} /*IgnorePatch*/

static boolean GetPatchCurves ( GHoleDomaind *domain, GHoptions *opt )
{
  void   *sp;
  double *acoeff;
  int    dimV0;

  sp = pkv_GetScratchMemTop ();
  if ( opt->order == 1 )
    dimV0 = g1h_V0SpaceDimd ( domain );
  else
    dimV0 = g2h_V0SpaceDimd ( domain );
  acoeff = pkv_GetScratchMemd ( 3*dimV0 );
  if ( !acoeff )
    goto failure;
  constropt = opt;
  nccurves = 0;
  if ( opt->order == 1 ) {
    if ( !g1h_Q2FillHoled ( domain, 3, &hole_cp[0].x, acoeff, NULL,
                            IgnorePatch ) )
      goto failure;
    if ( !g1h_GetFinalPatchCurvesd ( domain, 3, &hole_cp[0].x, acoeff,
                                     OutConstrCurve ) )
      goto failure;
  }
  else {
    if ( !g2h_FillHoled ( domain, 3, &hole_cp[0].x, acoeff, NULL,
                          IgnorePatch ) )
      goto failure;
    if ( !g2h_GetFinalPatchCurvesd ( domain, 3, &hole_cp[0].x, acoeff,
                                     OutConstrCurve ) )
      goto failure;
  }

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*GetPatchCurves*/

void InitConstraintFrame ( char surfno )
{
  void     *sp;
  GHoleDomaind *domain;
  GHoptions    *opt;
  int       i, j, ncp0, ncp1;
  point3d   bcp[16], *constrcp;
  vector3d  v1, v2;
  double    *knots, t;
  int       lkn;

  sp = pkv_GetScratchMemTop ();
  if ( surfno == 1 ) {
    domain = domain1;
    opt = &options1;
  }
  else {
    domain = domain2;
    opt = &options2;
  }
  if ( opt->order == 1 )
    ncp0 = ncp1 = 5;  /* degree is 4 */
  else
    ncp0 = ncp1 = 8;  /* degree is 7 */
  if ( opt->spline )
    ncp1 += opt->nk*opt->m1;
  opt->ncp = ncp1;
  lkn = ncp0+ncp1-1;
  if ( GetPatchCurves ( domain, opt ) ) {
    if ( ncp1 > ncp0 ) {
      pkv_Rearranged ( hole_k, 3*ncp0, 3*ncp0, 3*ncp1, opt->constrcp );
      if ( !(knots = pkv_GetScratchMemd ( lkn+1 )) )
        goto failure;
      for ( i = 0; i < ncp0; i++ )
        { knots[i] = 0.0;  knots[ncp0+i] = 1.0; }
      lkn = 2*ncp0-1;
      for ( i = 0; i < opt->nk; i++ ) {
        t = (double)(i+1)/(double)(opt->nk+1);
        for ( j = 0; j < opt->m1; j++ )
          mbs_multiKnotInsd ( ncp0-1, &lkn, knots, hole_k, 3, 3*ncp1, 3*ncp1,
                              &opt->constrcp[0].x, t );
      }
    }
  }
  else {
failure:
    constrcp = opt->constrcp;
    memset ( constrcp, 0, sizeof(point3d) );
    for ( i = 0; i < hole_k; i++ ) {
      GetSurfBezPatch ( i, 2, bcp );
      AddVector3d ( &constrcp[0], &bcp[0], &constrcp[0] );
      constrcp[(i+1)*ncp1-1] = bcp[0];
    }
    MultVector3d ( 1.0/(double)hole_k, &constrcp[0], &constrcp[0] );
    for ( i = 1; i < hole_k; i++ )
      constrcp[i*ncp1] = constrcp[0];
    for ( i = 0; i < hole_k; i++ )
      for ( j = 1; j < ncp1-1; j++ )
        InterPoint3d ( &constrcp[i*ncp1], &constrcp[(i+1)*ncp1-1],
                       (double)j/(double)(ncp1-1), &constrcp[i*ncp1+j] );
  }
        /* setup the "normal" vector for the constraints */
  SubtractPoints3d ( &opt->constrcp[1], &opt->constrcp[0], &v1 );
  SubtractPoints3d ( &opt->constrcp[ncp1+1], &opt->constrcp[0], &v2 );
  CrossProduct3d ( &v1, &v2, &opt->constrnv );
  NormalizeVector3d ( &opt->constrnv );

  constraints1st = constraints2nd = false;
  pkv_SetScratchMemTop ( sp );
} /*InitConstraintFrame*/

boolean GetCurrentConstraintFrame ( char surfno )
{
  GHoleDomaind *domain;
  GHoptions    *opt;
  double       *acoeff;
  vector3d     v1, v2;

  if ( surfno == 1 ) {
    domain = domain1;
    opt = &options1;
    acoeff = acoeff1;
  }
  else {
    domain = domain2;
    opt = &options2;
    acoeff = acoeff2;
  }
  if ( !acoeff )
    return false;

  constropt = opt;
  nccurves = 0;
  switch ( opt->order ) {
case 1:
    if ( opt->coons ) {
      if ( !g1h_GetFinalPatchCurvesd ( domain, 3, &hole_cp[0].x, acoeff,
                                       OutConstrCurve ) )
        return false;
    }
    else if ( opt->bezier ) {
      if ( !g1h_GetExtFinalPatchCurvesd ( domain, 3, &hole_cp[0].x, acoeff,
                                          OutConstrCurve ) )
        return false;
    }
    else if ( opt->spline ) {
      if ( !g1h_GetSplFinalPatchCurvesd ( domain, 3, &hole_cp[0].x, acoeff,
                                          OutSplConstrCurve ) )
        return false;
    }
    else
      return false;
    break;
case 2:
    if ( opt->coons ) {
      if ( !g2h_GetFinalPatchCurvesd ( domain, 3, &hole_cp[0].x, acoeff,
                                       OutConstrCurve ) )
        return false;
    }
    else if ( opt->bezier ) {
      if ( !g2h_GetExtFinalPatchCurvesd ( domain, 3, &hole_cp[0].x, acoeff,
                                          OutConstrCurve ) )
        return false;
    }
    else if ( opt->spline ) {
      if ( !g2h_GetSplFinalPatchCurvesd ( domain, 3, &hole_cp[0].x, acoeff,
                                          OutSplConstrCurve ) )
        return false;
    }
    else
      return false;
    break;
default:
    return false;
  }
        /* setup the "normal" vector for the constraints */
  SubtractPoints3d ( &opt->constrcp[1], &opt->constrcp[0], &v1 );
  SubtractPoints3d ( &opt->constrcp[opt->ncp+1], &opt->constrcp[0], &v2 );
  CrossProduct3d ( &v1, &v2, &opt->constrnv );
  NormalizeVector3d ( &opt->constrnv );

  return true;
} /*GetCurrentConstraintFrame*/

static GHoptions *GetNConstrPoints ( int *ncs, int *ncp )
{
  GHoptions *opt;

  switch ( constr_surfno ) {
case 1:
    opt = &options1;
    break;
case 2:
    opt = &options2;
    break;
default:
    return NULL;
  }
  if ( opt->order == 1 )
    *ncs = 3;
  else
    *ncs = 5;
  if ( opt->spline )
    *ncs += opt->nk*opt->m1;
  *ncp = opt->ncp;
  return opt;
} /*GetNConstrPoints*/

void DrawConstraintFrame ( int id, char surfno )
{
  void      *sp;
  GHoptions *opt;
  int       i, j, k, ncp, ncs;
  point3d   *constrcp, cp;
  point2s   *pp;

  id &= 0x03;
  ncs = ncp = 0;
  opt = GetNConstrPoints ( &ncs, &ncp );
  sp = pkv_GetScratchMemTop ();
  constrcp = opt->constrcp;
  pp = pkv_GetScratchMem ( (hole_k*ncp+1)*sizeof(point3d) );
  if ( !pp )
    goto way_out;

  for ( i = 0; i < hole_k*ncp; i++ ) {
    CameraProjectPoint3d ( &swind.CPos[id], &constrcp[i], &cp );
    pp[i].x = (short)(cp.x+0.5);  pp[i].y = (short)(cp.y+0.5);
  }
  if ( opt->constr_type == 2 ) {
    AddVector3d ( &constrcp[0], &opt->constrnv, &cp  );
    CameraProjectPoint3d ( &swind.CPos[id], &cp, &cp );
    pp[hole_k*ncp].x = (short)(cp.x+0.5);  pp[hole_k*ncp].y = (short)(cp.y+0.5);
    xgeSetForeground ( xgec_Grey1 );
    xgeDrawLine ( pp[0].x, pp[0].y, pp[hole_k*ncp].x, pp[hole_k*ncp].y );
  }
  xgeSetForeground ( xgec_Green3 );
  for ( i = k = 0;  i < hole_k;  i++, k++ )
    for ( j = 1;  j < ncp;  j++, k++ )
      xgeDrawLine ( pp[k].x, pp[k].y, pp[k+1].x, pp[k+1].y );
  xgeSetForeground ( xgec_Green2 );
  MyMarkPoint2s ( &pp[0] );
  for ( i = 0, k = 1;  i < hole_k;  i++, k++ )
    for ( j = 1;  j < ncp;  j++, k++ )
      MyMarkPoint2s ( &pp[i*ncp+j] );
  xgeSetForeground ( xgec_Yellow );
  if ( opt->constrsw[0] && !opt->mccp[0] )
    MyMarkPoint2s ( &pp[0] );
  for ( i = 0, k = 1; i < hole_k; i++ )
    for ( j = 1;  j < ncs;  j++, k++ )
      if ( opt->constrsw[k] && !opt->mccp[k] )
        MyMarkPoint2s ( &pp[i*ncp+j] );
  xgeSetForeground ( xgec_Red4 );
  if ( !opt->constrsw[0] && opt->mccp[0] )
    MyMarkPoint2s ( &pp[0] );
  for ( i = 0, k = 1; i < hole_k; i++ )
    for ( j = 1;  j < ncs;  j++, k++ )
      if ( !opt->constrsw[k] && opt->mccp[k] )
        MyMarkPoint2s ( &pp[i*ncp+j] );
  xgeSetForeground ( xgec_Red );
  if ( opt->constrsw[0] && opt->mccp[0] )
    MyMarkPoint2s ( &pp[0] );
  for ( i = 0, k = 1; i < hole_k; i++ )
    for ( j = 1;  j < ncs;  j++, k++ )
      if ( opt->constrsw[k] && opt->mccp[k] )
        MyMarkPoint2s ( &pp[i*ncp+j] );
  if ( opt->constr_type == 2 ) {
    xgeSetForeground ( xgec_Yellow );
    MyMarkPoint2s ( &pp[ncp*hole_k] );
  }

way_out:
  pkv_SetScratchMemTop ( sp );
} /*DrawConstraintFrame*/

static void ProjectConstrPoints ( GHoptions *opt, CameraRecd *CPos,
                                  int ncs, int ncp, point2d *pp )
{
  point3d q;
  int     i, j, k;

  CameraProjectPoint3d ( CPos, &opt->constrcp[0], &q );
  pp[0].x = q.x;  pp[0].y = q.y;
  for ( i = 0, k = 1;  i < hole_k;  i++ )
    for ( j = 1; j < ncs; j++, k++ ) {
      CameraProjectPoint3d ( CPos, &opt->constrcp[i*ncp+j], &q );
      pp[k].x = q.x;  pp[k].y = q.y;
    }
  AddVector3d ( &opt->constrcp[0], &opt->constrnv, &q );
  CameraProjectPoint3d ( CPos, &q, &q );
  pp[k].x = q.x;  pp[k].y = q.y;
} /*ProjectConstrPoints*/

boolean FindNearestConstraintCPoint ( int id, short x, short y, short mindist )
{
  void      *sp;
  int       i, j, ncs, ncp;
  point2d   *pp;
  GHoptions *opt;

  id &= 0x03;
  ncs = ncp = 0;
  opt = GetNConstrPoints ( &ncs, &ncp );
  sp = pkv_GetScratchMemTop ();
  pp = pkv_GetScratchMem ( ((ncs-1)*hole_k+2)*sizeof(point2d) );
  if ( !pp )
    goto failure;
  ProjectConstrPoints ( opt, &swind.CPos[id], ncs, ncp, pp );
  if ( FindNearestPoint ( hole_k*(ncs-1)+2, pp, x, y, mindist, &i ) ) {
    if ( i == (ncs-1)*hole_k+1 ) {
      if ( opt->constr_type == 2 )
        swind.current_point = -1;
      else
        goto failure;
    }
    else if ( i > 0 ) {
      i--;
      j = i % (ncs-1);
      i = i / (ncs-1);
      swind.current_point = i*ncp+j+1;
    }
    else swind.current_point = 0;
    pkv_SetScratchMemTop ( sp );
    return true;
  }
  
failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*FindNearestConstraintCPoint*/

void SetConstraintCPoint ( int id, short x, short y )
{
  GHoptions *opt;
  int       i, ncp;
  point3d   q;

  id &= 0x03;
  switch ( constr_surfno ) {
case 1:
    opt = &options1;
    break;
case 2:
    opt = &options2;
    break;
default:
    return;
  }
  swind_picture = false;
  i = swind.current_point;
  if ( i == -1 ) {  /* normal vector */
    if ( opt->constr_type == 2 ) {
      AddVector3d ( &opt->constrcp[0], &opt->constrnv, &q );
      CameraProjectPoint3d ( &swind.CPos[id], &q, &q );
      q.x = (double)x;
      q.y = (double)y;
      CameraUnProjectPoint3d ( &swind.CPos[id], &q, &q );
      SubtractPoints3d ( &q, &opt->constrcp[0], &opt->constrnv );
      NormalizeVector3d ( &opt->constrnv );
    }
  }
  else {
    CameraProjectPoint3d ( &swind.CPos[id], &opt->constrcp[i], &q );
    q.x = (double)x;
    q.y = (double)y;
    CameraUnProjectPoint3d ( &swind.CPos[id], &q, &opt->constrcp[i] );
    if ( i == 0 ) {   /* central point */
      ncp = opt->ncp;
      for ( i = 1; i < hole_k; i++ )
        opt->constrcp[i*ncp] = opt->constrcp[0];
    }
  }
} /*SetConstraintCPoint*/

void SelectConstraintCPoints ( int id )
{
  void      *sp;
  GHoptions *opt;
  int       ncs, ncp;
  point2d   *pp;

  id &= 0x03;
  ncs = ncp = 0;
  opt = GetNConstrPoints ( &ncs, &ncp );
  sp = pkv_GetScratchMemTop ();
  pp = pkv_GetScratchMem ( ((ncs-1)*hole_k+1)*sizeof(point2d) );
  if ( !pp )
    goto failure;
  ProjectConstrPoints ( opt, &swind.CPos[id], ncs, ncp, pp );
  SelectCPoints ( &swind.selection_rect, (ncs-1)*hole_k+1, pp,
                  opt->mccp, true );
failure:
  pkv_SetScratchMemTop ( sp );
} /*SelectConstraintCPoints*/

void UnselectConstraintCPoints ( int id )
{
  void      *sp;
  GHoptions *opt;
  int       ncs, ncp;
  point2d   *pp;

  id &= 0x03;
  ncs = ncp = 0;
  opt = GetNConstrPoints ( &ncs, &ncp );
  sp = pkv_GetScratchMemTop ();
  pp = pkv_GetScratchMem ( ((ncs-1)*hole_k+1)*sizeof(point2d) );
  if ( !pp )
    goto failure;
  ProjectConstrPoints ( opt, &swind.CPos[id], ncs, ncp, pp );
  SelectCPoints ( &swind.selection_rect, (ncs-1)*hole_k+1, pp,
                  opt->mccp, false );
failure:
  pkv_SetScratchMemTop ( sp );
} /*UnselectConstraintCPoints*/

void SaveConstraintCPoints ( void )
{
  GHoptions *opt;
  int       ncs, ncp;

  opt = NULL;  ncs = ncp = 0;
  opt = GetNConstrPoints ( &ncs, &ncp );
  memcpy ( saved_cp, opt->constrcp, hole_k*opt->ncp*sizeof(point3d) );
} /*SaveConstraintCPoints*/

void TransformConstraintCPoints ( void )
{
  GHoptions *opt;
  int       ncs, ncp;
  int       i, j, k;

  ncs = ncp = 0;
  opt = GetNConstrPoints ( &ncs, &ncp );
  if ( opt->mccp[0] ) {
    TransPoint3d ( &swind.gwtrans, &saved_cp[0], &opt->constrcp[0] );
    for ( i = 1; i < hole_k; i++ )
      opt->constrcp[i*ncp] = opt->constrcp[0];
  }
  for ( i = 0, k = 1;  i < hole_k;  i++ )
    for ( j = 1;  j < ncs;  j++, k++ )
      if ( opt->mccp[k] )
        TransPoint3d ( &swind.gwtrans, &saved_cp[i*ncp+j], &opt->constrcp[i*ncp+j] );
} /*TransformConstraintCPoints*/

/* ////////////////////////////////////////////////////////////////////////// */
void CountTheConstraints ( void )
{
  GHoptions *opt;
  int       i, cnt, ncs, ncp;
/*
printf ( "Constraints:\n" );
*/
  ncs = ncp = 0;
  opt = GetNConstrPoints ( &ncs, &ncp );
  for ( i = cnt = 0;  i < MAX_CONSTRAINTS;  i++ )
    if ( opt->constrsw[i] ) {
      if ( !i ) {
        opt->constr_spec[cnt][0] = opt->constr_spec[cnt][1] = 0;
      }
      else {
        opt->constr_spec[cnt][0] = (i-1) / (ncs-1);      /* curve number */
        opt->constr_spec[cnt][1] = (i-1) % (ncs-1) + 1;  /* control point number */
      }
/*
printf ( "%d %d\n", opt->constr_spec[cnt][0], opt->constr_spec[cnt][1] );
*/
      cnt ++;
    }
  opt->constr_no = (short)cnt;
/*
printf ( "constr_no: %d\n", cnt );
*/
} /*CountTheConstraints*/

boolean ComputeConstraintMatrices ( int surfno )
{
  void         *sp;
  double       *cmat, *fcp, *kncp;
  GHoleDomaind *domain;
  GHoptions    *opt;
  int          nfunc_a, nfunc_b, nfunc_c, nfunc_d, dimV0, size;
  int          constr_no, i, fn, nc, np, lkn;

  switch ( surfno ) {
case 1:
    opt    = &options1;
    domain = domain1;
    break;
case 2:
    opt    = &options2;
    domain = domain2;
    break;
default:
    return false;
  }
  sp = pkv_GetScratchMemTop ();
  CountTheConstraints ();
  if ( opt->bcmat ) {
    free ( opt->bcmat );
    opt->bcmat = NULL;
  }
  nfunc_b = 6*hole_k+1;
  constr_no = opt->constr_no;
  fcp = pkv_GetScratchMemd ( opt->ncp );
  if ( !fcp )
    goto failure;
  opt->bcmat = malloc ( nfunc_b*constr_no*sizeof(double) );
  if ( !opt->bcmat )
    goto failure;
  if ( opt->order == 1 ) {
    if ( opt->coons ) {
      dimV0 = nfunc_a = g1h_V0SpaceDimd ( domain );
      size = dimV0*constr_no;
      cmat = pkv_GetScratchMemd ( size );
      if ( !cmat )
        goto failure;
      memset ( cmat, 0, size*sizeof(double) );
      for ( i = 0;  i < constr_no;  i++ ) {
        nc = opt->constr_spec[i][0];
        np = opt->constr_spec[i][1];
        for ( fn = 0; fn < nfunc_a; fn++ ) {
          g1h_GetABasisFPatchCurved ( domain, fn, nc, fcp );
          cmat[i*dimV0+fn] = fcp[np];
        }
      }
      if ( g1h_SetConstraintMatrixd ( domain, constr_no, cmat ) )
        goto comp_poly_bcmat1;
      else
        goto failure;
    }
    else if ( opt->bezier ) {
      nfunc_a = g1h_V0SpaceDimd ( domain );
      dimV0 = g1h_ExtV0SpaceDimd ( domain );
      nfunc_c = dimV0-nfunc_a;
      size = dimV0*constr_no;
      cmat = pkv_GetScratchMemd ( size );
      if ( !cmat )
        goto failure;
      memset ( cmat, 0, size*sizeof(double) );
      for ( i = 0;  i < constr_no;  i++ ) {
        nc = opt->constr_spec[i][0];
        np = opt->constr_spec[i][1];
        for ( fn = 0; fn < nfunc_a; fn++ ) {
          g1h_GetABasisFPatchCurved ( domain, fn, nc, fcp );
          cmat[i*dimV0+nfunc_c+fn] = fcp[np];
        }
      }
      if ( !g1h_SetExtConstraintMatrixd ( domain, constr_no, cmat ) )
        goto failure;
comp_poly_bcmat1:
      for ( i = 0;  i < constr_no;  i++ ) {
        nc = opt->constr_spec[i][0];
        np = opt->constr_spec[i][1];
        for ( fn = 0; fn < nfunc_b; fn++ ) {
          g1h_GetBBasisFPatchCurved ( domain, fn, nc, fcp );
          opt->bcmat[i*nfunc_b+fn] = fcp[np];
        }
      }
    }
    else {
      g1h_DrawSplBasFuncNumd ( domain, &nfunc_a, &nfunc_b, &nfunc_c, &nfunc_d );
      dimV0 = nfunc_a+nfunc_c+nfunc_d;
      size = dimV0*constr_no;
      kncp = pkv_GetScratchMemd ( opt->ncp+G1H_OMCDEG+1 );
      cmat = pkv_GetScratchMemd ( size );
      if ( !kncp || !cmat )
        goto failure;
      memset ( cmat, 0, size*sizeof(double) );
      for ( i = 0;  i < constr_no;  i++ ) {
        nc = opt->constr_spec[i][0];
        np = opt->constr_spec[i][1];
        for ( fn = 0; fn < nfunc_a; fn++ ) {
          g1h_GetSplABasisFPatchCurved ( domain, fn, nc, &lkn, kncp, fcp );
          cmat[i*dimV0+nfunc_c+nfunc_d+fn] = fcp[np];
        }
        for ( fn = 0; fn < nfunc_d;  fn++ ) {
          g1h_GetSplDBasisFPatchCurved ( domain, fn, nc, &lkn, kncp, fcp );
          cmat[i*dimV0+nfunc_c+fn] = fcp[np];
        }
      }
      if ( !g1h_SetSplConstraintMatrixd ( domain, constr_no, cmat ) )
        goto failure;
      for ( i = 0;  i < constr_no;  i++ ) {
        nc = opt->constr_spec[i][0];
        np = opt->constr_spec[i][1];
        for ( fn = 0; fn < nfunc_b; fn++ ) {
          g1h_GetSplBBasisFPatchCurved ( domain, fn, nc, &lkn, kncp, fcp );
          opt->bcmat[i*nfunc_b+fn] = fcp[np];
        }
      }
    }
  }
  else {  /* order == 2 */
    if ( opt->coons ) {
      dimV0 = nfunc_a = g2h_V0SpaceDimd ( domain );
      size = dimV0*constr_no;
      cmat = pkv_GetScratchMemd ( size );
      if ( !cmat )
        goto failure;
      memset ( cmat, 0, size*sizeof(double) );
      for ( i = 0;  i < constr_no;  i++ ) {
        nc = opt->constr_spec[i][0];
        np = opt->constr_spec[i][1];
        for ( fn = 0; fn < nfunc_a; fn++ ) {
          g2h_GetABasisFPatchCurved ( domain, fn, nc, fcp );
          cmat[i*dimV0+fn] = fcp[np];
        }
      }
      if ( g2h_SetConstraintMatrixd ( domain, constr_no, cmat ) )
        goto comp_poly_bcmat2;
      else
        goto failure;
    }
    else if ( opt->bezier ) {
      nfunc_a = g2h_V0SpaceDimd ( domain );
      dimV0 = g2h_ExtV0SpaceDimd ( domain );
      nfunc_c = dimV0-nfunc_a;
      size = dimV0*constr_no;
      cmat = pkv_GetScratchMemd ( size );
      if ( !cmat )
        goto failure;
      memset ( cmat, 0, size*sizeof(double) );
      for ( i = 0;  i < constr_no;  i++ ) {
        nc = opt->constr_spec[i][0];
        np = opt->constr_spec[i][1];
        for ( fn = 0; fn < nfunc_a; fn++ ) {
          g2h_GetABasisFPatchCurved ( domain, fn, nc, fcp );
          cmat[i*dimV0+nfunc_c+fn] = fcp[np];
        }
      }
      if ( !g2h_SetExtConstraintMatrixd ( domain, constr_no, cmat ) )
        goto failure;
comp_poly_bcmat2:
      for ( i = 0;  i < constr_no;  i++ ) {
        nc = opt->constr_spec[i][0];
        np = opt->constr_spec[i][1];
        for ( fn = 0; fn < nfunc_b; fn++ ) {
          g2h_GetBBasisFPatchCurved ( domain, fn, nc, fcp );
          opt->bcmat[i*nfunc_b+fn] = fcp[np];
        }
      }
    }
    else {
      g2h_DrawSplBasFuncNumd ( domain, &nfunc_a, &nfunc_b, &nfunc_c, &nfunc_d );
      dimV0 = nfunc_a+nfunc_c+nfunc_d;
      size = dimV0*constr_no;
      kncp = pkv_GetScratchMemd ( opt->ncp+G2H_OMCDEG+1 );
      cmat = pkv_GetScratchMemd ( size );
      if ( !kncp || !cmat )
        goto failure;
      memset ( cmat, 0, size*sizeof(double) );
      for ( i = 0;  i < constr_no;  i++ ) {
        nc = opt->constr_spec[i][0];
        np = opt->constr_spec[i][1];
        for ( fn = 0; fn < nfunc_a; fn++ ) {
          g2h_GetSplABasisFPatchCurved ( domain, fn, nc, &lkn, kncp, fcp );
          cmat[i*dimV0+nfunc_c+nfunc_d+fn] = fcp[np];
        }
        for ( fn = 0; fn < nfunc_d;  fn++ ) {
          g2h_GetSplDBasisFPatchCurved ( domain, fn, nc, &lkn, kncp, fcp );
          cmat[i*dimV0+nfunc_c+fn] = fcp[np];
        }
      }
      if ( !g2h_SetSplConstraintMatrixd ( domain, constr_no, cmat ) )
        goto failure;
      for ( i = 0;  i < constr_no;  i++ ) {
        nc = opt->constr_spec[i][0];
        np = opt->constr_spec[i][1];
        for ( fn = 0; fn < nfunc_b; fn++ ) {
          g2h_GetSplBBasisFPatchCurved ( domain, fn, nc, &lkn, kncp, fcp );
          opt->bcmat[i*nfunc_b+fn] = fcp[np];
        }
      }
    }
  }
  opt->constr_matrix_valid = true;

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  opt->constr_matrix_valid = false;
  pkv_SetScratchMemTop ( sp );
  return false;
} /*ComputeConstraintMatrices*/

boolean ComputeConstraintRightSide ( int surfno, vector3d *rhs )
{
  void   *sp;
  GHoptions *opt;
  int    nfunc_b, constr_no;
  int    i, j, nc, np, ncp;
  unsigned char *bfcpn;

  sp = pkv_GetScratchMemTop ();
  switch ( surfno ) {
case 1:
    opt = &options1;
    break;
case 2:
    opt = &options2;
    break;
default:
    goto failure;
  }
  if ( !opt->bcmat )
    goto failure;

  nfunc_b = 6*hole_k+1;
  bfcpn = pkv_GetScratchMem ( (6*hole_k+1)*sizeof(unsigned char) );
  if ( !bfcpn )
    goto failure;
  if ( !gh_DrawBFcpn ( hole_k, bfcpn ) )
    goto failure;
  constr_no = opt->constr_no;
  ncp = opt->ncp;
  for ( i = 0; i < constr_no; i++ ) {
    nc = opt->constr_spec[i][0];
    np = opt->constr_spec[i][1];
    rhs[i] = opt->constrcp[ncp*nc+np];
    for ( j = 0; j < nfunc_b; j++ )
      AddVector3Md ( &rhs[i], &hole_cp[bfcpn[j]],
                     -opt->bcmat[nfunc_b*i+j], &rhs[i] );
  }

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*ComputeConstraintRightSide*/

boolean ComputeConstraintMatricesAlt ( int surfno )
{
  void         *sp;
  double       *cmat, *fcp, *kncp;
  GHoleDomaind *domain;
  GHoptions    *opt;
  int          nfunc_a, nfunc_b, nfunc_c, nfunc_d, dimV0, size;
  int          constr_no, i, fn, nc, np, lkn;

  switch ( surfno ) {
case 1:
    opt    = &options1;
    domain = domain1;
    break;
case 2:
    opt    = &options2;
    domain = domain2;
    break;
default:
    return false;
  }
  sp = pkv_GetScratchMemTop ();
  CountTheConstraints ();
  if ( opt->bcmat ) {
    free ( opt->bcmat );
    opt->bcmat = NULL;
  }
  nfunc_b = 6*hole_k+1;
  constr_no = opt->constr_no;
  fcp = pkv_GetScratchMemd ( opt->ncp );
  if ( !fcp )
    goto failure;
  opt->bcmat = malloc ( nfunc_b*constr_no*sizeof(double) );
  if ( !opt->bcmat )
    goto failure;
  if ( opt->order == 1 ) {
    if ( opt->coons ) {
      dimV0 = nfunc_a = g1h_V0SpaceDimd ( domain );
      size = 3*dimV0*constr_no;
      cmat = pkv_GetScratchMemd ( size );
      if ( !cmat )
        goto failure;
      memset ( cmat, 0, size*sizeof(double) );
      for ( i = 0;  i < constr_no; i++ ) {
        nc = opt->constr_spec[i][0];
        np = opt->constr_spec[i][1];
        for ( fn = 0; fn < nfunc_a; fn++ ) {
          g1h_GetABasisFPatchCurved ( domain, fn, nc, fcp );
          cmat[i*3*dimV0+fn] = fcp[np];
        }
      }
      pkv_Moved ( constr_no, dimV0, 3*dimV0, dimV0, cmat );
      pkv_Moved ( constr_no, dimV0, 3*dimV0, 2*dimV0, cmat );
      pkn_MultMatrixNumd ( constr_no, dimV0, 3*dimV0, cmat,
                           opt->constrnv.x, 3*dimV0, cmat );
      pkn_MultMatrixNumd ( constr_no, dimV0, 3*dimV0, &cmat[dimV0],
                           opt->constrnv.y, 3*dimV0, &cmat[dimV0] );
      pkn_MultMatrixNumd ( constr_no, dimV0, 3*dimV0, &cmat[2*dimV0],
                           opt->constrnv.z, 3*dimV0, &cmat[2*dimV0] );
      if ( g1h_SetAltConstraintMatrixd ( domain, 3, constr_no, cmat ) )
        goto comp_poly_bcmat1;
      else
        goto failure;
    }
    else if ( opt->bezier ) {
      nfunc_a = g1h_V0SpaceDimd ( domain );
      dimV0 = g1h_ExtV0SpaceDimd ( domain );
      nfunc_c = dimV0-nfunc_a;
      size = 3*dimV0*constr_no;
      cmat = pkv_GetScratchMemd ( size );
      if ( !cmat )
        goto failure;
      memset ( cmat, 0, size*sizeof(double) );
      for ( i = 0;  i < constr_no;  i++ ) {
        nc = opt->constr_spec[i][0];
        np = opt->constr_spec[i][1];
        for ( fn = 0; fn < nfunc_a; fn++ ) {
          g1h_GetABasisFPatchCurved ( domain, fn, nc, fcp );
          cmat[i*3*dimV0+nfunc_c+fn] = fcp[np];
        }
      }
      pkv_Moved ( constr_no, nfunc_a, 3*dimV0, dimV0, &cmat[nfunc_c] );
      pkv_Moved ( constr_no, nfunc_a, 3*dimV0, 2*dimV0, &cmat[nfunc_c] );
      pkn_MultMatrixNumd ( constr_no, nfunc_a, 3*dimV0, &cmat[nfunc_c],
                           opt->constrnv.x, 3*dimV0, &cmat[nfunc_c] );
      pkn_MultMatrixNumd ( constr_no, nfunc_a, 3*dimV0, &cmat[nfunc_c+dimV0],
                           opt->constrnv.y, 3*dimV0, &cmat[nfunc_c+dimV0] );
      pkn_MultMatrixNumd ( constr_no, nfunc_a, 3*dimV0, &cmat[nfunc_c+2*dimV0],
                           opt->constrnv.z, 3*dimV0, &cmat[nfunc_c+2*dimV0] );
      if ( !g1h_SetExtAltConstraintMatrixd ( domain, 3, constr_no, cmat ) )
        goto failure;
comp_poly_bcmat1:
      for ( i = 0;  i < constr_no;  i++ ) {
        nc = opt->constr_spec[i][0];
        np = opt->constr_spec[i][1];
        for ( fn = 0; fn < nfunc_b; fn++ ) {
          g1h_GetBBasisFPatchCurved ( domain, fn, nc, fcp );
          opt->bcmat[i*nfunc_b+fn] = fcp[np];
        }
      }
    }
    else {  /* spline G1 */
      g1h_DrawSplBasFuncNumd ( domain, &nfunc_a, &nfunc_b, &nfunc_c, &nfunc_d );
      dimV0 = nfunc_a+nfunc_c+nfunc_d;
      size = 3*dimV0*constr_no;
      kncp = pkv_GetScratchMemd ( opt->ncp+G1H_OMCDEG+1 );
      cmat = pkv_GetScratchMemd ( size );
      if ( !kncp || !cmat )
        goto failure;
      memset ( cmat, 0, size*sizeof(double) );
      for ( i = 0;  i < constr_no;  i++ ) {
        nc = opt->constr_spec[i][0];
        np = opt->constr_spec[i][1];
        for ( fn = 0; fn < nfunc_a; fn++ ) {
          g1h_GetSplABasisFPatchCurved ( domain, fn, nc, &lkn, kncp, fcp );
          cmat[i*3*dimV0+nfunc_c+nfunc_d+fn] = fcp[np];
        }
        for ( fn = 0; fn < nfunc_d; fn++ ) {
          g1h_GetSplDBasisFPatchCurved ( domain, fn, nc, &lkn, kncp, fcp );
          cmat[i*3*dimV0+nfunc_c+fn] = fcp[np];
        }
      }
      pkv_Moved ( constr_no, nfunc_a+nfunc_d, 3*dimV0, dimV0, &cmat[nfunc_c] );
      pkv_Moved ( constr_no, nfunc_a+nfunc_d, 3*dimV0, 2*dimV0, &cmat[nfunc_c] );
      pkn_MultMatrixNumd ( constr_no, nfunc_a+nfunc_d, 3*dimV0, &cmat[nfunc_c],
                           opt->constrnv.x, 3*dimV0, &cmat[nfunc_c] );
      pkn_MultMatrixNumd ( constr_no, nfunc_a+nfunc_d, 3*dimV0, &cmat[dimV0+nfunc_c],
                           opt->constrnv.y, 3*dimV0, &cmat[dimV0+nfunc_c] );
      pkn_MultMatrixNumd ( constr_no, nfunc_a+nfunc_d, 3*dimV0, &cmat[2*dimV0+nfunc_c],
                           opt->constrnv.z, 3*dimV0, &cmat[2*dimV0+nfunc_c] );
      if ( !g1h_SetSplAltConstraintMatrixd ( domain, 3, constr_no, cmat ) )
        goto failure;
      for ( i = 0;  i < constr_no;  i++ ) {
        nc = opt->constr_spec[i][0];
        np = opt->constr_spec[i][1];
        for ( fn = 0; fn < nfunc_b; fn++ ) {
          g1h_GetSplBBasisFPatchCurved ( domain, fn, nc, &lkn, kncp, fcp );
          opt->bcmat[i*nfunc_b+fn] = fcp[np];
        }
      }
    }
  }
  else { /* G2 */
    if ( opt->coons ) {
      dimV0 = nfunc_a = g2h_V0SpaceDimd ( domain );
      size = 3*dimV0*constr_no;
      cmat = pkv_GetScratchMemd ( size );
      if ( !cmat )
        goto failure;
      memset ( cmat, 0, size*sizeof(double) );
      for ( i = 0;  i < constr_no; i++ ) {
        nc = opt->constr_spec[i][0];
        np = opt->constr_spec[i][1];
        for ( fn = 0; fn < nfunc_a; fn++ ) {
          g2h_GetABasisFPatchCurved ( domain, fn, nc, fcp );
          cmat[i*3*dimV0+fn] = fcp[np];
        }
      }
      pkv_Moved ( constr_no, dimV0, 3*dimV0, dimV0, cmat );
      pkv_Moved ( constr_no, dimV0, 3*dimV0, 2*dimV0, cmat );
      pkn_MultMatrixNumd ( constr_no, dimV0, 3*dimV0, cmat,
                           opt->constrnv.x, 3*dimV0, cmat );
      pkn_MultMatrixNumd ( constr_no, dimV0, 3*dimV0, &cmat[dimV0],
                           opt->constrnv.y, 3*dimV0, &cmat[dimV0] );
      pkn_MultMatrixNumd ( constr_no, dimV0, 3*dimV0, &cmat[2*dimV0],
                           opt->constrnv.z, 3*dimV0, &cmat[2*dimV0] );
      if ( g2h_SetAltConstraintMatrixd ( domain, 3, constr_no, cmat ) )
        goto comp_poly_bcmat2;
      else
        goto failure;
    }
    else if ( opt->bezier ) {
      nfunc_a = g2h_V0SpaceDimd ( domain );
      dimV0 = g2h_ExtV0SpaceDimd ( domain );
      nfunc_c = dimV0-nfunc_a;
      size = 3*dimV0*constr_no;
      cmat = pkv_GetScratchMemd ( size );
      if ( !cmat )
        goto failure;
      memset ( cmat, 0, size*sizeof(double) );
      for ( i = 0;  i < constr_no;  i++ ) {
        nc = opt->constr_spec[i][0];
        np = opt->constr_spec[i][1];
        for ( fn = 0; fn < nfunc_a; fn++ ) {
          g2h_GetABasisFPatchCurved ( domain, fn, nc, fcp );
          cmat[i*3*dimV0+nfunc_c+fn] = fcp[np];
        }
      }
      pkv_Moved ( constr_no, nfunc_a, 3*dimV0, dimV0, &cmat[nfunc_c] );
      pkv_Moved ( constr_no, nfunc_a, 3*dimV0, 2*dimV0, &cmat[nfunc_c] );
      pkn_MultMatrixNumd ( constr_no, nfunc_a, 3*dimV0, &cmat[nfunc_c],
                           opt->constrnv.x, 3*dimV0, &cmat[nfunc_c] );
      pkn_MultMatrixNumd ( constr_no, nfunc_a, 3*dimV0, &cmat[nfunc_c+dimV0],
                           opt->constrnv.y, 3*dimV0, &cmat[nfunc_c+dimV0] );
      pkn_MultMatrixNumd ( constr_no, nfunc_a, 3*dimV0, &cmat[nfunc_c+2*dimV0],
                           opt->constrnv.z, 3*dimV0, &cmat[nfunc_c+2*dimV0] );
      if ( !g2h_SetExtAltConstraintMatrixd ( domain, 3, constr_no, cmat ) )
        goto failure;
comp_poly_bcmat2:
      for ( i = 0;  i < constr_no;  i++ ) {
        nc = opt->constr_spec[i][0];
        np = opt->constr_spec[i][1];
        for ( fn = 0; fn < nfunc_b; fn++ ) {
          g2h_GetBBasisFPatchCurved ( domain, fn, nc, fcp );
          opt->bcmat[i*nfunc_b+fn] = fcp[np];
        }
      }
    }
    else {  /* spline G2 */
      g2h_DrawSplBasFuncNumd ( domain, &nfunc_a, &nfunc_b, &nfunc_c, &nfunc_d );
      dimV0 = nfunc_a+nfunc_c+nfunc_d;
      size = 3*dimV0*constr_no;
      kncp = pkv_GetScratchMemd ( opt->ncp+G2H_OMCDEG+1 );
      cmat = pkv_GetScratchMemd ( size );
      if ( !kncp || !cmat )
        goto failure;
      memset ( cmat, 0, size*sizeof(double) );
      for ( i = 0;  i < constr_no;  i++ ) {
        nc = opt->constr_spec[i][0];
        np = opt->constr_spec[i][1];
        for ( fn = 0; fn < nfunc_a; fn++ ) {
          g2h_GetSplABasisFPatchCurved ( domain, fn, nc, &lkn, kncp, fcp );
          cmat[i*3*dimV0+nfunc_c+nfunc_d+fn] = fcp[np];
        }
        for ( fn = 0; fn < nfunc_d; fn++ ) {
          g2h_GetSplDBasisFPatchCurved ( domain, fn, nc, &lkn, kncp, fcp );
          cmat[i*3*dimV0+nfunc_c+fn] = fcp[np];
        }
      }
      pkv_Moved ( constr_no, nfunc_a+nfunc_d, 3*dimV0, dimV0, &cmat[nfunc_c] );
      pkv_Moved ( constr_no, nfunc_a+nfunc_d, 3*dimV0, 2*dimV0, &cmat[nfunc_c] );
      pkn_MultMatrixNumd ( constr_no, nfunc_a+nfunc_d, 3*dimV0, &cmat[nfunc_c],
                           opt->constrnv.x, 3*dimV0, &cmat[nfunc_c] );
      pkn_MultMatrixNumd ( constr_no, nfunc_a+nfunc_d, 3*dimV0, &cmat[dimV0+nfunc_c],
                           opt->constrnv.y, 3*dimV0, &cmat[dimV0+nfunc_c] );
      pkn_MultMatrixNumd ( constr_no, nfunc_a+nfunc_d, 3*dimV0, &cmat[2*dimV0+nfunc_c],
                           opt->constrnv.z, 3*dimV0, &cmat[2*dimV0+nfunc_c] );
      if ( !g2h_SetSplAltConstraintMatrixd ( domain, 3, constr_no, cmat ) )
        goto failure;
      for ( i = 0;  i < constr_no;  i++ ) {
        nc = opt->constr_spec[i][0];
        np = opt->constr_spec[i][1];
        for ( fn = 0; fn < nfunc_b; fn++ ) {
          g2h_GetSplBBasisFPatchCurved ( domain, fn, nc, &lkn, kncp, fcp );
          opt->bcmat[i*nfunc_b+fn] = fcp[np];
        }
      }
    }
  }
  opt->constr_matrix_valid = true;

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  opt->constr_matrix_valid = false;
  pkv_SetScratchMemTop ( sp );
  return false;
} /*ComputeConstraintMatricesAlt*/

boolean ComputeConstraintRightSideAlt ( int surfno, double *altrhs )
{
  void      *sp;
  GHoptions *opt;
  int       nfunc_b, constr_no;
  int       i, j, nc, np, ncp;
  unsigned char *bfcpn;
  vector3d  ccp;

  sp = pkv_GetScratchMemTop ();
  switch ( surfno ) {
case 1:
    opt = &options1;
    break;
case 2:
    opt = &options2;
    break;
default:
    goto failure;
  }
  if ( !opt->bcmat )
    goto failure;

  nfunc_b = 6*hole_k+1;
  bfcpn = pkv_GetScratchMem ( nfunc_b*sizeof(unsigned char) );
  if ( !bfcpn )
    goto failure;
  if ( !gh_DrawBFcpn ( hole_k, bfcpn ) )
    goto failure;
  constr_no = opt->constr_no;
  ncp = opt->ncp;
  for ( i = 0; i < constr_no; i++ ) {
    nc = opt->constr_spec[i][0];
    np = opt->constr_spec[i][1];
    ccp = opt->constrcp[ncp*nc+np];
    for ( j = 0; j < nfunc_b; j++ )
      AddVector3Md ( &ccp, &hole_cp[bfcpn[j]],
                     -opt->bcmat[nfunc_b*i+j], &ccp );
    altrhs[i] = DotProduct3d ( &opt->constrnv, &ccp );
  }

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*ComputeConstraintRightSideAlt*/

