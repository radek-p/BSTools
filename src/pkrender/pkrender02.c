
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2015                                  */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#include "pkvaria.h"
#include "pkvthreads.h"
#include "pknum.h"
#include "pkgeom.h"
#include "multibs.h"
#include "camera.h"
#include "raybez.h"
#include "pkrender.h"

#include "pkrenderprivate.h"

static boolean ReallocObjTab ( pkRenderer *rend )
{
  renderobj *newtab;
  int       newtablength;

  newtablength = 2*rend->obj_tab_length;
  PKV_MALLOC ( newtab, newtablength*sizeof(renderobj) );
  if ( newtab ) {
    memcpy ( newtab, rend->obj_tab, rend->nobjects*sizeof(renderobj) );
    PKV_FREE ( rend->obj_tab );
    rend->obj_tab = newtab;
    rend->obj_tab_length = newtablength;
    return true;
  }
  else
    return false;
} /*ReallocObjTab*/

boolean RendEnterTriangle3d ( pkRenderer *rend,
                              point3d *p0, point3d *p1, point3d *p2,
                              double *colour )
{
#define TOL 1.0e-14
  renderobj_triangle *tr;
  triangle_data      *trd;
  double             m11, m21, m22, l11, l21, l22;
  vector3d           v1, v2;

  trd = NULL;
  if ( rend->nobjects >= rend->obj_tab_length ) {
    if ( !ReallocObjTab ( rend ) )
      goto failure;
  }
  PKV_MALLOC ( trd, sizeof(triangle_data) );
  if ( trd ) {
    tr = &rend->obj_tab[rend->nobjects].triang;
    tr->type = obj_TRIANGLE;
    memcpy ( tr->colour, colour, 3*sizeof(double) );
    tr->trdata = trd;
        /* find the bounding box */
    rbez_InitBBox3d ( &trd->bbox, p0 );
    rbez_Extend2BBox3d ( &trd->bbox, p1, p2 );
        /* compute the unit normal vector */
    trd->p0 = *p0;
    SubtractPoints3d ( p1, p0, &v1 );
    SubtractPoints3d ( p2, p0, &v2 );
    CrossProduct3d ( &v1, &v2, &trd->n );
    NormalizeVector3d ( &trd->n );
        /* find the pseudo-inversion of the 3x2 matrix A=[v1,v2] */
          /* compute A^T*A */
    m11 = DotProduct3d ( &v1, &v1 );
    if ( m11 <= TOL )
      goto failure;
    m21 = DotProduct3d ( &v2, &v1 );
    m22 = DotProduct3d ( &v2, &v2 );
          /* this is the Cholesky's decomposition of a 2x2 matrix A^T*A */
    l11 = sqrt ( m11 );
    l21 = m21/l11;
    l22 = m22-m21*m21/m11;
    if ( l22 <= TOL )
      goto failure;
          /* compute the rows a1^T, a2^T, of (A^T*A)^{-1}A^T */
    MultVector3d ( 1.0/l11, &v1, &v1 );
    AddVector3Md ( &v2, &v1, -l21, &v2 );
    MultVector3d ( 1.0/l22, &v2, &trd->a2 );
    AddVector3Md ( &v1, &trd->a2, -l21, &v1 );
    MultVector3d ( 1.0/l11, &v1, &trd->a1 );
    rend->nobjects ++;
    return true;
  }
  else {
failure:
    if ( trd ) PKV_FREE ( trd );
    return false;
  }
#undef TOL
} /*RendEnterTriangle3d*/

boolean RendEnterBezPatch3d ( pkRenderer *rend,
                              int n, int m, const point3d *cp, double *colour )
{
  if ( rend->nobjects >= rend->obj_tab_length ) {
    if ( !ReallocObjTab ( rend ) )
      goto failure;
  }
  if ( rend->nobjects < rend->obj_tab_length ) {
    rend->obj_tab[rend->nobjects].type = obj_BSPATCH;
    rend->obj_tab[rend->nobjects].bsp.ptree =
          rbez_NewBezPatchTreed ( rend->nobjects, n, m, 0.0, 1.0, 0.0, 1.0, cp );
    if ( rend->obj_tab[rend->nobjects].bsp.ptree ) {
      memcpy ( &rend->obj_tab[rend->nobjects].bsp.colour[0],
               colour, 3*sizeof(double) );
      rend->nobjects ++;
      return true;
    }
    else
      goto failure;
  }
failure:
  return false;
} /*RendEnterBezPatch3d*/

boolean RendEnterBSPatch3d ( pkRenderer *rend,
                             int n, int lknu, const double *knu,  
                             int m, int lknv, const double *knv,  
                             const point3d *cp, double *colour )
{
  int pitch;

  if ( rend->nobjects >= rend->obj_tab_length ) {
    if ( !ReallocObjTab ( rend ) )
      goto failure;
  }
  if ( rend->nobjects < rend->obj_tab_length ) {
    pitch = lknv-m;
    rend->obj_tab[rend->nobjects].type = obj_BSPATCH;
    rend->obj_tab[rend->nobjects].bsp.ptree =
          rbez_NewBSPatchTreed ( rend->nobjects, n, lknu, knu, m, lknv, knv,
                                 pitch, cp );
    if ( rend->obj_tab[rend->nobjects].bsp.ptree ) {
      memcpy ( &rend->obj_tab[rend->nobjects].bsp.colour[0],
               colour, 3*sizeof(double) );
      rend->nobjects ++;
      return true;
    }
    else
      goto failure;
  }
failure:
  return false;
} /*RendEnterBSPatch3d*/

boolean RendEnterBezPatch3Rd ( pkRenderer *rend,
                               int n, int m, const point4d *cp, double *colour )
{
  if ( rend->nobjects >= rend->obj_tab_length ) {
    if ( !ReallocObjTab ( rend ) )
      goto failure;
  }
  if ( rend->nobjects < rend->obj_tab_length ) {
    rend->obj_tab[rend->nobjects].type = obj_RBSPATCH;
    rend->obj_tab[rend->nobjects].rbsp.ptree =
          rbez_NewRBezPatchTreed ( rend->nobjects, n, m, 0.0, 1.0, 0.0, 1.0, cp );
    if ( rend->obj_tab[rend->nobjects].rbsp.ptree ) {
      memcpy ( &rend->obj_tab[rend->nobjects].rbsp.colour[0],
               colour, 3*sizeof(double) );
      rend->nobjects ++;
      return true;
    }
    else
      goto failure;
  }
failure:
  return false;
} /*RendEnterBezPatch3Rd*/

boolean RendEnterBSPatch3Rd ( pkRenderer *rend,
                              int n, int lknu, const double *knu,  
                              int m, int lknv, const double *knv,  
                              const point4d *cp, double *colour )
{
  int pitch;

  if ( rend->nobjects >= rend->obj_tab_length ) {
    if ( !ReallocObjTab ( rend ) )
      goto failure;
  }
  if ( rend->nobjects < rend->obj_tab_length ) {
    pitch = lknv-m;
    rend->obj_tab[rend->nobjects].type = obj_RBSPATCH;
    rend->obj_tab[rend->nobjects].rbsp.ptree =
          rbez_NewRBSPatchTreed ( rend->nobjects, n, lknu, knu, m, lknv, knv,
                                  pitch, cp );
    if ( rend->obj_tab[rend->nobjects].rbsp.ptree ) {
      memcpy ( &rend->obj_tab[rend->nobjects].rbsp.colour[0],
               colour, 3*sizeof(double) );
      rend->nobjects ++;
      return true;
    }
    else
      goto failure;
  }
failure:
  return false;
} /*RendEnterBSPatch3Rd*/

boolean RendEnterBezCurve3d ( pkRenderer *rend,
                              int n, const point3d *cp, double r, double *colour )
{
  if ( rend->nobjects >= rend->obj_tab_length ) {
    if ( !ReallocObjTab ( rend ) )
      goto failure;
  }
  if ( rend->nobjects < rend->obj_tab_length ) {
    rend->obj_tab[rend->nobjects].type = obj_BEZCURVE;
    rend->obj_tab[rend->nobjects].bezc.ctree =
          rbez_NewBezCurveTreed ( rend->nobjects, n, 0.0, 1.0, cp, r );
    if ( rend->obj_tab[rend->nobjects].bezc.ctree ) {
      memcpy ( &rend->obj_tab[rend->nobjects].bezc.colour[0],
               colour, 3*sizeof(double) );
      rend->nobjects ++;
      return true;
    }
    else
      goto failure;
  }
failure:
  return false;
} /*RendEnterBezCurve3d*/

boolean RendEnterBSCurve3d ( pkRenderer *rend,
                             int n, int lkn, const double *kn,
                             const point3d *cp, double r, double *colour )
{
  if ( rend->nobjects >= rend->obj_tab_length ) {
    if ( !ReallocObjTab ( rend ) )
      goto failure;
  }
  if ( rend->nobjects < rend->obj_tab_length ) {
    rend->obj_tab[rend->nobjects].type = obj_BEZCURVE;
    rend->obj_tab[rend->nobjects].bezc.ctree =
          rbez_NewBSCurveTreed ( rend->nobjects, n, lkn, kn, cp, r );
    if ( rend->obj_tab[rend->nobjects].bezc.ctree ) {
      memcpy ( &rend->obj_tab[rend->nobjects].bezc.colour[0],
               colour, 3*sizeof(double) );
      rend->nobjects ++;
      return true;
    }
    else
      goto failure;
  }
failure:
  return false;
} /*RendEnterBSCurve3d*/

boolean RendEnterBezCurve3Rd ( pkRenderer *rend,
                               int n, const point4d *cp, double r, double *colour )
{
  if ( rend->nobjects >= rend->obj_tab_length ) {
    if ( !ReallocObjTab ( rend ) )
      goto failure;
  }
  if ( rend->nobjects < rend->obj_tab_length ) {
    rend->obj_tab[rend->nobjects].type = obj_RBEZCURVE;
    rend->obj_tab[rend->nobjects].rbezc.ctree =
          rbez_NewRBezCurveTreed ( rend->nobjects, n, 0.0, 1.0, cp, r );
    if ( rend->obj_tab[rend->nobjects].rbezc.ctree ) {
      memcpy ( &rend->obj_tab[rend->nobjects].rbezc.colour[0],
               colour, 3*sizeof(double) );
      rend->nobjects ++;
      return true;
    }
    else
      goto failure;
  }
failure:
  return false;
} /*RendEnterBezCurve3Rd*/

boolean RendEnterBSCurve3Rd ( pkRenderer *rend,
                              int n, int lkn, const double *kn,
                              const point4d *cp, double r, double *colour )
{
  if ( rend->nobjects >= rend->obj_tab_length ) {
    if ( !ReallocObjTab ( rend ) )
      goto failure;
  }
  if ( rend->nobjects < rend->obj_tab_length ) {
    rend->obj_tab[rend->nobjects].type = obj_RBEZCURVE;
    rend->obj_tab[rend->nobjects].rbezc.ctree =
          rbez_NewRBSCurveTreed ( rend->nobjects, n, lkn, kn, cp, r );
    if ( rend->obj_tab[rend->nobjects].rbezc.ctree ) {
      memcpy ( &rend->obj_tab[rend->nobjects].rbezc.colour[0],
               colour, 3*sizeof(double) );
      rend->nobjects ++;
      return true;
    }
    else
      goto failure;
  }
failure:
  return false;
} /*RendEnterBSCurve3Rd*/

