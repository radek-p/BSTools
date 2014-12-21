
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2012, 2014                            */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

#include <pthread.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

#include "pkvaria.h"
#include "pkgeom.h"
#include "pknum.h"
#include "multibs.h"
#include "raybez.h"

/* ////////////////////////////////////////////////////////////////////////// */
#define MAX_LEVEL 32

typedef struct {
    RBezPatchTreeVertexdp vert1, vert2;
    byte                  level;
    boolean               nvpassed;
  } rbiHyperDomaind;


/* ////////////////////////////////////////////////////////////////////////// */
vector3d *_rbi_GetRBezNVPatchd ( RBezPatchTreedp tree,
                                 RBezPatchTreeVertexdp vertex )
{
  RBezPatchTreeVertexdp up, lv, rv;
  int                   vn, vm, nvsize;
  vector3d              *nvcp;

  if ( !vertex->nvcpoints ) {
    if ( (up = vertex->up) ) {
      if ( (nvcp = _rbi_GetRBezNVPatchd ( tree, up )) ) {
        nvsize = tree->nvsize;
        lv = up->left;
        rv = up->right;
        PKV_MALLOC ( lv->nvcpoints, nvsize );
        PKV_MALLOC ( rv->nvcpoints, nvsize );
        if ( !lv->nvcpoints || !rv->nvcpoints ) {
          if ( lv->nvcpoints ) PKV_FREE ( lv->nvcpoints );
          if ( rv->nvcpoints ) PKV_FREE ( rv->nvcpoints );
          return NULL;
        }
        memcpy ( rv->nvcpoints, nvcp, nvsize );
        if ( !up->divdir )
          mbs_BisectBP3ud ( tree->vn, tree->vm, rv->nvcpoints, lv->nvcpoints );
        else
          mbs_BisectBP3vd ( tree->vn, tree->vm, rv->nvcpoints, lv->nvcpoints );
      }
    }
    else {    /* the vertex is the tree root */
      mbs_BezP3RNormalDeg ( tree->n, tree->m, &vn, &vm );
      tree->vn = (unsigned char)vn;
      tree->vm = (unsigned char)vm;
      tree->nvsize = nvsize = (vn+1)*(vm+1)*sizeof(vector3d);
      PKV_MALLOC ( vertex->nvcpoints, nvsize );
      if ( !vertex->nvcpoints )
        return NULL;
      mbs_BezP3RNormald ( tree->n, tree->m, vertex->ctlpoints,
                          &vn, &vm, vertex->nvcpoints );
    }
  }
  return vertex->nvcpoints;
} /*_rbi_GetRBezNVPatchd*/

/* ////////////////////////////////////////////////////////////////////////// */
double _rbi_AuxPCHdpRd ( point4d *v, vector3d *w )
{
  return (v->x*w->x + v->y*w->y + v->z*w->z)/v->w;
} /*_rbi_AuxPCHdpRd*/

boolean _rbi_AuxPConvexHullTestRd ( int k1, point4d *p1, int k2, point4d *p2,
                                    vector3d *w )
{
#define EPS 1.0e-11
  int    i, m;
  double a, b, c;

  if ( DotProduct3d ( w, w ) < EPS )
    return true;
  m = min ( k1, k2 );
  a = _rbi_AuxPCHdpRd ( &p1[0], w );
  b = _rbi_AuxPCHdpRd ( &p2[0], w );
  if ( a <b ) {
    a = -a;  b = -b;
    w->x = -w->x;  w->y = -w->y;  w->z = -w->z;
  }
  for ( i = 1; i < m; i++ ) {
    c = _rbi_AuxPCHdpRd ( &p1[i], w );  a = min ( a, c );
    c = _rbi_AuxPCHdpRd ( &p2[i], w );  b = max ( b, c );
    if ( a <= b )
      return true;
  }
  for ( ; i < k1; i++ ) {
    c = _rbi_AuxPCHdpRd ( &p1[i], w );
    if ( c <= b )
      return true;
  }
  for ( ; i < k2; i++ ) {
    c = _rbi_AuxPCHdpRd ( &p2[i], w );
    if ( a <= c )
      return true;
  }
  return false;
#undef EPS
} /*_rbi_AuxPConvexHullTestRd*/

boolean _rbi_PatchConvexHullTestRd ( int n1, int m1, point4d *p1,
                                     point3d *pcent1, vector3d *nvcent1,
                                     int n2, int m2, point4d *p2,
                                     point3d *pcent2, vector3d *nvcent2 )
{
  vector3d w, w1, w2;
  int      k1, k2;

  k1 = (n1+1)*(m1+1);
  k2 = (n2+1)*(m2+1);
  SubtractPoints3d ( pcent2, pcent1, &w );
  if ( _rbi_AuxPConvexHullTestRd ( k1, p1, k2, p2, &w ) ) {
    w1 = *nvcent1;  NormalizeVector3d ( &w1 );
    w2 = *nvcent2;  NormalizeVector3d ( &w2 );
    AddVector3d ( &w1, &w2, &w );
    if ( _rbi_AuxPConvexHullTestRd ( k1, p1, k2, p2, &w ) ) {
      SubtractPoints3d ( &w1, &w2, &w );
      return _rbi_AuxPConvexHullTestRd ( k1, p1, k2, p2, &w );
    }
  }
  return false;
} /*_rbi_PatchConvexHullTestRd*/

/* ////////////////////////////////////////////////////////////////////////// */
boolean _rbi_NormalTestRd ( int vn1, int vm1, vector3d *nvcp1, vector3d *nv1,
                            int vn2, int vm2, vector3d *nvcp2, vector3d *nv2 )
{
#define EPS_NT 1.0e-6
#define BALL_R 0.6
  vector3d    we, tv;
  double      m[9];
  int         P[3], Q[3];
  int         i, k;
  ConvexCone2d cone1, cone2, cone;

        /* find the direction of the wedge edge */
  CrossProduct3d ( nv1, nv2, &we );
  m[0] = nv1->x;  m[1] = nv2->x;  m[2] = we.x;
  m[3] = nv1->y;  m[4] = nv2->y;  m[5] = we.y;
  m[6] = nv1->z;  m[7] = nv2->z;  m[8] = we.z;
  if ( !pkn_GaussDecomposePLUQd ( 3, m, P, Q ) )
    return false;

        /* test the first normal vector patch */
  cone1.code = 0;
  k = (vn1+1)*(vm1+1);
  for ( i = 0; i < k; i++ ) {
    tv = nvcp1[i];
    pkn_multiSolvePLUQd ( 3, m, P, Q, 1, 1, &tv.x );
    if ( ExtendConvexCone2dv ( &cone1, (vector2d*)&tv ) == 2 )
      return false;
  }
        /* test the second normal vector patch */
  cone2.code = 0;
  k = (vn2+1)*(vm2+1);
  for ( i = 0; i < k; i++ ) {
    tv = nvcp2[i];
    pkn_multiSolvePLUQd ( 3, m, P, Q, 1, 1, &tv.x );
    if ( ExtendConvexCone2dv ( &cone2, (vector2d*)&tv ) == 2 )
      return false;
  }
        /* test if the wedges are disjoint */
  if ( InsideConvexCone2d ( &cone1, cone2.amin ) ) return false;
  if ( InsideConvexCone2d ( &cone1, cone2.amax ) ) return false;
  if ( InsideConvexCone2d ( &cone2, cone1.amin ) ) return false;
  if ( InsideConvexCone2d ( &cone2, cone1.amax ) ) return false;
        /* test if the wedges are contained in a convex cone */
  cone = cone1;
  if ( ExtendConvexCone2d ( &cone, cone2.amin ) == 2 )
    return false;
  if ( ExtendConvexCone2d ( &cone, cone2.amax ) == 2 )
    return false;
  cone = cone2;
  if ( ExtendConvexCone2d ( &cone, cone1.amin ) == 2 )
    return false;
  if ( ExtendConvexCone2d ( &cone, cone1.amax ) == 2 )
    return false;
  return true;
} /*_rbi_NormalTestRd*/

/* ////////////////////////////////////////////////////////////////////////// */
#define POINT_BUF_LENGTH 101

typedef struct {
    RBezPatchTreeVertexdp pvertex;
    RBezCurveTreeVertexdp cvertex;
    short                 nc;
    short                 level;
  } pc_domain;

typedef struct {
    RBezPatchTreedp       ptree;
    RBezPatchTreeVertexdp pvertex;
    RBezCurveTreedp       ctree;
    RBezCurveTreeVertexdp cvertex;
    int                   cdir;  /* which curve */
    double                t, s0, s1;
  } pci_key;

boolean _rbi_ArcPointsTooFarRd ( vector4d *p1, vector4d *p2,
                                 double epsilon2 )
{
  vector4d v;

  SubtractPoints4d ( p1, p2, &v );
  return v.x*v.x + v.y*v.y > epsilon2 ||
         v.z*v.z + v.w*v.w > epsilon2;
} /*_rbi_ArcPointsTooFarRd*/

void _rbi_ComputeAFDFRd ( RBezPatchTreedp tree1, RBezPatchTreedp tree2,
                          vector4d *p, point3d *F, double *DF )
{
  vector4d p1, p1u, p1v, p2, p2u, p2v;

  mbs_BCHornerDerP4d ( tree1->n, tree1->m, tree1->root->ctlpoints, p->x, p->y,
                       &p1, &p1u, &p1v );
  mbs_BCHornerDerP4d ( tree2->n, tree2->m, tree2->root->ctlpoints, p->z, p->w,
                       &p2, &p2u, &p2v );

  SetVector3d ( F, p1.x*p2.w - p1.w*p2.x, p1.y*p2.w - p1.w*p2.y, 
                p1.z*p2.w - p1.w*p2.z );

  DF[ 0] = p1u.x*p2.w - p1u.w*p2.x;
  DF[ 3] = p1v.x*p2.w - p1v.w*p2.x;
  DF[ 6] = p1.x*p2u.w - p1.w*p2u.x;
  DF[ 9] = p1.x*p2v.w - p1.w*p2v.x;
  DF[ 1] = p1u.y*p2.w - p1u.w*p2.y;
  DF[ 4] = p1v.y*p2.w - p1v.w*p2.y;
  DF[ 7] = p1.y*p2u.w - p1.w*p2u.y;
  DF[10] = p1.y*p2v.w - p1.w*p2v.y;
  DF[ 2] = p1u.z*p2.w - p1u.w*p2.z;
  DF[ 5] = p1v.z*p2.w - p1v.w*p2.z;
  DF[ 8] = p1.z*p2u.w - p1.w*p2u.z;
  DF[11] = p1.z*p2v.w - p1.w*p2v.z;
} /*_rbi_ComputeAFDFRd*/

boolean _rbi_ImproveArcPointRd ( RBezPatchTreedp tree1, RBezPatchTreedp tree2,
                                 vector4d *p )
{
#define MAXIT 10
#define EPS 1.0e-10
  vector3d F;
  double   DF[12], aa[6], d;
  vector4d dp;
  int      i;

  for ( i = 0; i < MAXIT; i++ ) {
        /* the Newton's method with pseudo-inversion */
    _rbi_ComputeAFDFRd ( tree1, tree2, p, &F, DF );
    if ( pkn_QRDecomposeMatrixd ( 4, 3, DF, aa ) ) {
      pkn_multiMultInvTrUTVectord ( 3, DF, 1, 1, &F.x, 1, &dp.x );
      dp.w = 0.0;
      d = DotProduct3d ( (vector3d*)&dp.x, (vector3d*)&dp.x );
      pkn_multiInvReflectVectord ( 4, 3, DF, aa, 1, 1, &dp.x );
      SubtractPoints4d ( p, &dp, p );
      if ( d < EPS )
        return true;
    }
    else break;
  }
  return false;
#undef EPS
#undef MAXIT
} /*_rbi_ImproveArcPointRd*/

int _rbi_ProcessTheArcRd ( RBezPatchTreedp tree1, RBezPatchTreedp tree2,
                           rbiHyperDomaind *hc, vector4d *ppar,
                           double epsilon2,
                           rbiArcOutd *outproc, void *usrptr )
{
#define AP_STACK_LENGTH 10
  void *sp;
  rbiIntersArcd   rbiarc;
  vector4d *stack, *buf, last, next;
  int             stp;

  sp = pkv_GetScratchMemTop ();
  stack = pkv_GetScratchMem ( AP_STACK_LENGTH*sizeof(vector4d) );
  buf = pkv_GetScratchMem ( POINT_BUF_LENGTH*sizeof(vector4d) );
  if ( !stack || !buf )
    goto failure;
  rbiarc.u0 = hc->vert1->u0;
  rbiarc.u1 = hc->vert1->u1;
  rbiarc.v0 = hc->vert1->v0;
  rbiarc.v1 = hc->vert1->v1;
  rbiarc.s0 = hc->vert2->u0;
  rbiarc.s1 = hc->vert2->u1;
  rbiarc.t0 = hc->vert2->v0;
  rbiarc.t1 = hc->vert2->v1;
        /* output the first arc end point */
  last = buf[0] = ppar[0];
  rbiarc.npoints = 1;
        /* push the second arc end point */
  stp = 0;
  stack[0] = ppar[1];
        /* now the arc subdivision */
  do {
          /* get a copy from the stack */
    next = stack[stp];
    if ( _rbi_ArcPointsTooFarRd ( &last, &next, epsilon2 ) ) {
      MidPoint4d ( &last, &next, &next );
      if ( !_rbi_ImproveArcPointRd ( tree1, tree2, &next ) )
        goto not_quite;
      if ( stp >= AP_STACK_LENGTH-1 || rbiarc.npoints >= POINT_BUF_LENGTH-1 )
        goto failure;
          /* push this point, it will be output later */
      stack[++stp] = next;
    }
    else {
      stp --;                              /* pop */
      last = buf[rbiarc.npoints++] = next; /* output */
    }
  } while ( stp >= 0 );
        /* output the arc */
  outproc ( usrptr, &rbiarc, buf );

  pkv_SetScratchMemTop ( sp );
  return 1;

not_quite:
  pkv_SetScratchMemTop ( sp );
  return 0;

failure:
  pkv_SetScratchMemTop ( sp );
  return -1;
#undef AP_STACK_LENGTH
} /*_rbi_ProcessTheArcRd*/

void _rbi_ExtractCurveRd ( int n, int m, point4d *pcp,
                           int nc, int *cdeg, point4d *ccp )
{
  switch ( nc ) {
case 0:  /* v = 0 */
    pkv_Selectc ( n+1, sizeof(point4d), (m+1)*sizeof(point4d),
                  sizeof(point4d), (char*)pcp, (char*)ccp );
    *cdeg = n;
    break;
case 1:  /* v = 1 */
    pkv_Selectc ( n+1, sizeof(point4d), (m+1)*sizeof(point4d),
                  sizeof(point4d), (char*)&pcp[m], (char*)ccp );
    *cdeg = n;
    break;
case 2:  /* u = 0 */
    memcpy ( ccp, pcp, (m+1)*sizeof(point4d) );
    *cdeg = m;
    break;
case 3:  /* u = 1 */
    memcpy ( ccp, &pcp[n*(m+1)], (m+1)*sizeof(point4d) );
    *cdeg = m;
    break;
  }
} /*_rbi_ExtractCurveRd*/

boolean _rbi_PCConvexHullTestRd ( RBezPatchTreedp ptree, RBezCurveTreedp ctree,
                                  pc_domain *ha )
{
  vector3d v;
  int      npp, ncp;

  npp = (ptree->n+1)*(ptree->m+1);
  ncp = ctree->degree+1;
  if ( _rbi_AuxPConvexHullTestRd ( npp, ha->pvertex->ctlpoints,
                                   ncp, ha->cvertex->ctlpoints,
                                   &ha->pvertex->nvcent ) ) {
    SubtractPoints3d ( &ha->pvertex->pcent, &ha->cvertex->ccent, &v );
    return _rbi_AuxPConvexHullTestRd ( npp, ha->pvertex->ctlpoints,
                                       ncp, ha->cvertex->ctlpoints, &v );
  }
  else
    return false;
} /*_rbi_PCConvexHullTestRd*/

int _rbi_PCHodogTestRd ( RBezPatchTreedp ptree, RBezCurveTreedp ctree,
                         pc_domain *ha, point3d *kk )
{
#define BALL 0.9 /* must be < 1.0 */
  void                  *sp;
  int                   n, m, d, size;
  vector4d              *pcp, *ccp, *pu, *pv, *ct, pm, pum, pvm, cm, ctm;
  RBezPatchTreeVertexdp pvertex;
  RBezCurveTreeVertexdp cvertex;
  double                DF[9];
  int                   P[3], Q[3];
  point3d               tv;
  double                k1, k2, k3, ks, rr;
  int                   i, j, l;

  sp = pkv_GetScratchMemTop ();
  n = ptree->n;
  m = ptree->m;
  d = ctree->degree;
  pvertex = ha->pvertex;
  pcp = pvertex->ctlpoints;
  cvertex = ha->cvertex;
  ccp = cvertex->ctlpoints;
        /* find the derivatives */
  size = (n*(m+1) + (n+1)*m + d)*sizeof(vector4d);
  pu = pkv_GetScratchMem ( size );
  if ( !pu )
    goto failure;
  pv = &pu[n*(m+1)];
  ct = &pv[(n+1)*m];
  mbs_multiFindBezDerivatived ( n, 1, 4*(m+1), 0, (double*)pcp, 0, (double*)pu );
  mbs_multiFindBezDerivatived ( m, n+1, 4, 4*(m+1), (double*)pcp, 4*m, (double*)pv );
  mbs_FindBezDerivativeC4d ( d, ccp, ct );
        /* find the partial derivatives at (0.5,0.5,0.5) */
  mbs_BCHornerDerP4d ( n, m, (double*)pcp, 0.5, 0.5, &pm, &pum, &pvm );
  mbs_BCHornerDerC4d ( d, (double*)ccp, 0.5, &cm, &ctm );
        /* compute the differential of F at (0.5,0.5,0.5) */
  DF[0] = pum.x*cm.w-pum.w*cm.x;
  DF[1] = pvm.x*cm.w-pvm.w*cm.x;
  DF[2] = pm.x*ctm.w-pm.w*ctm.x;
  DF[3] = pum.y*cm.w-pum.w*cm.y;
  DF[4] = pvm.y*cm.w-pvm.w*cm.y;
  DF[5] = pm.y*ctm.w-pm.w*ctm.y;
  DF[6] = pum.z*cm.w-pum.w*cm.z;
  DF[7] = pvm.z*cm.w-pvm.w*cm.z;
  DF[8] = pm.z*ctm.w-pm.w*ctm.z;
        /* find a triangular decomposition of DF */
  if ( !pkn_GaussDecomposePLUQd ( 3, DF, P, Q ) )
    goto nodiff;
        /* the actual hodograph test */
  k1 = 0.0;
  l = n*(m+1);
  for ( i = 0; i < l; i++ )
    for ( j = 0; j <= d; j++ ) {
      tv.x = pu[i].x*ccp[j].w - pu[i].w*ccp[j].x;
      tv.y = pu[i].y*ccp[j].w - pu[i].w*ccp[j].y;
      tv.z = pu[i].z*ccp[j].w - pu[i].w*ccp[j].z;
      pkn_multiSolvePLUQd ( 3, DF, P, Q, 1, 1, &tv.x );
      tv.x -= 1.0;
      rr = DotProduct3d ( &tv, &tv );
      if ( rr > k1 ) {
        if ( rr >= BALL ) goto nodiff;
                     else k1 = rr;
      }
    }
  ks = k1;
  k2 = 0.0;
  l = (n+1)*m;
  for ( i = 0; i < l; i++ )
    for ( j = 0; j <= d; j++ ) {
      tv.x = pv[i].x*ccp[j].w - pv[i].w*ccp[j].x;
      tv.y = pv[i].y*ccp[j].w - pv[i].w*ccp[j].y;
      tv.z = pv[i].z*ccp[j].w - pv[i].w*ccp[j].z;
      pkn_multiSolvePLUQd ( 3, DF, P, Q, 1, 1, &tv.x );
      tv.y -= 1.0;
      rr = DotProduct3d ( &tv, &tv );
      if ( rr > k2 ) {
        if ( ks+rr >= BALL ) goto nodiff;
                             else k2 = rr;
      }
    }
  ks += k2;
  k3 = 0.0;
  l = (n+1)*(m+1);
  for ( i = 0; i < l; i++ )
    for ( j = 0; j < d; j++ ) {
      tv.x = pcp[i].x*ct[j].w - pcp[i].w*ct[j].x;
      tv.y = pcp[i].y*ct[j].w - pcp[i].w*ct[j].y;
      tv.z = pcp[i].z*ct[j].w - pcp[i].w*ct[j].z;
      pkn_multiSolvePLUQd ( 3, DF, P, Q, 1, 1, &tv.x );
      tv.z -= 1.0;
      rr = DotProduct3d ( &tv, &tv );
      if ( rr > k3 ) {
        if ( ks+rr >= BALL ) goto nodiff;
                             else k3 = rr;
      }
    }
        /* passed */
  SetPoint3d ( kk, k1, k2, k3 );
  pkv_SetScratchMemTop ( sp );
  return 1;

nodiff:
  pkv_SetScratchMemTop ( sp );
  return 0;

failure:
  pkv_SetScratchMemTop ( sp );
  return -1;
#undef BALL
} /*_rbi_PCHodogTestRd*/

void _rbi_ComputeFDFRd ( int n, int m, point4d *pcp, int d, point4d *ccp,
                         point3d *x, point3d *F, double *DF )
{
  vector4d p, pu, pv, c, ct;

  mbs_BCHornerDerP4d ( n, m, pcp, x->x, x->y, &p, &pu, &pv );
  mbs_BCHornerDerC4d ( d, ccp, x->z, &c, &ct );
        /* compute the value of F */
  F->x = p.x*c.w - p.w*c.x;
  F->y = p.y*c.w - p.w*c.y;
  F->z = p.z*c.w - p.w*c.z;
        /* compute DF, the derivative of F */
  DF[0] = pu.x*c.w - pu.w*c.x;
  DF[1] = pv.x*c.w - pv.w*c.x;
  DF[2] = p.x*ct.w - p.w*ct.x;
  DF[3] = pu.y*c.w - pu.w*c.y;
  DF[4] = pv.y*c.w - pv.w*c.y;
  DF[5] = p.y*ct.w - p.w*ct.y;
  DF[6] = pu.z*c.w - pu.w*c.z;
  DF[7] = pv.z*c.w - pv.w*c.z;
  DF[8] = p.z*ct.w - p.w*ct.z;
} /*_rbi_ComputeFDFRd*/

boolean _rbi_PCNewtonRd ( RBezPatchTreedp ptree, RBezCurveTreedp ctree,
                          pc_domain *ha, point3d *x )
{
#define MAX_ITER 7
#define EPS      1.0e-10
  int     n, m, d;
  RBezPatchTreeVertexdp pvertex;
  RBezCurveTreeVertexdp cvertex;
  point4d *pcp, *ccp;
  int     i;
  point3d F;
  double  DF[9], fx, fy, fz, dist;
  int     P[3], Q[3];

  n = ptree->n;
  m = ptree->m;
  pvertex = ha->pvertex;
  pcp = pvertex->ctlpoints;
  fx = pvertex->u1 - pvertex->u0;
  fy = pvertex->v1 - pvertex->v0;
  d = ctree->degree;
  cvertex = ha->cvertex;
  ccp = cvertex->ctlpoints;
  fz = cvertex->t1 - cvertex->t0;

  SetPoint3d ( x, 0.5, 0.5, 0.5 );
  for ( i = 0; i < MAX_ITER; i++ ) {
    _rbi_ComputeFDFRd ( n, m, pcp, d, ccp, x, &F, DF );
    if ( !pkn_GaussDecomposePLUQd ( 3, DF, P, Q ) )
      return false;
    pkn_multiSolvePLUQd ( 3, DF, P, Q, 1, 1, &F.x );
    SubtractPoints3d ( x, &F, x );
    if ( x->x < -0.25 || x->x > 1.25 ||
         x->y < -0.25 || x->y > 1.25 ||
         x->z < -0.25 || x->z > 1.25 )
      return false;
    F.x *= fx;  F.y *= fy;  F.z *= fz;
    dist = DotProduct3d ( &F, &F );
    if ( dist < EPS )
      return true;
  }
  return false;
#undef EPS
#undef MAX_ITER
} /*_rbi_PCNewtonRd*/

int _rbi_AssignPCIntersParamRd ( RBezPatchTreedp ptree, RBezCurveTreedp ctree,
                                 pc_domain *ha, double t, point3d *x, point3d *kk,
                                 vector4d *ip )
{
#define TOL  5.0e-5
#define BALL 1.0  /* must be <= 1.0 */
  int     n, m, d;
  RBezPatchTreeVertexdp pvertex;
  RBezCurveTreeVertexdp cvertex;
  double  out, ax, ay, az, rx, ry, rz;

  n = ptree->n;
  m = ptree->m;
  pvertex = ha->pvertex;
  d = ctree->degree;
  cvertex = ha->cvertex;

  ax = fabs ( 2.0*x->x - 1.0 );
  ay = fabs ( 2.0*x->y - 1.0 );
  az = fabs ( 2.0*x->z - 1.0 );
  out = max ( ax, ay );
  out = 0.5*(max ( out, az ) - 1.0);
  if ( out >= TOL ) {  /* the point sticks out too far */
    if ( ax > 1.0 ) {
      rx = pkv_rpower ( ax, n-1 );
      ry = pkv_rpower ( ax, m );
      rz = pkv_rpower ( ax, d );
    }
    else
      rx = ry = rz = 1.0;
    if ( ay > 1.0 ) {
      rx *= pkv_rpower ( ay, n );
      ry *= pkv_rpower ( ay, m-1 );
      rz *= pkv_rpower ( ay, d );
    }
    if ( az > 1.0 ) {
      rx *= pkv_rpower ( az, n );
      ry *= pkv_rpower ( az, m );
      rz *= pkv_rpower ( az, d-1 );
    }
    if ( rz+ry+rz >= BALL )
      return 2;
    else
      return 0;
  }
  else {  /* the point sticks out just a little */
    if ( out > 0.0 ) {
      if ( x->x < 0.0 ) x->x = 0.0; else if ( x->x > 1.0 ) x->x = 1.0;
      if ( x->y < 0.0 ) x->y = 0.0; else if ( x->y > 1.0 ) x->y = 1.0;
      if ( x->z < 0.0 ) x->z = 0.0; else if ( x->z > 1.0 ) x->z = 1.0;
    }
    switch ( ha->nc ) {
  case 0:  case 1:
      ip->z = (1.0-x->z)*cvertex->t0 + x->z*cvertex->t1;
      ip->w = t;
      ip->x = (1.0-x->x)*pvertex->u0 + x->x*pvertex->u1;
      ip->y = (1.0-x->y)*pvertex->v0 + x->y*pvertex->v1;
      break;
  case 2:  case 3:
      ip->z = t;
      ip->w = (1.0-x->z)*cvertex->t0 + x->z*cvertex->t1;
      ip->x = (1.0-x->x)*pvertex->u0 + x->x*pvertex->u1;
      ip->y = (1.0-x->y)*pvertex->v0 + x->y*pvertex->v1;
      break;
  case 4:  case 5:
      ip->x = (1.0-x->z)*cvertex->t0 + x->z*cvertex->t1;
      ip->y = t;
      ip->z = (1.0-x->x)*pvertex->u0 + x->x*pvertex->u1;
      ip->w = (1.0-x->y)*pvertex->v0 + x->y*pvertex->v1;
      break;
  case 6:  case 7:
      ip->x = t;
      ip->y = (1.0-x->z)*cvertex->t0 + x->z*cvertex->t1;
      ip->z = (1.0-x->x)*pvertex->u0 + x->x*pvertex->u1;
      ip->w = (1.0-x->y)*pvertex->v0 + x->y*pvertex->v1;
      break;
    }
    return 1;
  }
#undef TOL
#undef BALL
} /*_rbi_AssignPCIntersParamRd*/

void _rbi_IdentifySolutionsRd ( vector4d *ppar, int *nppar )
{
#define EPS 1.0e-9
  int      i, np;
  vector4d v;
  double   d;

  np = *nppar;
  for ( i = 0; i < np; i++ ) {
    SubtractPoints4d ( (point4d*)&ppar[i], (point4d*)&ppar[np], &v );
    d = DotProduct4d ( &v, &v );
    if ( d <= EPS )
      return;
  }
  (*nppar)++;
#undef EPS
} /*_rbi_IdentifySolutionsRd*/

void _rbi_FindRBezCurvePCentd ( RBezCurveTreedp tree,
                                RBezCurveTreeVertexdp vertex )
{
  mbs_BCHornerC3Rd ( tree->degree, vertex->ctlpoints, 0.5, &vertex->ccent );
} /*_rbi_FindRBezCurvePCentd*/

RBezCurveTreeVertexdp _rbi_GetRBezCurveLeftVertexd ( RBezCurveTreedp tree,
                                                     RBezCurveTreeVertexdp vertex )
{
  boolean               leaf;
  RBezCurveTreeVertexdp lv;

  leaf = !vertex->left;
  lv = rbez_GetRBezCurveLeftVertexd ( tree, vertex );
  if ( leaf && lv ) {
    _rbi_FindRBezCurvePCentd ( tree, lv );
    _rbi_FindRBezCurvePCentd ( tree, vertex->right );
  }
  return lv;
} /*_rbi_GetRBezCurveLeftVertexd*/

RBezCurveTreeVertexdp _rbi_GetRBezCurveRightVertexd ( RBezCurveTreedp tree,
                                                      RBezCurveTreeVertexdp vertex )
{
  boolean               leaf;
  RBezCurveTreeVertexdp rv;

  leaf = !vertex->left;
  rv = rbez_GetRBezCurveRightVertexd ( tree, vertex );
  if ( leaf && rv ) {
    _rbi_FindRBezCurvePCentd ( tree, vertex->left );
    _rbi_FindRBezCurvePCentd ( tree, rv );
  }
  return rv;
} /*_rbi_GetRBezCurveRightVertexd*/

boolean _rbi_DividePCDomainRd ( RBezPatchTreedp ptree, RBezCurveTreedp ctree,
                                pc_domain *ha, pc_domain *hb )
{
  double                dp, dc;
  RBezPatchTreeVertexdp pvert;
  RBezCurveTreeVertexdp cvert;

  pvert = ha->pvertex;
  cvert = ha->cvertex;
  dp = pvert->maxder;
  dc = cvert->maxder;
  if ( dp > dc ) {
    hb->pvertex = rbez_GetRBezRightVertexd ( ptree, pvert );
    ha->pvertex = rbez_GetRBezLeftVertexd ( ptree, pvert );
    if ( !hb->pvertex || !ha->pvertex )
      goto failure;
    hb->cvertex = cvert;
  }
  else {
    hb->pvertex = pvert;
    hb->cvertex = _rbi_GetRBezCurveRightVertexd ( ctree, cvert );
    ha->cvertex = _rbi_GetRBezCurveLeftVertexd ( ctree, cvert );
    if ( !hb->cvertex || !ha->cvertex )
      goto failure;
  }
  hb->level = ++ ha->level;
  hb->nc = ha->nc;  /* probably unnecessary */
  return true;

failure:
  return false;
} /*_rbi_DividePCDomainRd*/

int _rbi_OtherSingularityRd ( pci_key *pcii,
                              RBezPatchTreedp ptree, RBezCurveTreedp ctree,
                              pc_domain *hc, double tsing )
{
#define STACK_LENGTH 32
#define EPS           0.01
  void      *sp;
  pc_domain ha, *pcdstack;
  int       pcdsp;
  RBezCurveTreeVertexdp cvertex;
  vector3d  kk, x;
  boolean   result;

  sp = pkv_GetScratchMemTop ();
  pcdstack = pkv_GetScratchMem ( STACK_LENGTH*sizeof(pc_domain) );
  if ( !pcdstack )
    goto failure;
        /* push the cube on the stack */
  pcdsp = 0;
  pcdstack[0] = *hc;
  result = false;
  do {
    ha = pcdstack[pcdsp--];
    if ( _rbi_PCConvexHullTestRd ( ptree, ctree, &ha ) ) {
      switch ( _rbi_PCHodogTestRd ( ptree, ctree, &ha, &kk ) ) {
    case -1:
        goto failure;
    case 1:
        if ( _rbi_PCNewtonRd ( ptree, ctree, &ha, &x ) )
          goto not_quite;
        else
          goto subdivide;
    case 0:
subdivide:
        if ( ha.level < STACK_LENGTH ) {
          pcdsp += 2;
          if ( !_rbi_DividePCDomainRd ( ptree, ctree,
                                        &pcdstack[pcdsp-1], &pcdstack[pcdsp] ) )
            goto failure;
        }
        else {
          cvertex = ha.cvertex;
          result = fabs ( tsing - 0.5*(cvertex->t1-cvertex->t0) ) >= EPS;
          break;
        }
      }
    }
  } while ( pcdsp >= 0 );
  if ( !result )
    goto not_quite;

  pkv_SetScratchMemTop ( sp );
  return 1;

not_quite:
  pkv_SetScratchMemTop ( sp );
  return 0;

failure:
  pkv_SetScratchMemTop ( sp );
  return -1;
#undef EPS
#undef STACK_LENGTH
} /*_rbi_OtherSingularityRd*/

int _rbi_FindPCIntersRd ( pci_key *pcii, int cdeg, point4d *ccp,
                          int nc, int level,
                          vector4d *ppar, int *nppar )
{
  void                  *sp;
  RBezPatchTreedp       ptree;
  RBezCurveTreedp       ctree;
  RBezCurveTreeVertexdp cvertex;
  int                   pcdsp;     /* stack pointer */
  pc_domain             *pcdstack; /* stack */
  pc_domain             ha;
  byte                  nsing;
  point3d               x, kk;

  sp = pkv_GetScratchMemTop ();
  ctree = NULL;
  ptree = pcii->ptree;
  pcdstack = pkv_GetScratchMem ( (MAX_LEVEL+2)*sizeof(pc_domain) );
  if ( !pcdstack )
    goto failure;
  ctree = rbez_NewRBezCurveTreed ( 0, cdeg, pcii->s0, pcii->s1, 0.0, ccp );
  if ( !ctree )
    goto failure;
  _rbi_FindRBezCurvePCentd ( ctree, ctree->root );
  nsing = 0;
        /* make the stack empty */
  pcdsp = 0;
        /* push the initial cube */
  pcdstack[pcdsp].pvertex = pcii->pvertex;
  pcdstack[pcdsp].cvertex = ctree->root;
  pcdstack[pcdsp].nc = nc;
  pcdstack[pcdsp].level = level;
  do {
        /* pop the cube */
    ha = pcdstack[pcdsp--];
    if ( _rbi_PCConvexHullTestRd ( ptree, ctree, &ha ) ) {
      switch ( _rbi_PCHodogTestRd ( ptree, ctree, &ha, &kk ) ) {
    case -1:
        goto failure;
    case 1:
        if ( !_rbi_PCNewtonRd ( ptree, ctree, &ha, &x ) )
          goto subdivide;
        switch ( _rbi_AssignPCIntersParamRd ( ptree, ctree, &ha, pcii->t,
                                            &x, &kk, &ppar[*nppar] ) ) {
      case 0:
          break;
      case 1:
          _rbi_IdentifySolutionsRd ( ppar, nppar );
          if ( *nppar > 2 )
            goto subdivide;
          break;
      case 2:
          goto subdivide;
        }
        break;
    case 0:
subdivide:
        if ( ha.level < MAX_LEVEL ) {
          pcdsp += 2;
          if ( !_rbi_DividePCDomainRd ( ptree, ctree,
                                        &pcdstack[pcdsp-1], &pcdstack[pcdsp] ) )
            goto failure;
        }
        else {  /* found a singular solution */
          if ( ++nsing > 2 ) {
            cvertex = ha.cvertex;
            switch ( _rbi_OtherSingularityRd ( pcii, ptree, ctree, &ha,
                                               0.5*(cvertex->t1-cvertex->t0) ) ) {
          case -1:
              goto failure;
          case 0:
              goto not_quite;
          default:
              break;
            }
          }
          SetPoint3d ( &x, 0.5, 0.5, 0.5 );
          SetPoint3d ( &kk, 0.0, 0.0, 0.0 );
          if ( _rbi_AssignPCIntersParamRd ( ptree, ctree, &ha, pcii->t,
                                       &x, &kk, &ppar[*nppar] ) == -1 )
            goto failure;
          _rbi_IdentifySolutionsRd ( ppar, nppar );
          if ( *nppar > 2 )
            goto not_quite;
        }
      }
    }
  } while ( pcdsp >= 0 );

  rbez_DestroyRBezCurveTreed ( ctree );
  pkv_SetScratchMemTop ( sp );
  return true;

not_quite:
  if ( ctree )
    rbez_DestroyRBezCurveTreed ( ctree );
  pkv_SetScratchMemTop ( sp );
  return 0;

failure:
  if ( ctree )
    rbez_DestroyRBezCurveTreed ( ctree );
  pkv_SetScratchMemTop ( sp );
  return -1;
} /*_rbi_FindPCIntersRd*/

int _rbi_TryIntersArcRd ( RBezPatchTreedp tree1, RBezPatchTreedp tree2,
                          rbiHyperDomaind *hc,  vector4d *ppar,
                          int *singcurve )
{
  void                  *sp;
  int                   nppar, i;
  RBezPatchTreeVertexdp vert1, vert2;
  pci_key               pcii;
  int                   cdeg;
  point4d               *ccp;

  sp = pkv_GetScratchMemTop ();
  cdeg = max ( tree1->n, tree1->m );
  cdeg = max ( cdeg, tree2->n );
  cdeg = max ( cdeg, tree2->m );
  ccp = pkv_GetScratchMem ( (cdeg+1)*sizeof(point4d) );
  if ( !ccp )
    return -1;
  vert1 = hc->vert1;
  vert2 = hc->vert2;
  *singcurve = 0;
  nppar = 0;
        /* intersect boundary curves of the second patch with the first patch */
  pcii.ptree = tree1;
  pcii.pvertex = vert1;
  pcii.cdir = 0;
  for ( i = 0; i < 4; i++ ) {
    _rbi_ExtractCurveRd ( tree2->n, tree2->m, vert2->ctlpoints, i, &cdeg, ccp );
    switch ( i ) {
  case 0:
      pcii.s0 = vert2->u0;  pcii.s1 = vert2->u1;  pcii.t  = vert2->v0;
      break;
  case 1:
      pcii.s0 = vert2->u0;  pcii.s1 = vert2->u1;  pcii.t  = vert2->v1;
      break;
  case 2:
      pcii.s0 = vert2->v0;  pcii.s1 = vert2->v1;  pcii.t  = vert2->u0;
      break;
  case 3:
      pcii.s0 = vert2->v0;  pcii.s1 = vert2->v1;  pcii.t  = vert2->u1;
      break;
    }
    switch ( _rbi_FindPCIntersRd ( &pcii, cdeg, ccp, i, hc->level,
                                   ppar, &nppar ) ) {
  case -1:
      goto failure;
  case 0:
      *singcurve = i+1;
      pkv_SetScratchMemTop ( sp );
      return 2;
  default:
      break;
    }
    if ( nppar > 2 ) return 2;
  }
        /* intersect boundary curves of the first patch with the second patch */
  pcii.ptree = tree2;
  pcii.pvertex = vert2;
  pcii.cdir = 0;
  for ( i = 0; i < 4; i++ ) {
    _rbi_ExtractCurveRd ( tree1->n, tree1->m, vert1->ctlpoints, i, &cdeg, ccp );
    switch ( i ) {
  case 0:
      pcii.s0 = vert1->u0;  pcii.s1 = vert1->u1;  pcii.t  = vert1->v0;
      break;
  case 1:
      pcii.s0 = vert1->u0;  pcii.s1 = vert1->u1;  pcii.t  = vert1->v1;
      break;
  case 2:
      pcii.s0 = vert1->v0;  pcii.s1 = vert1->v1;  pcii.t  = vert1->u0;
      break;
  case 3:
      pcii.s0 = vert1->v0;  pcii.s1 = vert1->v1;  pcii.t  = vert1->u1;
      break;
    }
    switch ( _rbi_FindPCIntersRd ( &pcii, cdeg, ccp, i, hc->level,
                                   ppar, &nppar ) ) {
  case -1:
      goto failure;
  case 0:
      *singcurve = i+5;
      pkv_SetScratchMemTop ( sp );
      return 2;
  default:
      break;
    }
    if ( nppar > 2 ) return 2;
  }
  pkv_SetScratchMemTop ( sp );
  return nppar == 2 ? 1 : 0;

failure:
  pkv_SetScratchMemTop ( sp );
  return -1;
} /*_rbi_TryIntersArcRd*/

int _rbi_FindIntersectionArcRd ( RBezPatchTreedp tree1, RBezPatchTreedp tree2,
                                 rbiHyperDomaind *hc, double epsilon2,
                                 rbiArcOutd *outproc, void *usrptr )
{
#define DELTA 0.25
  void                  *sp;
  int                   trial_res, singcurve;
  vector4d              ppar[4], *aux1, *aux2;
  vector3d              *nvcp1, *nvcp2;
  RBezPatchTreedp       etree;
  RBezPatchTreeVertexdp vert1, vert2, evert;
  rbiHyperDomaind       hq;

  sp = pkv_GetScratchMemTop ();
  etree = NULL;
  trial_res = _rbi_TryIntersArcRd ( tree1, tree2, hc, ppar, &singcurve );
  if ( trial_res < 0 ) {
    pkv_SetScratchMemTop ( sp );
    return -1;
  }
  if ( trial_res == 1 ) {
    switch ( _rbi_ProcessTheArcRd ( tree1, tree2, hc, ppar, epsilon2,
                                    outproc, usrptr ) ) {
  case -1:
      goto failure;
  case 0:
      trial_res = 2;
      break;
  default:
      break;
    }
    pkv_SetScratchMemTop ( sp );
    return trial_res;
  }

  if ( trial_res == 0 || singcurve == 0 ) {
    pkv_SetScratchMemTop ( sp );
    return trial_res;
  }
          /* there is a trouble with a singularity, the patch must be extended */
  hq = *hc;
  if ( singcurve < 5 ) {  /* patch 2 to extend */
    if ( !(aux1 = pkv_GetScratchMem ( tree2->cpsize )) )
      goto failure;
    vert2 = hc->vert2;
    etree = rbez_NewRBezPatchTreed ( 2, tree2->n, tree2->m,
                                     vert2->u0, vert2->u1, vert2->v0, vert2->v1,
                                     vert2->ctlpoints );
    if ( !etree )
      goto failure;
    hq.vert2 = evert = etree->root;
    aux2 = evert->ctlpoints;
    switch ( singcurve ) {
  case 1:  /* t = t0 of p2 */
      mbs_DivideBP4vd ( etree->n, etree->m, -DELTA, aux2, aux1 );
      evert->v0 = (1.0+DELTA)*vert2->v0 - DELTA*vert2->v1;
      break;
  case 2:  /* t = t1 of p2 */
      memcpy ( aux1, aux2, tree2->cpsize );
      mbs_DivideBP4vd ( etree->n, etree->m, 1.0+DELTA, aux1, aux2 );
      evert->v1 = -DELTA*vert2->v0 + (1.0+DELTA)*vert2->v1;
      break;
  case 3:  /* s = s0 of p2 */
      mbs_DivideBP4ud ( etree->n, etree->m, -DELTA, aux2, aux1 );
      evert->u0 = (1.0+DELTA)*vert2->u0 - DELTA*vert2->u1;
      break;
  case 4:  /* s = s1 of p2 */
      memcpy ( aux1, aux2, tree2->cpsize );
      mbs_DivideBP4ud ( etree->n, etree->m, 1.0+DELTA, aux1, aux2 );
      evert->u1 = -DELTA*vert2->u0 + (1.0+DELTA)*vert2->u1;
      break;
    }
    nvcp1 = _rbi_GetRBezNVPatchd ( tree1, hq.vert1 );
    nvcp2 = _rbi_GetRBezNVPatchd ( etree, evert );
    if ( !nvcp1 || !nvcp2 )
      goto failure;
    hq.nvpassed = _rbi_NormalTestRd ( tree1->vn, tree1->vm,
                                      nvcp1, &hq.vert1->nvcent,
                                      etree->vn, etree->vm,
                                      nvcp2, &evert->nvcent );
    if ( hq.nvpassed )
      trial_res = _rbi_FindIntersectionArcRd ( tree1, etree, &hq,
                                               epsilon2, outproc, usrptr );
    else
      trial_res = 2;
  }
  else {  /* patch 1 to extend */
    if ( !(aux1 = pkv_GetScratchMem ( tree1->cpsize )) )
      goto failure;
    vert1 = hc->vert1;
    etree = rbez_NewRBezPatchTreed ( 2, tree1->n, tree1->m,
                                     vert1->u0, vert1->u1, vert1->v0, vert1->v1,
                                     vert1->ctlpoints );
    if ( !etree )
      goto failure;
    hq.vert1 = evert = etree->root;
    aux2 = evert->ctlpoints;
    switch ( singcurve ) {
  case 5:  /* v = v0 of p1 */
      mbs_DivideBP4vd ( etree->n, etree->m, -DELTA, aux2, aux1 );
      evert->v0 = (1.0+DELTA)*vert1->v0 - DELTA*vert1->v1;
      break;
  case 6:  /* v = v1 of p1 */
      memcpy ( aux1, aux2, tree1->cpsize );
      mbs_DivideBP4vd ( etree->n, etree->m, 1.0+DELTA, aux1, aux2 );
      evert->v1 = -DELTA*vert1->v0 + (1.0+DELTA)*vert1->v1;
      break;
  case 7:  /* u = u0 of p1 */
      mbs_DivideBP4ud ( etree->n, etree->m, -DELTA, aux2, aux1 );
      evert->u0 = (1.0+DELTA)*vert1->u0 - DELTA*vert1->u1;
      break;
  case 8:  /* u = v1 of p1 */
      memcpy ( aux1, aux2, tree1->cpsize );
      mbs_DivideBP4ud ( etree->n, etree->m, 1.0+DELTA, aux1, aux2 );
      evert->u1 = -DELTA*vert1->u0 + (1.0+DELTA)*vert1->u1;
      break;
    }
    nvcp1 = _rbi_GetRBezNVPatchd ( etree, evert );
    nvcp2 = _rbi_GetRBezNVPatchd ( tree2, hq.vert2 );
    if ( !nvcp1 || !nvcp2 )
      goto failure;
    hq.nvpassed = _rbi_NormalTestRd ( etree->vn, etree->vm,
                                      nvcp1, &evert->nvcent,
                                      tree2->vn, tree2->vm,
                                      nvcp2, &hq.vert2->nvcent );
    if ( hq.nvpassed )
      trial_res = _rbi_FindIntersectionArcRd ( etree, tree2, &hq,
                                               epsilon2, outproc, usrptr );
    else
      trial_res = 2;
  }
  rbez_DestroyRBezPatchTreed ( etree );
  pkv_SetScratchMemTop ( sp );
  return trial_res;

failure:
  if ( etree )
    rbez_DestroyRBezPatchTreed ( etree );
  pkv_SetScratchMemTop ( sp );
  return -1;
#undef DELTA
} /*_rbi_FindIntersectionArcRd*/

/* ////////////////////////////////////////////////////////////////////////// */
boolean _rbi_DivideDomainRd ( RBezPatchTreedp tree1, RBezPatchTreedp tree2,
                              rbiHyperDomaind *hca, rbiHyperDomaind *hcb )
{
  RBezPatchTreeVertexdp vert1, vert2;
  boolean               divide_first;

        /* choose the patch to subdivide, spline first */
  vert1 = hca->vert1;
  vert2 = hca->vert2;
  if ( vert1->ctlpoints ) {
    if ( vert2->ctlpoints )
      goto by_length;
    else
      divide_first = false;
  }
  else if ( vert2->ctlpoints )
    divide_first = true;
  else {
by_length:
    divide_first = vert1->maxder > vert2->maxder;
  }
        /* divide the hypercube */
  hcb->level = ++hca->level;
  if ( divide_first ) {
    hcb->vert1 = rbez_GetRBezLeftVertexd ( tree1, hca->vert1 );
    hca->vert1 = rbez_GetRBezRightVertexd ( tree1, hca->vert1 );
    if ( !hcb->vert1 || !hca->vert1 )
      goto failure;
    hcb->vert2 = hca->vert2;
  }
  else {
    hcb->vert2 = rbez_GetRBezLeftVertexd ( tree2, hca->vert2 );
    hca->vert2 = rbez_GetRBezRightVertexd ( tree2, hca->vert2 );
    if ( !hcb->vert2 || !hca->vert2 )
      goto failure;
    hcb->vert1 = hca->vert1;
  }
  hcb->nvpassed = hca->nvpassed;
  return true;

failure:
  return false;
} /*_rbi_DivideDomainRd*/

/* ////////////////////////////////////////////////////////////////////////// */
boolean rbi_FindRBezIntersectiond ( int n1, int m1, point4d *p1,
                                    int n2, int m2, point4d *p2,
                                    double epsilon, byte maxlevel,
                                    rbiArcOutd *outproc, void *usrptr )
{
  void *sp;
  RBezPatchTreedp tree1, tree2;  /* trees of patch subdivision */
  rbiHyperDomaind *hcstack, hc;
  int             hcsp;          /* hypercube stack poinnter */
  vector3d        *nvcp1, *nvcp2;

  sp = pkv_GetScratchMemTop ();
  tree1 = tree2 = NULL;
  tree1 = rbez_NewRBezPatchTreed ( 1, n1, m1, 0.0, 1.0, 0.0, 1.0, p1 );
  tree2 = rbez_NewRBezPatchTreed ( 2, n2, m2, 0.0, 1.0, 0.0, 1.0, p2 );
  if ( !tree1 || !tree2 )
    goto failure;
  maxlevel = min ( maxlevel, MAX_LEVEL );
  epsilon *= epsilon;

        /* prepare the hypercube subdivision stack */
        /* and push the entire unit hypercube */
  hcstack = pkv_GetScratchMem ( (MAX_LEVEL+2)*sizeof(rbiHyperDomaind) );
  if ( !hcstack )
    goto failure;
  hcsp = 0;
  hcstack[hcsp].vert1 = tree1->root;
  hcstack[hcsp].vert2 = tree2->root;
  hcstack[hcsp].level = 0;
  hcstack[hcsp].nvpassed = false;
        /* now the numerical computation, with the hypercube subdivision */
  do {
          /* pop the hypercube */
    hc = hcstack[hcsp--];
    if ( hc.vert1->ctlpoints == NULL || hc.vert2->ctlpoints == NULL )
      goto subdivide;
    if ( _rbi_PatchConvexHullTestRd ( n1, m1, hc.vert1->ctlpoints,
                                      &hc.vert1->pcent, &hc.vert1->nvcent,
                                      n2, m2, hc.vert2->ctlpoints,
                                      &hc.vert2->pcent, &hc.vert2->nvcent ) ) {
      if ( !hc.nvpassed ) {
        nvcp1 = _rbi_GetRBezNVPatchd ( tree1, hc.vert1 );
        nvcp2 = _rbi_GetRBezNVPatchd ( tree2, hc.vert2 );
        if ( !nvcp1 || !nvcp2 )
          goto failure;
        hc.nvpassed = _rbi_NormalTestRd ( tree1->vn, tree1->vm,
                                          nvcp1, &hc.vert1->nvcent,
                                          tree2->vn, tree2->vm,
                                          nvcp2, &hc.vert2->nvcent );
      }
      if ( hc.nvpassed ) {
        switch ( _rbi_FindIntersectionArcRd ( tree1, tree2, &hc,
                                              epsilon, outproc, usrptr ) ) {
 case  0:
 case  1: break;
 case  2: goto subdivide;
 case -1: goto failure;
        }
      }
      else {
subdivide:
        if ( hc.level < MAX_LEVEL ) {
          hcsp += 2;
          if ( !_rbi_DivideDomainRd ( tree1, tree2,
                                      &hcstack[hcsp-1], &hcstack[hcsp] ) )
            goto failure;
        }
      }
    }
  } while ( hcsp >= 0 );

  rbez_DestroyRBezPatchTreed ( tree1 );
  rbez_DestroyRBezPatchTreed ( tree2 );
  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  if ( tree1 ) rbez_DestroyRBezPatchTreed ( tree1 );
  if ( tree2 ) rbez_DestroyRBezPatchTreed ( tree2 );
  pkv_SetScratchMemTop ( sp );
  return false;
} /*rbi_FindRBezIntersectiond*/

