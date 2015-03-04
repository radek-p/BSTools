
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2015                                  */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "pkvaria.h"
#include "pkgeom.h"
#include "multibs.h"


#define EPS 1.0e-6

typedef struct {
    short int deg;     /* degree */
    float     *coeff;  /* data array */
  } _mbs_BezPolyf;

typedef struct {
    float t;              /* clipping line parameter */
    int   ein, eout;      /* indices of polycurve edges */
    int   opv;            /* index of the other vertex making a pair */
    char  icase, icase1;  /* 0 if on the clipping line, -1 or +1 if not */
    char  tag;            /* 0 if valid, 1 else */
    char  pad;
  } _mbs_polycurve_vertexf;

typedef struct {
    int   v0, v1;   /* indices of the vertices */
    byte  cpdimen;  /* 2, polynomial, or 3, rational */
    byte  degree;   /* Bezier curve degree, from 1 to 61 */
    char  icase;    /* 0 if on the clipping line, */
                    /*-1 or +1 if on one side, 2 else */
    char  tag;      /* 0 if valid, 1 else */
    float *cp;      /* control points: cpdimen*(degree+1) numbers */
  } _mbs_polycurve_edgef;

typedef struct {
    float t;        /* clipping line parameter */
    int   v;        /* vertex index */
  } _mbs_div_vertexf;

static boolean _mbs_ReallocPolycurveTab ( int *tablength, int amount,
                                          _mbs_polycurve_vertexf **cv,
                                          _mbs_polycurve_edgef **ce )
{
  _mbs_polycurve_vertexf *ncv;
  _mbs_polycurve_edgef   *nce;
  int                    tl;

  tl = *tablength + amount;
  PKV_MALLOC ( ncv, tl*sizeof(_mbs_polycurve_vertexf) );
  PKV_MALLOC ( nce, tl*sizeof(_mbs_polycurve_edgef) );
  if ( ncv && nce ) {
    memcpy ( ncv, *cv, *tablength*sizeof(_mbs_polycurve_vertexf) );
    memcpy ( nce, *ce, *tablength*sizeof(_mbs_polycurve_edgef) );
    *tablength = tl;
    PKV_FREE ( *cv );
    PKV_FREE ( *ce );
    *cv = ncv;
    *ce = nce;
    return true;
  }
  else {
    if ( ncv ) PKV_FREE ( ncv );
    if ( nce ) PKV_FREE ( nce );
    return false;
  }
} /*_mbs_ReallocPolycurveTab*/

static boolean _mbs_InsertLineSegmentf ( int *tablength, int *narcs,
                                         _mbs_polycurve_vertexf **cv,
                                         _mbs_polycurve_edgef **ce,
                                         point2f *a, point2f *b,
                                         int fvind, boolean closing )
{
  _mbs_polycurve_edgef *edg;
  int                  n;

  n = *narcs;
  if ( a->x == b->x && a->y == b->y ) { /* no need to insert anything */
    if ( closing ) {         /* but the polycurve may have to be closed */
      edg = &(*ce)[n-1];
      edg->v1 = fvind;
      (*cv)[fvind].ein = n-1;
    }
    return true;
  }
  if ( n == *tablength ) {
    if ( !_mbs_ReallocPolycurveTab ( tablength, 16, cv, ce ) )
      return false;
  }
  edg = &(*ce)[n];
  PKV_MALLOC ( edg->cp, 4*sizeof(float) );
  if ( !edg->cp )
    return false;
  edg->v0 = n;
  edg->v1 = closing ? fvind : n+1;
  edg->cpdimen = 2;
  edg->degree = 1;
  memcpy ( edg->cp, a, 2*sizeof(float) );
  memcpy ( &edg->cp[2], b, 2*sizeof(float) );
  (*cv)[n].eout = (*cv)[edg->v1].ein = n;
  (*narcs) ++;
  return true;
} /*_mbs_InsertLineSegmentf*/

static boolean _mbs_InsertBezierArcf ( int *tablength, int *narcs,
                                       _mbs_polycurve_vertexf **cv,
                                       _mbs_polycurve_edgef **ce,
                                       byte dim, byte deg, const float *cp )
{
  _mbs_polycurve_edgef *edg;
  int                  n, size;

  n = *narcs;
  if ( n == *tablength ) {
    if ( !_mbs_ReallocPolycurveTab ( tablength, 16, cv, ce ) )
      return false;
  }
  edg = &(*ce)[n];
  size = (deg+1)*dim*sizeof(float);
  PKV_MALLOC ( edg->cp, size );
  if ( !edg->cp )
    return false;
  edg->v0 = n;
  edg->v1 = n+1;
  edg->cpdimen = dim;
  edg->degree = deg;
  memcpy ( edg->cp, cp, size );
  (*cv)[n].eout = (*cv)[edg->v1].ein = n;
  (*narcs) ++;
  return true;
} /*_mbs_InsertBezierArcf*/

static boolean _mbs_DivideBezierArcf ( int *tablength, int *narcs,
                                       _mbs_polycurve_vertexf **cv,
                                       _mbs_polycurve_edgef **ce,
                                       int e, float t )
{
  _mbs_polycurve_edgef *ed0, *ed1;
  int                  n, size;
  short int            dim, deg;

  n = *narcs;
  if ( n == *tablength ) {
    if ( !_mbs_ReallocPolycurveTab ( tablength, 16, cv, ce ) )
      return false;
  }
  ed0 = &(*ce)[e];
  ed1 = &(*ce)[n];
  dim = ed0->cpdimen;
  deg = ed0->degree;
  size = (deg+1)*dim*sizeof(float);
  PKV_MALLOC ( ed1->cp, size );
  if ( !ed1->cp )
    return false;
  ed1->cpdimen = dim;
  ed1->degree = deg;
  ed1->v1 = ed0->v1;
  ed1->v0 = ed0->v1 = n;
  memcpy ( ed1->cp, ed0->cp, size );
  mbs_multiDivideBezCurvesf ( deg, 1, dim, 0, t, ed1->cp, ed0->cp );
  (*cv)[n].ein = e;
  (*cv)[n].eout = n;
  (*cv)[n].icase = 0;
  (*cv)[ed1->v1].ein = n;
  (*narcs) ++;
  return true;
} /*_mbs_DivideBezierArcf*/

static void _mbs_RemovePolycurveVertexf ( _mbs_polycurve_vertexf *cv,
                                          _mbs_polycurve_edgef *ce, int v )
{
  int   e0, e1, v0, v1;
  float *cp0, *cp1;

  e0 = cv[v].ein;
  e1 = cv[v].eout;
  v0 = ce[e0].v0;
  v1 = ce[e1].v1;
  if ( v0 != v1 ) {
        /* find the edge end points */
    cp0 = ce[e0].cp;
    if ( ce[e0].cpdimen == 3 ) {
      cp0[0] /= cp0[2];
      cp0[1] /= cp0[2];
    }
    cp1 = &ce[e1].cp[ce[e1].cpdimen*ce[e1].degree];
    if ( ce[e1].cpdimen == 2 )
      memcpy ( &cp0[2], cp1, 2*sizeof(float) );
    else {
      cp0[2] = cp1[0]/cp1[2];
      cp0[3] = cp1[1]/cp1[2];
    }
        /* the edge e0 is a line segment, with icase unchanged */
    ce[e0].cpdimen = 2;
    ce[e0].degree = 1;
    ce[e0].v1 = v1;
    cv[v1].ein = e0;
        /* mark the vertex and the outcoming edge invalid */
    ce[e1].tag = 1;
    cv[v].tag = 1;
  }
  else { /* remove the cycle of length 2 */
    ce[e0].tag = ce[e1].tag = 1;
    cv[v0].tag = cv[v].tag = 1;
  }
} /*_mbs_RemovePolycurveVertexf*/

static void _mbs_CleanupCyclef ( int narcs,
                                 _mbs_polycurve_vertexf *cv,
                                 _mbs_polycurve_edgef *ce )
{
  int i;

  for ( i = 0; i < narcs; i++ ) {
    cv[i].tag = ce[i].tag = 0;
    cv[i].icase1 = cv[i].icase;
  }
  for ( i = 0; i < narcs; i++ )
    if ( cv[i].tag == 0 && cv[i].icase == 0 ) {
      if ( ce[cv[i].ein].icase == 0 && ce[cv[i].eout].icase == 0 )
        _mbs_RemovePolycurveVertexf ( cv, ce, i );
    }
  for ( i = 0; i < narcs; i++ )
    if ( cv[i].tag == 0 && cv[i].icase == 0 ) {
      if ( ce[cv[i].eout].icase == 0 )
        cv[i].icase1 = ce[cv[i].eout].icase = ce[cv[i].ein].icase;
    }
} /*_mbs_CleanupCyclef*/

static float _mbs_apolyf ( void *usrptr, float t )
{
  _mbs_BezPolyf *poly;
  float         value;

  poly = (_mbs_BezPolyf*)usrptr;
  mbs_BCHornerC1f ( poly->deg, poly->coeff, t, &value );
  return value;
} /*_mbs_apolyf*/

static boolean _mbs_OutputPolycurvef ( int narcs,
                                       _mbs_polycurve_vertexf *cv,
                                       _mbs_polycurve_edgef *ce,
                                       char side,
                                       int *nel, mbs_polycurvef **cont )
{
  int            ncarcs, i, e, v0, v1;
  mbs_polycurvef *_cont;
  int            _nel;
  int            size;

  ncarcs = 0;
  for ( i = 0; i < narcs; i++ ) {
    if ( !ce[i].tag && ce[i].icase == side )
      ncarcs ++;
    if ( !cv[i].tag && cv[i].opv > i )
      ncarcs ++;
  }
  if ( !ncarcs ) {  /* empty polycurve */
    *nel = 0;
    *cont = NULL;
    return true;
  }

  PKV_MALLOC ( _cont, ncarcs*sizeof(mbs_polycurvef) );
  if ( !_cont )
    return false;
  _nel = 0;
  for ( i = 0; i < narcs; i++ ) {
    if ( !ce[i].tag && ce[i].icase == side ) {
      v0 = ce[i].v0;
      e = i;
      do {
            /* output the edge of the original polycurve */
        _cont[_nel].ident = _nel;
        _cont[_nel].cpdimen = ce[e].cpdimen;
        _cont[_nel].degree = ce[e].degree;
        _cont[_nel].lastknot = -1;
        _cont[_nel].knots = NULL;
        size = (ce[e].degree+1)*ce[e].cpdimen*sizeof(float);
        PKV_MALLOC ( _cont[_nel].points, size );
        if ( !_cont[_nel].points )
          goto failure;
        memcpy ( _cont[_nel].points, ce[e].cp, size );
        ce[e].tag = 2;
        v1 = ce[e].v1;
        _cont[_nel].closing = v1 == v0;
        _nel ++;
        e = cv[v1].eout;
        if ( v1 != v0 && cv[v1].opv >= 0 ) {
            /* output the line segment on the clipping line */
          _cont[_nel].ident = _nel;
          _cont[_nel].cpdimen = 2;
          _cont[_nel].degree = 1;
          _cont[_nel].lastknot = -1;
          _cont[_nel].knots = NULL;
          PKV_MALLOC ( _cont[_nel].points, 4*sizeof(float) );
          if ( !_cont[_nel].points )
            goto failure;
          if ( ce[e].cpdimen == 2 )
            memcpy ( _cont[_nel].points, ce[e].cp, 2*sizeof(float) );
          else {
            _cont[_nel].points[0] = ce[e].cp[0]/ce[e].cp[2];
            _cont[_nel].points[1] = ce[e].cp[1]/ce[e].cp[2];
          }
          v1 = cv[v1].opv;
          e = cv[v1].eout;
          if ( ce[e].cpdimen == 2 )
            memcpy ( &_cont[_nel].points[2], ce[e].cp, 2*sizeof(float) );
          else {
            _cont[_nel].points[2] = ce[e].cp[0]/ce[e].cp[2];
            _cont[_nel].points[3] = ce[e].cp[1]/ce[e].cp[2];
          }
          _cont[_nel].closing = v1 == v0;
          _nel ++;
        }
      } while ( v1 != v0 );
    }
  }
  *nel = _nel;
  *cont = _cont;
  return true;

failure:
  mbs_DeallocatePolycurvef ( _nel, _cont );
  return false;
} /*_mbs_OutputPolycurvef*/

boolean mbs_ClipPolycurvef ( int nelem, const mbs_polycurvef *poly,
                             point2f *p0, vector2f *v,
                             int *nelem1, mbs_polycurvef **poly1,
                             int *nelem2, mbs_polycurvef **poly2 )
{
#define CONVERTP(p,q) \
  { if ( dim == 2 ) memcpy ( &q, p, 2*sizeof(float) ); \
    else Point3to2f ( (point3f*)p, &q ); }

#define INSERTLS(a,b,c) \
  { if ( !_mbs_InsertLineSegmentf ( &tablength, &narcs, &cv, &ce, \
                                    &a, &b, fvind, c ) ) \
      goto failure; }

#define INSERTARC(cp) \
  { if ( !_mbs_InsertBezierArcf ( &tablength, &narcs, &cv, &ce, \
                                  dim, deg, cp ) ) \
      goto failure; }

  void                   *sp;
  short int              deg, dim, maxdeg;
  int                    lkn, narcs, ku, ncv;
  int                    tablength;
  _mbs_polycurve_vertexf *cv;
  _mbs_polycurve_edgef   *ce, *edg;
  _mbs_div_vertexf       *dv;
  float                  *points, *bsbuf;
  boolean                first;
  point2f                firstp, lastp, nextp;
  vector2f               w;
  int                    fvind;
  int                    i, j, vi, vj;
  float                  t, t0, t1, *p, *q;
  _mbs_BezPolyf          bpoly;
  boolean                nonzero, error;
  int                    ic0, ic1, veind;
  unsigned int           *vperm, *vstk, vsp;

  sp = pkv_GetScratchMemTop ();
  cv = NULL;
  ce = NULL;

        /* compute the number of arcs */
  maxdeg = 1;
  narcs = 0;
  for ( i = 0; i < nelem; i++ ) {
    dim = poly[i].cpdimen;
    if ( dim != 2 && dim != 3 )
      goto failure;
    deg = poly[i].degree;
    if ( deg > maxdeg ) maxdeg = deg;
    lkn = poly[i].lastknot;
    if ( poly[i].knots ) {  /* a B-spline curve - count its arcs */
      ku = mbs_NumKnotIntervalsf ( deg, lkn, poly[i].knots );
      narcs += ku;
    }
    else if ( lkn == -1 ) /* a Bezier curve - one arc */
      narcs ++;
    else if ( deg == 1 && lkn > 0 ) /* a polyline - count its line segments */
      narcs += lkn;
    else
      goto failure;
  }
  if ( !narcs )
    goto way_out;

        /* allocate arrays, if necessary, may be reallocated */
  tablength = narcs + nelem;
  PKV_MALLOC ( cv, tablength*sizeof(_mbs_polycurve_vertexf) );
  PKV_MALLOC ( ce, tablength*sizeof(_mbs_polycurve_edgef) );
  if ( !cv || !ce )
    goto failure;

        /* construct the internal polycurve representation */
  narcs = fvind = 0;
  first = true;
  for ( i = 0; i < nelem; i++ ) {
    dim = poly[i].cpdimen;
    deg = poly[i].degree;
    lkn = poly[i].lastknot;
    points = &poly[i].points[0];
    if ( poly[i].knots ) { /* a B-spline curve */
      ku = mbs_NumKnotIntervalsf ( deg, lkn, poly[i].knots );
      bsbuf = pkv_GetScratchMemf ( ku*dim*(deg+1) );
      if ( !bsbuf )
        goto failure;
      if ( !mbs_multiBSCurvesToBezf ( dim, 1, deg, lkn, poly[i].knots,
                                      0, points, &ku, NULL, NULL, 0, bsbuf ) )
        goto failure;
      if ( first ) CONVERTP ( bsbuf, firstp )
      else {
        CONVERTP ( bsbuf, nextp )
        INSERTLS ( lastp, nextp, false )
      }
      for ( j = 0; j < ku; j++ )
        INSERTARC ( &bsbuf[j*dim*(deg+1)] )
      CONVERTP ( &bsbuf[ku*dim*(deg+1)-dim], lastp )
      pkv_SetScratchMemTop ( bsbuf );
    }
    else if ( poly[i].lastknot == -1 ) { /* a Bezier curve */
      if ( first ) CONVERTP ( points, firstp )
      else {
        CONVERTP ( points, nextp )
        INSERTLS ( lastp, nextp, false )
      }
      INSERTARC ( points )
      CONVERTP ( &points[deg*dim], lastp )
    }
    else {  /* a polyline */
      if ( first ) CONVERTP ( points, firstp )
      else {
        CONVERTP ( points, nextp )
        INSERTLS ( lastp, nextp, false )
      }
      for ( j = 0; j < lkn-1; j++ )
        INSERTARC ( &points[j*dim] )
      CONVERTP ( &points[(lkn-1)*dim], lastp )
    }

    if ( poly[i].closing ) {
      INSERTLS ( lastp, firstp, true )
      fvind = narcs;
      first = true;
    }
    else
      first = false;
  }

        /* determine the vertices located on the clipping line */
  for ( i = 0; i < narcs; i++ ) {
    edg = &ce[cv[i].eout];
    points = &edg->cp[0];
    if ( edg->cpdimen == 2 )
      SubtractPoints2f ( (point2f*)points, p0, &w );
    else
      SetPoint2f ( &w, points[0] - p0->x*points[2], points[1] - p0->y*points[2] );
    cv[i].icase = pkv_signf ( DotProduct2f ( &w, v ) );
  }
        /* intersect the arcs with the clipping line */
  p = pkv_GetScratchMemf ( 2*(maxdeg+1) );
  if ( !p )
    goto failure;
  q = &p[maxdeg+1];
  bpoly.coeff = p;
  for ( i = 0; i < narcs; i++ ) {
        /* compute the polynomial coefficients */
    edg = &ce[i];
    bpoly.deg = deg = edg->degree;
    points = &edg->cp[0];
    if ( edg->cpdimen == 2 ) {
      for ( j = 0; j <= deg; j++ ) {
        SubtractPoints2f ( (point2f*)&points[2*j], p0, &w );
        p[j] = DotProduct2f ( &w, v );
      }
    }
    else {
      for ( j = 0; j <= deg; j++ ) {
        SetPoint2f ( &w, points[3*j] - p0->x*points[3*j+2],
                     points[3*j+1] - p0->y*points[3*j+2] );
        p[j] = DotProduct2f ( &w, v );
      }
    }
        /* compensate rounding errors */
    if ( (ic0 = cv[edg->v0].icase) == 0 ) p[0] = 0.0;
    if ( (ic1 = cv[edg->v1].icase) == 0 ) p[deg] = 0.0;
        /* find the smallest polynomial zero in (0,1) */
        /* other intersection points will be found later */
    if ( pkn_VarSignf ( deg+1, p ) ) {
      t0 = 0.0;  t1 = 1.0;
      for (;;) {
        if ( pkn_OneSignChangef ( deg+1, p, &nonzero ) ) {
          if ( nonzero && ic0 && ic1 ) {
            t = pkn_Illinoisf ( _mbs_apolyf, (void*)&bpoly, 0.0, 1.0, EPS, &error );
              /* split the arc */
            if ( !_mbs_DivideBezierArcf ( &tablength, &narcs, &cv, &ce,
                                          i, t0 + t*(t1-t0) ) )
              goto failure;
            for ( j = 0; j <= deg; j++ )
              if ( (edg->icase = pkv_signf ( p[j] )) )
                break;
          }
          else
            edg->icase = cv[edg->v0].icase + cv[edg->v1].icase;
          break;
        }
        else if ( nonzero ) {
          mbs_BisectBC1f ( deg, p, q );
          if ( pkn_VarSignf ( deg+1, q ) && ic0 ) {
            memcpy ( p, q, (deg+1)*sizeof(float) );
            t1 = 0.5*(t0+t1);
            ic1 = 1;
          }
          else if ( pkn_VarSignf ( deg+1, p ) && ic1 ) {
            t0 = 0.5*(t0+t1);
            ic0 = 1;
          }
          else {
            for ( j = 0; j <= deg; j++ )
              if ( (edg->icase = pkv_signf ( q[j] )) )
                break;
            break;
          }
        }
        else {  /* the edge is on the clipping line */
          edg->icase = 0;
          break;
        }
      } /* for */
    }
    else        /* the edge is on one side or on the line */
      edg->icase = ic0;
  }
  _mbs_CleanupCyclef ( narcs, cv, ce );

        /* order the vertices on the clipping line */
  ncv = 0;
  for ( i = 0; i < narcs; i++ ) {
    if ( cv[i].tag == 0 ) {
      ic0 = ce[cv[i].ein].icase;
      ic1 = ce[cv[i].eout].icase;
      if ( ic0 == ic1 )
        cv[i].icase1 = ic0;
      if ( !cv[i].icase1 )
        ncv ++;
    }
    cv[i].opv = -1;
  }
  if ( ncv > 0 ) {
    dv = pkv_GetScratchMem ( ncv*(sizeof(_mbs_div_vertexf)+2*sizeof(int)) );
    if ( !dv )
      goto failure;
    vperm = (unsigned int*)&dv[ncv];
    vstk = &vperm[ncv];
    ncv = 0;
    for ( i = 0; i < narcs; i++ )
      if ( cv[i].tag == 0 && !cv[i].icase1 ) {
        edg = &ce[cv[i].eout];
        points = &edg->cp[0];
        if ( edg->cpdimen == 2 )
          SubtractPoints2f ( (point2f*)points, p0, &w );
        else
          SetVector2f ( &w, points[0]/points[2]-p0->x, points[1]/points[2]-p0->y );
        dv[ncv].v = i;
        dv[ncv++].t = cv[i].t = det2f ( &w, v );
      }
    if ( ncv > 1 ) {
      for ( i = 0; i < ncv; i++ ) vperm[i] = i;
      if ( pkv_SortKernel ( sizeof(float), ID_IEEE754_FLOAT,
                            sizeof(_mbs_div_vertexf), 0,
                            ncv, &dv[0].t, vperm ) != SORT_OK )
        goto failure;
        /* make pairs of the vertices on the clipping line */
      vsp = veind = 0;
      ic0 = 0;
      for ( i = 0; i < ncv; i++ ) {
        vi = dv[vperm[i]].v;
        ic1 = ce[cv[vi].eout].icase;
        if ( ic1 ) {
          if ( vsp == 0 || ic1*veind > 0 ) { /* push vertex index */
            vstk[vsp++] = vi;
            ic0 = ic1;
          }
          else if ( vsp > 0 ) {  /* pop and make a pair */
            vj = vstk[--vsp];
            cv[vi].opv = vj;
            cv[vj].opv = vi;
          }
          veind += ic1;
        }
      }
    }
  }
        /* output the two polycurves */
  if ( nelem1 && poly1 ) {
    if ( !_mbs_OutputPolycurvef ( narcs, cv, ce, +1, nelem1, poly1 ) )
      goto failure;
  }
  if ( poly2 && poly2 ) {
    if ( !_mbs_OutputPolycurvef ( narcs, cv, ce, -1, nelem2, poly2 ) ) {
      if ( *nelem1 && *poly1 )
        mbs_DeallocatePolycurvef ( *nelem1, *poly1 );
      goto failure;
    }
  }

way_out:
        /* deallocate the internal polycurve */
  if ( cv ) PKV_FREE ( cv );
  if ( ce ) {
    for ( i = 0; i < narcs; i++ )
      PKV_FREE ( ce[i].cp );
    PKV_FREE ( ce );
  }
  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  if ( cv ) PKV_FREE ( cv );
  if ( ce ) {
    for ( i = 0; i < narcs; i++ )
      if ( ce[i].cp ) PKV_FREE ( ce[i].cp );
    PKV_FREE ( ce );
  }
  pkv_SetScratchMemTop ( sp );
  return false;
#undef CONVERTP
#undef REALLOCTAB
} /*mbs_ClipPolycurvef*/

