
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2005, 2013                            */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdio.h>  /* for testing */

#include "pkvaria.h"
#include "pkgeom.h"
#include "multibs.h"

/* ////////////////////////////////////////////////// */
/* Making a "grid" picture of a trimmed patch domain, */
/* or the patch itself                                */
static int _mbs_densd ( int n0, double a, double b, double h )
{
  double c;

  c = h/n0;
  if ( c < a ) c = a; else if ( c > b ) c = b;
  n0 = (int)(h/c+0.5);
  return max ( n0, 1 );
} /*_mbs_densd*/

void mbs_DrawTrimBSPatchDomd ( int degu, int lastuknot, const double *uknots,
                               int degv, int lastvknot, const double *vknots,
                               int nelem, const polycurved *bound,
                               int nu, double au, double bu,
                               int nv, double av, double bv,
                               int maxinters,
                               void (*NotifyLine)(char,int,point2d*,point2d*),
                               void (*DrawLine)(point2d*,point2d*,int),
                               void (*DrawCurve)(int,int,const double*) )
{
  void    *buf, *sp;
  char    *bufp;
  int     i, j, k, l, s, s0, N, ninters;
  char    dim;
  point2d p, q, c[2];
  double   *cp, h, u, v;
  signpoint1d *inters;

  sp = pkv_GetScratchMemTop ();
  buf = mbs_CompileTrimPatchBoundd ( nelem, bound, NULL );
  inters = (signpoint1d*)pkv_GetScratchMem ( maxinters*sizeof(signpoint1d) );
  if ( !buf || !inters ) {
    PKV_SIGNALERROR ( LIB_MULTIBS, ERRCODE_2, ERRMSG_2 );
    exit ( 0 );
  }

  if ( DrawLine != NULL ) {
    for ( i = degu+1; i <= lastuknot-degu; i++ ) {
      if ( (h = uknots[i] - uknots[i-1]) > 0.0 ) {
        l = _mbs_densd ( nu, au, bu, h );
        h /= l;
        if ( i == lastuknot-degu ) s0 = 1; else s0 = 0;
        for ( s = l-1; s >= s0; s-- ) {
          u = uknots[i] - s*h;
          SetPoint2d ( &p, u, vknots[degv] );
          SetPoint2d ( &q, u, vknots[lastvknot-degv] );
          ninters = maxinters;
          mbs_FindBoundLineIntersectionsd ( buf, &p, p.y, &q, q.y,
                                            inters, &ninters );
          if ( ninters ) {
            if ( NotifyLine != NULL )
              NotifyLine ( 1, i, &p, &q );
            InterPoint2d ( &p, &q, (inters[0].t-p.y)/(q.y-p.y), &c[0] );
            for ( k = j = 0; k < ninters-1; k++ ) {
              j += inters[k].sign1;
              InterPoint2d ( &p, &q, (inters[k+1].t-p.y)/(q.y-p.y), &c[1] );
              DrawLine ( &c[0], &c[1], j );
              c[0] = c[1];
            }
          }
        }
      }
    }
    for ( i = degv+1; i <= lastvknot-degv; i++ ) {
      if ( (h = vknots[i] - vknots[i-1]) > 0.0 ) {
        l = _mbs_densd ( nv, av, bv, h );
        h /= l;
        if ( i == lastvknot-degv ) s0 = 1; else s0 = 0;
        for ( s = l-1; s >= s0; s-- ) {
          v = vknots[i] - s*h;
          SetPoint2d ( &p, uknots[degu], v );
          SetPoint2d ( &q, uknots[lastvknot-degv], v );
          ninters = maxinters;
          mbs_FindBoundLineIntersectionsd ( buf, &p, p.x, &q, q.x,
                                            inters, &ninters );
          if ( ninters ) {
            if ( NotifyLine != NULL )
              NotifyLine ( 2, i, &p, &q );
            InterPoint2d ( &p, &q, (inters[0].t-p.x)/(q.x-p.x), &c[0] );
            for ( k = j = 0; k < ninters-1; k++ ) {
              j += inters[k].sign1;
              InterPoint2d ( &p, &q, (inters[k+1].t-p.x)/(q.x-p.x), &c[1] );
              DrawLine ( &c[0], &c[1], j );
              c[0] = c[1];
            }
          }
        }
      }
    }
  }

  if ( DrawCurve != NULL ) {
    bufp = (char*)buf;
    while ( bufp[0] != 4 ) {
      dim = bufp[1];
      N   = *((short*)(&bufp[2]));
      cp  = (double*)(&bufp[2+sizeof(short)]);
      switch ( bufp[0] ) {
      case 0:
      case 1:
        for ( i = 0; i < N; i++ )
          DrawCurve ( dim, 1, &cp[dim*i] );
        break;

      case 2:
      case 3:
        DrawCurve ( dim, N, cp );
        break;
 
      default:
        goto out;
      }
      bufp = &bufp[2+sizeof(short)+(N+1)*dim*sizeof(double)];
    }
  }

out:
  pkv_SetScratchMemTop ( sp );
} /*mbs_DrawTrimBSPatchDomd*/

