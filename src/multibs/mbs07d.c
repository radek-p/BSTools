
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2005, 2009                            */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

#include <math.h>
#include <string.h>
#include <stdio.h>  /* for testing */

#include "pkvaria.h"
#include "pkgeom.h"
#include "multibs.h"


/* /////////////////////////////////////////// */
/* knot removal                                */

int mbs_multiKnotRemoved ( int degree, int *lastknot,
                           double *knots,
                           int ncurves, int spdimen, int inpitch, int outpitch,
                           double *ctlpoints, int knotnum )
{
  int i, j, k, l, r;
  int scratchsize, nkn, mp, neq, neqs;
  double  t0;
  long double e, gc, gs;
  long double *a, *b;
  double  *c, *d;

  /* find the last instance of the knot to remove */
  nkn = (*lastknot)-1;
  k = knotnum;
  t0 = knots[k];
  while ( k <= nkn && knots[k+1] == t0 )
    k++;
  k--;
  knotnum = k;

  memmove ( &knots[k+1], &knots[k+2], (nkn-k)*sizeof(double) );
  i = knotnum;
  r = 0;
  while ( i > 0 && knots[i] == t0 ) {
    r++;
    i--;
  }

  if ( r < degree ) {
                      /* solve the system of equations with Givens rotations */
                        /* allocate scratch memory */

    neq  = degree-r+2;     /* number of equations */
    neqs = neq*spdimen;
    scratchsize = neq*(2*sizeof(long double)+ncurves*spdimen*sizeof(double));
    a = (long double*)pkv_GetScratchMem ( scratchsize );
    b = &a[neq];
    c = (double*)&b[neq];

                        /* copy data */
    pkv_Selectd ( ncurves, neq*spdimen, inpitch, neq*spdimen,
                  &ctlpoints[(knotnum-degree)*spdimen], c );

                        /* compute matrix coefficients */
    b[0] = 1.0;
    a[0] = 0.0;
    for ( i = 1; i < neq-1; i++ ) {
      b[i] = (t0-knots[k-degree+i])/(knots[k+i]-knots[k-degree+i]);
      a[i] = 1.0-b[i];
    }
    b[neq-1] = 0.0;
    a[neq-1] = 1.0;
                        /* solve the equations */
    for ( i = 0; i < neq-1; i++ ) {
      e      = sqrt( (double)(b[i]*b[i]+a[i+1]*a[i+1]) );
      gc     = b[i]/e;        /* cosine */
      gs     = -a[i+1]/e;     /* sine */
      b[i]   =  gc*b[i] - gs*a[i+1];
      a[i]   = -gs*b[i+1];
      b[i+1] =  gc*b[i+1];
      for ( l = 0, d = &c[i*spdimen]; l < ncurves; l++, d += neqs )
        for ( j = 0; j < spdimen; j++ ) {
          e            = gc*d[j] - gs*d[j+spdimen];
          d[j+spdimen] = (double)(gs*d[j] + gc*d[j+spdimen]);
          d[j]         = (double)e;
        }
    }

    for ( l = 0, d = &c[(neq-2)*spdimen]; l < ncurves; l++, d += neqs )
      for ( j = 0; j < spdimen; j++ )
        d[j] = (double)(d[j]/b[neq-2]);
    for ( i = neq-3; i >= 0; i-- )
      for ( l = 0, d = &c[i*spdimen]; l < ncurves; l++, d+= neqs )
        for ( j = 0; j < spdimen; j++ )
          d[j] = (double)((d[j]-a[i]*d[j+spdimen])/b[i]);

                        /* solve the equations */
                        /* insert result and rearrange data */
    pkv_Selectd ( ncurves, (neq-1)*spdimen, neq*spdimen, inpitch,
                  c, &ctlpoints[(knotnum-degree)*spdimen] );
    mp = (nkn-degree-1-(knotnum-r))*spdimen;
    pkv_Moved ( ncurves, mp, inpitch, -spdimen,
                &ctlpoints[(knotnum-r+2)*spdimen] );

                        /* cleanup */
    pkv_FreeScratchMem ( scratchsize );
  }
  else {
    if ( r == degree ) {
                       /* here the least squares problem is solved */
                       /* by taking averages */
      for ( l = 0, d = &ctlpoints[(k-r)*spdimen];
            l < ncurves;
            l++, d += inpitch )
        for ( j = 0; j < spdimen; j++ )
          d[j] = 0.5*(d[j]+d[j+spdimen]);
    }

                        /* move the control points */
    mp = (nkn-degree-(k-r)-1)*spdimen;
    pkv_Moved ( ncurves, mp, inpitch, -spdimen, &ctlpoints[(k-r+2)*spdimen] );
  }
  *lastknot = nkn;
  if ( outpitch != inpitch )
    pkv_Rearranged ( ncurves, (*lastknot-degree)*spdimen,
                     inpitch, outpitch, ctlpoints );   
  return knotnum;
} /*mbs_multiKnotRemoved*/

