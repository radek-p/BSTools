
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2010, 2011                            */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

#include "pkvaria.h"
#include "pknum.h"
#include "pkgeom.h"
#include "multibs.h"
#include "bsmesh.h"
#include "egholed.h"
#include "g2blendingd.h"
#include "g2mblendingd.h"

#include "g2blprivated.h"
#include "g2mblprivated.h"

/* ///////////////////////////////////////////////////////////////////////// */
double g2mbl_UFuncd ( int nkn, double *qcoeff, double **Nitabs, double **Jac,
                      int nv, point3d *mvcp,
                      int ndomelems, meshdom_elem *domelem, int *domelcpind,
                      boolean force, double *ftab )
{
  double f;
  int    i, k, ind;

        /* compute the integrals */
  if ( force )
    for ( i = 0; i < ndomelems; i++ ) {
      k = domelem[i].type;
      ind = domelem[i].firstcpi;
      if ( k == 4 )  /* integrate over one square only */
        g2mbl_UFuncRSQd ( nkn, qcoeff, Nitabs[1], &domelcpind[ind], mvcp,
                          domelem[i].C, &ftab[i] );
      else           /* integrate over k squares */
        g2mbl_UFuncSSQd ( nkn, qcoeff, k, Nitabs[k-3], Jac[k-3],
                          &domelcpind[ind], mvcp,
                          domelem[i].C, &ftab[i] );
    }
        /* sum the integrals over the elements */
  f = 0.0;
  for ( i = 0; i < ndomelems; i++ )
    f += ftab[i];
  return f;
} /*g2mbl_UFuncd*/

/* ///////////////////////////////////////////////////////////////////////// */
boolean g2mbl_UFuncGradd ( int nkn, double *qcoeff, double **Nitabs, double **Jac,
                           int nv, point3d *mvcp,
                           int ndomelems, meshdom_elem *domelem,
                           int *domelcpind, int *nncpi,
                           boolean force, double *ftab, double *gtab,
                           double *func, int nvars, double *grad )
{
  double f;
  int i, j, k, ind, fp, ncp;

        /* compute the integrals */
  if ( force )
    for ( i = 0; i < ndomelems; i++ ) {
      k = domelem[i].type;
      ind = domelem[i].firstcpi;
      if ( k == 4 )
        g2mbl_UFuncGradRSQd ( nkn, qcoeff, Nitabs[1], &domelcpind[ind],
                              nncpi, mvcp, domelem[i].C,
                              &ftab[i], &gtab[3*ind] );
      else {
        if ( !g2mbl_UFuncGradSSQd ( nkn, qcoeff, k, Nitabs[k-3], Jac[k-3],
                          &domelcpind[ind],
                          nncpi, mvcp, domelem[i].C,
                          &ftab[i], &gtab[3*ind] ) )
          return false;
      }
    }
        /* sum the integrals over the elements */
  f = 0.0;
  memset ( grad, 0, nvars*sizeof(double) );
  for ( i = 0; i < ndomelems; i++ ) {
    f += ftab[i];
    fp = domelem[i].firstcpi;
    ncp = domelem[i].ncp;
    for ( j = 0; j < ncp; j++ ) {
      k = nncpi[domelcpind[fp+j]];
      if ( k >= 0 )
        AddVector3d ( (point3d*)&grad[3*k], (vector3d*)&gtab[3*(fp+j)],
                      (point3d*)&grad[3*k] );
    }
  }
  *func = f;
  return true;
} /*g2mbl_UFuncGradd*/

/* ///////////////////////////////////////////////////////////////////////// */
boolean g2mbl_UFuncGradHessiand ( int nkn, double *qcoeff,
                        double **Nitabs, double **Nijtabs, double **Mijtabs,
                        double **Jac,
                        int nv, point3d *mvcp,
                        int ndomelems, meshdom_elem *domelem,
                        int *domelcpind, int *nncpi,
                        boolean force,
                        double *ftab, double *gtab, double *htab,
                        double *func, int nvars, double *grad,
                        int hsize, int *hprof, double **hrows )
{
  double f;
  int i, j, k, l, m, ind, fp, ncp, hti, s, t, b;

        /* compute the integrals */
  if ( force )
    for ( i = 0; i < ndomelems; i++ ) {

#ifdef __DEBUG
printf ( "%5d\b\b\b\b\b", i );
#endif

      k = domelem[i].type;
      ind = domelem[i].firstcpi;
      hti = domelem[i].hti;
      if ( k == 4 )
        g2mbl_UFuncGradHessRSQd ( nkn, qcoeff,
                                  Nitabs[1], Nijtabs[1], Mijtabs[1],
                                  &domelcpind[ind], nncpi, mvcp,
                                  domelem[i].C,
                                  &ftab[i], &gtab[3*ind], &htab[hti] );
      else {
        if ( !g2mbl_UFuncGradHessSSQd ( nkn, qcoeff, k,
                          Nitabs[k-3], Nijtabs[k-3], Mijtabs[k-3], Jac[k-3],
                          &domelcpind[ind], nncpi, mvcp,
                          domelem[i].C,
                          &ftab[i], &gtab[3*ind], &htab[hti] ) )
          return false;
      }
    }
        /* sum the integrals over the elements */
  f = 0.0;
  memset ( grad, 0, nvars*sizeof(double) );
  memset ( hrows[0], 0, hsize*sizeof(double) );
  for ( i = 0; i < ndomelems; i++ ) {
    f += ftab[i];
    fp = domelem[i].firstcpi;
    ncp = domelem[i].ncp;
    hti = domelem[i].hti;
    for ( j = 0; j < ncp; j++ ) {
      k = nncpi[domelcpind[fp+j]];
      if ( k >= 0 ) {
        AddVector3d ( (point3d*)&grad[3*k], (vector3d*)&gtab[3*(fp+j)],
                      (point3d*)&grad[3*k] );
        for ( l = 0; l <= j; l++ ) {
          m = nncpi[domelcpind[fp+l]];
          if ( m >= 0 ) {
            b = hti + 9*pkn_SymMatIndex ( j, l );
            if ( k > m ) {
              for ( s = 0; s < 3; s++ )
                for ( t = 0; t < 3; t++ )
                  hrows[3*k+s][3*m+t] += htab[b+3*s+t];
            }
            else if ( k < m ) {
              for ( s = 0; s < 3; s++ )
                for ( t = 0; t < 3; t++ )
                  hrows[3*m+s][3*k+t] += htab[b+3*s+t];
            }
            else {
              for ( s = 0; s < 3; s++ )
                for ( t = 0; t <= s; t++ )
                  hrows[3*k+s][3*k+t] += htab[b+3*s+t];
            }
          }
        }
      }
    }
  }
  *func = f;
  return true;
} /*g2mbl_UFuncGradHessiand*/

