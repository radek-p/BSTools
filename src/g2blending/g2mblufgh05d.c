
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2010                                  */
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
double g2mbl_AFuncd ( int nkn, double *qcoeff, double **Nitabs, double **Jac,
                      int nv, point3d *mvcp, int bdomelems, int *domelind,
                      meshdom_elem *domelem, int *domelcpind,
                      boolean force, double *ftab )
{
  double f;
  int    i, j, k, ind;

        /* compute the integrals in dirty elements */
  if ( force )
    for ( i = 0; i < bdomelems; i++ ) {
      j = domelind[i];
      k = domelem[j].type;
      ind = domelem[j].firstcpi;
      if ( k == 4 )  /* integrate over one square only */
        g2mbl_UFuncRSQd ( nkn, qcoeff, Nitabs[1], &domelcpind[ind], mvcp,
                          domelem[j].C, &ftab[j] );
      else           /* integrate over k squares */
        g2mbl_UFuncSSQd ( nkn, qcoeff, k, Nitabs[k-3], Jac[k-3],
                          &domelcpind[ind], mvcp,
                          domelem[j].C, &ftab[j] );
    }
        /* sum the integrals over the elements */
  f = 0.0;
  for ( i = 0; i < bdomelems; i++ )
    f += ftab[domelind[i]];
  return f;
} /*g2mbl_AFuncd*/

/* ///////////////////////////////////////////////////////////////////////// */
boolean g2mbl_AFuncGradd ( int nkn, double *qcoeff, double **Nitabs, double **Jac,
                           int nv, point3d *mvcp, int bdomelems, int *domelind,
                           meshdom_elem *domelem, int *domelcpind, int *bncpi,
                           boolean force, double *ftab, double *gtab,
                           double *func, int nvars, double *grad )
{
  double f;
  int i, j, k, l, ind, fp, ncp;

        /* compute the integrals in dirty elements */
  if ( force )
    for ( i = 0; i < bdomelems; i++ ) {
      l = domelind[i];
      k = domelem[l].type;
      ind = domelem[l].firstcpi;
      if ( k == 4 )
        g2mbl_UFuncGradRSQd ( nkn, qcoeff, Nitabs[1], &domelcpind[ind],
                              bncpi, mvcp, domelem[l].C,
                              &ftab[l], &gtab[3*ind] );
      else {
        if ( !g2mbl_UFuncGradSSQd ( nkn, qcoeff, k, Nitabs[k-3], Jac[k-3],
                          &domelcpind[ind],
                          bncpi, mvcp, domelem[l].C,
                          &ftab[l], &gtab[3*ind] ) )
          return false;
      }
    }
        /* sum the integrals over the elements */
  f = 0.0;
  memset ( grad, 0, nvars*sizeof(double) );
  for ( i = 0; i < bdomelems; i++ ) {
    l = domelind[i];
    f += ftab[l];
    fp = domelem[l].firstcpi;
    ncp = domelem[l].ncp;
    for ( j = 0; j < ncp; j++ ) {
      k = bncpi[domelcpind[fp+j]];
      if ( k >= 0 )
        AddVector3d ( (point3d*)&grad[3*k], (vector3d*)&gtab[3*(fp+j)],
                      (point3d*)&grad[3*k] );
    }
  }
  *func = f;
  return true;
} /*g2mbl_AFuncGradd*/

/* ///////////////////////////////////////////////////////////////////////// */
boolean g2mbl_AFuncGradHessiand ( int nkn, double *qcoeff,
                        double **Nitabs, double **Nijtabs, double **Mijtabs,
                        double **Jac, int nv, point3d *mvcp,
                        int bdomelems, int *domelind, meshdom_elem *domelem,
                        int *domelcpind, int *nncpi, int *bncpi,
                        boolean force,
                        double *ftab, double *gtab, double *htab,
                        double *func, int nvars, double *grad,
                        int hsize, int *hprof, double **hrows )
{
  double f;
  int i, j, k, l, m, kk, mm, p, ind, fp, ncp, hti, s, t, b;

        /* compute the integrals in dirty elements */
  if ( force )
    for ( i = 0; i < bdomelems; i++ ) {
      p = domelind[i];

#ifdef __DEBUG
printf ( "%5d\b\b\b\b\b", i );
#endif

      k = domelem[p].type;
      ind = domelem[p].firstcpi;
      hti = domelem[p].hti;
      if ( k == 4 )
        g2mbl_UFuncGradHessRSQd ( nkn, qcoeff,
                                  Nitabs[1], Nijtabs[1], Mijtabs[1],
                                  &domelcpind[ind], nncpi, mvcp,
                                  domelem[p].C,
                                  &ftab[p], &gtab[3*ind], &htab[hti] );
      else {
        if ( !g2mbl_UFuncGradHessSSQd ( nkn, qcoeff, k,
                          Nitabs[k-3], Nijtabs[k-3], Mijtabs[k-3], Jac[k-3],
                          &domelcpind[ind], nncpi, mvcp,
                          domelem[p].C,
                          &ftab[p], &gtab[3*ind], &htab[hti] ) )
          return false;
      }
    }
        /* sum the integrals over the elements */
  f = 0.0;
  memset ( grad, 0, nvars*sizeof(double) );
  memset ( hrows[0], 0, hsize*sizeof(double) );
  for ( i = 0; i < bdomelems; i++ ) {
    p = domelind[i];
    f += ftab[p];
    fp = domelem[p].firstcpi;
    ncp = domelem[p].ncp;
    hti = domelem[p].hti;
    for ( j = 0; j < ncp; j++ ) {
      k = bncpi[domelcpind[fp+j]];
      if ( k >= 0 ) {
        kk = nncpi[domelcpind[fp+j]];
        AddVector3d ( (point3d*)&grad[3*k], (vector3d*)&gtab[3*(fp+j)],
                      (point3d*)&grad[3*k] );
        for ( l = 0; l <= j; l++ ) {
          m = bncpi[domelcpind[fp+l]];
          if ( m >= 0 ) {
            mm = nncpi[domelcpind[fp+l]];
            b = hti + 9*pkn_SymMatIndex ( j, l );
            if ( k > m ) {
              if ( kk > mm ) {
                for ( s = 0; s < 3; s++ )
                  for ( t = 0; t < 3; t++ )
                    hrows[3*k+s][3*m+t] += htab[b+3*s+t];
              }
              else {
                for ( s = 0; s < 3; s++ )
                  for ( t = 0; t < 3; t++ )
                    hrows[3*k+s][3*m+t] += htab[b+3*t+s];
              }
            }
            else if ( k < m ) {
              if ( kk > mm ) {
                for ( s = 0; s < 3; s++ )
                  for ( t = 0; t < 3; t++ )
                    hrows[3*m+s][3*k+t] += htab[b+3*t+s];
              }
              else {
                for ( s = 0; s < 3; s++ )
                  for ( t = 0; t < 3; t++ )
                    hrows[3*m+s][3*k+t] += htab[b+3*s+t];
              }
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
} /*g2mbl_AFuncGradHessiand*/

