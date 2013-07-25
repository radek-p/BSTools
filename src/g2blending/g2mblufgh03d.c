
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
double g2mbl_BFuncd ( int nkn, double *qcoeff, double **Nitabs, double **Jac,
                      int nv, point3d *mvcp, int bdomelems, int *domelind,
                      meshdom_elem *domelem, int *domelcpind,
                      boolean force, double *ftab )
{
  double f;
  int    i, k, el, ind;

        /* compute the integrals */
  if ( force )
    for ( i = 0; i < bdomelems; i++ ) {
      el = domelind[i];
      k = domelem[el].type;
      ind  = domelem[el].firstcpi;
      if ( k == 4 )
        g2mbl_UFuncRSQd ( nkn, qcoeff, Nitabs[1], &domelcpind[ind], mvcp,
                          domelem[el].C, &ftab[el] );
      else
        g2mbl_UFuncSSQd ( nkn, qcoeff, k, Nitabs[k-3], Jac[k-3],
                          &domelcpind[ind], mvcp,
                          domelem[el].C, &ftab[el] );
    }
        /* sum the integrals over the elements */
  f = 0.0;
  for ( i = 0; i < bdomelems; i++ )
    f += ftab[domelind[i]];
  return f;
} /*g2mbl_BFuncd*/

/* ///////////////////////////////////////////////////////////////////////// */
boolean g2mbl_BFuncGradd ( int nkn, double *qcoeff, double **Nitabs, double **Jac,
                        int nv, point3d *mvcp, int bdomelems, int *domelind,
                        meshdom_elem *domelem, int *domelcpind, int *nncpi,
                        int cp0, int cp1, int *vpermut,
                        boolean force, double *ftab, double *gtab,
                        double *func, int nvars, double *grad )
{
  double f;
  int i, j, k, el, ind, fp, ncp;

        /* compute the integrals in dirty elements */
  if ( force )
    for ( i = 0; i < bdomelems; i++ ) {
      el = domelind[i];
      k = domelem[el].type;
      ind = domelem[el].firstcpi;
      if ( k == 4 )
        g2mbl_UFuncGradRSQd ( nkn, qcoeff, Nitabs[1], &domelcpind[ind],
                              nncpi, mvcp, domelem[el].C,
                              &ftab[el], &gtab[3*ind] );
      else {
        if ( !g2mbl_UFuncGradSSQd ( nkn, qcoeff, k, Nitabs[k-3], Jac[k-3],
                          &domelcpind[ind],
                          nncpi, mvcp, domelem[el].C,
                          &ftab[el], &gtab[3*ind] ) )
          return false;
      }
    }
        /* sum the integrals over the elements */
  f = 0.0;
  memset ( grad, 0, nvars*sizeof(double) );
  for ( i = 0; i < bdomelems; i++ ) {
    el = domelind[i];
    f += ftab[el];
    fp = domelem[el].firstcpi;
    ncp = domelem[el].ncp;
    for ( j = 0; j < ncp; j++ ) {
      k = nncpi[domelcpind[fp+j]];
      if ( k >= cp0 && k < cp1 ) {
        ind = 3*vpermut[k-cp0];
        AddVector3d ( (point3d*)&grad[ind], (vector3d*)&gtab[3*(fp+j)],
                      (point3d*)&grad[ind] );
      }
    }
  }
  *func = f;
  return true;
} /*g2mbl_BFuncGradd*/

/* ///////////////////////////////////////////////////////////////////////// */
boolean g2mbl_BFuncGradHessiand ( int nkn, double *qcoeff,
                        double **Nitabs, double **Nijtabs, double **Mijtabs,
                        double **Jac,
                        int nv, point3d *mvcp, int bdomelems, int *domelind,
                        meshdom_elem *domelem, int *domelcpind, int *nncpi,
                        int cp0, int cp1, int *vpermut,
                        boolean force,
                        double *ftab, double *gtab, double *htab,
                        double *func, int nvars, double *grad,
                        int hsize, int *hprof, double **hrows )
{
  double f;
  int i, j, k, l, m, el, ind, inm, fp, ncp, hti, s, t, b;

        /* compute the integrals in dirty elements */
  if ( force )
    for ( i = 0; i < bdomelems; i++ ) {
      el = domelind[i];

#ifdef __DEBUG
printf ( "%5d\b\b\b\b\b", i );
#endif

      k = domelem[el].type;
      ind = domelem[el].firstcpi;
      hti = domelem[el].hti;
      if ( k == 4 )
        g2mbl_UFuncGradHessRSQd ( nkn, qcoeff,
                                  Nitabs[1], Nijtabs[1], Mijtabs[1],
                                  &domelcpind[ind], nncpi, mvcp,
                                  domelem[el].C,
                                  &ftab[el], &gtab[3*ind], &htab[hti] );
      else {
        if ( !g2mbl_UFuncGradHessSSQd ( nkn, qcoeff, k,
                          Nitabs[k-3], Nijtabs[k-3], Mijtabs[k-3], Jac[k-3],
                          &domelcpind[ind], nncpi, mvcp,
                          domelem[el].C,
                          &ftab[el], &gtab[3*ind], &htab[hti] ) )
          return false;
      }
    }
        /* sum the integrals over the elements */
  f = 0.0;
  memset ( grad, 0, nvars*sizeof(double) );
  memset ( hrows[0], 0, hsize*sizeof(double) );
  for ( i = 0; i < bdomelems; i++ ) {
    el = domelind[i];
    f += ftab[el];
    fp = domelem[el].firstcpi;
    ncp = domelem[el].ncp;
    hti = domelem[el].hti;
    for ( j = 0; j < ncp; j++ ) {
      k = nncpi[domelcpind[fp+j]];
      if ( k >= cp0 && k < cp1 ) {
        ind = 3*vpermut[k-cp0];
        AddVector3d ( (point3d*)&grad[ind], (vector3d*)&gtab[3*(fp+j)],
                      (point3d*)&grad[ind] );
        for ( l = 0; l <= j; l++ ) {
          m = nncpi[domelcpind[fp+l]];
          if ( m >= cp0 && m < cp1 ) {
            inm = 3*vpermut[m-cp0];
            b = hti+9*pkn_SymMatIndex ( j, l );
            if ( ind > inm ) {
              if ( k > m ) {
                for ( s = 0; s < 3; s++ )
                  for ( t = 0; t < 3; t++ )
                    hrows[ind+s][inm+t] += htab[b+3*s+t];
              }
              else {
                for ( s = 0; s < 3; s++ )
                  for ( t = 0; t < 3; t++ )
                    hrows[ind+s][inm+t] += htab[b+3*t+s];
              }
            }
            else if ( ind < inm ) {
              if ( k > m ) {
                for ( s = 0; s < 3; s++ )
                  for ( t = 0; t < 3; t++ )
                    hrows[inm+s][ind+t] += htab[b+3*t+s];
              }
              else {
                for ( s = 0; s < 3; s++ )
                  for ( t = 0; t < 3; t++ )
                    hrows[inm+s][ind+t] += htab[b+3*s+t];
              }
            }
            else {
              for ( s = 0; s < 3; s++ )
                for ( t = 0; t <= s; t++ )
                  hrows[ind+s][inm+t] += htab[b+3*s+t];
            }
          }
        }
      }
    }
  }
  *func = f;
  return true;
} /*g2mbl_BFuncGradHessiand*/

