
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2010, 2012                            */
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
boolean g2mbl_B3x3FuncGradHessiand ( int nkn, double *qcoeff,
                        double **Nitabs, double **Nijtabs, double **Mijtabs,
                        double **Jac,
                        int nv, point3d *mvcp,
                        int ndomelems, meshdom_elem *domelem,
                        int *domelcpind, int *nncpi,
                        boolean force,
                        double *ftab, double *gtab, double *htab,
                        double *func, int nvars, double *grad,
                        int Hblsize, nzHbl_rowdesc *iHbl, int *cHbl, double *Hbl )
{
  double   f;
  int      i, j, k, l, m, ind, fp, ncp, hti, s, t, b, hi;

        /* compute the integrals in dirty elements */
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
  memset ( Hbl, 0, Hblsize*9*sizeof(double) );
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
            hi = hti + 9*pkn_SymMatIndex ( j, l );
            if ( k > m ) {       /* (k,m) is below diagonal */
                  /* binary search */
              s = iHbl[k].firsthbl;
              t = s + iHbl[k].nhbl-1;
              do {
                b = (s+t) / 2;
                if ( cHbl[b] > m )
                  t = b;
                else
                  s = b;
              } while ( t-s > 1 );
              b = 9*s;
              pkn_AddMatrixd ( 1, 9, 0, &Hbl[b], 0, &htab[hi], 0, &Hbl[b] );
            }
            else if ( k < m ) {  /* (m,k) is below diagonal  */
                  /* binary search */
              s = iHbl[m].firsthbl;
              t = s + iHbl[m].nhbl-1;
              do {
                b = (s+t) / 2;
                if ( cHbl[b] > k )
                  t = b;
                else
                  s = b;
              } while ( t-s > 1 );
              b = 9*s;
              pkn_AddMatrixd ( 1, 9, 0, &Hbl[b], 0, &htab[hi], 0, &Hbl[b] );
            }
            else {               /* diagonal block */
              b = 9*(iHbl[k].firsthbl+iHbl[k].nhbl-1);
              pkn_AddMatrixd ( 1, 9, 0, &Hbl[b], 0, &htab[hi], 0, &Hbl[b] );
            }
          }
        }
      }
    }
  }

  *func = f;
  return true;
} /*g2mbl_B3x3FuncGradHessiand*/

