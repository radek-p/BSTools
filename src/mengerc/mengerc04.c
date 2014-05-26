
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2014                                  */
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
#include "mengerc.h"

#include "mengercprivate.h"

/* ///////////////////////////////////////////////////////////////////////// */
boolean mengerc_intD ( mengerc_data *md,
                       int lkn, double *knots, point3d *cpoints,
                       double *dl, double *acp )
{
  double   *qk, *qc;  /* wezly i wspolczynniki kwadratury */
  double   L, M, A, C, C2, f;
  int      n, nqkn, i, j;
  point3d  p;
  vector3d d;

  nqkn = md->nqkn;
  qk = md->qkn;
  qc = md->qc;
  n = md->n;
  C = md->L/(lkn-2*n);
  C2 = C*C;
        /* calkujemy i sumujemy dlugosci lukow krzywej */
        /* to jest o.k. jesli wezly sa rownoodlegle i odleglosc */
        /* wezlow sasiednich jest 1.0 */
  L = A = 0.0;
  for ( j = 0; j < nqkn; j++ ) {
    for ( i = 0; i < lkn-2*n; i++ ) {
      mbs_deBoorDerC3d ( n, 2*n+1, knots, &cpoints[i], (double)n+qk[j], &p, &d );
      M = DotProduct3d ( &d, &d );
      f = sqrt ( M );
      L += f*qc[j];
      A += (0.5*(f+C2/f)-C)*qc[j];
    }
  }
  *dl = L;
  *acp = A;
  return true;
} /*mengerc_intD*/

boolean mengerc_gradIntD ( mengerc_data *md,
                           int lkn, double *knots, point3d *cpoints,
                           double *dl, double *grdl, double *acp, double *gracp )
{
  void     *sp;
  double   *qk, *qc;  /* wezly i wspolczynniki kwadratury */
  double   L, M, f, C, C2, *dM;
  double   A, Nl, Nm, a, b;
  int      n, ncp, nvcp, N, nqkn;
  int      i, j, l, m, ii, jj;
  point3d  p;
  vector3d d, *dcp;

  sp = pkv_GetScratchMemTop ();
  nqkn = md->nqkn;
  qk = md->qkn;
  qc = md->qc;
  n = md->n;
  ncp = lkn - n;  /* liczba wszystkich punktow kontrolnych */
  nvcp = ncp - n; /* liczba niezaleznych p. k., a takze lukow wielomianowych */
  N = 3*nvcp;     /* liczba zmiennych */
  dcp = pkv_GetScratchMem ( (ncp-1)*sizeof(vector3d) );
  C = md->L/(lkn-2*n);
  C2 = C*C;
  dM  = pkv_GetScratchMemd ( 3*(n+1) );
  if ( !dcp || !dM )
    goto failure;

        /* obliczamy roznice punktow kontrolnych */
  for ( i = 0; i < ncp-1; i++ )
    SubtractPoints3d ( &cpoints[i+1], &cpoints[i], &dcp[i] );

        /* calkujemy i sumujemy dlugosci lukow krzywej */
        /* to jest o.k. jesli wezly sa rownoodlegle i odleglosc */
        /* wezlow sasiednich jest 1.0 */
  L = A = 0.0;
  memset ( grdl, 0, N*sizeof(double) );
  memset ( gracp, 0, N*sizeof(double) );
  for ( j = 0; j < nqkn; j++ ) {
    for ( i = 0; i < lkn-2*n; i++ ) {
      mbs_deBoorDerC3d ( n, 2*n+1, knots, &cpoints[i], (double)n+qk[j], &p, &d );
      M = DotProduct3d ( &d, &d );
      f = sqrt ( M );
      L += f*qc[j];
      A += (0.5*(f+C2/f)-C)*qc[j];
          /* te dwie petle przebiegaja po punktach kontrolnych, */
          /* od ktorych zalezy i-ty luk krzywej sklejanej */
          /* obliczamy pochodne czastkowe M wzgledem wspolrzednych p.k. */
      memset ( dM, 0, 3*(n+1)*sizeof(double) );
      for ( l = 0; l < n ; l++ ) {
        Nl = md->bsf1[n*j+l];
        for ( m = 0; m < n; m++ ) {
          Nm = md->bsf1[n*j+m];
          a = 2.0*Nl*Nm;
          ii = 3*l;
          jj = i+m;
          dM[ii]   -= a*dcp[jj].x;
          dM[ii+1] -= a*dcp[jj].y;
          dM[ii+2] -= a*dcp[jj].z;
          ii += 3;
          dM[ii]   += a*dcp[jj].x;
          dM[ii+1] += a*dcp[jj].y;
          dM[ii+2] += a*dcp[jj].z;
        }
      }
          /* w tej petli jest numeryczne calkowanie */
      a = qc[j]/(f+f);
      b = 0.5*a*(1.0-C2/M);
      for ( l = 0; l <= n; l++ ) {
        ii = i+l < nvcp ? 3*(i+l) : 3*(i+l-nvcp);
        jj = 3*l;
        grdl[ii]    += a*dM[jj];
        grdl[ii+1]  += a*dM[jj+1];
        grdl[ii+2]  += a*dM[jj+2];
        gracp[ii]   += b*dM[jj];
        gracp[ii+1] += b*dM[jj+1];
        gracp[ii+2] += b*dM[jj+2];
      }
    }
  }
  *dl = L;
  *acp = A;
  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*mengerc_gradIntD*/

boolean mengerc_hessIntD ( mengerc_data *md,
                           int lkn, double *knots, point3d *cpoints,
                           double *dl, double *grdl, double *hesdl,
                           double *acp, double *gracp, double *hesacp )
{
  void *sp;
  double   *qk, *qc;  /* wezly i wspolczynniki kwadratury */
  double   L, M, f, C, C2, *dM, *ddM;
  double   A;
  double   Nl, Nm, a, b, c, e;
  int      n, ncp, nvcp, N, nqkn;
  int      i, j, l, m, ii, jj, ik, jk;
  point3d  p;
  vector3d d, *dcp;

  sp = pkv_GetScratchMemTop ();
  nqkn = md->nqkn;
  qk = md->qkn;
  qc = md->qc;
  n = md->n;
  ncp = lkn - n;  /* liczba wszystkich punktow kontrolnych */
  nvcp = ncp - n; /* liczba niezaleznych p. k., a takze lukow wielomianowych */
  N = 3*nvcp;     /* liczba zmiennych */
  dcp = pkv_GetScratchMem ( (ncp-1)*sizeof(vector3d) );
  C = md->L/(lkn-2*n);
  C2 = C*C;
  dM = pkv_GetScratchMemd ( 3*(n+1) );
  ddM = pkv_GetScratchMemd ( ((n+1)*(n+2))/2 );
  if ( !dcp || !dM || !ddM )
    goto failure;
        /* obliczamy roznice punktow kontrolnych */
  for ( i = 0; i < ncp-1; i++ )
    SubtractPoints3d ( &cpoints[i+1], &cpoints[i], &dcp[i] );

        /* calkujemy i sumujemy dlugosci lukow krzywej */
        /* to jest o.k. jesli wezly sa rownoodlegle i odleglosc */
        /* wezlow sasiednich jest 1.0 */
  L = A = 0.0;
  memset ( grdl, 0, N*sizeof(double) );
  memset ( hesdl, 0, (N*(N+1))/2*sizeof(double) );
  memset ( gracp, 0, N*sizeof(double) );
  memset ( hesacp, 0, (N*(N+1))/2*sizeof(double) );
  for ( j = 0; j < nqkn; j++ ) {
    for ( i = 0; i < lkn-2*n; i++ ) {
      mbs_deBoorDerC3d ( n, 2*n+1, knots, &cpoints[i], (double)n+qk[j], &p, &d );
      M = DotProduct3d ( &d, &d );
      f = sqrt ( M );
      L += f*qc[j];
      A += (0.5*(f+C2/f)-C)*qc[j];
          /* te dwie petle przebiegaja po punktach kontrolnych, */
          /* od ktorych zalezy i-ty luk krzywej sklejanej */
          /* obliczamy pochodne czastkowe M wzgledem wspolrzednych p.k. */
      memset ( dM, 0, 3*(n+1)*sizeof(double) );
      memset ( ddM, 0, ((n+1)*(n+2))/2*sizeof(double) );
      for ( l = 0; l < n ; l++ ) {
        Nl = md->bsf1[n*j+l];
        for ( m = 0; m < n; m++ ) {
          Nm = md->bsf1[n*j+m];
          a = 2.0*Nl*Nm;
               /* pochodne pierwszego rzedu wyrazenia M */
          ii = 3*l;
          jj = i+m;
          dM[ii]   -= a*dcp[jj].x;
          dM[ii+1] -= a*dcp[jj].y;
          dM[ii+2] -= a*dcp[jj].z;
          ii += 3;
          dM[ii]   += a*dcp[jj].x;
          dM[ii+1] += a*dcp[jj].y;
          dM[ii+2] += a*dcp[jj].z;
               /* pochodne drugiego rzedu wyrazenia M */
          if ( l >= m ) {
            ddM[pkn_LowerTrMatIndex(l,m)]     += a;
            ddM[pkn_LowerTrMatIndex(l+1,m+1)] += a;
          }
          if ( l+1 >= m )
            ddM[pkn_LowerTrMatIndex(l+1,m)]   -= a;
          if ( l >= m+1 )
            ddM[pkn_LowerTrMatIndex(l,m+1)]   -= a;
        }
      }
          /* w tej petli jest numeryczne calkowanie */
      a = qc[j]/(f+f);
      b = 0.5*a*(1.0-C2/M);
      for ( l = 0; l <= n; l++ ) {
        ii = i+l < nvcp ? 3*(i+l) : 3*(i+l-nvcp);
        jj = 3*l;
        grdl[ii]    += a*dM[jj];
        grdl[ii+1]  += a*dM[jj+1];
        grdl[ii+2]  += a*dM[jj+2];
        gracp[ii]   += b*dM[jj];
        gracp[ii+1] += b*dM[jj+1];
        gracp[ii+2] += b*dM[jj+2];
      }
          /* w tej petli jest numeryczne calkowanie wspolrzednych hesjanow */
      b = qc[j]*(1.0-C2/M)/(4.0*f);
      c = qc[j]*(3.0*C2/M-1.0)/(8.0*M*f);
      for ( l = 0; l <= n; l++ ) {
        ii = i+l < nvcp ? 3*(i+l) : 3*(i+l-nvcp);
        ik = 3*l;
        for ( m = 0; m < l; m++ ) {
          jj = i+m < nvcp ? 3*(i+m) : 3*(i+m-nvcp);
          jk = 3*m;
          e = ddM[pkn_SymMatIndex(l,m)];
          hesdl[pkn_SymMatIndex(ii,jj)]     += a*(e-dM[ik]*dM[jk]/(M+M));
          hesdl[pkn_SymMatIndex(ii,jj+1)]   -= a*(dM[ik]*dM[jk+1]/(M+M));
          hesdl[pkn_SymMatIndex(ii,jj+2)]   -= a*(dM[ik]*dM[jk+2]/(M+M));
          hesdl[pkn_SymMatIndex(ii+1,jj)]   -= a*(dM[ik+1]*dM[jk]/(M+M));
          hesdl[pkn_SymMatIndex(ii+1,jj+1)] += a*(e-dM[ik+1]*dM[jk+1]/(M+M));
          hesdl[pkn_SymMatIndex(ii+1,jj+2)] -= a*(dM[ik+1]*dM[jk+2]/(M+M));
          hesdl[pkn_SymMatIndex(ii+2,jj)]   -= a*(dM[ik+2]*dM[jk]/(M+M));
          hesdl[pkn_SymMatIndex(ii+2,jj+1)] -= a*(dM[ik+2]*dM[jk+1]/(M+M));
          hesdl[pkn_SymMatIndex(ii+2,jj+2)] += a*(e-dM[ik+2]*dM[jk+2]/(M+M));

          hesacp[pkn_SymMatIndex(ii,jj)]     += b*e + c*dM[ik]*dM[jk];
          hesacp[pkn_SymMatIndex(ii,jj+1)]   += c*dM[ik]*dM[jk+1];
          hesacp[pkn_SymMatIndex(ii,jj+2)]   += c*dM[ik]*dM[jk+2];
          hesacp[pkn_SymMatIndex(ii+1,jj)]   += c*dM[ik+1]*dM[jk];
          hesacp[pkn_SymMatIndex(ii+1,jj+1)] += b*e + c*dM[ik+1]*dM[jk+1];
          hesacp[pkn_SymMatIndex(ii+1,jj+2)] += c*dM[ik+1]*dM[jk+2];
          hesacp[pkn_SymMatIndex(ii+2,jj)]   += c*dM[ik+2]*dM[jk];
          hesacp[pkn_SymMatIndex(ii+2,jj+1)] += c*dM[ik+2]*dM[jk+1];
          hesacp[pkn_SymMatIndex(ii+2,jj+2)] += b*e + c*dM[ik+2]*dM[jk+2];
        }
            /* obliczenie dla m == l */
        e = ddM[pkn_LowerTrMatIndex(l,l)];
        hesdl[pkn_LowerTrMatIndex(ii,ii)]     += a*(e-dM[ik]*dM[ik]/(M+M));
        hesdl[pkn_LowerTrMatIndex(ii+1,ii)]   -= a*(dM[ik+1]*dM[ik]/(M+M));
        hesdl[pkn_LowerTrMatIndex(ii+1,ii+1)] += a*(e-dM[ik+1]*dM[ik+1]/(M+M));
        hesdl[pkn_LowerTrMatIndex(ii+2,ii)]   -= a*(dM[ik+2]*dM[ik]/(M+M));
        hesdl[pkn_LowerTrMatIndex(ii+2,ii+1)] -= a*(dM[ik+2]*dM[ik+1]/(M+M));
        hesdl[pkn_LowerTrMatIndex(ii+2,ii+2)] += a*(e-dM[ik+2]*dM[ik+2]/(M+M));

        hesacp[pkn_LowerTrMatIndex(ii,ii)]     += b*e + c*dM[ik]*dM[ik];
        hesacp[pkn_LowerTrMatIndex(ii+1,ii)]   += c*dM[ik+1]*dM[ik];
        hesacp[pkn_LowerTrMatIndex(ii+1,ii+1)] += b*e + c*dM[ik+1]*dM[ik+1];
        hesacp[pkn_LowerTrMatIndex(ii+2,ii)]   += c*dM[ik+2]*dM[ik];
        hesacp[pkn_LowerTrMatIndex(ii+2,ii+1)] += c*dM[ik+2]*dM[ik+1];
        hesacp[pkn_LowerTrMatIndex(ii+2,ii+2)] += b*e + c*dM[ik+2]*dM[ik+2];
      }
    }
  }
  *dl = L;
  *acp = A;
  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*mengerc_hessIntD*/

