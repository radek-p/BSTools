
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2014                                  */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */ 

#include <pthread.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h> 
#include <math.h>

#include "pkvaria.h"
#include "pkvthreads.h"
#include "pknum.h"
#include "pkgeom.h"
#include "multibs.h"
#include "mengerc.h"

#include "mengercprivate.h"

typedef struct {
    mengerc_data *md;
    int          npoints;
    vector3d     *pt, *dpt, *ddpt, *dcp;
    double       *ldpt, *dldpt, *ddldpt;
    double       *_f, *_g, *_h;
  } intKM_job_desc;

/* ///////////////////////////////////////////////////////////////////////// */
static boolean _intKM_F ( void *usrdata, int3 *jobnum )
{
  intKM_job_desc *data;
  mengerc_data   *md;
  int            n, lkn, nqkn;
  vector3d       *pt, *dpt, *ddpt;
  double         *ldpt;
  int            q1, q2, q3, i, j, k;
  point3d        pt1, pt2, pt3;
  vector3d       v12, v23, v31, auxv1;
  double         fv12, fv23, fv31, fl12, fl23, fl31,
                 half, K2num, K2den, K2, KtoP, jac, KtoPj;
  double         intf3, intf4, intf5;

  data = (intKM_job_desc*)usrdata;
  md   = data->md;
  pt   = data->pt;
  dpt  = data->dpt;
  ddpt = data->ddpt;
  ldpt = data->ldpt;

  n    = md->deg;
  lkn  = md->lkn;
  nqkn = md->nqkn;

  q1 = jobnum->x;
  q2 = jobnum->y;
  q3 = jobnum->z;

  intf3 = 0.0;
          /* the three loops iterate over the cubes (i,j,k) dla i>j>k */
  for ( i = 2; i < lkn-2*n; i++ ) {
    pt1 = pt[i*nqkn+q1];
    for ( j = 1; j < i; j++ ) {
      pt2 = pt[j*nqkn+q2];
      SubtractPoints3d ( &pt1, &pt2, &v12 );
      for ( k = 0; k < j; k++ ) {
        pt3 = pt[k*nqkn+q3];
        SubtractPoints3d ( &pt2, &pt3, &v23 );
        SubtractPoints3d ( &pt3, &pt1, &v31 );
        fv12 = DotProduct3d ( &v12, &v12 );  fl12 = sqrt ( fv12 );
        fv23 = DotProduct3d ( &v23, &v23 );  fl23 = sqrt ( fv23 );
        fv31 = DotProduct3d ( &v31, &v31 );  fl31 = sqrt ( fv31 );

        half = 0.5*(fl12+fl23+fl31);
        K2num = 16.0*(half*(half-fl12))*((half-fl23)*(half-fl31));
        K2den = fv12*fv23*fv31;
          /* kwadrat krzywizny */
        K2 = K2num/K2den;
          /* krzywizna w potedze p */
        KtoP = pow ( K2, 0.5*md->w );
          /* jakobian */
        jac = ldpt[i*nqkn+q1]*ldpt[j*nqkn+q2]*ldpt[k*nqkn+q3];
        KtoPj = KtoP*jac;

        intf3 += KtoPj;
      }
    }
  }
  intf4 = 0.0;
          /* the two loops iterate over the cubes (i,j,j) for i > j */
  for ( i = 1; i < lkn-2*n; i++ ) {
    pt1 = pt[i*nqkn+q1];
    for ( j = 0; j < i; j++ )
      if ( i != j ) {
        pt2 = pt[j*nqkn+q2];
        SubtractPoints3d ( &pt1, &pt2, &v12 );
        pt3 = pt[j*nqkn+q3];
        SubtractPoints3d ( &pt2, &pt3, &v23 );
        SubtractPoints3d ( &pt3, &pt1, &v31 );
        fv12 = DotProduct3d ( &v12, &v12 );  fl12 = sqrt ( fv12 );
        fv23 = DotProduct3d ( &v23, &v23 );  fl23 = sqrt ( fv23 );
        fv31 = DotProduct3d ( &v31, &v31 );  fl31 = sqrt ( fv31 );
        if ( q2 == q3 ) {
            /* licznik wyrazenia dla krzywizny */
          CrossProduct3d ( &dpt[j*nqkn+q2], &v12, &auxv1 );
          K2num = 4.0*DotProduct3d ( &auxv1, &auxv1 );
            /* mianownik wyrazenia dla krzywizny */
          K2den = DotProduct3d ( &dpt[j*nqkn+q2], &dpt[j*nqkn+q2] )*(fv12*fv12);
        }
        else {
          half = 0.5*(fl12+fl23+fl31);
            /* licznik wyrazenia dla krzywizny */
          K2num = 16.0*(half*(half-fl12))*((half-fl23)*(half-fl31));
            /* mianownik wyrazenia dla krzywizny */
          K2den = fv12*fv23*fv31;
        }
            /* kwadrat krzywizny */
        K2 = K2num/K2den;
            /* krzywizna w potedze p */
        KtoP = pow ( K2, 0.5*md->w );
            /* jakobian */
        jac = ldpt[i*nqkn+q1]*ldpt[j*nqkn+q2]*ldpt[j*nqkn+q3];
        KtoPj = KtoP*jac;

        intf4 += KtoPj;
      }
  }
          /* the loops iterate over the cubes (i,j,j) for i < j */
  for ( i = 0; i < lkn-2*n-1; i++ ) {
    pt1 = pt[i*nqkn+q1];
    for ( j = i+1; j < lkn-2*n; j++ )
      if ( i != j ) {
        pt2 = pt[j*nqkn+q2];
        SubtractPoints3d ( &pt1, &pt2, &v12 );
        pt3 = pt[j*nqkn+q3];
        SubtractPoints3d ( &pt2, &pt3, &v23 );
        SubtractPoints3d ( &pt3, &pt1, &v31 );
        fv12 = DotProduct3d ( &v12, &v12 );  fl12 = sqrt ( fv12 );
        fv23 = DotProduct3d ( &v23, &v23 );  fl23 = sqrt ( fv23 );
        fv31 = DotProduct3d ( &v31, &v31 );  fl31 = sqrt ( fv31 );
        if ( q2 == q3 ) {
            /* licznik wyrazenia dla krzywizny */
          CrossProduct3d ( &dpt[j*nqkn+q2], &v12, &auxv1 );
          K2num = 4.0*DotProduct3d ( &auxv1, &auxv1 );
            /* mianownik wyrazenia dla krzywizny */
          K2den = DotProduct3d ( &dpt[j*nqkn+q2], &dpt[j*nqkn+q2] )*(fv12*fv12);
        }
        else {
          half = 0.5*(fl12+fl23+fl31);
            /* licznik wyrazenia dla krzywizny */
          K2num = 16.0*(half*(half-fl12))*((half-fl23)*(half-fl31));
            /* mianownik wyrazenia dla krzywizny */
          K2den = fv12*fv23*fv31;
        }
            /* kwadrat krzywizny */
        K2 = K2num/K2den;
            /* krzywizna w potedze p */
        KtoP = pow ( K2, 0.5*md->w );
            /* jakobian */
        jac = ldpt[i*nqkn+q1]*ldpt[j*nqkn+q2]*ldpt[j*nqkn+q3];
        KtoPj = KtoP*jac;
        intf4 += KtoPj;
      }
  }
          /* this loop iterates over the cubes (i,i,i) */
  intf5 = 0.0;
  for ( i = 0; i < lkn-2*n; i++ ) {
    pt1 = pt[i*nqkn+q1];
    pt2 = pt[i*nqkn+q2];
    SubtractPoints3d ( &pt1, &pt2, &v12 );
    pt3 = pt[i*nqkn+q3];
    SubtractPoints3d ( &pt2, &pt3, &v23 );
    SubtractPoints3d ( &pt3, &pt1, &v31 );
    fv12 = DotProduct3d ( &v12, &v12 );  fl12 = sqrt ( fv12 );
    fv23 = DotProduct3d ( &v23, &v23 );  fl23 = sqrt ( fv23 );
    fv31 = DotProduct3d ( &v31, &v31 );  fl31 = sqrt ( fv31 );
    if ( q1 == q2 ) {
      if ( q2 == q3 ) {
          /* przypadek specjalny - trzeba obliczyc krzywizne okregu */
          /* scisle stycznego do krzywej */
            /* licznik wyrazenia dla krzywizny */
        CrossProduct3d ( &ddpt[i*nqkn+q1], &dpt[i*nqkn+q1], &auxv1 );
        K2num = DotProduct3d ( &auxv1, &auxv1 );
            /* mianownik wyrazenia dla krzywizny */
        K2den = DotProduct3d ( &dpt[i*nqkn+q1], &dpt[i*nqkn+q1] );
        K2den = K2den*K2den*K2den;
      }
      else {
            /* licznik wyrazenia dla krzywizny */
        CrossProduct3d ( &dpt[i*nqkn+q1], &v31, &auxv1 );
        K2num = 4.0*DotProduct3d ( &auxv1, &auxv1 );
            /* mianownik wyrazenia dla krzywizny */
        K2den = DotProduct3d ( &dpt[i*nqkn+q1], &dpt[i*nqkn+q1] )*(fv31*fv31);
      }
    }
    else if ( q2 == q3 ) {
            /* licznik wyrazenia dla krzywizny */
      CrossProduct3d ( &dpt[i*nqkn+q2], &v12, &auxv1 );
      K2num = 4.0*DotProduct3d ( &auxv1, &auxv1 );
            /* mianownik wyrazenia dla krzywizny */
      K2den = DotProduct3d ( &dpt[i*nqkn+q2], &dpt[i*nqkn+q2] )*(fv12*fv12);
    }
    else if ( q1 == q3 ) {
            /* licznik wyrazenia dla krzywizny */
      CrossProduct3d ( &dpt[i*nqkn+q3], &v23, &auxv1 );
      K2num = 4.0*DotProduct3d ( &auxv1, &auxv1 );
            /* mianownik wyrazenia dla krzywizny */
      K2den = DotProduct3d ( &dpt[i*nqkn+q3], &dpt[i*nqkn+q3] )*(fv23*fv23);
    }
    else {
      half = 0.5*(fl12+fl23+fl31);
            /* licznik wyrazenia dla krzywizny */
      K2num = 16.0*(half*(half-fl12))*((half-fl23)*(half-fl31));
            /* mianownik wyrazenia dla krzywizny */
      K2den = fv12*fv23*fv31;
    }
          /* kwadrat krzywizny */
    K2 = K2num/K2den;
          /* krzywizna w potedze p */
    KtoP = pow ( K2, 0.5*md->w );
          /* jakobian */
    jac = ldpt[i*nqkn+q1]*ldpt[i*nqkn+q2]*ldpt[i*nqkn+q3];
    KtoPj = KtoP*jac;
    intf5 += KtoPj;
  }
  data->_f[(q1*nqkn+q2)*nqkn+q3] = 6.0*intf3 + 3.0*intf4 + intf5;
  return true;
} /*_intKM_F*/

boolean _mengerc_intF ( mengerc_data *md, double *func )
{
  void           *sp;
  intKM_job_desc data;
  int            lkn, n, npoints, nqkn;
  double         *knots, *qkn, *qc;
  point3d        *cpoints;
  vector3d       *pt, *dpt, *ddpt;
  double         *ldpt;
  double         *_f, ff1, ff2, ff3;
  int            i, j, k, l, q1;
  int3           jobsize;
  boolean        success;

  sp = pkv_GetScratchMemTop ();
  memset ( &data, 0, sizeof(intKM_job_desc) );
  data.md = md;
  n       = md->deg;
  lkn     = md->lkn;
  knots   = md->knots;
  cpoints = md->cpoints;
  nqkn    = md->nqkn;
  qkn     = md->qkn;
        /* liczba punktow krzywej */
  npoints = data.npoints = nqkn*(lkn-2*n);
  pt = data.pt = pkv_GetScratchMem ( 3*npoints*sizeof(point3d) +
                                     npoints*sizeof(double) );
  _f = pkv_GetScratchMemd ( nqkn*nqkn*nqkn );
  if ( !pt || !_f )
    goto failure;
  data._f = _f;
  dpt = data.dpt = &pt[npoints];
  ddpt = data.ddpt = &dpt[npoints];
  ldpt = data.ldpt = (double*)&ddpt[npoints];
        /* znajdz punkty krzywej i pochodne pierwszego i drugiego rzedu */
  for ( i = j = 0;  i < lkn-2*n;  i++ )
    for ( q1 = 0;  q1 < nqkn;  q1 ++, j++ ) {
      mbs_deBoorDer2C3d ( n, 2*n+1, knots, &cpoints[i], (double)n+qkn[q1],
                          &pt[j], &dpt[j], &ddpt[j] );
      ldpt[j] = sqrt ( DotProduct3d ( &dpt[j], &dpt[j] ) );
    }
        /* oblicz calki */
  if ( md->npthr > 1 ) {
    jobsize.x = jobsize.y = jobsize.z = nqkn;
    pkv_SetPThreadsToWork ( &jobsize, md->npthr, 1048576, 1048576,
                            (void*)&data, _intKM_F, NULL, NULL, &success );
    if ( !success )
      goto failure;
  }
  else {
    for ( jobsize.x = 0; jobsize.x < nqkn; jobsize.x++ )
      for ( jobsize.y = 0; jobsize.y < nqkn; jobsize.y++ )
        for ( jobsize.z = 0; jobsize.z < nqkn; jobsize.z++ )
          _intKM_F ( (void*)&data, &jobsize );
  }
        /* sumuj calki */
  qc = md->qc;
  ff1 = 0.0;
  for ( i = l = 0;  i < nqkn;  i++ ) {
    ff2 = 0.0;
    for ( j = 0; j < nqkn; j++ ) {
      ff3 = 0.0;
      for ( k = 0;  k < nqkn;  k++, l++ )
        ff3 += qc[k]*_f[l];
      ff2 += qc[j]*ff3;
    }
    ff1 += qc[i]*ff2;
  }

  *func = ff1;
  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*_mengerc_intF*/

/* ///////////////////////////////////////////////////////////////////////// */
static int _imc_AdjustGrad ( int n, int N3, int i, int k, int ci3, int nia3,
                             double *df )
{
  int z, t, l;

  z = 3*(i-k+n+1) - N3;
  if ( z > 0 ) {
    t = nia3-z;
    for ( l = 0; l < z; l++ )
      df[l] += df[l+t];
    return t;
  }
  else
    return nia3;
} /*_imc_AdjustGrad*/

static void _imc_CopyGrad ( int N3, int i, int j, int k,
                            int ci3, int cj3, int nia3,
                            double *df, double *g )
{
  int c3, ic3;

  for ( c3 = 0; c3 < cj3; c3 += 3 ) {
    ic3 = 3*k+c3;  if ( ic3 >= N3 ) ic3 -= N3;
    g[ic3]   += df[c3];
    g[ic3+1] += df[c3+1];
    g[ic3+2] += df[c3+2];
  }
  for ( ; c3 < ci3; c3 += 3 ) {
    ic3 = 3*j+c3-cj3;  if ( ic3 >= N3 ) ic3 -= N3;
    g[ic3]   += df[c3];
    g[ic3+1] += df[c3+1];
    g[ic3+2] += df[c3+2];
  }
  for ( ; c3 < nia3; c3 += 3 ) {
    ic3 = 3*i+c3-ci3;  if ( ic3 >= N3 ) ic3 -= N3;
    g[ic3]   += df[c3];
    g[ic3+1] += df[c3+1];
    g[ic3+2] += df[c3+2];
  }
} /*_imc_CopyGrad*/

static void _imc_dProduct ( int n, double a, double b, double *ab,
                            double *da, double *db, double *dab )
{
  *ab = a*b;
  pkn_MatrixLinCombd ( 1, n, 0, da, b, 0, db, a, 0, dab );
} /*_imc_dProduct*/

static void _imc_dRatio ( int n, double g, double h, double *f,
                          double *dg, double *dh, double *df )
{
  double _f;

  *f = _f = g/h;
  pkn_AddMatrixMd ( 1, n, 0, dg, 0, dh, -_f, 0, df );
  pkn_MultMatrixNumd ( 1, n, 0, df, 1.0/h, 0, df );
} /*_imc_dRatio*/

static void _imc_dPower ( int n, double a, double *da, double p,
                          double *atop, double *datop )
{
  double _atop;

  *atop = _atop = pow ( a, p );
  pkn_MultMatrixNumd ( 1, n, 0, da, p*_atop/a, 0, datop );
} /*_imc_dPower*/

static void _imc_dK2num ( int nia,
                          double fl12, double fl23, double fl31, double half,
                          double *dfl12, double *dfl23, double *dfl31, double *dhalf,
                          double *daux1, double *daux2, double *daux3, double *daux4,
                          double *K2num, double *dK2num )
{
  double a1, a2, a3, a4, a5;

  a1 = half-fl12;  a2 = half-fl23;  a3 = half-fl31;
  a4 = half*a1;  a5 = a2*a3;
  *K2num = 16.0*a4*a5;
  pkn_SubtractMatrixd ( 1, nia, 0, dhalf, 0, dfl12, 0, daux1 );
  pkn_MatrixLinCombd ( 1, nia, 0, dhalf, a1, 0, daux1, half, 0, daux2 );
  pkn_SubtractMatrixd ( 1, nia, 0, dhalf, 0, dfl23, 0, daux1 );
  pkn_SubtractMatrixd ( 1, nia, 0, dhalf, 0, dfl31, 0, daux3 );
  pkn_MatrixLinCombd ( 1, nia, 0, daux1, a3, 0, daux3, a2, 0, daux4 );
  pkn_MatrixLinCombd ( 1, nia, 0, daux2, 16.0*a5, 0, daux4, 16.0*a4, 0, dK2num );
} /*_imc_dK2num*/

static void _imc_dK2den ( int nia, double fv12, double fv23, double fv31,
                          double *dfv12, double *dfv23, double *dfv31,
                          double *K2den, double *dK2den )
{
  double a1, a2, a3;

  a1 = fv12*fv23;  a2 = fv23*fv31;  a3 = fv31*fv12;
  *K2den = a1*fv31;
  pkn_MatrixLinCombd ( 1, nia, 0, dfv12, a2, 0, dfv23, a3, 0, dK2den );
  pkn_AddMatrixMd ( 1, nia, 0, dK2den, 0, dfv31, a1, 0, dK2den );
} /*_imc_dK2den*/

static void _imc_dSpecialK2num ( int n, int N, double *bsf, double *bsf1,
                                 int q1, int q2, int i, int j, int ci, int nia,
                                 vector3d *v12, vector3d *dpt,
                                 double *K2num, double *dK2num )
{
  double   a1, a2;
  vector3d auxv1;
  int      ii, ic, jc;

  memset ( dK2num, 0, nia*sizeof(double) );
  CrossProduct3d ( dpt, v12, &auxv1 );
  *K2num = 4.0*DotProduct3d ( &auxv1, &auxv1 );
  if ( i >= j ) {
    for ( ii = ic = 0;  ii < n;  ii++, ic += 3 ) {
      a1 = bsf1[q2*n+ii];
      a2 = (v12->y*auxv1.z - v12->z*auxv1.y)*a1;
      dK2num[ic  ] -= a2;  dK2num[ic+3] += a2;
      a2 = (v12->z*auxv1.x - v12->x*auxv1.z)*a1;
      dK2num[ic+1] -= a2;  dK2num[ic+4] += a2;
      a2 = (v12->x*auxv1.y - v12->y*auxv1.x)*a1;
      dK2num[ic+2] -= a2;  dK2num[ic+5] += a2;
    }
    for ( ii = ic = 0, jc = ci;  ii <= n;  ii++, ic += 3, jc += 3 ) {
      a1 = bsf[q2*(n+1)+ii];
      dK2num[ic  ] -= (dpt->z*auxv1.y - dpt->y*auxv1.z)*a1;
      dK2num[ic+1] -= (dpt->x*auxv1.z - dpt->z*auxv1.x)*a1;
      dK2num[ic+2] -= (dpt->y*auxv1.x - dpt->x*auxv1.y)*a1;
      a1 = bsf[q1*(n+1)+ii];
      dK2num[jc  ] += (dpt->z*auxv1.y - dpt->y*auxv1.z)*a1;
      dK2num[jc+1] += (dpt->x*auxv1.z - dpt->z*auxv1.x)*a1;
      dK2num[jc+2] += (dpt->y*auxv1.x - dpt->x*auxv1.y)*a1;
    }
    if ( i > j )
      nia = _imc_AdjustGrad ( n, N, i, j, ci, nia, dK2num );
  }
  else {
    for ( ii = 0, ic = ci;  ii < n;  ii++, ic += 3 ) {
      a1 = bsf1[q2*n+ii];
      a2 = (v12->y*auxv1.z - v12->z*auxv1.y)*a1;
      dK2num[ic  ] -= a2;  dK2num[ic+3] += a2;
      a2 = (v12->z*auxv1.x - v12->x*auxv1.z)*a1;
      dK2num[ic+1] -= a2;  dK2num[ic+4] += a2;
      a2 = (v12->x*auxv1.y - v12->y*auxv1.x)*a1;
      dK2num[ic+2] -= a2;  dK2num[ic+5] += a2;
    }
    for ( ii = jc = 0, ic = ci;  ii <= n;  ii++, ic += 3, jc += 3 ) {
      a1 = bsf[q2*(n+1)+ii];
      dK2num[ic  ] -= (dpt->z*auxv1.y - dpt->y*auxv1.z)*a1;
      dK2num[ic+1] -= (dpt->x*auxv1.z - dpt->z*auxv1.x)*a1;
      dK2num[ic+2] -= (dpt->y*auxv1.x - dpt->x*auxv1.y)*a1;
      a1 = bsf[q1*(n+1)+ii];
      dK2num[jc  ] += (dpt->z*auxv1.y - dpt->y*auxv1.z)*a1;
      dK2num[jc+1] += (dpt->x*auxv1.z - dpt->z*auxv1.x)*a1;
      dK2num[jc+2] += (dpt->y*auxv1.x - dpt->x*auxv1.y)*a1;
    }
    nia = _imc_AdjustGrad ( n, N, j, i, ci, nia, dK2num );
  }
  pkn_MultMatrixNumd ( 1, nia, 0, dK2num, 8.0, 0, dK2num );
} /*_imc_dSpecialK2num*/

static void _imc_dSpecialK2den ( int n, int N, double *bsf1, int i, int j,
                                 int ci, int nia,
                                 vector3d *dcp, vector3d *dpt, double fv12,
                                 double *dfv12, double *daux1, double *daux2,
                                 double *K2den, double *dK2den )
{
  double a1, a2, a3, a4;
  int    ii, jj, ic;

  a1 = DotProduct3d ( dpt, dpt );
  memset ( daux1, 0, nia*sizeof(double) );
  for ( ii = 0, ic = ci;  ii < n;  ii++, ic += 3 ) {
    a2 = 2.0*bsf1[ii];
    for ( jj = 0;  jj < n;  jj++ ) {
      a3 = a2*bsf1[jj];
      a4 = dcp[jj].x*a3;  daux1[ic  ] -= a4;  daux1[ic+3] += a4;
      a4 = dcp[jj].y*a3;  daux1[ic+1] -= a4;  daux1[ic+4] += a4;
      a4 = dcp[jj].z*a3;  daux1[ic+2] -= a4;  daux1[ic+5] += a4;
    }
  }
  if ( i < j )
    nia = _imc_AdjustGrad ( n, N, j, i, ci, nia, daux1 );
  pkn_MultMatrixNumd ( 1, nia, 0, dfv12, 2.0*fv12, 0, daux2 );
  _imc_dProduct ( nia, a1, fv12*fv12, K2den, daux1, daux2, dK2den );
} /*_imc_dSpecialK2den*/

static boolean _intKM_FG ( void *usrdata, int3 *jobnum )
{
  void           *sp;
  intKM_job_desc *data;
  mengerc_data       *md;
  int            n, lkn, nqkn, ncp, nvcp, N;
  vector3d       *dcp, *pt, *dpt, *ddpt;
  double         *ldpt, *dldpt;
  int            q1, q2, q3, i, j, k, c, ic, ii, jj, kk;
  double         *bsf, *bsf1, *dbsf1;
  point3d        pt1, pt2, pt3;
  vector3d       v12, v23, v31, auxv1;
  int            ti, tj, tk, nia, nia1, maxnia, ci, cj;
  double         a, a1, a2, a3, a1a2, a2a3, a3a1;
  double         fv12, fv23, fv31, fl12, fl23, fl31,
                 half, K2num, K2den, K2, KtoP, jac, KtoPj;
  double         *dfv12, *dfv23, *dfv31, *dfl12, *dfl23, *dfl31,
                 *dhalf, *dK2num, *dK2den, *dK2, *dKtoP, *djac, *dKtoPj,
                 *daux1, *daux2, *daux3, *daux4;
  double         intf3, intf4, intf5;
  double         *g3, *g4, *g5;

  sp = pkv_GetScratchMemTop ();
  data = (intKM_job_desc*)usrdata;
  md    = data->md;
  pt    = data->pt;
  dpt   = data->dpt;
  ddpt  = data->ddpt;
  ldpt  = data->ldpt;
  dldpt = data->dldpt;
  dcp   = data->dcp;

  n     = md->deg;
  lkn   = md->lkn;
  ncp   = lkn-n;
  nvcp  = ncp-n;
  N     = 3*nvcp;
  nqkn  = md->nqkn;
  bsf   = md->bsf;
  bsf1  = md->bsf1;
  dbsf1 = md->dbsf1;

  q1 = jobnum->x;
  q2 = jobnum->y;
  q3 = jobnum->z;

  g3 = pkv_GetScratchMemd ( 3*N );
  if ( !g3 )
    goto failure;
  g4 = &g3[N];  g5 = &g4[N];
        /* maximal number of relevant variables in any unit cube */
  maxnia = 9*(n+1);
  dfv12 = pkv_GetScratchMemd ( 13*maxnia );
  if ( !dfv12 )
    goto failure;
  dfv23 = &dfv12[maxnia];  dfv31 = &dfv23[maxnia];
  dfl12 = &dfv31[maxnia];  dfl23 = &dfl12[maxnia];   dfl31 = &dfl23[maxnia];
  dhalf = &dfl31[maxnia];  dK2num = &dhalf[maxnia];  dK2den = &dK2num[maxnia];
  dK2 = &dK2den[maxnia];  dKtoP = &dK2[maxnia];  djac = &dKtoP[maxnia];
  dKtoPj = &djac[maxnia];
  daux1 = dKtoP, daux2 = dKtoPj;  daux3 = dK2den;  daux4 = djac;
          /* the three loops iterate over the cubes (i,j,k) for i>j>k */
  intf3 = 0.0;
  memset ( g3, 0, N*sizeof(double) );
  for ( i = 2; i < lkn-2*n; i++ ) {
    ti = 6*(n+1);
    pt1 = pt[i*nqkn+q1];
    for ( j = 1; j < i; j++ ) {
      tj = ti - 3*min(i-j,n+1);
      pt2 = pt[j*nqkn+q2];
      SubtractPoints3d ( &pt1, &pt2, &v12 );
      for ( k = 0; k < j; k++ ) {
        tk = tj - 3*min(j-k,n+1);
        nia = maxnia-tk;  ci = ti-tk;  cj = tj-tk;
        memset ( dfv12, 0, nia*sizeof(double) );
        memset ( dfv23, 0, nia*sizeof(double) );
        memset ( dfv31, 0, nia*sizeof(double) );

        pt3 = pt[k*nqkn+q3];
        SubtractPoints3d ( &pt2, &pt3, &v23 );
        SubtractPoints3d ( &pt3, &pt1, &v31 );
        fv12 = DotProduct3d ( &v12, &v12 );  fl12 = sqrt ( fv12 );
        fv23 = DotProduct3d ( &v23, &v23 );  fl23 = sqrt ( fv23 );
        fv31 = DotProduct3d ( &v31, &v31 );  fl31 = sqrt ( fv31 );

        for ( c = ic = 0;  c <= n;  c++, ic += 3 ) {
          a1 = 2.0*bsf[(n+1)*q1+c];
          a2 = 2.0*bsf[(n+1)*q2+c];
          a3 = 2.0*bsf[(n+1)*q3+c];
          dfv12[ci+ic]   += a1*v12.x;
          dfv12[ci+ic+1] += a1*v12.y;
          dfv12[ci+ic+2] += a1*v12.z;
          dfv12[cj+ic]   -= a2*v12.x;
          dfv12[cj+ic+1] -= a2*v12.y;
          dfv12[cj+ic+2] -= a2*v12.z;
          dfv23[cj+ic]   += a2*v23.x;
          dfv23[cj+ic+1] += a2*v23.y;
          dfv23[cj+ic+2] += a2*v23.z;
          dfv23[ic]   -= a3*v23.x;
          dfv23[ic+1] -= a3*v23.y;
          dfv23[ic+2] -= a3*v23.z;
          dfv31[ic]   += a3*v31.x;
          dfv31[ic+1] += a3*v31.y;
          dfv31[ic+2] += a3*v31.z;
          dfv31[ci+ic]   -= a1*v31.x;
          dfv31[ci+ic+1] -= a1*v31.y;
          dfv31[ci+ic+2] -= a1*v31.z;
        }
        nia1 = _imc_AdjustGrad ( n, N, i, k, ci, nia, dfv12 );
        _imc_AdjustGrad ( n, N, i, k, ci, nia, dfv23 );
        _imc_AdjustGrad ( n, N, i, k, ci, nia, dfv31 );

        pkn_MultMatrixNumd ( 1, nia1, 0, dfv12, 0.5/fl12, 0, dfl12 );
        pkn_MultMatrixNumd ( 1, nia1, 0, dfv23, 0.5/fl23, 0, dfl23 );
        pkn_MultMatrixNumd ( 1, nia1, 0, dfv31, 0.5/fl31, 0, dfl31 );

        half = 0.5*(fl12+fl23+fl31);
        pkn_AddMatrixd ( 1, nia1, 0, dfl12, 0, dfl23, 0, dhalf );
        pkn_AddMatrixd ( 1, nia1, 0, dhalf, 0, dfl31, 0, dhalf );
        pkn_MultMatrixNumd ( 1, nia1, 0, dhalf, 0.5, 0, dhalf );
          /* licznik wyrazenia dla krzywizny */
        _imc_dK2num ( nia1, fl12, fl23, fl31, half,
                      dfl12, dfl23, dfl31, dhalf, daux1, daux2, daux3, daux4,
                      &K2num, dK2num );
          /* mianownik wyrazenia dla krzywizny */
        _imc_dK2den ( nia1, fv12, fv23, fv31, dfv12, dfv23, dfv31,
                      &K2den, dK2den );
          /* kwadrat krzywizny */
        _imc_dRatio ( nia1, K2num, K2den, &K2, dK2num, dK2den, dK2 );
          /* krzywizna w potedze p */
        _imc_dPower ( nia1, K2, dK2, 0.5*md->w, &KtoP, dKtoP );
          /* jakobian */
        memset ( djac, 0, nia*sizeof(double) );
        a1 = ldpt[i*nqkn+q1];
        a2 = ldpt[j*nqkn+q2];
        a3 = ldpt[k*nqkn+q3];
        a1a2 = a1*a2;
        a2a3 = a2*a3;
        a3a1 = a3*a1;
        jac = a1a2*a3;
        ii = 3*(n+1)*(nqkn*i+q1);
        jj = 3*(n+1)*(nqkn*j+q2);
        kk = 3*(n+1)*(nqkn*k+q3);
        for ( c = 0; c < 3*(n+1); c++ ) {
          djac[ci+c] += dldpt[ii+c]*a2a3;
          djac[cj+c] += dldpt[jj+c]*a3a1;
          djac[c]    += dldpt[kk+c]*a1a2;
        }
        _imc_AdjustGrad ( n, N, i, k, ci, nia, djac );
         /* pelna funkcja podcalkowa */
        _imc_dProduct ( nia1, KtoP, jac, &KtoPj, dKtoP, djac, dKtoPj );
        intf3 += KtoPj;
        _imc_CopyGrad ( N, i, j, k, ci, cj, nia1, dKtoPj, g3 );
      }
    }
  }
  intf4 = 0.0;
  memset ( g4, 0, N*sizeof(double) );
          /* the two loops iterate over the cubes (i,j,j) for i > j */
  for ( i = 1; i < lkn-2*n; i++ ) {
    ti = 6*(n+1);
    pt1 = pt[i*nqkn+q1];
    for ( j = 0; j < i; j++ ) {
      tj = ti - 3*min(i-j,n+1);
      nia = maxnia-tj;  ci = ti-tj;
      memset ( dfv12, 0, nia*sizeof(double) );
      memset ( dfv23, 0, nia*sizeof(double) );
      memset ( dfv31, 0, nia*sizeof(double) );

      pt2 = pt[j*nqkn+q2];
      SubtractPoints3d ( &pt1, &pt2, &v12 );
      pt3 = pt[j*nqkn+q3];
      SubtractPoints3d ( &pt2, &pt3, &v23 );
      SubtractPoints3d ( &pt3, &pt1, &v31 );
      fv12 = DotProduct3d ( &v12, &v12 );  fl12 = sqrt ( fv12 );
      fv23 = DotProduct3d ( &v23, &v23 );  fl23 = sqrt ( fv23 );
      fv31 = DotProduct3d ( &v31, &v31 );  fl31 = sqrt ( fv31 );

      for ( c = ic = 0;  c <= n;  c++, ic += 3 ) {
        a1 = 2.0*bsf[(n+1)*q1+c];
        a2 = 2.0*bsf[(n+1)*q2+c];
        a3 = 2.0*bsf[(n+1)*q3+c];
        dfv12[ci+ic]   += a1*v12.x;
        dfv12[ci+ic+1] += a1*v12.y;
        dfv12[ci+ic+2] += a1*v12.z;
        dfv12[ic]   -= a2*v12.x;
        dfv12[ic+1] -= a2*v12.y;
        dfv12[ic+2] -= a2*v12.z;
        dfv23[ic]   += a2*v23.x;
        dfv23[ic+1] += a2*v23.y;
        dfv23[ic+2] += a2*v23.z;
        dfv23[ic]   -= a3*v23.x;
        dfv23[ic+1] -= a3*v23.y;
        dfv23[ic+2] -= a3*v23.z;
        dfv31[ic]   += a3*v31.x;
        dfv31[ic+1] += a3*v31.y;
        dfv31[ic+2] += a3*v31.z;
        dfv31[ci+ic]   -= a1*v31.x;
        dfv31[ci+ic+1] -= a1*v31.y;
        dfv31[ci+ic+2] -= a1*v31.z;
      }
      nia1 = _imc_AdjustGrad ( n, N, i, j, ci, nia, dfv12 );
      _imc_AdjustGrad ( n, N, i, j, ci, nia, dfv23 );
      _imc_AdjustGrad ( n, N, i, j, ci, nia, dfv31 );
      pkn_MultMatrixNumd ( 1, nia1, 0, dfv12, 0.5/fl12, 0, dfl12 );
      memset ( dK2num, 0, nia1*sizeof(double) );
      memset ( dK2den, 0, nia1*sizeof(double) );
      if ( q2 == q3 ) {
          /* licznik wyrazenia dla krzywizny */
        _imc_dSpecialK2num ( n, N, bsf, bsf1, q1, q2, i, j, ci, nia,
                             &v12, &dpt[j*nqkn+q2], &K2num, dK2num );
          /* mianownik wyrazenia dla krzywizny */
        _imc_dSpecialK2den ( n, N, &bsf1[q2*n], i, j, 0, nia, &dcp[j], &dpt[j*nqkn+q2],
                             fv12, dfv12, daux1, daux2, &K2den, dK2den );
      }
      else {
        pkn_MultMatrixNumd ( 1, nia1, 0, dfv23, 0.5/fl23, 0, dfl23 );
        pkn_MultMatrixNumd ( 1, nia1, 0, dfv31, 0.5/fl31, 0, dfl31 );
        half = 0.5*(fl12+fl23+fl31);
        pkn_AddMatrixd ( 1, nia1, 0, dfl12, 0, dfl23, 0, dhalf );
        pkn_AddMatrixd ( 1, nia1, 0, dhalf, 0, dfl31, 0, dhalf );
        pkn_MultMatrixNumd ( 1, nia1, 0, dhalf, 0.5, 0, dhalf );

          /* licznik wyrazenia dla krzywizny */
        _imc_dK2num ( nia1, fl12, fl23, fl31, half,
                      dfl12, dfl23, dfl31, dhalf, daux1, daux2, daux3, daux4,
                      &K2num, dK2num );
          /* mianownik wyrazenia dla krzywizny */
        _imc_dK2den ( nia1, fv12, fv23, fv31, dfv12, dfv23, dfv31,
                      &K2den, dK2den );
      }
          /* kwadrat krzywizny */
      _imc_dRatio ( nia1, K2num, K2den, &K2, dK2num, dK2den, dK2 );
          /* krzywizna w potedze p */
      _imc_dPower ( nia1, K2, dK2, 0.5*md->w, &KtoP, dKtoP );
          /* jakobian */
      memset ( djac, 0, nia*sizeof(double) );
      a1 = ldpt[i*nqkn+q1];
      a2 = ldpt[j*nqkn+q2];
      a3 = ldpt[j*nqkn+q3];
      a1a2 = a1*a2;
      a2a3 = a2*a3;
      a3a1 = a3*a1;
      jac = a1a2*a3;
      ii = 3*(n+1)*(nqkn*i+q1);
      jj = 3*(n+1)*(nqkn*j+q2);
      kk = 3*(n+1)*(nqkn*j+q3);
      for ( c = 0; c < 3*(n+1); c++ ) {
        djac[ci+c] += dldpt[ii+c]*a2a3;
        djac[c]    += dldpt[jj+c]*a3a1;
        djac[c]    += dldpt[kk+c]*a1a2;
      }
      _imc_AdjustGrad ( n, N, i, j, ci, nia, djac );
    /* pelna funkcja podcalkowa */
      _imc_dProduct ( nia1, KtoP, jac, &KtoPj, dKtoP, djac, dKtoPj );
      intf4 += KtoPj;
      _imc_CopyGrad ( N, i, j, j, ci, 0, nia1, dKtoPj, g4 );
    }
  }
          /* the two loops iterate over the cubes (i,j,j) for i < j */
  tj = 6*(n+1);
  for ( i = 0; i < lkn-2*n-1; i++ ) {
    pt1 = pt[i*nqkn+q1];
    for ( j = i+1; j < lkn-2*n; j++ ) {
      ti = tj - 3*min(j-i,n+1);
      nia = maxnia-ti;  cj = tj-ti;
      memset ( dfv12, 0, nia*sizeof(double) );
      memset ( dfv23, 0, nia*sizeof(double) );
      memset ( dfv31, 0, nia*sizeof(double) );

      pt2 = pt[j*nqkn+q2];
      SubtractPoints3d ( &pt1, &pt2, &v12 );
      pt3 = pt[j*nqkn+q3];
      SubtractPoints3d ( &pt2, &pt3, &v23 );
      SubtractPoints3d ( &pt3, &pt1, &v31 );
      fv12 = DotProduct3d ( &v12, &v12 );  fl12 = sqrt ( fv12 );
      fv23 = DotProduct3d ( &v23, &v23 );  fl23 = sqrt ( fv23 );
      fv31 = DotProduct3d ( &v31, &v31 );  fl31 = sqrt ( fv31 );

      for ( c = ic = 0;  c <= n;  c++, ic += 3 ) {
        a1 = 2.0*bsf[(n+1)*q1+c];
        a2 = 2.0*bsf[(n+1)*q2+c];
        a3 = 2.0*bsf[(n+1)*q3+c];
        dfv12[ic]   += a1*v12.x;
        dfv12[ic+1] += a1*v12.y;
        dfv12[ic+2] += a1*v12.z;
        dfv12[cj+ic]   -= a2*v12.x;
        dfv12[cj+ic+1] -= a2*v12.y;
        dfv12[cj+ic+2] -= a2*v12.z;
        dfv23[cj+ic]   += a2*v23.x;
        dfv23[cj+ic+1] += a2*v23.y;
        dfv23[cj+ic+2] += a2*v23.z;
        dfv23[cj+ic]   -= a3*v23.x;
        dfv23[cj+ic+1] -= a3*v23.y;
        dfv23[cj+ic+2] -= a3*v23.z;
        dfv31[cj+ic]   += a3*v31.x;
        dfv31[cj+ic+1] += a3*v31.y;
        dfv31[cj+ic+2] += a3*v31.z;
        dfv31[ic]   -= a1*v31.x;
        dfv31[ic+1] -= a1*v31.y;
        dfv31[ic+2] -= a1*v31.z;
      }
      nia1 = _imc_AdjustGrad ( n, N, j, i, cj, nia, dfv12 );
      _imc_AdjustGrad ( n, N, j, i, cj, nia, dfv23 );
      _imc_AdjustGrad ( n, N, j, i, cj, nia, dfv31 );
      pkn_MultMatrixNumd ( 1, nia1, 0, dfv12, 0.5/fl12, 0, dfl12 );

      memset ( dK2num, 0, nia1*sizeof(double) );
      memset ( dK2den, 0, nia1*sizeof(double) );
      if ( q2 == q3 ) {
          /* licznik wyrazenia dla krzywizny */
        _imc_dSpecialK2num ( n, N, bsf, bsf1, q1, q2, i, j, cj, nia,
                             &v12, &dpt[j*nqkn+q2], &K2num, dK2num );
          /* mianownik wyrazenia dla krzywizny */
        _imc_dSpecialK2den ( n, N, &bsf1[q2*n], i, j, cj, nia, &dcp[j], &dpt[j*nqkn+q2],
                             fv12, dfv12, daux1, daux2, &K2den, dK2den );
      }
      else {
        pkn_MultMatrixNumd ( 1, nia1, 0, dfv23, 0.5/fl23, 0, dfl23 );
        pkn_MultMatrixNumd ( 1, nia1, 0, dfv31, 0.5/fl31, 0, dfl31 );
        half = 0.5*(fl12+fl23+fl31);
        pkn_AddMatrixd ( 1, nia1, 0, dfl12, 0, dfl23, 0, dhalf );
        pkn_AddMatrixd ( 1, nia1, 0, dhalf, 0, dfl31, 0, dhalf );
        pkn_MultMatrixNumd ( 1, nia1, 0, dhalf, 0.5, 0, dhalf );

          /* licznik wyrazenia dla krzywizny */
        _imc_dK2num ( nia1, fl12, fl23, fl31, half,
                      dfl12, dfl23, dfl31, dhalf, daux1, daux2, daux3, daux4,
                      &K2num, dK2num );
          /* mianownik wyrazenia dla krzywizny */
        _imc_dK2den ( nia1, fv12, fv23, fv31, dfv12, dfv23, dfv31,
                      &K2den, dK2den );
      }
          /* kwadrat krzywizny */
      _imc_dRatio ( nia1, K2num, K2den, &K2, dK2num, dK2den, dK2 );
          /* krzywizna w potedze p */
      _imc_dPower ( nia1, K2, dK2, 0.5*md->w, &KtoP, dKtoP );
          /* jakobian */
      memset ( djac, 0, nia*sizeof(double) );
      a1 = ldpt[i*nqkn+q1];
      a2 = ldpt[j*nqkn+q2];
      a3 = ldpt[j*nqkn+q3];
      a1a2 = a1*a2;
      a2a3 = a2*a3;
      a3a1 = a3*a1;
      jac = a1a2*a3;
      ii = 3*(n+1)*(nqkn*i+q1);
      jj = 3*(n+1)*(nqkn*j+q2);
      kk = 3*(n+1)*(nqkn*j+q3);
      for ( c = 0; c < 3*(n+1); c++ ) {
        djac[c]    += dldpt[ii+c]*a2a3;
        djac[cj+c] += dldpt[jj+c]*a3a1;
        djac[cj+c] += dldpt[kk+c]*a1a2;
      }
      _imc_AdjustGrad ( n, N, j, i, cj, nia, djac );
         /* pelna funkcja podcalkowa */
      _imc_dProduct ( nia1, KtoP, jac, &KtoPj, dKtoP, djac, dKtoPj );
      intf4 += KtoPj;
      _imc_CopyGrad ( N, j, j, i, cj, cj, nia1, dKtoPj, g4 );
    }
  }
          /* this loop iterates over the cubes (i,i,i) */
  intf5 = 0.0;
  memset ( g5, 0, N*sizeof(double) );
  nia = nia1 = 3*(n+1);
  for ( i = 0; i < lkn-2*n; i++ ) {
    memset ( dfv12, 0, nia*sizeof(double) );
    memset ( dfv23, 0, nia*sizeof(double) );
    memset ( dfv31, 0, nia*sizeof(double) );

    pt1 = pt[i*nqkn+q1];
    pt2 = pt[i*nqkn+q2];
    SubtractPoints3d ( &pt1, &pt2, &v12 );
    pt3 = pt[i*nqkn+q3];
    SubtractPoints3d ( &pt2, &pt3, &v23 );
    SubtractPoints3d ( &pt3, &pt1, &v31 );
    fv12 = DotProduct3d ( &v12, &v12 );  fl12 = sqrt ( fv12 );
    fv23 = DotProduct3d ( &v23, &v23 );  fl23 = sqrt ( fv23 );
    fv31 = DotProduct3d ( &v31, &v31 );  fl31 = sqrt ( fv31 );

    for ( c = ic = 0;  c <= n;  c++, ic += 3 ) {
      a1 = 2.0*bsf[(n+1)*q1+c];
      a2 = 2.0*bsf[(n+1)*q2+c];
      a3 = 2.0*bsf[(n+1)*q3+c];
      dfv12[ic]   += a1*v12.x;
      dfv12[ic+1] += a1*v12.y;
      dfv12[ic+2] += a1*v12.z;
      dfv12[ic]   -= a2*v12.x;
      dfv12[ic+1] -= a2*v12.y;
      dfv12[ic+2] -= a2*v12.z;
      dfv23[ic]   += a2*v23.x;
      dfv23[ic+1] += a2*v23.y;
      dfv23[ic+2] += a2*v23.z;
      dfv23[ic]   -= a3*v23.x;
      dfv23[ic+1] -= a3*v23.y;
      dfv23[ic+2] -= a3*v23.z;
      dfv31[ic]   += a3*v31.x;
      dfv31[ic+1] += a3*v31.y;
      dfv31[ic+2] += a3*v31.z;
      dfv31[ic]   -= a1*v31.x;
      dfv31[ic+1] -= a1*v31.y;
      dfv31[ic+2] -= a1*v31.z;
    }
    memset ( dK2num, 0, nia1*sizeof(double) );
    memset ( dK2den, 0, nia1*sizeof(double) );
    if ( q1 == q2 ) {
      if ( q2 == q3 ) {
          /* przypadek specjalny - trzeba obliczyc krzywizne okregu */
          /* scisle stycznego do krzywej */
            /* licznik wyrazenia dla krzywizny */
        CrossProduct3d ( &ddpt[i*nqkn+q1], &dpt[i*nqkn+q1], &auxv1 );
        K2num = DotProduct3d ( &auxv1, &auxv1 );
        for ( ii = 1, ci = 3;  ii < n;  ii++, ci += 3 )
          for ( jj = cj = 0;  jj < ii;  jj++, cj += 3 ) {
            a1 = dbsf1[q1*n+ii]*bsf1[q1*n+jj] - dbsf1[q1*n+jj]*bsf1[q1*n+ii];
            a2 = (dcp[i+jj].y*auxv1.z - dcp[i+jj].z*auxv1.y)*a1;
            dK2num[ci  ] -= a2;  dK2num[ci+3] += a2;
            a2 = (dcp[i+jj].z*auxv1.x - dcp[i+jj].x*auxv1.z)*a1;
            dK2num[ci+1] -= a2;  dK2num[ci+4] += a2;
            a2 = (dcp[i+jj].x*auxv1.y - dcp[i+jj].y*auxv1.x)*a1;
            dK2num[ci+2] -= a2;  dK2num[ci+5] += a2;
            a2 = (dcp[i+ii].y*auxv1.z - dcp[i+ii].z*auxv1.y)*a1;
            dK2num[cj  ] += a2;  dK2num[cj+3] -= a2;
            a2 = (dcp[i+ii].z*auxv1.x - dcp[i+ii].x*auxv1.z)*a1;
            dK2num[cj+1] += a2;  dK2num[cj+4] -= a2;
            a2 = (dcp[i+ii].x*auxv1.y - dcp[i+ii].y*auxv1.x)*a1;
            dK2num[cj+2] += a2;  dK2num[cj+5] -= a2;
          }
        pkn_AddMatrixd ( 1, nia1, 0, dK2num, 0, dK2num, 0, dK2num );
            /* mianownik wyrazenia dla krzywizny */
        a2 = DotProduct3d ( &dpt[i*nqkn+q1], &dpt[i*nqkn+q1] );
        a1 = a2*a2;
        K2den = a1*a2;
        for ( ii = ci = 0;  ii < n;  ii++, ci += 3 ) {
          a2 = 2.0*bsf1[q1*n+ii];
          for ( jj = 0;  jj < n;  jj++ ) {
            a3 = a2*bsf1[q1*n+jj];
            a = dcp[i+jj].x*a3;  dK2den[ci  ] -= a;  dK2den[ci+3] += a;
            a = dcp[i+jj].y*a3;  dK2den[ci+1] -= a;  dK2den[ci+4] += a;
            a = dcp[i+jj].z*a3;  dK2den[ci+2] -= a;  dK2den[ci+5] += a;
          }
        }
        pkn_MultMatrixNumd ( 1, nia1, 0, dK2den, 3.0*a1, 0, dK2den );
      }
      else {
            /* licznik wyrazenia dla krzywizny */
        _imc_dSpecialK2num ( n, N, bsf, bsf1, q3, q1, i, i, 0, nia,
                             &v31, &dpt[i*nqkn+q1], &K2num, dK2num );
            /* mianownik wyrazenia dla krzywizny */
        _imc_dSpecialK2den ( n, N, &bsf1[q1*n], i, i, 0, nia, &dcp[i], &dpt[i*nqkn+q1],
                             fv31, dfv31, daux1, daux2, &K2den, dK2den );
      }
    }
    else if ( q2 == q3 ) {
            /* licznik wyrazenia dla krzywizny */
      _imc_dSpecialK2num ( n, N, bsf, bsf1, q1, q2, i, i, 0, nia,
                           &v12, &dpt[i*nqkn+q2], &K2num, dK2num );
            /* mianownik wyrazenia dla krzywizny */
      _imc_dSpecialK2den ( n, N, &bsf1[q2*n], i, i, 0, nia, &dcp[i], &dpt[i*nqkn+q2],
                           fv12, dfv12, daux1, daux2, &K2den, dK2den );
    }
    else if ( q1 == q3 ) {
            /* licznik wyrazenia dla krzywizny */
      _imc_dSpecialK2num ( n, N, bsf, bsf1, q2, q3, i, i, 0, nia,
                           &v23, &dpt[i*nqkn+q3], &K2num, dK2num );
            /* mianownik wyrazenia dla krzywizny */
      _imc_dSpecialK2den ( n, N, &bsf1[q3*n], i, i, 0, nia, &dcp[i], &dpt[i*nqkn+q3],
                           fv23, dfv23, daux1, daux2, &K2den, dK2den );
    }
    else {
      pkn_MultMatrixNumd ( 1, nia1, 0, dfv12, 0.5/fl12, 0, dfl12 );
      pkn_MultMatrixNumd ( 1, nia1, 0, dfv23, 0.5/fl23, 0, dfl23 );
      pkn_MultMatrixNumd ( 1, nia1, 0, dfv31, 0.5/fl31, 0, dfl31 );
      half = 0.5*(fl12+fl23+fl31);
      pkn_AddMatrixd ( 1, nia1, 0, dfl12, 0, dfl23, 0, dhalf );
      pkn_AddMatrixd ( 1, nia1, 0, dhalf, 0, dfl31, 0, dhalf );
      pkn_MultMatrixNumd ( 1, nia1, 0, dhalf, 0.5, 0, dhalf );
            /* licznik wyrazenia dla krzywizny */
      _imc_dK2num ( nia1, fl12, fl23, fl31, half,
                    dfl12, dfl23, dfl31, dhalf, daux1, daux2, daux3, daux4,
                    &K2num, dK2num );
            /* mianownik wyrazenia dla krzywizny */
      _imc_dK2den ( nia1, fv12, fv23, fv31, dfv12, dfv23, dfv31,
                    &K2den, dK2den );
    }
          /* kwadrat krzywizny */
    _imc_dRatio ( nia1, K2num, K2den, &K2, dK2num, dK2den, dK2 );
          /* krzywizna w potedze p */
    _imc_dPower ( nia1, K2, dK2, 0.5*md->w, &KtoP, dKtoP );
          /* jakobian */
    memset ( djac, 0, nia*sizeof(double) );
    a1 = ldpt[i*nqkn+q1];
    a2 = ldpt[i*nqkn+q2];
    a3 = ldpt[i*nqkn+q3];
    a1a2 = a1*a2;
    a2a3 = a2*a3;
    a3a1 = a3*a1;
    jac = a1a2*a3;
    ii = 3*(n+1)*(nqkn*i+q1);
    jj = 3*(n+1)*(nqkn*i+q2);
    kk = 3*(n+1)*(nqkn*i+q3);
    for ( c = 0; c < 3*(n+1); c++ ) {
      djac[c] += dldpt[ii+c]*a2a3;
      djac[c] += dldpt[jj+c]*a3a1;
      djac[c] += dldpt[kk+c]*a1a2;
    }
         /* pelna funkcja podcalkowa */
    _imc_dProduct ( nia1, KtoP, jac, &KtoPj, dKtoP, djac, dKtoPj );
    intf5 += KtoPj;
    _imc_CopyGrad ( N, i, i, i, 0, 0, nia1, dKtoPj, g5 );
  }
  i = (q1*nqkn+q2)*nqkn+q3;
  data->_f[i] = 6.0*intf3 + 3.0*intf4 + intf5;
  pkn_AddMatrixMd ( 1, N, 0, g4, 0, g3, 2.0, 0, g3 );
  pkn_AddMatrixMd ( 1, N, 0, g5, 0, g3, 3.0, 0, &data->_g[i*N] );
  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*_intKM_FG*/

boolean _mengerc_gradIntF ( mengerc_data *md, double *func, double *grad )
{
  void           *sp;
  intKM_job_desc data;
  int            lkn, n, ncp, nvcp, N, npoints, nqkn;
  double         *knots, *qkn, *qc;
  point3d        *cpoints;
  vector3d       *dcp, *pt, *dpt, *ddpt;
  double         *ldpt, *dldpt, Nl, Nm, a, b;
  double         *_f, *_g, ff1, ff2, ff3, *gf2, *gf3;
  int            i, j, k, l, ll, m, im, q1;
  int3           jobsize;
  boolean        success;

  sp = pkv_GetScratchMemTop ();
  memset ( &data, 0, sizeof(intKM_job_desc) );
  data.md = md;
  n       = md->deg;
  lkn     = md->lkn;
  knots   = md->knots;
  cpoints = md->cpoints;
  nqkn    = md->nqkn;
  qkn     = md->qkn;
  ncp     = lkn - n;
  nvcp    = ncp - n;
  N = 3*nvcp;
  dcp = data.dcp = pkv_GetScratchMem ( (ncp-1)*sizeof(vector3d) );
  _f = pkv_GetScratchMemd ( nqkn*nqkn*nqkn*(1+N) );
  if ( !dcp || !_f )
    goto failure;
  data._f = _f;
  data._g = _g = &_f[nqkn*nqkn*nqkn];
        /* liczba punktow krzywej */
  npoints = data.npoints = nqkn*(lkn-2*n);
  pt = data.pt = pkv_GetScratchMem ( 3*npoints*sizeof(point3d) +
                                     (3*n+4)*npoints*sizeof(double) );
  if ( !pt )
    goto failure;
  dpt = data.dpt = &pt[npoints];
  ddpt = data.ddpt = &dpt[npoints];
  ldpt = data.ldpt = (double*)&ddpt[npoints];
  dldpt = data.dldpt = &ldpt[npoints];
        /* znajdz roznice kolejnych punktow kontrolnych */
  for ( i = 0; i < ncp-1; i++ )
    SubtractPoints3d ( &cpoints[i+1], &cpoints[i], &dcp[i] );
        /* znajdz punkty krzywej i pochodne pierwszego i drugiego rzedu */
        /* oraz pochodne dlugosci pochodnych parametrryzacji (do jakobianu) */
  memset ( dldpt, 0, 3*(n+1)*npoints*sizeof(double) );
  for ( i = j = 0;  i < lkn-2*n;  i++ )
    for ( q1 = 0;  q1 < nqkn;  q1 ++, j++ ) {
      mbs_deBoorDer2C3d ( n, 2*n+1, knots, &cpoints[i], (double)n+qkn[q1],
                          &pt[j], &dpt[j], &ddpt[j] );
          /* obliczenia pomocnicze do jakobianu: dlugosc wektora pochodnej */
      ldpt[j] = sqrt ( DotProduct3d ( &dpt[j], &dpt[j] ) );
          /* pochodne dlugosci wektora pochodnej wzgledem wspolrzednych */
          /* punktow kontrolnych */
      for ( l = 0; l < n; l++ ) {
        Nl = 2.0*md->bsf1[n*q1+l];
        ll = 3*((n+1)*j+l);
        for ( m = 0; m < n; m++ ) {
          Nm = md->bsf1[n*q1+m];
          a = Nl*Nm;
          im = i+m;
          dldpt[ll  ] -= b = a*dcp[im].x;
          dldpt[ll+3] += b;
          dldpt[ll+1] -= b = a*dcp[im].y;
          dldpt[ll+4] += b;
          dldpt[ll+2] -= b = a*dcp[im].z;
          dldpt[ll+5] += b;
        }
      }
      ll = 3*(n+1)*j;
      pkn_MultMatrixNumd ( 1, 3*(n+1), 0, &dldpt[ll], 0.5/ldpt[j], 0, &dldpt[ll] );
    }
        /* oblicz calki */
  if ( md->npthr > 1 ) {
    jobsize.x = jobsize.y = jobsize.z = nqkn;
    pkv_SetPThreadsToWork ( &jobsize, md->npthr, 1048576, 1048576,
                            (void*)&data, _intKM_FG, NULL, NULL, &success );
    if ( !success )
      goto failure;
  }
  else {
    for ( jobsize.x = 0; jobsize.x < nqkn; jobsize.x++ )
      for ( jobsize.y = 0; jobsize.y < nqkn; jobsize.y++ )
        for ( jobsize.z = 0; jobsize.z < nqkn; jobsize.z++ )
          _intKM_FG ( (void*)&data, &jobsize );
  }
        /* sumuj calki */
  gf2 = pkv_GetScratchMemd ( 2*N );
  if ( !gf2 )
    goto failure;
  gf3 = &gf2[N];
  qc = md->qc;
  ff1 = 0.0;
  memset ( grad, 0, N*sizeof(double) );
  for ( i = l = m = 0;  i < nqkn;  i++ ) {
    ff2 = 0.0;
    memset ( gf2, 0, N*sizeof(double) );
    for ( j = 0; j < nqkn; j++ ) {
      ff3 = 0.0;
      memset ( gf3, 0, N*sizeof(double) );
      for ( k = 0;  k < nqkn;  k++, l++, m += N ) {
        ff3 += qc[k]*_f[l];
        pkn_AddMatrixMd ( 1, N, 0, gf3, 0, &_g[m], qc[k], 0, gf3 );
      }
      ff2 += qc[j]*ff3;
      pkn_AddMatrixMd ( 1, N, 0, gf2, 0, gf3, qc[j], 0, gf2 );
    }
    ff1 += qc[i]*ff2;
    pkn_AddMatrixMd ( 1, N, 0, grad, 0, gf2, qc[i], 0, grad );
  }
  *func = ff1;

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*_mengerc_gradIntF*/

/* ///////////////////////////////////////////////////////////////////////// */
static int _imc_AdjustHes1 ( int n, int N, int i, int k,
                             int ci, int nia, double *ddf )
{
  int z, t, l, m;

  z = i-k+n+1-N;
  if ( z > 0 ) {
    t = nia-z;
    for ( l = 0; l < z; l++ ) {
      for ( m = 0; m < t; m++ )
        ddf[pkn_SymMatIndex(l,m)] += ddf[pkn_SymMatIndex(l+t,m)];
      for ( m = 0; m <= l; m++ )
        ddf[pkn_SymMatIndex(l,m)] += ddf[pkn_SymMatIndex(l+t,m+t)];
      ddf[pkn_SymMatIndex(l,l)] += ddf[pkn_SymMatIndex(l+t,l)];
    }
    return t;
  }
  else
    return nia;
} /*_imc_AdjustHes1*/

static int _imc_AdjustHes2 ( int n, int N, int i, int k,
                             int ci3, int nia3, double *ddf )
{
  int z, t, l, m;

  z = 3*(i-k+n+1-N);
  if ( z > 0 ) {
    t = nia3-z;
    for ( l = 0; l < z; l++ ) {
      for ( m = 0; m < t; m++ )
        ddf[pkn_SymMatIndex(l,m)] += ddf[pkn_SymMatIndex(l+t,m)];
      for ( m = 0; m <= l; m++ )
        ddf[pkn_SymMatIndex(l,m)] += ddf[pkn_SymMatIndex(l+t,m+t)];
      ddf[pkn_SymMatIndex(l,l)] += ddf[pkn_SymMatIndex(l+t,l)];
    }
    return t;
  }
  else
    return nia3;
} /*_imc_AdjustHes2*/

static void _imc_CopyHes2 ( int N, int i, int j, int k,
                            int ci, int cj, int nia,
                            double *ddf, double *h )
{
  int c, ic, d, id;

  i *= 3;  j *= 3;  k *= 3;
  for ( c = 0; c < cj; c++ ) {
    ic = k+c;  if ( ic >= N ) ic -= N;
    for ( d = 0; d <= c; d++ ) {
      id = k+d;  if ( id >= N ) id -= N;
      h[pkn_SymMatIndex(ic,id)] += ddf[pkn_LowerTrMatIndex(c,d)];
    }
  }
  for ( ; c < ci; c++ ) {
    ic = j+c-cj;  if ( ic >= N ) ic -= N;
    for ( d = 0; d < cj; d++ ) {
      id = k+d;  if ( id >= N ) id -= N;
      h[pkn_SymMatIndex(ic,id)] += ddf[pkn_LowerTrMatIndex(c,d)];
    }
    for ( ; d <= c; d++ ) {
      id = j+d-cj;  if ( id >= N ) id -= N;
      h[pkn_SymMatIndex(ic,id)] += ddf[pkn_LowerTrMatIndex(c,d)];
    }
  }
  for ( ; c < nia; c++ ) {
    ic = i+c-ci;  if ( ic >= N ) ic -= N;
    for ( d = 0; d < cj; d++ ) {
      id = k+d;  if ( id >= N ) id -= N;
      h[pkn_SymMatIndex(ic,id)] += ddf[pkn_LowerTrMatIndex(c,d)];
    }
    for ( ; d < ci; d++ ) {
      id = j+d-cj;  if ( id >= N ) id -= N;
      h[pkn_SymMatIndex(ic,id)] += ddf[pkn_LowerTrMatIndex(c,d)];
    }
    for ( ; d <= c; d++ ) {
      id = i+d-ci;  if ( id >= N ) id -= N;
      h[pkn_SymMatIndex(ic,id)] += ddf[pkn_LowerTrMatIndex(c,d)];
    }
  }
} /*_imc_CopyHes2*/

static void _imc_HLab ( int nia, double fvab, double flab, double a1,
                        double *dfvab, double *ddfvab,
                        double *dflab, double *ddflab )
{
  double a2, a3;
  int    c, c3, d, d3, nia3;

  nia3 = 3*nia;
  a2 = -0.25/(flab*fvab);
  for ( c = 0; c < nia3; c++ ) {
    a3 = a2*dfvab[c];
    for ( d = 0; d <= c; d++ )
      ddflab[pkn_LowerTrMatIndex(c,d)] = a3*dfvab[d];
  }
  for ( c = c3 = 0;  c < nia;  c++, c3 += 3 )
    for ( d = d3 = 0;  d <= c;  d++, d3 += 3 ) {
      a3 = a1*ddfvab[pkn_LowerTrMatIndex(c,d)];
      ddflab[pkn_LowerTrMatIndex(c3,d3)]     += a3;
      ddflab[pkn_LowerTrMatIndex(c3+1,d3+1)] += a3;
      ddflab[pkn_LowerTrMatIndex(c3+2,d3+2)] += a3;
    }
} /*_imc_HLab*/

static void _imc_ddProduct ( int n, double a, double b, double *ab,
                             double *da, double *db, double *dab,
                             double *dda, double *ddb, double *ddab )
{
  double ca, cb;
  int    i, j, k, nn;

  *ab = a*b;
  pkn_MatrixLinCombd ( 1, n, 0, da, b, 0, db, a, 0, dab );
  nn = (n*(n+1))/2;
  pkn_MatrixLinCombd ( 1, nn, 0, dda, b, 0, ddb, a, 0, ddab );
  for ( i = k = 0;  i < n;  i++ ) {
    ca = da[i];  cb = db[i];
    for ( j = 0;  j <= i;  j++, k++ )
      ddab[k] += ca*db[j] + cb*da[j];
  }
} /*_imc_ddProduct*/

static void _imc_ddRatio ( int n, double g, double h, double *f,
                           double *dg, double *dh, double *df,
                           double *ddg, double *ddh, double *ddf )
{
  double _f;
  int    np, ii, jj, kk;

  np = (n*(n+1))/2;
  *f = _f = g/h;
  pkn_AddMatrixMd ( 1, n, 0, dg, 0, dh, -_f, 0, df );
  pkn_MultMatrixNumd ( 1, n, 0, df, 1.0/h, 0, df );
  pkn_AddMatrixMd ( 1, np, 0, ddg, 0, ddh, -_f, 0, ddf );
  for ( ii = kk = 0;  ii < n;  ii++ )
    for ( jj = 0;  jj <= ii; jj++, kk++ )
      ddf[kk] -= df[ii]*dh[jj] + df[jj]*dh[ii];
  pkn_MultMatrixNumd ( 1, np, 0, ddf, 1.0/h, 0, ddf );
} /*_imc_ddRatio*/

static void _imc_ddPower ( int n, double a, double *da, double *dda, double p,
                           double *atop, double *datop, double *ddatop )
{
  double _atop, a1;
  int    nn, ii, jj, kk;

  nn = (n*(n+1))/2;
  *atop = _atop = pow ( a, p );
  a1 = p*_atop/a;
  pkn_MultMatrixNumd ( 1, n, 0, da, a1, 0, datop );
  pkn_MultMatrixNumd ( 1, nn, 0, dda, a1, 0, ddatop );
  a1 = (p-1.0)/(p*_atop);
  for ( ii = kk = 0;  ii < n;  ii++ )
    for ( jj = 0;  jj <= ii;  jj++, kk++ )
      ddatop[kk] += a1*datop[ii]*datop[jj];
} /*_imc_ddPower*/

static void _imc_ddK2num ( int nia,
                    double fl12, double fl23, double fl31, double half,
                    double *dfl12, double *dfl23, double *dfl31, double *dhalf,
                    double *ddfl12, double *ddfl23, double *ddfl31, double *ddhalf,
                    double *daux1, double *daux2, double *daux3,
                    double *daux4, double *daux5,
                    double *ddaux1, double *ddaux2, double *ddaux3, double *ddaux4,
                    double *K2num, double *dK2num, double *ddK2num )
{
  double a1, a2, a3, a4, a5;
  int    npia;

  npia = (nia*(nia+1))/2;
  a1 = half-fl12;  a2 = half-fl23;  a3 = half-fl31;
  pkn_SubtractMatrixd ( 1, nia, 0, dhalf, 0, dfl12, 0, daux1 );
  pkn_SubtractMatrixd ( 1, npia, 0, ddhalf, 0, ddfl12, 0, ddaux1 );
  _imc_ddProduct ( nia, half, a1, &a4, dhalf, daux1, daux2, ddhalf, ddaux1, ddaux2 );

  pkn_SubtractMatrixd ( 1, nia, 0, dhalf, 0, dfl23, 0, daux5 );
  pkn_SubtractMatrixd ( 1, npia, 0, ddhalf, 0, ddfl23, 0, ddaux1 );
  pkn_SubtractMatrixd ( 1, nia, 0, dhalf, 0, dfl31, 0, daux3 );
  pkn_SubtractMatrixd ( 1, npia, 0, ddhalf, 0, ddfl31, 0, ddaux3 );
  _imc_ddProduct ( nia, a2, a3, &a5, daux5, daux3, daux4, ddaux1, ddaux3, ddaux4 );

  pkn_MatrixLinCombd ( 1, nia, 0, daux5, a3, 0, daux3, a2, 0, daux4 );
  pkn_MatrixLinCombd ( 1, nia, 0, daux2, a5, 0, daux4, a4, 0, dK2num );
  _imc_ddProduct ( nia, a4, a5, K2num, daux2, daux4, dK2num, ddaux2, ddaux4, ddK2num );
  *K2num *= 16.0;
  pkn_MultMatrixNumd ( 1, nia, 0, dK2num, 16.0, 0, dK2num );
  pkn_MultMatrixNumd ( 1, npia, 0, ddK2num, 16.0, 0, ddK2num );
} /*_imc_ddK2num*/

static void _imc_ddK2den ( int nia3, double fv12, double fv23, double fv31,
                    double *dfv12, double *dfv23, double *dfv31, double *daux,
                    double *ddfv12, double *ddfv23, double *ddfv31, double *ddaux,
                    double *ddaux1, double *ddaux2,
                    double *K2den, double *dK2den, double *ddK2den )
{
  double aux;
  int    nia, npia;
  int    i, j, k;

  nia = nia3/3;
  npia = (nia3*(nia3+1))/2;
  memset ( ddaux1, 0, npia*sizeof(double) );
  for ( i = k = 0;  i < nia;  i++ )
    for ( j = 0;  j <= i;  j++, k++ )
      ddaux1[pkn_LowerTrMatIndex(3*i,3*j)] =
      ddaux1[pkn_LowerTrMatIndex(3*i+1,3*j+1)] =
      ddaux1[pkn_LowerTrMatIndex(3*i+2,3*j+2)] = ddfv12[k];
  memset ( ddaux2, 0, npia*sizeof(double) );
  for ( i = k = 0;  i < nia;  i++ )
    for ( j = 0;  j <= i;  j++, k++ )
      ddaux2[pkn_LowerTrMatIndex(3*i,3*j)] =
      ddaux2[pkn_LowerTrMatIndex(3*i+1,3*j+1)] =
      ddaux2[pkn_LowerTrMatIndex(3*i+2,3*j+2)] = ddfv23[k];
  _imc_ddProduct ( nia3, fv12, fv23, &aux, dfv12, dfv23, daux, ddaux1, ddaux2, ddaux );
  memset ( ddaux1, 0, npia*sizeof(double) );
  for ( i = k = 0;  i < nia;  i++ )
    for ( j = 0;  j <= i;  j++, k++ )
      ddaux1[pkn_LowerTrMatIndex(3*i,3*j)] =
      ddaux1[pkn_LowerTrMatIndex(3*i+1,3*j+1)] =
      ddaux1[pkn_LowerTrMatIndex(3*i+2,3*j+2)] = ddfv31[k];
  _imc_ddProduct ( nia3, aux, fv31, K2den, daux, dfv31, dK2den, ddaux, ddaux1, ddK2den );
} /*_imc_ddK2den*/

static void _imc_ddSpecialK2num ( int n, int N, int nvcp,
                           double *bsf, double *dbsf, double *bsf1,
                           int q1, int q2, int i, int j, int ci, int nia, int npia,
                           vector3d *v12, vector3d *dpt, vector3d *auxv,
                           double *K2num, double *dK2num, double *ddK2num )
{
  double   a1, a2;
  vector3d auxv1;
  int      ii, jj, kk, ic, jc, id, jd;

  memset ( dK2num, 0, nia*sizeof(double) );
  memset ( ddK2num, 0, npia*sizeof(double) );
  memset ( auxv, 0, nia*sizeof(vector3d) );
  CrossProduct3d ( dpt, v12, &auxv1 );
  *K2num = 4.0*DotProduct3d ( &auxv1, &auxv1 );
  if ( i >= j ) {
    for ( ii = ic = 0;  ii < n;  ii++, ic += 3 ) {
      a1 = bsf1[q2*n+ii];
      a2 = v12->z*a1;  auxv[ic  ].y += a2;  auxv[ic+3].y -= a2;
      a2 = v12->y*a1;  auxv[ic  ].z -= a2;  auxv[ic+3].z += a2;
      a2 = (v12->y*auxv1.z - v12->z*auxv1.y)*a1;
      dK2num[ic  ] -= a2;  dK2num[ic+3] += a2;
      a2 = v12->x*a1;  auxv[ic+1].z += a2;  auxv[ic+4].z -= a2;
      a2 = v12->z*a1;  auxv[ic+1].x -= a2;  auxv[ic+4].x += a2;
      a2 = (v12->z*auxv1.x - v12->x*auxv1.z)*a1;
      dK2num[ic+1] -= a2;  dK2num[ic+4] += a2;
      a2 = v12->y*a1;  auxv[ic+2].x += a2;  auxv[ic+5].x -= a2;
      a2 = v12->x*a1;  auxv[ic+2].y -= a2;  auxv[ic+5].y += a2;
      a2 = (v12->x*auxv1.y - v12->y*auxv1.x)*a1;
      dK2num[ic+2] -= a2;  dK2num[ic+5] += a2;
    }
    for ( ii = ic = 0, jc = ci;  ii <= n;  ii++, ic += 3, jc += 3 ) {
      a1 = bsf[q2*(n+1)+ii];
      auxv[ic  ].y -= dpt->z*a1;  auxv[ic  ].z += dpt->y*a1;
      dK2num[ic  ] -= (dpt->z*auxv1.y - dpt->y*auxv1.z)*a1;
      auxv[ic+1].z -= dpt->x*a1;  auxv[ic+1].x += dpt->z*a1;
      dK2num[ic+1] -= (dpt->x*auxv1.z - dpt->z*auxv1.x)*a1;
      auxv[ic+2].x -= dpt->y*a1;  auxv[ic+2].y += dpt->x*a1;
      dK2num[ic+2] -= (dpt->y*auxv1.x - dpt->x*auxv1.y)*a1;
      a1 = bsf[q1*(n+1)+ii];
      auxv[jc  ].y += dpt->z*a1;  auxv[jc  ].z -= dpt->y*a1;
      dK2num[jc  ] += (dpt->z*auxv1.y - dpt->y*auxv1.z)*a1;
      auxv[jc+1].z += dpt->x*a1;  auxv[jc+1].x -= dpt->z*a1;
      dK2num[jc+1] += (dpt->x*auxv1.z - dpt->z*auxv1.x)*a1;
      auxv[jc+2].x += dpt->y*a1;  auxv[jc+2].y -= dpt->x*a1;
      dK2num[jc+2] += (dpt->y*auxv1.x - dpt->x*auxv1.y)*a1;
    }
    for ( ii = kk = 0;  ii < nia;  ii++ ) {
      for ( jj = 0;  jj <= ii;  jj++, kk++ )
        ddK2num[kk] += DotProduct3d ( &auxv[ii], &auxv[jj] );
    }
    for ( ii = ic = 0, id = ci;  ii <= n;  ii++, ic += 3, id += 3 )
      for ( jj = jc = 0, jd = ci;  jj <= n;  jj++, jc += 3, jd += 3 ) {
        if ( ic > jc ) {
          a1 = dbsf[q2*(n+1)+ii]*bsf[q2*(n+1)+jj] - dbsf[q2*(n+1)+jj]*bsf[q2*(n+1)+ii];
          a2 = a1*auxv1.z;
          ddK2num[pkn_LowerTrMatIndex(ic,  jc+1)] -= a2;
          ddK2num[pkn_LowerTrMatIndex(ic+1,jc  )] += a2;
          a2 = a1*auxv1.y;
          ddK2num[pkn_LowerTrMatIndex(ic,  jc+2)] += a2;
          ddK2num[pkn_LowerTrMatIndex(ic+2,jc  )] -= a2;
          a2 = a1*auxv1.x;
          ddK2num[pkn_LowerTrMatIndex(ic+1,jc+2)] -= a2;
          ddK2num[pkn_LowerTrMatIndex(ic+2,jc+1)] += a2;
        }
        if ( ic > jd ) {
          a1 = dbsf[q2*(n+1)+ii]*bsf[q1*(n+1)+jj];
          a2 = a1*auxv1.z;
          ddK2num[pkn_LowerTrMatIndex(ic,  jd+1)] += a2;
          ddK2num[pkn_LowerTrMatIndex(ic+1,jd  )] -= a2;
          a2 = a1*auxv1.y;
          ddK2num[pkn_LowerTrMatIndex(ic,  jd+2)] -= a2;
          ddK2num[pkn_LowerTrMatIndex(ic+2,jd  )] += a2;
          a2 = a1*auxv1.x;
          ddK2num[pkn_LowerTrMatIndex(ic+1,jd+2)] += a2;
          ddK2num[pkn_LowerTrMatIndex(ic+2,jd+1)] -= a2;
        }
        if ( id > jc ) {
          a1 = dbsf[q2*(n+1)+jj]*bsf[q1*(n+1)+ii];
          a2 = a1*auxv1.z;
          ddK2num[pkn_LowerTrMatIndex(id,  jc+1)] -= a2;
          ddK2num[pkn_LowerTrMatIndex(id+1,jc  )] += a2;
          a2 = a1*auxv1.y;
          ddK2num[pkn_LowerTrMatIndex(id,  jc+2)] += a2;
          ddK2num[pkn_LowerTrMatIndex(id+2,jc  )] -= a2;
          a2 = a1*auxv1.x;
          ddK2num[pkn_LowerTrMatIndex(id+1,jc+2)] -= a2;
          ddK2num[pkn_LowerTrMatIndex(id+2,jc+1)] += a2;
        }
      }
    if ( i > j ) {
      _imc_AdjustHes2 ( n, nvcp, i, j, ci, nia, ddK2num );
      nia = _imc_AdjustGrad ( n, N, i, j, ci, nia, dK2num );
      npia = (nia*(nia+1))/2;
    }
  }
  else {
    for ( ii = 0, ic = ci;  ii < n;  ii++, ic += 3 ) {
      a1 = bsf1[q2*n+ii];
      a2 = v12->z*a1;  auxv[ic  ].y += a2;  auxv[ic+3].y -= a2;
      a2 = v12->y*a1;  auxv[ic  ].z -= a2;  auxv[ic+3].z += a2;
      a2 = (v12->y*auxv1.z - v12->z*auxv1.y)*a1;
      dK2num[ic  ] -= a2;  dK2num[ic+3] += a2;
      a2 = v12->x*a1;  auxv[ic+1].z += a2;  auxv[ic+4].z -= a2;
      a2 = v12->z*a1;  auxv[ic+1].x -= a2;  auxv[ic+4].x += a2;
      a2 = (v12->z*auxv1.x - v12->x*auxv1.z)*a1;
      dK2num[ic+1] -= a2;  dK2num[ic+4] += a2;
      a2 = v12->y*a1;  auxv[ic+2].x += a2;  auxv[ic+5].x -= a2;
      a2 = v12->x*a1;  auxv[ic+2].y -= a2;  auxv[ic+5].y += a2;
      a2 = (v12->x*auxv1.y - v12->y*auxv1.x)*a1;
      dK2num[ic+2] -= a2;  dK2num[ic+5] += a2;
    }
    for ( ii = jc = 0, ic = ci;  ii <= n;  ii++, ic += 3, jc += 3 ) {
      a1 = bsf[q2*(n+1)+ii];
      auxv[ic  ].y -= dpt->z*a1;  auxv[ic  ].z += dpt->y*a1;
      dK2num[ic  ] -= (dpt->z*auxv1.y - dpt->y*auxv1.z)*a1;
      auxv[ic+1].z -= dpt->x*a1;  auxv[ic+1].x += dpt->z*a1;
      dK2num[ic+1] -= (dpt->x*auxv1.z - dpt->z*auxv1.x)*a1;
      auxv[ic+2].x -= dpt->y*a1;  auxv[ic+2].y += dpt->x*a1;
      dK2num[ic+2] -= (dpt->y*auxv1.x - dpt->x*auxv1.y)*a1;
      a1 = bsf[q1*(n+1)+ii];
      auxv[jc  ].y += dpt->z*a1;  auxv[jc  ].z -= dpt->y*a1;
      dK2num[jc  ] += (dpt->z*auxv1.y - dpt->y*auxv1.z)*a1;
      auxv[jc+1].z += dpt->x*a1;  auxv[jc+1].x -= dpt->z*a1;
      dK2num[jc+1] += (dpt->x*auxv1.z - dpt->z*auxv1.x)*a1;
      auxv[jc+2].x += dpt->y*a1;  auxv[jc+2].y -= dpt->x*a1;
      dK2num[jc+2] += (dpt->y*auxv1.x - dpt->x*auxv1.y)*a1;
    }
    for ( ii = kk = 0;  ii < nia;  ii++ ) {
      for ( jj = 0;  jj <= ii;  jj++, kk++ )
        ddK2num[kk] += DotProduct3d ( &auxv[ii], &auxv[jj] );
    }
    for ( ii = id = 0, ic = ci;  ii <= n;  ii++, ic += 3, id += 3 )
      for ( jj = jd = 0, jc = ci;  jj <= n;  jj++, jc += 3, jd += 3 ) {
        if ( ic > jc ) {
          a1 = dbsf[q2*(n+1)+ii]*bsf[q2*(n+1)+jj] - dbsf[q2*(n+1)+jj]*bsf[q2*(n+1)+ii];
          a2 = a1*auxv1.z;
          ddK2num[pkn_LowerTrMatIndex(ic,  jc+1)] -= a2;
          ddK2num[pkn_LowerTrMatIndex(ic+1,jc  )] += a2;
          a2 = a1*auxv1.y;
          ddK2num[pkn_LowerTrMatIndex(ic,  jc+2)] += a2;
          ddK2num[pkn_LowerTrMatIndex(ic+2,jc  )] -= a2;
          a2 = a1*auxv1.x;
          ddK2num[pkn_LowerTrMatIndex(ic+1,jc+2)] -= a2;
          ddK2num[pkn_LowerTrMatIndex(ic+2,jc+1)] += a2;
        }
        if ( ic > jd ) {
          a1 = dbsf[q2*(n+1)+ii]*bsf[q1*(n+1)+jj];
          a2 = a1*auxv1.z;
          ddK2num[pkn_LowerTrMatIndex(ic,  jd+1)] += a2;
          ddK2num[pkn_LowerTrMatIndex(ic+1,jd  )] -= a2;
          a2 = a1*auxv1.y;
          ddK2num[pkn_LowerTrMatIndex(ic,  jd+2)] -= a2;
          ddK2num[pkn_LowerTrMatIndex(ic+2,jd  )] += a2;
          a2 = a1*auxv1.x;
          ddK2num[pkn_LowerTrMatIndex(ic+1,jd+2)] += a2;
          ddK2num[pkn_LowerTrMatIndex(ic+2,jd+1)] -= a2;
        }
        if ( id > jc ) {
          a1 = dbsf[q2*(n+1)+jj]*bsf[q1*(n+1)+ii];
          a2 = a1*auxv1.z;
          ddK2num[pkn_LowerTrMatIndex(id,  jc+1)] -= a2;
          ddK2num[pkn_LowerTrMatIndex(id+1,jc  )] += a2;
          a2 = a1*auxv1.y;
          ddK2num[pkn_LowerTrMatIndex(id,  jc+2)] += a2;
          ddK2num[pkn_LowerTrMatIndex(id+2,jc  )] -= a2;
          a2 = a1*auxv1.x;
          ddK2num[pkn_LowerTrMatIndex(id+1,jc+2)] -= a2;
          ddK2num[pkn_LowerTrMatIndex(id+2,jc+1)] += a2;
        }
      }
    _imc_AdjustHes2 ( n, nvcp, j, i, ci, nia, ddK2num );
    nia = _imc_AdjustGrad ( n, N, j, i, ci, nia, dK2num );
    npia = (nia*(nia+1))/2;
  }
  pkn_MultMatrixNumd ( 1, nia, 0, dK2num, 8.0, 0, dK2num );
  pkn_MultMatrixNumd ( 1, npia, 0, ddK2num, 8.0, 0, ddK2num );
} /*_imc_ddSpecialK2num*/

static void _imc_ddSpecialK2den ( int n, int N, int nvcp, double *dbsf, double *bsf1,
                           int i, int j, int ci, int nia3,
                           vector3d *dcp, vector3d *dpt, double fv12,
                           double *dfv12, double *daux1, double *daux2,
                           double *ddfv12, double *ddaux1, double *ddaux2,
                           double *K2den, double *dK2den, double *ddK2den )
{
  double a1, a2, a3, a4;
  int    nia, npia3, ii, jj, kk, ic, jc;

  nia = nia3/3;
  npia3 = (nia3*(nia3+1))/2;
  memset ( daux1, 0, nia3*sizeof(double) );
  a1 = DotProduct3d ( dpt, dpt );
  for ( ii = 0, ic = ci;  ii < n;  ii++, ic += 3 ) {
    a2 = 2.0*bsf1[ii];
    for ( jj = 0;  jj < n;  jj++ ) {
      a3 = a2*bsf1[jj];
      a4 = dcp[jj].x*a3;  daux1[ic  ] -= a4;  daux1[ic+3] += a4;
      a4 = dcp[jj].y*a3;  daux1[ic+1] -= a4;  daux1[ic+4] += a4;
      a4 = dcp[jj].z*a3;  daux1[ic+2] -= a4;  daux1[ic+5] += a4;
    }
  }
  memset ( ddaux1, 0, npia3*sizeof(double) );
  for ( ii = 0, ic = ci;  ii <= n;  ii++, ic += 3 ) {
    a2 = 2.0*dbsf[ii];
    for ( jj = 0, jc = ci;  jj <= ii;  jj++, jc += 3 ) {
      ddaux1[pkn_LowerTrMatIndex(ic,jc)] =
      ddaux1[pkn_LowerTrMatIndex(ic+1,jc+1)] =
      ddaux1[pkn_LowerTrMatIndex(ic+2,jc+2)] = a2*dbsf[jj];
    }
  }
  memset ( ddaux2, 0, npia3*sizeof(double) );
  a3 = 2.0*fv12;
  for ( ii = kk = 0, ic = 0;  ii <= nia;  ii++, ic += 3 ) {
    for ( jj = 0, jc = 0;  jj <= ii;  jj++, jc += 3, kk++ ) {
      ddaux2[pkn_LowerTrMatIndex(ic,jc)] =
      ddaux2[pkn_LowerTrMatIndex(ic+1,jc+1)] =
      ddaux2[pkn_LowerTrMatIndex(ic+2,jc+2)] = a3*ddfv12[kk];
    }
  }
  if ( i < j ) {
    _imc_AdjustHes2 ( n, nvcp, j, i, ci, nia3, ddaux1 );
/* nie trzeba adjustowac ddaux2, bo tablica ddfv12 juz byla zadjustowana */
    nia3 = _imc_AdjustGrad ( n, N, j, i, ci, nia3, daux1 );
  }
  pkn_MultMatrixNumd ( 1, nia3, 0, dfv12, a3, 0, daux2 );
  a3 = fv12*fv12;
  a2 = 0.5/a3;
  for ( ii = kk = 0;  ii < nia3;  ii++ )
    for ( jj = 0;  jj <= ii;  jj++, kk++ )
      ddaux2[kk] += a2*daux2[ii]*daux2[jj];
  _imc_ddProduct ( nia3, a1, a3, K2den, daux1, daux2, dK2den,
                   ddaux1, ddaux2, ddK2den );
} /*_imc_ddSpecialK2den*/

static boolean _intKM_FGH ( void *usrdata, int3 *jobnum )
{
  void           *sp;
  intKM_job_desc *data;
  mengerc_data       *md;
  int            n, lkn, nqkn, ncp, nvcp, N, M;
  vector3d       *dcp, *pt, *dpt, *ddpt;
  double         *ldpt, *dldpt, *ddldpt;
  int            q1, q2, q3, i, j, k, c, c3, d, ic, jc,
                 ii, jj, kk, iii, jjj, kkk;
  int            maxnia, maxnia3, maxnpia, maxnpia3;
  double         *bsf, *dbsf, *bsf1, *dbsf1;
  point3d        pt1, pt2, pt3;
  vector3d       v12, v23, v31, auxv1, *auxv;
  int            ti, tj, tk, nia, nia1, nia3, nia31, npia, npia3, npia31,
                 ci, ci3, cj, cj3, cc, dd;
  double         a, a1, a2, a3, a1a2, a2a3, a3a1;
  double         fv12, fv23, fv31, fl12, fl23, fl31,
                 half, K2num, K2den, K2, KtoP, jac, KtoPj;
  double         *dfv12, *dfv23, *dfv31, *dfl12, *dfl23, *dfl31,
                 *dhalf, *dK2num, *dK2den, *dK2, *dKtoP, *djac, *dKtoPj,
                 *daux1, *daux2, *daux3, *daux4, *daux5;
  double         *ddfv12, *ddfv23, *ddfv31, *ddfl12, *ddfl23, *ddfl31,
                 *ddhalf, *ddK2num, *ddK2den, *ddK2, *ddKtoP, *ddjac, *ddKtoPj,
                 *ddaux1, *ddaux2, *ddaux3, *ddaux4;
  double         intf3, intf4, intf5;
  double         *g3, *g4, *g5, *h3, *h4, *h5;

  sp = pkv_GetScratchMemTop ();
  data = (intKM_job_desc*)usrdata;
  md     = data->md;
  pt     = data->pt;
  dpt    = data->dpt;
  ddpt   = data->ddpt;
  ldpt   = data->ldpt;
  dldpt  = data->dldpt;
  ddldpt = data->ddldpt;
  dcp    = data->dcp;

  n     = md->deg;
  lkn   = md->lkn;
  ncp   = lkn-n;
  nvcp  = ncp-n;
  N     = 3*nvcp;
  M     = (N*(N+1))/2;
  nqkn  = md->nqkn;
  bsf   = md->bsf;
  dbsf  = md->dbsf;
  bsf1  = md->bsf1;
  dbsf1 = md->dbsf1;

  q1 = jobnum->x;
  q2 = jobnum->y;
  q3 = jobnum->z;

  g3 = pkv_GetScratchMemd ( 3*(N+M) );
  if ( !g3 )
    goto failure;
  g4 = &g3[N];  g5 = &g4[N];
  h3 = &g5[N];  h4 = &h3[M];  h5 = &h4[M];
        /* maximal number of relevant variables in any unit cube */
  maxnia = 3*(n+1);
  maxnia3 = 3*maxnia;
  maxnpia = (maxnia*(maxnia+1))/2;
  maxnpia3 = (maxnia3*(maxnia3+1))/2;
  dfv12 = pkv_GetScratchMemd ( 16*maxnia3 + 3*maxnpia + 11*maxnpia3 );
  if ( !dfv12 )
    goto failure;
  dfv23 = &dfv12[maxnia3];  dfv31 = &dfv23[maxnia3];
  dfl12 = &dfv31[maxnia3];  dfl23 = &dfl12[maxnia3];  dfl31 = &dfl23[maxnia3];
  dhalf = &dfl31[maxnia3];  dK2num = &dhalf[maxnia3];  dK2den = &dK2num[maxnia3];
  dK2 = &dK2den[maxnia3];  dKtoP = &dK2[maxnia3];  djac = &dKtoP[maxnia3];
  dKtoPj = &djac[maxnia3];
  daux1 = dKtoPj;  daux2 = &daux1[maxnia3];  daux3 = &daux2[maxnia3];
  daux4 = &daux3[maxnia3];  daux5 = djac;
  ddfv12 = &daux4[maxnia3];  ddfv23 = &ddfv12[maxnpia];  ddfv31 = &ddfv23[maxnpia];
  ddfl12 = &ddfv31[maxnpia];  ddfl23 = &ddfl12[maxnpia3];  ddfl31 = &ddfl23[maxnpia3];
  ddhalf = &ddfl31[maxnpia3];  ddK2num = &ddhalf[maxnpia3];
  ddK2den = &ddK2num[maxnpia3];  ddK2 = &ddK2den[maxnpia3];
  ddKtoP = &ddK2[maxnpia3];  ddjac = &ddKtoP[maxnpia3];
  ddKtoPj = &ddjac[maxnpia3];
  ddaux1 = ddKtoPj;  ddaux2 = &ddaux1[maxnpia3];  ddaux3 = &daux2[maxnpia3];
  ddaux4 = ddjac;  auxv = (vector3d*)ddjac;
        /* compute the integrals */
  intf3 = 0.0;
  memset ( g3, 0, N*sizeof(double) );
  memset ( h3, 0, M*sizeof(double) );
          /* the three loops iterate over the cubes (i,j,k) for i>j>k */
  for ( i = 2; i < lkn-2*n; i++ ) {
    ti = 6*(n+1);
    pt1 = pt[i*nqkn+q1];
    for ( j = 1; j < i; j++ ) {
      tj = ti - 3*min(i-j,n+1);
      pt2 = pt[j*nqkn+q2];
      SubtractPoints3d ( &pt1, &pt2, &v12 );
      for ( k = 0; k < j; k++ ) {
        tk = tj - 3*min(j-k,n+1);
        nia3 = maxnia3-tk;  ci3 = ti-tk;  cj3 = tj-tk;
        nia = nia3/3;  ci = ci3/3;  cj = cj3/3;
        npia = (nia*(nia+1))/2;
        npia3 = (nia3*(nia3+1))/2;
        memset ( dfv12, 0, nia3*sizeof(double) );
        memset ( dfv23, 0, nia3*sizeof(double) );
        memset ( dfv31, 0, nia3*sizeof(double) );
        memset ( ddfv12, 0, npia*sizeof(double) );
        memset ( ddfv23, 0, npia*sizeof(double) );
        memset ( ddfv31, 0, npia*sizeof(double) );

        pt3 = pt[k*nqkn+q3];
        SubtractPoints3d ( &pt2, &pt3, &v23 );
        SubtractPoints3d ( &pt3, &pt1, &v31 );
        fv12 = DotProduct3d ( &v12, &v12 );  fl12 = sqrt ( fv12 );
        fv23 = DotProduct3d ( &v23, &v23 );  fl23 = sqrt ( fv23 );
        fv31 = DotProduct3d ( &v31, &v31 );  fl31 = sqrt ( fv31 );

        for ( c = c3 = 0;  c <= n;  c++, c3 += 3 ) {
          a1 = 2.0*bsf[(n+1)*q1+c];
          a2 = 2.0*bsf[(n+1)*q2+c];
          a3 = 2.0*bsf[(n+1)*q3+c];
          dfv12[ci3+c3]   += a1*v12.x;
          dfv12[ci3+c3+1] += a1*v12.y;
          dfv12[ci3+c3+2] += a1*v12.z;
          dfv12[cj3+c3]   -= a2*v12.x;
          dfv12[cj3+c3+1] -= a2*v12.y;
          dfv12[cj3+c3+2] -= a2*v12.z;
          dfv23[cj3+c3]   += a2*v23.x;
          dfv23[cj3+c3+1] += a2*v23.y;
          dfv23[cj3+c3+2] += a2*v23.z;
          dfv23[c3]   -= a3*v23.x;
          dfv23[c3+1] -= a3*v23.y;
          dfv23[c3+2] -= a3*v23.z;
          dfv31[c3]   += a3*v31.x;
          dfv31[c3+1] += a3*v31.y;
          dfv31[c3+2] += a3*v31.z;
          dfv31[ci3+c3]   -= a1*v31.x;
          dfv31[ci3+c3+1] -= a1*v31.y;
          dfv31[ci3+c3+2] -= a1*v31.z;
          for ( d = 0; d <= n; d++ ) {
            if ( c >= d ) {
              ddfv12[pkn_LowerTrMatIndex(ci+c,ci+d)] += a1*bsf[(n+1)*q1+d];
              ddfv12[pkn_LowerTrMatIndex(cj+c,cj+d)] += a2*bsf[(n+1)*q2+d];
              ddfv23[pkn_LowerTrMatIndex(cj+c,cj+d)] += a2*bsf[(n+1)*q2+d];
              ddfv23[pkn_LowerTrMatIndex(c,d)] += a3*bsf[(n+1)*q3+d];
              ddfv31[pkn_LowerTrMatIndex(c,d)] += a3*bsf[(n+1)*q3+d];
              ddfv31[pkn_LowerTrMatIndex(ci+c,ci+d)] += a1*bsf[(n+1)*q1+d];
            }
            if ( ci+c >= cj+d )
              ddfv12[pkn_LowerTrMatIndex(ci+c,cj+d)] -= a1*bsf[(n+1)*q2+d];
            if ( ci+c >= d )
              ddfv31[pkn_LowerTrMatIndex(ci+c,d)] -= a1*bsf[(n+1)*q3+d];
            if ( cj+c >= ci+d )
              ddfv12[pkn_LowerTrMatIndex(cj+c,ci+d)] -= a2*bsf[(n+1)*q1+d];
            if ( cj+c >= d )
              ddfv23[pkn_LowerTrMatIndex(cj+c,d)] -= a2*bsf[(n+1)*q3+d];
            if ( c >= cj+d )
              ddfv23[pkn_LowerTrMatIndex(c,cj+d)] -= a3*bsf[(n+1)*q2+d];
            if ( c >= ci+d )
              ddfv31[pkn_LowerTrMatIndex(c,ci+d)] -= a3*bsf[(n+1)*q1+d];
          }
        }
        nia31 = _imc_AdjustGrad ( n, N, i, k, ci3, nia3, dfv12 );
        _imc_AdjustGrad ( n, N, i, k, ci3, nia3, dfv23 );
        _imc_AdjustGrad ( n, N, i, k, ci3, nia3, dfv31 );
        nia1 = _imc_AdjustHes1 ( n, nvcp, i, k, ci, nia, ddfv12 );
        _imc_AdjustHes1 ( n, nvcp, i, k, ci, nia, ddfv23 );
        _imc_AdjustHes1 ( n, nvcp, i, k, ci, nia, ddfv31 );
        npia31 = (nia31*(nia31+1))/2;

        a1 = 0.5/fl12;
        pkn_MultMatrixNumd ( 1, nia31, 0, dfv12, a1, 0, dfl12 );
        _imc_HLab ( nia1, fv12, fl12, a1, dfv12, ddfv12, dfl12, ddfl12 );
        a1 = 0.5/fl23;
        pkn_MultMatrixNumd ( 1, nia31, 0, dfv23, a1, 0, dfl23 );
        _imc_HLab ( nia1, fv23, fl23, a1, dfv23, ddfv23, dfl23, ddfl23 );
        a1 = 0.5/fl31;
        pkn_MultMatrixNumd ( 1, nia31, 0, dfv31, a1, 0, dfl31 );
        _imc_HLab ( nia1, fv31, fl31, a1, dfv31, ddfv31, dfl31, ddfl31 );

        half = 0.5*(fl12+fl23+fl31);
        pkn_AddMatrixd ( 1, nia31, 0, dfl12, 0, dfl23, 0, dhalf );
        pkn_AddMatrixd ( 1, nia31, 0, dhalf, 0, dfl31, 0, dhalf );
        pkn_MultMatrixNumd ( 1, nia31, 0, dhalf, 0.5, 0, dhalf );
        pkn_AddMatrixd ( 1, npia31, 0, ddfl12, 0, ddfl23, 0, ddhalf );
        pkn_AddMatrixd ( 1, npia31, 0, ddhalf, 0, ddfl31, 0, ddhalf );
        pkn_MultMatrixNumd ( 1, npia31, 0, ddhalf, 0.5, 0, ddhalf );
            /* licznik wyrazenia dla krzywizny */
        _imc_ddK2num ( nia31, fl12, fl23, fl31, half,
                       dfl12, dfl23, dfl31, dhalf,
                       ddfl12, ddfl23, ddfl31, ddhalf,
                       daux1, daux2, daux3, daux4, daux5,
                       ddaux1, ddaux2, ddaux3, ddaux4,
                       &K2num, dK2num, ddK2num );
            /* mianownik wyrazenia dla krzywizny */
        _imc_ddK2den ( nia31, fv12, fv23, fv31, dfv12, dfv23, dfv31, daux1,
                       ddfv12, ddfv23, ddfv31, ddaux1, ddaux2, ddaux3,
                       &K2den, dK2den, ddK2den );
            /* kwadrat krzywizny */
        _imc_ddRatio ( nia31, K2num, K2den, &K2, dK2num, dK2den, dK2,
                       ddK2num, ddK2den, ddK2 );
            /* krzywizna w potedze p */
        _imc_ddPower ( nia31, K2, dK2, ddK2, 0.5*md->w, &KtoP, dKtoP, ddKtoP );
            /* jakobian */
        memset ( djac, 0, nia3*sizeof(double) );
        memset ( ddjac, 0, npia3*sizeof(double) );
        a1 = ldpt[i*nqkn+q1];  a2 = ldpt[j*nqkn+q2];  a3 = ldpt[k*nqkn+q3];
        a1a2 = a1*a2;  a2a3 = a2*a3;  a3a1 = a3*a1;
        jac = a1a2*a3;
        kk = 3*(n+1);
        ii = kk*(nqkn*i+q1);  jj = kk*(nqkn*j+q2);  kk *= nqkn*k+q3;
        kkk = ((3*n+3)*(3*n+4))/2;
        iii = kkk*(nqkn*i+q1);  jjj = kkk*(nqkn*j+q2);  kkk *= nqkn*k+q3;
        for ( c = 0; c < 3*(n+1); c++ ) {
          djac[ci3+c] += dldpt[ii+c]*a2a3;
          djac[cj3+c] += dldpt[jj+c]*a3a1;
          djac[c]     += dldpt[kk+c]*a1a2;
          for ( d = 0; d < 3*(n+1); d++ ) {
            if ( c >= d ) {
              ddjac[pkn_LowerTrMatIndex(ci3+c,ci3+d)] +=
                  ddldpt[iii+pkn_LowerTrMatIndex(c,d)]*a2a3;
              ddjac[pkn_LowerTrMatIndex(cj3+c,cj3+d)] +=
                  ddldpt[jjj+pkn_LowerTrMatIndex(c,d)]*a3a1;
              ddjac[pkn_LowerTrMatIndex(c,d)] +=
                  ddldpt[kkk+pkn_LowerTrMatIndex(c,d)]*a1a2;
            }
            if ( ci3+c >= cj3+d )
              ddjac[pkn_LowerTrMatIndex(ci3+c,cj3+d)] += dldpt[ii+c]*dldpt[jj+d]*a3;
            cc = c + 3*(i-j);
            dd = d - 3*(i-j);
            if ( dd >= 0 && cc < 3*(n+1) && cj3+cc >= ci3+dd )
              ddjac[pkn_LowerTrMatIndex(cj3+cc,ci3+dd)] += dldpt[ii+dd]*dldpt[jj+cc]*a3;
            if ( cj3+c >= d )
              ddjac[pkn_LowerTrMatIndex(cj3+c,d)] += dldpt[jj+c]*dldpt[kk+d]*a1;
            cc = c + 3*(j-k);
            dd = d - 3*(j-k);
            if ( dd >= 0 && cc < 3*(n+1) && cc >= cj3+dd )
              ddjac[pkn_LowerTrMatIndex(cc,cj3+dd)] += dldpt[jj+dd]*dldpt[kk+cc]*a1;

            if ( ci3+c >= d )
              ddjac[pkn_LowerTrMatIndex(ci3+c,d)] += dldpt[ii+c]*dldpt[kk+d]*a2;
            cc = c + 3*(i-k);
            dd = d - 3*(i-k);
            if ( dd >= 0 && cc < 3*(n+1) && cc >= ci3+dd )
              ddjac[pkn_LowerTrMatIndex(cc,ci3+dd)] += dldpt[ii+dd]*dldpt[kk+cc]*a2;
          }
        }
        _imc_AdjustGrad ( n, N, i, k, ci3, nia3, djac );
        _imc_AdjustHes2 ( n, nvcp, i, k, ci3, nia3, ddjac );
         /* pelna funkcja podcalkowa */
        _imc_ddProduct ( nia31, KtoP, jac, &KtoPj,
                         dKtoP, djac, dKtoPj, ddKtoP, ddjac, ddKtoPj );
        intf3 += KtoPj;
        _imc_CopyGrad ( N, i, j, k, ci3, cj3, nia31, dKtoPj, g3 );
        _imc_CopyHes2 ( N, i, j, k, ci3, cj3, nia31, ddKtoPj, h3 );
      }
    }
  }
          /* the two loops iterate over the cubes (i,j,j) for i  >j */
  intf4 = 0.0;
  memset ( g4, 0, N*sizeof(double) );
  memset ( h4, 0, M*sizeof(double) );
  for ( i = 1; i < lkn-2*n; i++ ) {
    ti = 6*(n+1);
    pt1 = pt[i*nqkn+q1];
    for ( j = 0; j < i; j++ ) {
      tj = ti - 3*min(i-j,n+1);
      nia3 = maxnia3-tj;  ci3 = ti-tj;
      nia = nia3/3;  ci = ci3/3;
      npia = (nia*(nia+1))/2;
      npia3 = (nia3*(nia3+1))/2;
      memset ( dfv12, 0, nia3*sizeof(double) );
      memset ( dfv23, 0, nia3*sizeof(double) );
      memset ( dfv31, 0, nia3*sizeof(double) );
      memset ( ddfv12, 0, npia*sizeof(double) );
      memset ( ddfv23, 0, npia*sizeof(double) );
      memset ( ddfv31, 0, npia*sizeof(double) );

      pt2 = pt[j*nqkn+q2];
      SubtractPoints3d ( &pt1, &pt2, &v12 );
      pt3 = pt[j*nqkn+q3];
      SubtractPoints3d ( &pt2, &pt3, &v23 );
      SubtractPoints3d ( &pt3, &pt1, &v31 );
      fv12 = DotProduct3d ( &v12, &v12 );  fl12 = sqrt ( fv12 );
      fv23 = DotProduct3d ( &v23, &v23 );  fl23 = sqrt ( fv23 );
      fv31 = DotProduct3d ( &v31, &v31 );  fl31 = sqrt ( fv31 );

      for ( c = c3 = 0;  c <= n;  c++, c3 += 3 ) {
        a1 = 2.0*bsf[(n+1)*q1+c];
        a2 = 2.0*bsf[(n+1)*q2+c];
        a3 = 2.0*bsf[(n+1)*q3+c];
        dfv12[ci3+c3]   += a1*v12.x;
        dfv12[ci3+c3+1] += a1*v12.y;
        dfv12[ci3+c3+2] += a1*v12.z;
        dfv12[c3]   -= a2*v12.x;
        dfv12[c3+1] -= a2*v12.y;
        dfv12[c3+2] -= a2*v12.z;
        dfv23[c3]   += a2*v23.x;
        dfv23[c3+1] += a2*v23.y;
        dfv23[c3+2] += a2*v23.z;
        dfv23[c3]   -= a3*v23.x;
        dfv23[c3+1] -= a3*v23.y;
        dfv23[c3+2] -= a3*v23.z;
        dfv31[c3]   += a3*v31.x;
        dfv31[c3+1] += a3*v31.y;
        dfv31[c3+2] += a3*v31.z;
        dfv31[ci3+c3]   -= a1*v31.x;
        dfv31[ci3+c3+1] -= a1*v31.y;
        dfv31[ci3+c3+2] -= a1*v31.z;
        for ( d = 0; d <= n; d++ ) {
          if ( c >= d ) {
            ddfv12[pkn_LowerTrMatIndex(ci+c,ci+d)] += a1*bsf[(n+1)*q1+d];
            ddfv12[pkn_LowerTrMatIndex(c,d)] += a2*bsf[(n+1)*q2+d];
            ddfv23[pkn_LowerTrMatIndex(c,d)] +=
                    a2*(bsf[(n+1)*q2+d]-bsf[(n+1)*q3+d])
                    -a3*(bsf[(n+1)*q2+d]-bsf[(n+1)*q3+d]);
            ddfv31[pkn_LowerTrMatIndex(c,d)] += a3*bsf[(n+1)*q3+d];
            ddfv31[pkn_LowerTrMatIndex(ci+c,ci+d)] += a1*bsf[(n+1)*q1+d];
          }
          if ( ci+c >= d ) {
            ddfv12[pkn_LowerTrMatIndex(ci+c,d)] -= a1*bsf[(n+1)*q2+d];
            ddfv31[pkn_LowerTrMatIndex(ci+c,d)] -= a1*bsf[(n+1)*q3+d];
          }
          if ( c >= ci+d ) {
            ddfv12[pkn_LowerTrMatIndex(c,ci+d)] -= a2*bsf[(n+1)*q1+d];
            ddfv31[pkn_LowerTrMatIndex(c,ci+d)] -= a3*bsf[(n+1)*q1+d];
          }
        }
      }
      nia31 = _imc_AdjustGrad ( n, N, i, j, ci3, nia3, dfv12 );
      _imc_AdjustGrad ( n, N, i, j, ci3, nia3, dfv23 );
      _imc_AdjustGrad ( n, N, i, j, ci3, nia3, dfv31 );
      nia1 = _imc_AdjustHes1 ( n, nvcp, i, j, ci, nia, ddfv12 );
      _imc_AdjustHes1 ( n, nvcp, i, j, ci, nia, ddfv23 );
      _imc_AdjustHes1 ( n, nvcp, i, j, ci, nia, ddfv31 );
      npia31 = (nia31*(nia31+1))/2;

      a1 = 0.5/fl12;
      pkn_MultMatrixNumd ( 1, nia31, 0, dfv12, a1, 0, dfl12 );
      _imc_HLab ( nia1, fv12, fl12, a1, dfv12, ddfv12, dfl12, ddfl12 );
      if ( q2 == q3 ) {
          /* licznik wyrazenia dla krzywizny */
        _imc_ddSpecialK2num ( n, N, nvcp, bsf, dbsf, bsf1, q1, q2,
                              i, j, ci3, nia3, npia3, &v12, &dpt[j*nqkn+q2],
                              auxv, &K2num, dK2num, ddK2num );
          /* mianownik wyrazenia dla krzywizny */
        _imc_ddSpecialK2den ( n, N, nvcp, &dbsf[q2*(n+1)], &bsf1[q2*n],
                              i, j, 0, nia3, &dcp[j], &dpt[j*nqkn+q2], fv12,
                              dfv12, daux1, daux2, ddfv12, ddaux1, ddaux2,
                              &K2den, dK2den, ddK2den );
      }
      else {
        a1 = 0.5/fl23;
        pkn_MultMatrixNumd ( 1, nia31, 0, dfv23, a1, 0, dfl23 );
        _imc_HLab ( nia1, fv23, fl23, a1, dfv23, ddfv23, dfl23, ddfl23 );
        a1 = 0.5/fl31;
        pkn_MultMatrixNumd ( 1, nia31, 0, dfv31, a1, 0, dfl31 );
        _imc_HLab ( nia1, fv31, fl31, a1, dfv31, ddfv31, dfl31, ddfl31 );
        pkn_MultMatrixNumd ( 1, nia31, 0, dfv23, 0.5/fl23, 0, dfl23 );
        pkn_MultMatrixNumd ( 1, nia31, 0, dfv31, 0.5/fl31, 0, dfl31 );
        half = 0.5*(fl12+fl23+fl31);
        pkn_AddMatrixd ( 1, nia31, 0, dfl12, 0, dfl23, 0, dhalf );
        pkn_AddMatrixd ( 1, nia31, 0, dhalf, 0, dfl31, 0, dhalf );
        pkn_MultMatrixNumd ( 1, nia31, 0, dhalf, 0.5, 0, dhalf );
        pkn_AddMatrixd ( 1, npia31, 0, ddfl12, 0, ddfl23, 0, ddhalf );
        pkn_AddMatrixd ( 1, npia31, 0, ddhalf, 0, ddfl31, 0, ddhalf );
        pkn_MultMatrixNumd ( 1, npia31, 0, ddhalf, 0.5, 0, ddhalf );
            /* licznik wyrazenia dla krzywizny */
        _imc_ddK2num ( nia31, fl12, fl23, fl31, half,
                       dfl12, dfl23, dfl31, dhalf,
                       ddfl12, ddfl23, ddfl31, ddhalf,
                       daux1, daux2, daux3, daux4, daux5,
                       ddaux1, ddaux2, ddaux3, ddaux4,
                       &K2num, dK2num, ddK2num );
            /* mianownik wyrazenia dla krzywizny */
        _imc_ddK2den ( nia31, fv12, fv23, fv31, dfv12, dfv23, dfv31, daux1,
                       ddfv12, ddfv23, ddfv31, ddaux1, ddaux2, ddaux3,
                       &K2den, dK2den, ddK2den );
      }
          /* kwadrat krzywizny */
      _imc_ddRatio ( nia31, K2num, K2den, &K2, dK2num, dK2den, dK2,
                     ddK2num, ddK2den, ddK2 );
          /* krzywizna w potedze p */
      _imc_ddPower ( nia31, K2, dK2, ddK2, 0.5*md->w, &KtoP, dKtoP, ddKtoP );
          /* jakobian */
      memset ( djac, 0, nia3*sizeof(double) );
      memset ( ddjac, 0, npia3*sizeof(double) );
      a1 = ldpt[i*nqkn+q1];  a2 = ldpt[j*nqkn+q2];  a3 = ldpt[j*nqkn+q3];
      a1a2 = a1*a2;  a2a3 = a2*a3;  a3a1 = a3*a1;
      jac = a1a2*a3;
      kk = 3*(n+1);
      ii = kk*(nqkn*i+q1);  jj = kk*(nqkn*j+q2);  kk *= nqkn*j+q3;
      kkk = ((3*n+3)*(3*n+4))/2;
      iii = kkk*(nqkn*i+q1);  jjj = kkk*(nqkn*j+q2);  kkk *= nqkn*j+q3;
      for ( c = 0; c < 3*(n+1); c++ ) {
        djac[ci3+c] += dldpt[ii+c]*a2a3;
        djac[c]     += dldpt[jj+c]*a3a1;
        djac[c]     += dldpt[kk+c]*a1a2;
        for ( d = 0; d < 3*(n+1); d++ ) {
          if ( c >= d ) {
            ddjac[pkn_LowerTrMatIndex(ci3+c,ci3+d)] +=
                ddldpt[iii+pkn_LowerTrMatIndex(c,d)]*a2a3;
            ddjac[pkn_LowerTrMatIndex(c,d)] +=
                ddldpt[jjj+pkn_LowerTrMatIndex(c,d)]*a3a1;
            ddjac[pkn_LowerTrMatIndex(c,d)] +=
                ddldpt[kkk+pkn_LowerTrMatIndex(c,d)]*a1a2;
            ddjac[pkn_LowerTrMatIndex(c,d)] +=
                (dldpt[jj+c]*dldpt[kk+d] + dldpt[jj+d]*dldpt[kk+c])*a1;
          }
          if ( ci3+c >= d )
            ddjac[pkn_LowerTrMatIndex(ci3+c,d)] +=
                  dldpt[ii+c]*dldpt[jj+d]*a3 + dldpt[ii+c]*dldpt[kk+d]*a2;
          cc = c + 3*(i-j);
          dd = d - 3*(i-j);
          if ( dd >= 0 && cc < 3*(n+1) && cc >= ci3+dd )
            ddjac[pkn_LowerTrMatIndex(cc,ci3+dd)] +=
                  dldpt[ii+dd]*dldpt[jj+cc]*a3 + dldpt[ii+dd]*dldpt[kk+cc]*a2;
        }
      }
      _imc_AdjustGrad ( n, N, i, j, ci3, nia3, djac );
      _imc_AdjustHes2 ( n, nvcp, i, j, ci3, nia3, ddjac );
         /* pelna funkcja podcalkowa */
      _imc_ddProduct ( nia31, KtoP, jac, &KtoPj,
                       dKtoP, djac, dKtoPj, ddKtoP, ddjac, ddKtoPj );
      intf4 += KtoPj;
      _imc_CopyGrad ( N, i, j, j, ci3, 0, nia31, dKtoPj, g4 );
      _imc_CopyHes2 ( N, i, j, j, ci3, 0, nia31, ddKtoPj, h4 );
    }
  }
          /* the two loops iterate over the cubes (i,j,j) for i < j */
  tj = 6*(n+1);
  for ( i = 0; i < lkn-2*n-1; i++ ) {
    pt1 = pt[i*nqkn+q1];
    for ( j = i+1; j < lkn-2*n; j++ ) {
      ti = tj - 3*min(j-i,n+1);
      nia3 = maxnia3-ti;  cj3 = tj-ti;
      nia = nia3/3;  cj = cj3/3;
      npia = (nia*(nia+1))/2;
      npia3 = (nia3*(nia3+1))/2;
      memset ( dfv12, 0, nia3*sizeof(double) );
      memset ( dfv23, 0, nia3*sizeof(double) );
      memset ( dfv31, 0, nia3*sizeof(double) );
      memset ( ddfv12, 0, npia*sizeof(double) );
      memset ( ddfv23, 0, npia*sizeof(double) );
      memset ( ddfv31, 0, npia*sizeof(double) );

      pt2 = pt[j*nqkn+q2];
      SubtractPoints3d ( &pt1, &pt2, &v12 );
      pt3 = pt[j*nqkn+q3];
      SubtractPoints3d ( &pt2, &pt3, &v23 );
      SubtractPoints3d ( &pt3, &pt1, &v31 );
      fv12 = DotProduct3d ( &v12, &v12 );  fl12 = sqrt ( fv12 );
      fv23 = DotProduct3d ( &v23, &v23 );  fl23 = sqrt ( fv23 );
      fv31 = DotProduct3d ( &v31, &v31 );  fl31 = sqrt ( fv31 );

      for ( c = c3 = 0;  c <= n;  c++, c3 += 3 ) {
        a1 = 2.0*bsf[(n+1)*q1+c];
        a2 = 2.0*bsf[(n+1)*q2+c];
        a3 = 2.0*bsf[(n+1)*q3+c];
        dfv12[c3]   += a1*v12.x;
        dfv12[c3+1] += a1*v12.y;
        dfv12[c3+2] += a1*v12.z;
        dfv12[cj3+c3]   -= a2*v12.x;
        dfv12[cj3+c3+1] -= a2*v12.y;
        dfv12[cj3+c3+2] -= a2*v12.z;
        dfv23[cj3+c3]   += a2*v23.x;
        dfv23[cj3+c3+1] += a2*v23.y;
        dfv23[cj3+c3+2] += a2*v23.z;
        dfv23[cj3+c3]   -= a3*v23.x;
        dfv23[cj3+c3+1] -= a3*v23.y;
        dfv23[cj3+c3+2] -= a3*v23.z;
        dfv31[cj3+c3]   += a3*v31.x;
        dfv31[cj3+c3+1] += a3*v31.y;
        dfv31[cj3+c3+2] += a3*v31.z;
        dfv31[c3]   -= a1*v31.x;
        dfv31[c3+1] -= a1*v31.y;
        dfv31[c3+2] -= a1*v31.z;
        for ( d = 0; d <= n; d++ ) {
          if ( c >= d ) {
            ddfv12[pkn_LowerTrMatIndex(c,d)] += a1*bsf[(n+1)*q1+d];
            ddfv12[pkn_LowerTrMatIndex(cj+c,cj+d)] += a2*bsf[(n+1)*q2+d];
            ddfv23[pkn_LowerTrMatIndex(cj+c,cj+d)] +=
                    a2*(bsf[(n+1)*q2+d]-bsf[(n+1)*q3+d])
                    -a3*(bsf[(n+1)*q2+d]-bsf[(n+1)*q3+d]);
            ddfv31[pkn_LowerTrMatIndex(cj+c,cj+d)] += a3*bsf[(n+1)*q3+d];
            ddfv31[pkn_LowerTrMatIndex(c,d)] += a1*bsf[(n+1)*q1+d];
          }
          if ( c >= cj+d ) {
            ddfv12[pkn_LowerTrMatIndex(c,cj+d)] -= a1*bsf[(n+1)*q2+d];
            ddfv31[pkn_LowerTrMatIndex(c,cj+d)] -= a1*bsf[(n+1)*q3+d];
          }
          if ( cj+c >= d ) {
            ddfv12[pkn_LowerTrMatIndex(cj+c,d)] -= a2*bsf[(n+1)*q1+d];
            ddfv31[pkn_LowerTrMatIndex(cj+c,d)] -= a3*bsf[(n+1)*q1+d];
          }
        }
      }
      nia31 = _imc_AdjustGrad ( n, N, j, i, cj3, nia3, dfv12 );
      _imc_AdjustGrad ( n, N, j, i, cj3, nia3, dfv23 );
      _imc_AdjustGrad ( n, N, j, i, cj3, nia3, dfv31 );
      nia1 = _imc_AdjustHes1 ( n, nvcp, j, i, cj, nia, ddfv12 );
      _imc_AdjustHes1 ( n, nvcp, j, i, cj, nia, ddfv23 );
      _imc_AdjustHes1 ( n, nvcp, j, i, cj, nia, ddfv31 );
      npia31 = (nia31*(nia31+1))/2;

      a1 = 0.5/fl12;
      pkn_MultMatrixNumd ( 1, nia31, 0, dfv12, a1, 0, dfl12 );
      _imc_HLab ( nia1, fv12, fl12, a1, dfv12, ddfv12, dfl12, ddfl12 );
      if ( q2 == q3 ) {
          /* licznik wyrazenia dla krzywizny */
        _imc_ddSpecialK2num ( n, N, nvcp, bsf, dbsf, bsf1, q1, q2,
                              i, j, cj3, nia3, npia3, &v12, &dpt[j*nqkn+q2],
                              auxv, &K2num, dK2num, ddK2num );
          /* mianownik wyrazenia dla krzywizny */
        _imc_ddSpecialK2den ( n, N, nvcp, &dbsf[q2*(n+1)], &bsf1[q2*n],
                              i, j, cj3, nia3, &dcp[j], &dpt[j*nqkn+q2], fv12,
                              dfv12, daux1, daux2, ddfv12, ddaux1, ddaux2,
                              &K2den, dK2den, ddK2den );
      }
      else {
        a1 = 0.5/fl23;
        pkn_MultMatrixNumd ( 1, nia31, 0, dfv23, a1, 0, dfl23 );
        _imc_HLab ( nia1, fv23, fl23, a1, dfv23, ddfv23, dfl23, ddfl23 );
        a1 = 0.5/fl31;
        pkn_MultMatrixNumd ( 1, nia31, 0, dfv31, a1, 0, dfl31 );
        _imc_HLab ( nia1, fv31, fl31, a1, dfv31, ddfv31, dfl31, ddfl31 );
        pkn_MultMatrixNumd ( 1, nia31, 0, dfv23, 0.5/fl23, 0, dfl23 );
        pkn_MultMatrixNumd ( 1, nia31, 0, dfv31, 0.5/fl31, 0, dfl31 );
        half = 0.5*(fl12+fl23+fl31);
        pkn_AddMatrixd ( 1, nia31, 0, dfl12, 0, dfl23, 0, dhalf );
        pkn_AddMatrixd ( 1, nia31, 0, dhalf, 0, dfl31, 0, dhalf );
        pkn_MultMatrixNumd ( 1, nia31, 0, dhalf, 0.5, 0, dhalf );
        pkn_AddMatrixd ( 1, npia31, 0, ddfl12, 0, ddfl23, 0, ddhalf );
        pkn_AddMatrixd ( 1, npia31, 0, ddhalf, 0, ddfl31, 0, ddhalf );
        pkn_MultMatrixNumd ( 1, npia31, 0, ddhalf, 0.5, 0, ddhalf );
            /* licznik wyrazenia dla krzywizny */
        _imc_ddK2num ( nia31, fl12, fl23, fl31, half,
                       dfl12, dfl23, dfl31, dhalf,
                       ddfl12, ddfl23, ddfl31, ddhalf,
                       daux1, daux2, daux3, daux4, daux5,
                       ddaux1, ddaux2, ddaux3, ddaux4,
                       &K2num, dK2num, ddK2num );
            /* mianownik wyrazenia dla krzywizny */
        _imc_ddK2den ( nia31, fv12, fv23, fv31, dfv12, dfv23, dfv31, daux1,
                       ddfv12, ddfv23, ddfv31, ddaux1, ddaux2, ddaux3,
                       &K2den, dK2den, ddK2den );
      }
          /* kwadrat krzywizny */
      _imc_ddRatio ( nia31, K2num, K2den, &K2, dK2num, dK2den, dK2,
                     ddK2num, ddK2den, ddK2 );
          /* krzywizna w potedze p */
      _imc_ddPower ( nia31, K2, dK2, ddK2, 0.5*md->w, &KtoP, dKtoP, ddKtoP );
          /* jakobian */
      memset ( djac, 0, nia3*sizeof(double) );
      memset ( ddjac, 0, npia3*sizeof(double) );
      a1 = ldpt[i*nqkn+q1];  a2 = ldpt[j*nqkn+q2];  a3 = ldpt[j*nqkn+q3];
      a1a2 = a1*a2;  a2a3 = a2*a3;  a3a1 = a3*a1;
      jac = a1a2*a3;
      kk = 3*(n+1);
      ii = kk*(nqkn*i+q1);  jj = kk*(nqkn*j+q2);  kk *= nqkn*j+q3;
      kkk = ((3*n+3)*(3*n+4))/2;
      iii = kkk*(nqkn*i+q1);  jjj = kkk*(nqkn*j+q2);  kkk *= nqkn*j+q3;
      for ( c = 0; c < 3*(n+1); c++ ) {
        djac[c]     += dldpt[ii+c]*a2a3;
        djac[cj3+c] += dldpt[jj+c]*a3a1;
        djac[cj3+c] += dldpt[kk+c]*a1a2;
      }
      for ( c = 0; c < 3*(n+1); c++ )
        for ( d = 0; d <= c; d++ )
          ddjac[pkn_LowerTrMatIndex(c,d)] +=
                ddldpt[iii+pkn_LowerTrMatIndex(c,d)]*a2a3;
      for ( c = 0; c < 3*(n+1); c++ )
        for ( d = 0; d <= c; d++ )
            ddjac[pkn_LowerTrMatIndex(cj3+c,cj3+d)] +=
                ddldpt[jjj+pkn_LowerTrMatIndex(c,d)]*a3a1 +
                ddldpt[kkk+pkn_LowerTrMatIndex(c,d)]*a1a2 +
                (dldpt[jj+c]*dldpt[kk+d] + dldpt[jj+d]*dldpt[kk+c])*a1;
      for ( c = cj3; c < 3*(n+1); c++ )
        for ( d = cj3; d <= c; d++ )
          ddjac[pkn_SymMatIndex(c,d)] +=
              dldpt[ii+c]*dldpt[jj+d-cj3]*a3 +
              dldpt[ii+c]*dldpt[kk+d-cj3]*a2;
      for ( c = cj3; c < cj3+3*(n+1); c++ )
        for ( d = 0; d <= c && d < 3*(n+1); d++ )
          ddjac[pkn_SymMatIndex(c,d)] +=
              dldpt[jj+c-cj3]*dldpt[ii+d]*a3 +
              dldpt[kk+c-cj3]*dldpt[ii+d]*a2;
              
      _imc_AdjustGrad ( n, N, j, i, cj3, nia3, djac );
      _imc_AdjustHes2 ( n, nvcp, j, i, cj3, nia3, ddjac );
         /* pelna funkcja podcalkowa */
      _imc_ddProduct ( nia31, KtoP, jac, &KtoPj,
                       dKtoP, djac, dKtoPj, ddKtoP, ddjac, ddKtoPj );
      intf4 += KtoPj;
      _imc_CopyGrad ( N, j, j, i, cj3, cj3, nia31, dKtoPj, g4 );
      _imc_CopyHes2 ( N, j, j, i, cj3, cj3, nia31, ddKtoPj, h4 );
    }
  }
          /* this loop iterates over the cubes (i,i,i) */
  intf5 = 0.0;
  memset ( g5, 0, N*sizeof(double) );
  memset ( h5, 0, M*sizeof(double) );
  nia = nia1 = n+1;
  nia3 = nia31 = 3*nia;
  npia = (nia*(nia+1))/2;
  npia3 = npia31 = (nia3*(nia3+1))/2;
  for ( i = 0; i < lkn-2*n; i++ ) {
    memset ( dfv12, 0, nia3*sizeof(double) );
    memset ( dfv23, 0, nia3*sizeof(double) );
    memset ( dfv31, 0, nia3*sizeof(double) );
    memset ( ddfv12, 0, npia*sizeof(double) );
    memset ( ddfv23, 0, npia*sizeof(double) );
    memset ( ddfv31, 0, npia*sizeof(double) );

    pt1 = pt[i*nqkn+q1];
    pt2 = pt[i*nqkn+q2];
    SubtractPoints3d ( &pt1, &pt2, &v12 );
    pt3 = pt[i*nqkn+q3];
    SubtractPoints3d ( &pt2, &pt3, &v23 );
    SubtractPoints3d ( &pt3, &pt1, &v31 );
    fv12 = DotProduct3d ( &v12, &v12 );  fl12 = sqrt ( fv12 );
    fv23 = DotProduct3d ( &v23, &v23 );  fl23 = sqrt ( fv23 );
    fv31 = DotProduct3d ( &v31, &v31 );  fl31 = sqrt ( fv31 );

    for ( c = c3 = 0;  c <= n;  c++, c3 += 3 ) {
      a1 = 2.0*bsf[(n+1)*q1+c];
      a2 = 2.0*bsf[(n+1)*q2+c];
      a3 = 2.0*bsf[(n+1)*q3+c];
      dfv12[c3]   += a1*v12.x;
      dfv12[c3+1] += a1*v12.y;
      dfv12[c3+2] += a1*v12.z;
      dfv12[c3]   -= a2*v12.x;
      dfv12[c3+1] -= a2*v12.y;
      dfv12[c3+2] -= a2*v12.z;
      dfv23[c3]   += a2*v23.x;
      dfv23[c3+1] += a2*v23.y;
      dfv23[c3+2] += a2*v23.z;
      dfv23[c3]   -= a3*v23.x;
      dfv23[c3+1] -= a3*v23.y;
      dfv23[c3+2] -= a3*v23.z;
      dfv31[c3]   += a3*v31.x;
      dfv31[c3+1] += a3*v31.y;
      dfv31[c3+2] += a3*v31.z;
      dfv31[c3]   -= a1*v31.x;
      dfv31[c3+1] -= a1*v31.y;
      dfv31[c3+2] -= a1*v31.z;
      for ( d = 0; d <= n; d++ )
        if ( c >= d ) {
          ddfv12[pkn_LowerTrMatIndex(c,d)] +=
                  a1*(bsf[(n+1)*q1+d]-bsf[(n+1)*q2+d])
                  -a2*(bsf[(n+1)*q1+d]-bsf[(n+1)*q2+d]);
          ddfv23[pkn_LowerTrMatIndex(c,d)] +=
                  a2*(bsf[(n+1)*q2+d]-bsf[(n+1)*q3+d])
                  -a3*(bsf[(n+1)*q2+d]-bsf[(n+1)*q3+d]);
          ddfv31[pkn_LowerTrMatIndex(c,d)] +=
                  a3*(bsf[(n+1)*q3+d]-bsf[(n+1)*q1+d])
                  -a1*(bsf[(n+1)*q3+d]-bsf[(n+1)*q1+d]);
        }
    }
    memset ( dK2num, 0, nia31*sizeof(double) );
    memset ( ddK2num, 0, npia31*sizeof(double) );
    memset ( dK2den, 0, nia31*sizeof(double) );
    memset ( ddK2den, 0, npia31*sizeof(double) );
    if ( q1 == q2 ) {
      if ( q2 == q3 ) {
          /* przypadek specjalny - trzeba obliczyc krzywizne okregu */
          /* scisle stycznego do krzywej */
            /* licznik wyrazenia dla krzywizny */
        CrossProduct3d ( &ddpt[i*nqkn+q1], &dpt[i*nqkn+q1], &auxv1 );
        K2num = DotProduct3d ( &auxv1, &auxv1 );
        memset ( auxv, 0, nia31*sizeof(vector3d) );
        for ( ii = 1, ic = 3;  ii < n;  ii++, ic += 3 )
          for ( jj = jc = 0;  jj < ii;  jj++, jc += 3 ) {
            a1 = dbsf1[q1*n+ii]*bsf1[q1*n+jj] - dbsf1[q1*n+jj]*bsf1[q1*n+ii];
            a2 = dcp[i+jj].z*a1;  auxv[ic  ].y -= a2;  auxv[ic+3].y += a2;
            a2 = dcp[i+jj].y*a1;  auxv[ic  ].z += a2;  auxv[ic+3].z -= a2;
            a2 = (dcp[i+jj].y*auxv1.z - dcp[i+jj].z*auxv1.y)*a1;
            dK2num[ic  ] -= a2;  dK2num[ic+3] += a2;
            a2 = dcp[i+jj].y*a1;  auxv[ic+2].x -= a2;  auxv[ic+5].x += a2;
            a2 = dcp[i+jj].x*a1;  auxv[ic+2].y += a2;  auxv[ic+5].y -= a2;
            a2 = (dcp[i+jj].z*auxv1.x - dcp[i+jj].x*auxv1.z)*a1;
            dK2num[ic+1] -= a2;  dK2num[ic+4] += a2;
            a2 = dcp[i+jj].x*a1;  auxv[ic+1].z -= a2;  auxv[ic+4].z += a2;
            a2 = dcp[i+jj].z*a1;  auxv[ic+1].x += a2;  auxv[ic+4].x -= a2;
            a2 = (dcp[i+jj].x*auxv1.y - dcp[i+jj].y*auxv1.x)*a1;
            dK2num[ic+2] -= a2;  dK2num[ic+5] += a2;
            a2 = dcp[i+ii].z*a1;  auxv[jc  ].y += a2;  auxv[jc+3].y -= a2;
            a2 = dcp[i+ii].y*a1;  auxv[jc  ].z -= a2;  auxv[jc+3].z += a2;
            a2 = (dcp[i+ii].y*auxv1.z - dcp[i+ii].z*auxv1.y)*a1;
            dK2num[jc  ] += a2;  dK2num[jc+3] -= a2;
            a2 = dcp[i+ii].y*a1;  auxv[jc+2].x += a2;  auxv[jc+5].x -= a2;
            a2 = dcp[i+ii].x*a1;  auxv[jc+2].y -= a2;  auxv[jc+5].y += a2;
            a2 = (dcp[i+ii].z*auxv1.x - dcp[i+ii].x*auxv1.z)*a1;
            dK2num[jc+1] += a2;  dK2num[jc+4] -= a2;
            a2 = dcp[i+ii].x*a1;  auxv[jc+1].z += a2;  auxv[jc+4].z -= a2;
            a2 = dcp[i+ii].z*a1;  auxv[jc+1].x -= a2;  auxv[jc+4].x += a2;
            a2 = (dcp[i+ii].x*auxv1.y - dcp[i+ii].y*auxv1.x)*a1;
            dK2num[jc+2] += a2;  dK2num[jc+5] -= a2;
          }
        for ( ii = 0; ii < nia31; ii++ )
          for ( jj = 0; jj <= ii; jj++ )
            ddK2num[pkn_LowerTrMatIndex(ii,jj)] += DotProduct3d ( &auxv[ii], &auxv[jj] );
        for ( ii = 1, ic = 3;  ii < n;  ii++, ic += 3 )
          for ( jj = jc = 0;  jj < ii;  jj++, jc += 3 ) {
            a1 = dbsf1[q1*n+ii]*bsf1[q1*n+jj] - dbsf1[q1*n+jj]*bsf1[q1*n+ii];
            a2 = auxv1.z*a1;
            ddK2num[pkn_LowerTrMatIndex(ic,  jc+1)] += a2;
            ddK2num[pkn_LowerTrMatIndex(ic+1,jc  )] -= a2;
            ddK2num[pkn_LowerTrMatIndex(ic+3,jc+1)] -= a2;
            ddK2num[pkn_LowerTrMatIndex(ic+4,jc  )] += a2;
            ddK2num[pkn_LowerTrMatIndex(ic+3,jc+4)] += a2;
            ddK2num[pkn_LowerTrMatIndex(ic+4,jc+3)] -= a2;
            if ( ii > jj+1 ) {
              ddK2num[pkn_LowerTrMatIndex(ic  ,jc+4)] -= a2;
              ddK2num[pkn_LowerTrMatIndex(ic+1,jc+3)] += a2;
            }
            a2 = auxv1.y*a1;
            ddK2num[pkn_LowerTrMatIndex(ic,  jc+2)] -= a2;
            ddK2num[pkn_LowerTrMatIndex(ic+2,jc  )] += a2;
            ddK2num[pkn_LowerTrMatIndex(ic+3,jc+2)] += a2;
            ddK2num[pkn_LowerTrMatIndex(ic+5,jc  )] -= a2;
            ddK2num[pkn_LowerTrMatIndex(ic+3,jc+5)] -= a2;
            ddK2num[pkn_LowerTrMatIndex(ic+5,jc+3)] += a2;
            if ( ii > jj+1 ) {
              ddK2num[pkn_LowerTrMatIndex(ic,  jc+5)] += a2;
              ddK2num[pkn_LowerTrMatIndex(ic+2,jc+3)] -= a2;
            }
            a2 = auxv1.x*a1;
            ddK2num[pkn_LowerTrMatIndex(ic+1,jc+2)] += a2;
            ddK2num[pkn_LowerTrMatIndex(ic+2,jc+1)] -= a2;
            ddK2num[pkn_LowerTrMatIndex(ic+4,jc+2)] -= a2;
            ddK2num[pkn_LowerTrMatIndex(ic+5,jc+1)] += a2;
            ddK2num[pkn_LowerTrMatIndex(ic+4,jc+5)] += a2;
            ddK2num[pkn_LowerTrMatIndex(ic+5,jc+4)] -= a2;
            if ( ii > jj+1 ) {
              ddK2num[pkn_LowerTrMatIndex(ic+1,jc+5)] -= a2;
              ddK2num[pkn_LowerTrMatIndex(ic+2,jc+4)] += a2;
            }
          }
        pkn_AddMatrixd ( 1, nia31, 0, dK2num, 0, dK2num, 0, dK2num );
        pkn_AddMatrixd ( 1, npia31, 0, ddK2num, 0, ddK2num, 0, ddK2num );
            /* mianownik wyrazenia dla krzywizny */
        a2 = DotProduct3d ( &dpt[i*nqkn+q1], &dpt[i*nqkn+q1] );
        a1 = a2*a2;
        K2den = a1*a2;
        for ( ii = ic = 0;  ii < n;  ii++, ic += 3 ) {
          a2 = 2.0*bsf1[q1*n+ii];
          for ( jj = 0;  jj < n;  jj++ ) {
            a3 = a2*bsf1[q1*n+jj];
            a = dcp[i+jj].x*a3;  dK2den[ic  ] -= a;  dK2den[ic+3] += a;
            a = dcp[i+jj].y*a3;  dK2den[ic+1] -= a;  dK2den[ic+4] += a;
            a = dcp[i+jj].z*a3;  dK2den[ic+2] -= a;  dK2den[ic+5] += a;
          }
        }
        for ( ii = ic = 0;  ii <= n;  ii++, ic += 3 ) {
          a2 = 2.0*dbsf[q1*(n+1)+ii];
          for ( jj = jc = 0;  jj <= ii;  jj++, jc += 3 )
            ddK2den[pkn_LowerTrMatIndex(ic,jc)] += a2*dbsf[q1*(n+1)+jj];
        }
        pkn_MultMatrixNumd ( 1, nia31, 0, dK2den, 3.0*a1, 0, dK2den );
        a3 = 2.0/(3.0*K2den);
        for ( ii = ic = 0;  ii <= n;  ii++, ic += 3 ) {
          for ( jj = jc = 0;  jj < ii;  jj++, jc += 3 ) {
            a2 = 3.0*a1*ddK2den[pkn_LowerTrMatIndex(ic,jc)];
            ddK2den[pkn_LowerTrMatIndex(ic,jc)] =
            ddK2den[pkn_LowerTrMatIndex(ic+1,jc+1)] =
            ddK2den[pkn_LowerTrMatIndex(ic+2,jc+2)] = a2;
            ddK2den[pkn_LowerTrMatIndex(ic  ,jc  )] += a3*dK2den[ic]*dK2den[jc];
            ddK2den[pkn_LowerTrMatIndex(ic  ,jc+1)] += a3*dK2den[ic]*dK2den[jc+1];
            ddK2den[pkn_LowerTrMatIndex(ic  ,jc+2)] += a3*dK2den[ic]*dK2den[jc+2];
            ddK2den[pkn_LowerTrMatIndex(ic+1,jc  )] += a3*dK2den[ic+1]*dK2den[jc];
            ddK2den[pkn_LowerTrMatIndex(ic+1,jc+1)] += a3*dK2den[ic+1]*dK2den[jc+1];
            ddK2den[pkn_LowerTrMatIndex(ic+1,jc+2)] += a3*dK2den[ic+1]*dK2den[jc+2];
            ddK2den[pkn_LowerTrMatIndex(ic+2,jc  )] += a3*dK2den[ic+2]*dK2den[jc];
            ddK2den[pkn_LowerTrMatIndex(ic+2,jc+1)] += a3*dK2den[ic+2]*dK2den[jc+1];
            ddK2den[pkn_LowerTrMatIndex(ic+2,jc+2)] += a3*dK2den[ic+2]*dK2den[jc+2];
          }
          a2 = 3.0*a1*ddK2den[pkn_LowerTrMatIndex(ic,ic)];
          ddK2den[pkn_LowerTrMatIndex(ic,ic)] =
          ddK2den[pkn_LowerTrMatIndex(ic+1,ic+1)] =
          ddK2den[pkn_LowerTrMatIndex(ic+2,ic+2)] = a2;
          ddK2den[pkn_LowerTrMatIndex(ic  ,ic  )] += a3*dK2den[ic]*dK2den[ic];
          ddK2den[pkn_LowerTrMatIndex(ic+1,ic  )] += a3*dK2den[ic+1]*dK2den[ic];
          ddK2den[pkn_LowerTrMatIndex(ic+1,ic+1)] += a3*dK2den[ic+1]*dK2den[ic+1];
          ddK2den[pkn_LowerTrMatIndex(ic+2,ic  )] += a3*dK2den[ic+2]*dK2den[ic];
          ddK2den[pkn_LowerTrMatIndex(ic+2,ic+1)] += a3*dK2den[ic+2]*dK2den[ic+1];
          ddK2den[pkn_LowerTrMatIndex(ic+2,ic+2)] += a3*dK2den[ic+2]*dK2den[ic+2];
        }
      }
      else {
            /* licznik wyrazenia dla krzywizny */
        _imc_ddSpecialK2num ( n, N, nvcp, bsf, dbsf, bsf1, q3, q1,
                              i, i, 0, nia3, npia3, &v31, &dpt[i*nqkn+q1],
                              auxv, &K2num, dK2num, ddK2num );
            /* mianownik wyrazenia dla krzywizny */
        _imc_ddSpecialK2den ( n, N, nvcp, &dbsf[q1*(n+1)], &bsf1[q1*n],
                              i, i, 0, nia3, &dcp[i], &dpt[i*nqkn+q1], fv31,
                              dfv31, daux1, daux2, ddfv31, ddaux1, ddaux2,
                              &K2den, dK2den, ddK2den );
      }
    }
    else if ( q2 == q3 ) {
            /* licznik wyrazenia dla krzywizny */
      _imc_ddSpecialK2num ( n, N, nvcp, bsf, dbsf, bsf1, q1, q2,
                            i, i, 0, nia3, npia3, &v12, &dpt[i*nqkn+q2],
                            auxv, &K2num, dK2num, ddK2num );
            /* mianownik wyrazenia dla krzywizny */
      _imc_ddSpecialK2den ( n, N, nvcp, &dbsf[q2*(n+1)], &bsf1[q2*n],
                            i, i, 0, nia3, &dcp[i], &dpt[i*nqkn+q2], fv12,
                            dfv12, daux1, daux2, ddfv12, ddaux1, ddaux2,
                            &K2den, dK2den, ddK2den );
    }
    else if ( q1 == q3 ) {
            /* licznik wyrazenia dla krzywizny */
      _imc_ddSpecialK2num ( n, N, nvcp, bsf, dbsf, bsf1, q2, q3,
                            i, i, 0, nia3, npia3, &v23, &dpt[i*nqkn+q3],
                            auxv, &K2num, dK2num, ddK2num );
            /* mianownik wyrazenia dla krzywizny */
      _imc_ddSpecialK2den ( n, N, nvcp, &dbsf[q3*(n+1)], &bsf1[q3*n],
                            i, i, 0, nia3, &dcp[i], &dpt[i*nqkn+q3], fv23,
                            dfv23, daux1, daux2, ddfv23, ddaux1, ddaux2,
                            &K2den, dK2den, ddK2den );
    }
    else {
      a1 = 0.5/fl12;
      pkn_MultMatrixNumd ( 1, nia31, 0, dfv12, a1, 0, dfl12 );
      _imc_HLab ( nia, fv12, fl12, a1, dfv12, ddfv12, dfl12, ddfl12 );
      a1 = 0.5/fl23;
      pkn_MultMatrixNumd ( 1, nia31, 0, dfv23, a1, 0, dfl23 );
      _imc_HLab ( nia, fv23, fl23, a1, dfv23, ddfv23, dfl23, ddfl23 );
      a1 = 0.5/fl31;
      pkn_MultMatrixNumd ( 1, nia31, 0, dfv31, a1, 0, dfl31 );
      _imc_HLab ( nia, fv31, fl31, a1, dfv31, ddfv31, dfl31, ddfl31 );
      half = 0.5*(fl12+fl23+fl31);
      pkn_AddMatrixd ( 1, nia31, 0, dfl12, 0, dfl23, 0, dhalf );
      pkn_AddMatrixd ( 1, nia31, 0, dhalf, 0, dfl31, 0, dhalf );
      pkn_MultMatrixNumd ( 1, nia31, 0, dhalf, 0.5, 0, dhalf );
      pkn_AddMatrixd ( 1, npia31, 0, ddfl12, 0, ddfl23, 0, ddhalf );
      pkn_AddMatrixd ( 1, npia31, 0, ddhalf, 0, ddfl31, 0, ddhalf );
      pkn_MultMatrixNumd ( 1, npia31, 0, ddhalf, 0.5, 0, ddhalf );
            /* licznik wyrazenia dla krzywizny */
      _imc_ddK2num ( nia31, fl12, fl23, fl31, half,
                     dfl12, dfl23, dfl31, dhalf,
                     ddfl12, ddfl23, ddfl31, ddhalf,
                     daux1, daux2, daux3, daux4, daux5,
                     ddaux1, ddaux2, ddaux3, ddaux4,
                     &K2num, dK2num, ddK2num );
            /* mianownik wyrazenia dla krzywizny */
      _imc_ddK2den ( nia31, fv12, fv23, fv31, dfv12, dfv23, dfv31, daux1,
                     ddfv12, ddfv23, ddfv31, ddaux1, ddaux2, ddaux3,
                     &K2den, dK2den, ddK2den );
    }
          /* kwadrat krzywizny */
    _imc_ddRatio ( nia31, K2num, K2den, &K2, dK2num, dK2den, dK2,
                   ddK2num, ddK2den, ddK2 );
          /* krzywizna w potedze p */
    _imc_ddPower ( nia31, K2, dK2, ddK2, 0.5*md->w, &KtoP, dKtoP, ddKtoP );
          /* jakobian */
    memset ( djac, 0, nia3*sizeof(double) );
    memset ( ddjac, 0, npia3*sizeof(double) );
    a1 = ldpt[i*nqkn+q1];  a2 = ldpt[i*nqkn+q2];  a3 = ldpt[i*nqkn+q3];
    a1a2 = a1*a2;  a2a3 = a2*a3;  a3a1 = a3*a1;
    jac = a1a2*a3;
    kk = 3*(n+1);
    ii = kk*(nqkn*i+q1);  jj = kk*(nqkn*i+q2);  kk *= nqkn*i+q3;
    kkk = ((3*n+3)*(3*n+4))/2;
    iii = kkk*(nqkn*i+q1);  jjj = kkk*(nqkn*i+q2);  kkk *= nqkn*i+q3;
    for ( c = 0; c < 3*(n+1); c++ ) {
      djac[c] += dldpt[ii+c]*a2a3;
      djac[c] += dldpt[jj+c]*a3a1;
      djac[c] += dldpt[kk+c]*a1a2;
      for ( d = 0; d <= c; d++ )
        ddjac[pkn_LowerTrMatIndex(c,d)] +=
            ddldpt[iii+pkn_LowerTrMatIndex(c,d)]*a2a3 +
            ddldpt[jjj+pkn_LowerTrMatIndex(c,d)]*a3a1 +
            ddldpt[kkk+pkn_LowerTrMatIndex(c,d)]*a1a2 +
            (dldpt[ii+c]*dldpt[jj+d] + dldpt[ii+d]*dldpt[jj+c])*a3 +
            (dldpt[jj+c]*dldpt[kk+d] + dldpt[jj+d]*dldpt[kk+c])*a1 +
            (dldpt[ii+c]*dldpt[kk+d] + dldpt[ii+d]*dldpt[kk+c])*a2;
    }
         /* pelna funkcja podcalkowa */
    _imc_ddProduct ( nia31, KtoP, jac, &KtoPj,
                     dKtoP, djac, dKtoPj, ddKtoP, ddjac, ddKtoPj );
    intf5 += KtoPj;
    _imc_CopyGrad ( N, i, i, i, 0, 0, nia31, dKtoPj, g5 );
    _imc_CopyHes2 ( N, i, i, i, 0, 0, nia31, ddKtoPj, h5 );
  }

  i = (q1*nqkn+q2)*nqkn+q3;
  data->_f[i] = 6.0*intf3 + 3.0*intf4 + intf5;
  pkn_AddMatrixMd ( 1, N, 0, g4, 0, g3, 2.0, 0, g3 );
  pkn_AddMatrixMd ( 1, N, 0, g5, 0, g3, 3.0, 0, &data->_g[i*N] );
  pkn_AddMatrixMd ( 1, M, 0, h4, 0, h3, 2.0, 0, h3 );
  pkn_AddMatrixMd ( 1, M, 0, h5, 0, h3, 3.0, 0, &data->_h[i*M] );
  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*_intKM_FGH*/

boolean _mengerc_hessIntF ( mengerc_data *md,
                            double *func, double *grad, double *hess )
{
  void           *sp;
  intKM_job_desc data;
  int            lkn, n, ncp, nvcp, N, M, npoints, nqkn;
  double         *knots, *qkn, *qc;
  point3d        *cpoints;
  vector3d       *dcp, *pt, *dpt, *ddpt;
  double         *ldpt, *dldpt, *ddldpt;
  int            i, ii, im, j, jj, k, kk, l, ll, m, q1, iii, jjj, p;
  double         Nl, Nm, a, b;
  double         *ddM;
  double         *_f, *_g, *_h, ff1, ff2, ff3, *gf2, *gf3, *hf2, *hf3;
  int3           jobsize;
  boolean        success;

  sp = pkv_GetScratchMemTop ();
  memset ( &data, 0, sizeof(intKM_job_desc) );
  data.md = md;
  n       = md->deg;
  lkn     = md->lkn;
  knots   = md->knots;
  cpoints = md->cpoints;
  nqkn    = md->nqkn;
  qkn     = md->qkn;
  ncp     = lkn - n;
  nvcp    = ncp - n;
  N = 3*nvcp;
  M = (N*(N+1))/2;
  dcp = data.dcp = pkv_GetScratchMem ( (ncp-1)*sizeof(vector3d) );
  _f = pkv_GetScratchMemd ( nqkn*nqkn*nqkn*(1+N+M) );
  if ( !dcp || !_f )
    goto failure;
  data._f = _f;
  data._g = _g = &_f[nqkn*nqkn*nqkn];
  data._h = _h = &_g[N*nqkn*nqkn*nqkn];
        /* liczba punktow krzywej */
  npoints = nqkn*(lkn-2*n);
  pt = data.pt = pkv_GetScratchMem ( 3*npoints*sizeof(point3d) +
                           (3*n+4+((3*n+3)*(3*n+4))/2)*npoints*sizeof(double) );
  if ( !pt )
    goto failure;
  dpt = data.dpt = &pt[npoints];
  ddpt = data.ddpt = &dpt[npoints];
  ldpt = data.ldpt = (double*)&ddpt[npoints];
  dldpt = data.dldpt = &ldpt[npoints];
  ddldpt = data.ddldpt = &dldpt[3*(n+1)*npoints];
  if ( M >= nqkn*(((n+1)*(n+2))/2) )
    ddM = hess;  /* ta tablica jest potrzebna tylko na chwile */
  else
    ddM = pkv_GetScratchMemd ( nqkn*(((n+1)*(n+2))/2) );
  if ( !ddM )
    goto failure;
        /* oblicz pochodne drugiego rzedu kwadratu dlugosci pochodnej */
        /* parametryzacji - one nie zaleza od krzywej (tylko od funkcji */
        /* B-sklejanych) */
  memset ( ddM, 0, nqkn*(((n+1)*(n+2))/2)*sizeof(double) );
  for ( q1 = 0; q1 < nqkn; q1++ ) {
    iii = q1*(((n+1)*(n+2))/2);
    for ( l = 0; l < n; l++ ) {
      Nl = 2.0*md->bsf1[n*q1+l];
      for ( m = 0; m < n; m++ ) {
        a = Nl*md->bsf1[n*q1+m];
        if ( l >= m ) {
          ddM[iii+pkn_LowerTrMatIndex(l,m)] += a;
          ddM[iii+pkn_LowerTrMatIndex(l+1,m+1)] += a;
        }
        if ( l+1 >= m )
          ddM[iii+pkn_LowerTrMatIndex(l+1,m)] -= a;
        if ( l >= m+1 )
          ddM[iii+pkn_LowerTrMatIndex(l,m+1)] -= a;
      }
    }
  }
        /* znajdz roznice kolejnych punktow kontrolnych */
  for ( i = 0; i < ncp-1; i++ )
    SubtractPoints3d ( &cpoints[i+1], &cpoints[i], &dcp[i] );
        /* znajdz punkty krzywej i pochodne pierwszego i drugiego rzedu */
  memset ( dldpt, 0, 3*(n+1)*npoints*sizeof(double) );
  memset ( ddldpt, 0, (((3*n+3)*(3*n+4))/2)*npoints*sizeof(double) );
  for ( i = j = 0;  i < lkn-2*n;  i++ )
    for ( q1 = 0;  q1 < nqkn;  q1 ++, j++ ) {
      mbs_deBoorDer2C3d ( n, 2*n+1, knots, &cpoints[i], (double)n+qkn[q1],
                          &pt[j], &dpt[j], &ddpt[j] );
      ldpt[j] = sqrt ( DotProduct3d ( &dpt[j], &dpt[j] ) );
          /* pochodne dlugosci wektora pochodnej wzgledem wspolrzednych */
          /* punktow kontrolnych */
      iii = q1*(((n+1)*(n+2))/2);
      for ( l = 0; l < n; l++ ) {
        Nl = 2.0*md->bsf1[n*q1+l];
        ll = 3*((n+1)*j+l);
        for ( m = 0; m < n; m++ ) {
          Nm = md->bsf1[n*q1+m];
          a = Nl*Nm;
          im = i+m;
          dldpt[ll  ] -= b = a*dcp[im].x;
          dldpt[ll+3] += b;
          dldpt[ll+1] -= b = a*dcp[im].y;
          dldpt[ll+4] += b;
          dldpt[ll+2] -= b = a*dcp[im].z;
          dldpt[ll+5] += b;
        }
      }
      a = 0.5/ldpt[j];
      ll = 3*(n+1)*j;
      pkn_MultMatrixNumd ( 1, 3*(n+1), 0, &dldpt[ll], a, 0, &dldpt[ll] );
      jj = j*((3*n+3)*(3*n+4))/2;
      for ( l = 0; l <= n; l++ ) {
        ii = 3*l;
        for ( m = 0; m < l; m++ ) {
          kk = 3*m;
          b = 0.5*ddM[iii+pkn_LowerTrMatIndex(l,m)];
          jjj = jj+pkn_LowerTrMatIndex(ii,kk);
          ddldpt[jjj]     = b - dldpt[ll+ii]*dldpt[ll+kk];
          ddldpt[jjj+1]   = -dldpt[ll+ii]*dldpt[ll+kk+1];
          ddldpt[jjj+2]   = -dldpt[ll+ii]*dldpt[ll+kk+2];
          jjj = jj+pkn_LowerTrMatIndex(ii+1,kk);
          ddldpt[jjj]   = -dldpt[ll+ii+1]*dldpt[ll+kk];
          ddldpt[jjj+1] = b - dldpt[ll+ii+1]*dldpt[ll+kk+1];
          ddldpt[jjj+2] = -dldpt[ll+ii+1]*dldpt[ll+kk+2];
          jjj = jj+pkn_LowerTrMatIndex(ii+2,kk);
          ddldpt[jjj]   = -dldpt[ll+ii+2]*dldpt[ll+kk];
          ddldpt[jjj+1] = -dldpt[ll+ii+2]*dldpt[ll+kk+1];
          ddldpt[jjj+2] = b - dldpt[ll+ii+2]*dldpt[ll+kk+2];
        }
        b = 0.5*ddM[iii+pkn_LowerTrMatIndex(l,l)];
        ddldpt[jj+pkn_LowerTrMatIndex(ii,ii)] = b - dldpt[ll+ii]*dldpt[ll+ii];
        jjj = jj+pkn_LowerTrMatIndex(ii+1,ii);
        ddldpt[jjj]   = -dldpt[ll+ii+1]*dldpt[ll+ii];
        ddldpt[jjj+1] = b - dldpt[ll+ii+1]*dldpt[ll+ii+1];
        jjj = jj+pkn_LowerTrMatIndex(ii+2,ii);
        ddldpt[jjj]   = -dldpt[ll+ii+2]*dldpt[ll+ii];
        ddldpt[jjj+1] = -dldpt[ll+ii+2]*dldpt[ll+ii+1];
        ddldpt[jjj+2] = b - dldpt[ll+ii+2]*dldpt[ll+ii+2];
      }
      pkn_MultMatrixNumd ( 1, ((3*n+3)*(3*n+4))/2, 0, &ddldpt[jj], a+a,
                           0, &ddldpt[jj] );
    }
        /* oblicz calki */
  if ( md->npthr > 1 ) {
    jobsize.x = jobsize.y = jobsize.z = nqkn;
    pkv_SetPThreadsToWork ( &jobsize, md->npthr, 1048576, 16*1048576,
                            (void*)&data, _intKM_FGH, NULL, NULL, &success );
    if ( !success )
      goto failure;
  }
  else {
  for ( jobsize.x = 0; jobsize.x < nqkn; jobsize.x++ )
    for ( jobsize.y = 0; jobsize.y < nqkn; jobsize.y++ )
      for ( jobsize.z = 0; jobsize.z < nqkn; jobsize.z++ )
        _intKM_FGH ( (void*)&data, &jobsize );
  }
        /* sumuj calki - hesjan juz jest zsumowany */
  gf2 = pkv_GetScratchMemd ( 2*(N+M) );
  if ( !gf2 )
    goto failure;
  gf3 = &gf2[N];
  hf2 = &gf3[N];
  hf3 = &hf2[M];
  qc = md->qc;
  ff1 = 0.0;
  memset ( grad, 0, N*sizeof(double) );  
  memset ( hess, 0, M*sizeof(double) );
  for ( i = l = m = p = 0;  i < nqkn;  i++ ) {
    ff2 = 0.0;
    memset ( gf2, 0, N*sizeof(double) );
    memset ( hf2, 0, M*sizeof(double) );
    for ( j = 0; j < nqkn; j++ ) {
      ff3 = 0.0;
      memset ( gf3, 0, N*sizeof(double) );
      memset ( hf3, 0, M*sizeof(double) );
      for ( k = 0;  k < nqkn;  k++, l++, m += N, p += M ) {
        ff3 += qc[k]*_f[l];
        pkn_AddMatrixMd ( 1, N, 0, gf3, 0, &_g[m], qc[k], 0, gf3 );
        pkn_AddMatrixMd ( 1, M, 0, hf3, 0, &_h[p], qc[k], 0, hf3 );
      }
      ff2 += qc[j]*ff3;
      pkn_AddMatrixMd ( 1, N, 0, gf2, 0, gf3, qc[j], 0, gf2 );
      pkn_AddMatrixMd ( 1, M, 0, hf2, 0, hf3, qc[j], 0, hf2 );
    }
    ff1 += qc[i]*ff2;
    pkn_AddMatrixMd ( 1, N, 0, grad, 0, gf2, qc[i], 0, grad );
    pkn_AddMatrixMd ( 1, M, 0, hess, 0, hf2, qc[i], 0, hess );
  }
  *func = ff1;

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*_mengerc_hessIntF*/

