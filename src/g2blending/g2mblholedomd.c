
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2010, 2013                            */
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
#include "eg2holed.h"
#include "g2blendingd.h"
#include "g2mblendingd.h"

#include "g2blprivated.h"
#include "g2mblprivated.h"

/*#define NOEXTDEB*/  /* for now */

GHoleDomaind *g2mbl_domaind[GH_MAX_K-3] =
  {NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL};
point2d *g2mbl_domainpatchd[GH_MAX_K-3] =
  {NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL};
double *g2mbl_patchmatrixd[GH_MAX_K-3] =
  {NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL};


void g2mbl_CleanupHoleDomainsd ( void )
{
  int i;

  for ( i = 0; i < GH_MAX_K-3; i++ ) {
    if ( g2mbl_domaind[i] ) {
      gh_DestroyDomaind ( g2mbl_domaind[i] );
      g2mbl_domaind[i] = NULL;
    }
    if ( g2mbl_domainpatchd[i] ) PKV_FREE ( g2mbl_domainpatchd[i] );
    if ( g2mbl_patchmatrixd[i] ) PKV_FREE ( g2mbl_patchmatrixd[i] );
  }
} /*g2mbl_CleanupHoleDomainsd*/

static int _g2mbl_HoleOptionProc ( GHoleDomaind *domain,
                                   int query, int qn,
                                   int *ndata, int **idata, double **fdata )
{
  switch ( query ) {
case G2HQUERY_QUADRATURE:
    return G2H_QUADRATURE_GAUSS_LEGENDRE;
default:
    return G2H_DEFAULT;
  }
} /*_g2mbl_HoleOptionProc*/

static void _g2mbl_DrawDiPatchd ( int n, int m, const point2d *cp )
{
    /* this is tricky - we do not need a hole with 4 sides, so we use */
    /* this pointer as a temporary storage */
  if ( !g2mbl_domainpatchd[1] ) {
    PKV_MALLOC ( g2mbl_domainpatchd[1], (n+1)*(m+2)*sizeof(point2d) );
    if ( g2mbl_domainpatchd[1] )
      memcpy ( g2mbl_domainpatchd[1], cp, (n+1)*(m+2)*sizeof(point2d) );
  }
} /*_g2mbl_DrawDiPatchd*/

boolean g2mbl_SetupHolePatchMatrixd ( int k )
{
  int    i, j;
  double s;

  if ( k < 2 || k == 4 || k > GH_MAX_K )
    return false;
  i = k-3;
  if ( g2mbl_patchmatrixd[i] )  /* it already has been computed */
    return true;
  if ( !g2mbl_domaind[i] ) {
    if ( i == 0 ) j = i;
    else j = i-1;
    g2mbl_domaind[i] = gh_CreateDomaind ( k, NULL, egh_eigendomcpd[j] );
    if ( !g2mbl_domaind[i] )
      return false;
    g2h_SetOptionProcd ( g2mbl_domaind[i], _g2mbl_HoleOptionProc );
  }
  PKV_MALLOC ( g2mbl_patchmatrixd[i],
               g2h_SymPatchMatrixSize ( k )*sizeof(double) );
  if ( !g2mbl_patchmatrixd[i] )
    return false;
#ifdef NOEXTDEB
  if ( !g2h_GetSymPatchMatrixd ( g2mbl_domaind[i],
                                 g2mbl_patchmatrixd[i] ) ) {
    PKV_FREE ( g2mbl_patchmatrixd[i] );
    return false;
  }
#else
  if ( !g2h_GetExtSymPatchMatrixd ( g2mbl_domaind[i],
                                    g2mbl_patchmatrixd[i] ) ) {
    PKV_FREE ( g2mbl_patchmatrixd[i] );
    return false;
  }
#endif
  g2h_DrawDiPatchesd ( g2mbl_domaind[i], _g2mbl_DrawDiPatchd );
  if ( !g2mbl_domainpatchd[1] )
    return false;
  g2mbl_domainpatchd[i] = g2mbl_domainpatchd[1];
  g2mbl_domainpatchd[1] = NULL;
  s = gh_HoleDomainAread ( g2mbl_domaind[i], true );
#ifdef __DEBUG
printf ( "%d %f\n", k, s );
#endif
  s = sqrt((double)k/s);
  for ( j = 0; j < (G2H_FINALDEG+1)*(G2H_FINALDEG+1); j++ )
    MultVector2d ( s, &g2mbl_domainpatchd[i][j], &g2mbl_domainpatchd[i][j] );
  return true;
} /*g2mbl_SetupHolePatchMatrix*/

/* ////////////////////////////////////////////////////////////////////////// */
boolean g2mbl_TabNid ( int hole_k, int nkn, double *qknots,
                       double *Nitab, double *Jac, boolean reparam )
{
  void     *sp;
  int      i, j, l, kn, nf;
  double   *fcp, p, dp[9], *df;
  point2d  *domcp, c;
  vector2d *dc, *adc;

  sp = pkv_GetScratchMemTop ();
  dc = pkv_GetScratchMem ( nkn*nkn*9*sizeof(vector2d) );
  if ( !dc )
    goto failure;
  domcp = g2mbl_domainpatchd[hole_k-3];
  if ( !domcp )
    goto failure;
  nf = 6*hole_k+1;
  if ( reparam ) {
    for ( i = kn = 0, adc = dc;  i < nkn;  i++ )
      for ( j = 0;  j < nkn;  j++, adc += 9, kn++ ) {
        if ( !mbs_BCHornerDer3P2d ( G2H_FINALDEG, G2H_FINALDEG, domcp,
                              qknots[i], qknots[j], &c, &adc[0], &adc[1], &adc[2],
                              &adc[3], &adc[4], &adc[5], &adc[6], &adc[7], &adc[8] ) )
          goto failure;
        Jac[kn] = fabs(adc[0].x*adc[1].y-adc[0].y*adc[1].x);
      }
    for ( l = 0, df = Nitab;  l < nf;  l++ ) {
      fcp = &g2mbl_patchmatrixd[hole_k-3][(G2H_FINALDEG+1)*(G2H_FINALDEG+1)*l];
      for ( i = 0, adc = dc;  i < nkn;  i++ )
        for ( j = 0;  j < nkn;  j++, df += 9, adc += 9 ) {
          if ( !mbs_BCHornerDer3P1d ( G2H_FINALDEG, G2H_FINALDEG, fcp,
                                qknots[i], qknots[j], &p, &dp[0], &dp[1], &dp[2],
                                &dp[3], &dp[4], &dp[5], &dp[6], &dp[7], &dp[8] ) )
            goto failure;
          if ( !pkn_Comp2iDerivatives3d ( adc[0].x, adc[0].y, adc[1].x, adc[1].y,
                  adc[2].x, adc[2].y, adc[3].x, adc[3].y, adc[4].x, adc[4].y,
                  adc[5].x, adc[5].y, adc[6].x, adc[6].y, adc[7].x, adc[7].y,
                  adc[8].x, adc[8].y, 1, &dp[0], &dp[1], &dp[2], &dp[3], &dp[4],
                  &dp[5], &dp[6], &dp[7], &dp[8], &df[0], &df[1], &df[2], &df[3],
                  &df[4], &df[5], &df[6], &df[7], &df[8] ) )
          goto failure;
        }
    }
  }
  else {
    for ( i = 0; i < nkn*nkn; i++ )
      Jac[i] = 1.0;
    for ( l = 0, df = Nitab;  l < nf;  l++ ) {
      fcp = &g2mbl_patchmatrixd[hole_k-3][(G2H_FINALDEG+1)*(G2H_FINALDEG+1)*l];
      for ( i = 0;  i < nkn;  i++ )
        for ( j = 0;  j < nkn;  j++, df += 9 ) {
          if ( !mbs_BCHornerDer3P1d ( G2H_FINALDEG, G2H_FINALDEG, fcp,
                                qknots[i], qknots[j], &p, &df[0], &df[1], &df[2],
                                &df[3], &df[4], &df[5], &df[6], &df[7], &df[8] ) )
            goto failure;
      }
    }
  }
  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*g2mbl_TabNid*/

void g2mbl_TabNijd ( int nf, int nkn, double *Nitab, double *Nijtab )
{
  int    ii, ij, i, l;
  double *ni, *nj, *nij;

  nkn *= nkn;
  for ( ii = 0; ii < nf; ii++ ) {
    ni = &Nitab[ii*nkn*9];
    for ( ij = 0; ij <= ii; ij++ ) {
      nj = &Nitab[ij*nkn*9];
      nij = &Nijtab[pkn_SymMatIndex(ii,ij)*nkn*9];
      for ( i = l = 0;  i < nkn;  i++, l += 9 ) {
        nij[l]   = 2.0*ni[l]*nj[l];
        nij[l+1] = ni[l]*nj[l+1]+ni[l+1]*nj[l];
        nij[l+2] = 2.0*ni[l+1]*nj[l+1];
        nij[l+3] = 2.0*(ni[l+2]*nj[l]+ni[l]*nj[l+2]);
        nij[l+4] = ni[l+2]*nj[l+1]+ni[l+1]*nj[l+2]+ni[l+3]*nj[l]+ni[l]*nj[l+3];
        nij[l+5] = 2.0*(ni[l+3]*nj[l+1]+ni[l+1]*nj[l+3]);
        nij[l+6] = 2.0*(ni[l+3]*nj[l]+ni[l]*nj[l+3]);
        nij[l+7] = ni[l+3]*nj[l+1]+ni[l+1]*nj[l+3]+ni[l+4]*nj[l]+ni[l]*nj[l+4];
        nij[l+8] = 2.0*(ni[l+4]*nj[l+1]+ni[l+1]*nj[l+4]);
      }
    }
  }
} /*g2mbl_TabNijd*/

void g2mbl_TabMijd ( int nf, int nkn, double *Nitab, double *Mijtab )
{
  int    ii, ij, i, j, l;
  double *ni, *nj, *mij;

  nkn *= nkn;
  for ( ii = 0; ii < nf-1; ii++ ) {
    ni = &Nitab[ii*nkn*9];
    for ( ij = ii+1; ij < nf; ij++ ) {
      nj = &Nitab[ij*nkn*9];
      mij = &Mijtab[pkn_SymMatIndex(ii,ij-1)*nkn*18];
      for ( i = j = l = 0;  i < nkn;  i++, j += 9, l += 18 ) {
        mij[l]    = ni[j+1]*nj[j+4]-ni[j+4]*nj[j+1]; /*0102*/
        mij[l+1]  = ni[j+1]*nj[j+8]-ni[j+8]*nj[j+1]; /*0103*/
        mij[l+2]  = ni[j+1]*nj[j+3]-ni[j+3]*nj[j+1]; /*0111*/
        mij[l+3]  = ni[j+1]*nj[j+7]-ni[j+7]*nj[j+1]; /*0112*/
        mij[l+4]  = ni[j+1]*nj[j+2]-ni[j+2]*nj[j+1]; /*0120*/
        mij[l+5]  = ni[j+1]*nj[j+6]-ni[j+6]*nj[j+1]; /*0121*/
        mij[l+6]  = ni[j+1]*nj[j+5]-ni[j+5]*nj[j+1]; /*0130*/
        mij[l+7]  = ni[j]*nj[j+1]-ni[j+1]*nj[j];     /*1001*/
        mij[l+8]  = ni[j]*nj[j+4]-ni[j+4]*nj[j];     /*1002*/
        mij[l+9]  = ni[j]*nj[j+8]-ni[j+8]*nj[j];     /*1003*/
        mij[l+10] = ni[j]*nj[j+3]-ni[j+3]*nj[j];     /*1011*/
        mij[l+11] = ni[j]*nj[j+7]-ni[j+7]*nj[j];     /*1012*/
        mij[l+12] = ni[j]*nj[j+2]-ni[j+2]*nj[j];     /*1020*/
        mij[l+13] = ni[j]*nj[j+6]-ni[j+6]*nj[j];     /*1021*/
        mij[l+14] = ni[j]*nj[j+5]-ni[j+5]*nj[j];     /*1030*/
        mij[l+15] = ni[j+3]*nj[j+4]-ni[j+4]*nj[j+3]; /*1102*/
        mij[l+16] = ni[j+2]*nj[j+4]-ni[j+4]*nj[j+2]; /*2002*/
        mij[l+17] = ni[j+2]*nj[j+3]-ni[j+3]*nj[j+2]; /*2011*/
      }
    }
  }
} /*g2mbl_TabMijd*/

