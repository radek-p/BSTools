
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2013                                  */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

#include "pkvaria.h"
#include "pknum.h"

boolean pkn_SymMatFindEigenvaluesd ( int n, double *a, double *eigenval )
{
#define NU 1.0e-16
  void    *sp;
  double  *cd, *sd, *dd1, *dd2;
  int     i, j, k, p, q;
  double  c, s, mu1, mu, sigma;
  boolean newblock;

  sp = pkv_GetScratchMemTop ();
  if ( n < 2 )
    goto failure;
  cd = pkv_GetScratchMemd ( 4*n-5 );
  if ( !cd )
    goto failure;
  sd = &cd[n-1];
  dd1 = &sd[n-1];
  dd2 = &dd1[n-1];
        /* transforming the matrix to the three-diagonal form, */
        /* using Givens rotations & Stewart algorithm */
  for ( j = 0; j < n-1; j++ )
    for ( i = j+2; i < n; i++ ) {
          /* find the rotation sine and cosine */
      pkn_FindGivensRotationd ( a[pkn_LowerTrMatIndex(j+1,j)],
                                a[pkn_LowerTrMatIndex(i,j)], &c, &s );
          /* apply the rotation */
      a[pkn_LowerTrMatIndex(j+1,j)] =
          c*a[pkn_LowerTrMatIndex(j+1,j)] - s*a[pkn_LowerTrMatIndex(i,j)];
      for ( k = j+2; k < i; k++ )
        pkn_ApplyGivensRotationd ( c, s, &a[pkn_LowerTrMatIndex(k,j+1)],
                                   &a[pkn_LowerTrMatIndex(i,k)] );
      pkn_ApplySymGivensRotationd ( c, s, &a[pkn_LowerTrMatIndex(j+1,j+1)],
                                    &a[pkn_LowerTrMatIndex(i,j+1)],
                                    &a[pkn_LowerTrMatIndex(i,i)] );
      for ( k = i+1; k < n; k++ )
        pkn_ApplyGivensRotationd ( c, s, &a[pkn_LowerTrMatIndex(k,j+1)],
                                   &a[pkn_LowerTrMatIndex(k,i)] );
    }
        /* copy the coefficients of the three-diagonal matrix */
  eigenval[0] = a[0];
  for ( i = 1, j = 2;  i < n;  i++, j += i+1 ) {
    cd[i-1] = a[j-1];
    eigenval[i] = a[j];
  }
        /* zero the small co-diagonal coefficients */
  for ( i = 0; i < n-1; i++ )
    if ( fabs(cd[i]) <= NU*(fabs(eigenval[i])+fabs(eigenval[i+1])) )
      cd[i] = 0.0;
        /* the QR algorithm with deflation */
  newblock = true;
  p = q = 0;
  sigma = 0.0;
  do {
          /* find a block to iterate */
    if ( newblock ) {
      while ( p < n-2 && cd[p] == 0.0 )
        p++;
      if ( p == n-2 && !cd[n-1] )
        break;    /* the end */
      q = p+2;
      while ( q < n && cd[q-1] != 0.0 )
        q++;
    }
    if ( q-p > 2 ) {
          /* find the Wilkinson's translation */
      pkn_SolveSqEqd ( -0.5*(eigenval[q-2]+eigenval[q-1]),
                       eigenval[q-2]*eigenval[q-1]-cd[q-2]*cd[q-2], &mu, &mu1 );
      if ( fabs(mu-eigenval[q-1]) > fabs(mu1-eigenval[q-1]) )
        mu = mu1;
      sigma += mu;
          /* the QR iteration */
      for ( i = p; i < q; i++ )
        eigenval[i] -= mu;
      memcpy ( &dd1[p], &cd[p], (q-p-1)*sizeof(double) );
            /* A -> QR decomposition */
      for ( i = p; i < q-2; i++ ) {
        pkn_FindGivensRotationd ( eigenval[i], cd[i], &c, &s );
        eigenval[i] = c*eigenval[i] - s*cd[i];
        pkn_ApplyGivensRotationd ( c, s, &dd1[i], &eigenval[i+1] );
        dd2[i] = -s*dd1[i+1];  dd1[i+1] *= c;
        cd[i] = c;  sd[i] = s;
      }
      pkn_FindGivensRotationd ( eigenval[q-2], cd[q-2], &c, &s );
      eigenval[q-2] = c*eigenval[q-2] - s*cd[q-2];
      pkn_ApplyGivensRotationd ( c, s, &dd1[q-2], &eigenval[q-1] );
      cd[q-2] = c;  sd[q-2] = s;
            /* computation A = RQ */
      c = cd[p];  s = sd[p];
      for ( i = p; i < q-2; i++ ) {
        pkn_ApplyGivensRotationd ( c, s, &eigenval[i], &dd1[i] );
        cd[i] = -s*eigenval[i+1];  eigenval[i+1] *= c;
        c = cd[i+1];  s = sd[i+1];
        cd[i] = c*dd1[i] - s*dd2[i];
      }
      c = cd[q-2];  s = sd[q-2];
      pkn_ApplyGivensRotationd ( c, s, &eigenval[q-2], &dd1[q-2] );
      cd[q-2] = 0.0;
      pkn_ApplyGivensRotationd ( c, s, &cd[q-2], &eigenval[q-1] );
            /* deflation */
      newblock = false;
      for ( i = p; i < q-1; i++ )
        if ( fabs(cd[i]) <= NU*(fabs(eigenval[i]+sigma)+
                                fabs(eigenval[i+1]+sigma)) ) {
          cd[i] = 0.0;
          newblock = true;
        }
      if ( newblock ) {
        for ( i = p; i < q; i++ )
          eigenval[i] += sigma;
        sigma = 0.0;
      }
    }
    else {
          /* find the eigenvalues of a 2 x 2 block */
      pkn_SolveSqEqd ( -0.5*(eigenval[p]+eigenval[p+1]),
                       eigenval[p]*eigenval[p+1]-cd[p]*cd[p],
                       &eigenval[p], &eigenval[p+1] );
      cd[p] = 0.0;
      p += 2;
      newblock = true;
      sigma = 0.0;
    }
  } while ( p < n );

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
#undef NU
} /*pkn_SymMatFindEigenvaluesd*/

