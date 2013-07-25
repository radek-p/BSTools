
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2010                                  */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "pkvaria.h"
#include "pknum.h"
#include "bsmesh.h"
#include "bsmprivate.h"

/* ////////////////////////////////////////////////////////////////////////// */
boolean bsm_FindRegularSubnets ( int nv, BSMvertex *mv, int *mvhei,
                                 int nhe, BSMhalfedge *mhe,
                                 int nfac, BSMfacet *mfac, int *mfhei,
                                 int d, void *usrptr,
                                 void (*output)( int d, int *vertnum, int *mtab,
                                                 void *usrptr ) )
{
#define MI(i,j) ((i)*(d+d-1)+(j))
  void *sp;
  int  *vertnum, *mtab;
  int  i, j, k, l, v, he, f, fhe;

  sp = pkv_GetScratchMemTop ();
  if ( d < 2 )
    goto failure;
  vertnum = pkv_GetScratchMemi ( d*d );
  mtab = pkv_GetScratchMemi ( (2*d-1)*(2*d-1) );
  if ( !vertnum || !mtab )
    goto failure;

  if ( d % 2 ) {  /* d is odd - begin from vertices */
    for ( i = 0; i < nv; i++ ) {
      if ( mv[i].degree == 4 ) {
        fhe = mv[i].firsthalfedge;
        if ( mhe[mvhei[fhe+3]].otherhalf < 0 )  /* must be inner */
          goto nextvert;
                  /* get the central vertex */
        mtab[MI(d-1,d-1)] = i;
        mtab[MI(d-1,d)] = he = mvhei[fhe];
        mtab[MI(d,d)] = f = mhe[he].facetnum;
        if ( mfac[f].degree != 4 )
          goto nextvert;
                  /* get the edges and vertices of the facet */
        fhe = mfac[f].firsthalfedge;
        for ( l = 0; l < 4; l++ )
          if ( mfhei[fhe+l] == he )
            break;
        if ( l == 4 )
          goto failure;
        mtab[MI(d,d+1)] = he = mfhei[fhe+(l+1) % 4];
        mtab[MI(d-1,d+1)] = mhe[he].v0;
        mtab[MI(d+1,d)] = he = mfhei[fhe+(l+2) % 4];
        mtab[MI(d+1,d+1)] = mhe[he].v0;
        mtab[MI(d,d-1)] = he = mfhei[fhe+(l+3) % 4];
        mtab[MI(d+1,d-1)] = mhe[he].v0;
                  /* go north */
        for ( k = d+2; k < 2*d-1; k += 2 ) {
          he = mhe[mtab[MI(d,k-1)]].otherhalf;
          if ( mv[mtab[MI(d-1,k-1)]].degree != 4 || he < 0 )
            goto nextvert;
          f = mhe[he].facetnum;
          if ( mfac[f].degree != 4 )
            goto nextvert;
          mtab[MI(d,k)] = f;
          fhe = mfac[f].firsthalfedge;
          for ( l = 0; l < 4; l++ )
            if ( mfhei[fhe+l] == he )
              break;
          if ( l == 4 )
            goto failure;
          mtab[MI(d-1,k)] = mfhei[fhe+(l+1) % 4];
          mtab[MI(d,k+1)] = he = mfhei[fhe+(l+2) % 4];
          mtab[MI(d-1,k+1)] = mhe[he].v0;
          mtab[MI(d+1,k)] = he = mfhei[fhe+(l+3) % 4];
          mtab[MI(d+1,k+1)] = mhe[he].v0;
        }
                  /* go south */
        for ( k = d-2; k > 0; k -= 2 ) {
          he = mhe[mtab[MI(d,k+1)]].otherhalf;
          if ( mv[mtab[MI(d-1,k+1)]].degree != 4 || he < 0 )
            goto nextvert;
          f = mhe[he].facetnum;
          if ( mfac[f].degree != 4 )
            goto nextvert;
          mtab[MI(d,k)] = f;
          fhe = mfac[f].firsthalfedge;
          for ( l = 0; l < 4; l++ )
            if ( mfhei[fhe+l] == he )
              break;
          if ( l == 4 )
            goto failure;
          mtab[MI(d+1,k)] = mfhei[fhe+(l+1) % 4];
          mtab[MI(d,k-1)] = he = mfhei[fhe+(l+2) % 4];
          mtab[MI(d+1,k-1)] = mhe[he].v0;
          mtab[MI(d-1,k)] = he = mfhei[fhe+(l+3) % 4];
          mtab[MI(d-1,k-1)] = mhe[he].v0;
        }
        for ( k = 1; k < 2*d-2; k += 2 ) {
                  /* go west */
          for ( j = d-2; j > 0; j -= 2 ) {
            he = mhe[mtab[MI(j+1,k)]].otherhalf;
            if ( he < 0 )
              goto nextvert;
            f = mhe[he].facetnum;
            if ( mfac[f].degree != 4 )
              goto nextvert;
            mtab[MI(j,k)] = f;
            fhe = mfac[f].firsthalfedge;
            for ( l = 0; l < 4; l++ )
              if ( mfhei[fhe+l] == he )
                break;
            if ( l == 4 )
              goto failure;
            mtab[MI(j,k-1)] = mfhei[fhe+(l+1) % 4];
            mtab[MI(j-1,k)] = he = mfhei[fhe+(l+2) % 4];
            mtab[MI(j-1,k-1)] = mhe[he].v0;
            mtab[MI(j,k+1)] = he = mfhei[fhe+(l+3) % 4];
            mtab[MI(j-1,k+1)] = mhe[he].v0;
          }
                  /* go east */
          for ( j = d+2; j < 2*d-1; j += 2 ) {
            he = mhe[mtab[MI(j-1,k)]].otherhalf;
            if ( he < 0 )
              goto nextvert;
            f = mhe[he].facetnum;
            if ( mfac[f].degree != 4 )
              goto nextvert;
            mtab[MI(j,k)] = f;
            fhe = mfac[f].firsthalfedge;
            for ( l = 0; l < 4; l++ )
              if ( mfhei[fhe+l] == he )
                break;
            if ( l == 4 )
              goto failure;
            mtab[MI(j,k-1)] = mfhei[fhe+(l+1) % 4];
            mtab[MI(j+1,k)] = he = mfhei[fhe+(l+2) % 4];
            mtab[MI(j+1,k+1)] = mhe[he].v0;
            mtab[MI(j,k-1)] = he = mfhei[fhe+(l+3) % 4];
            mtab[MI(j+1,k-1)] = mhe[he].v0;
          }
        }
                  /* verify the inner vertices */
        for ( j = 2; j < 2*d-3; j += 2 )
          for ( k = 2; k < 2*d-3; k += 2 ) {
            v = mtab[MI(j,k)];
            if ( mv[v].degree != 4 )
              goto nextvert;
            fhe = mv[v].firsthalfedge;
            if ( mhe[mvhei[fhe+3]].otherhalf < 0 )  /* must be inner */
              goto nextvert;
          }
                  /* output the subnet */
        for ( j = l = 0;  j < 2*d-1;  j += 2 )
          for ( k = 0;  k < 2*d-1;  k += 2, l++ )
            vertnum[l] = mtab[MI(j,k)];
        output ( d, vertnum, mtab, usrptr );
      }
nextvert:
      ;
    }
  }
  else {          /* d is even - begin from facets */
    for ( i = 0; i < nfac; i++ ) {
      if ( mfac[i].degree == 4 ) {
        fhe = mfac[i].firsthalfedge;
                  /* get the central facet */
        mtab[MI(d-1,d-1)] = i;
        mtab[MI(d-1,d)] = he = mfhei[fhe];
        mtab[MI(d-2,d)] = mhe[he].v0;
        mtab[MI(d,d-1)] = he = mfhei[fhe+1];
        mtab[MI(d,d)] = mhe[he].v0;
        mtab[MI(d-1,d-2)] = he = mfhei[fhe+2];
        mtab[MI(d,d-2)] = mhe[he].v0;
        mtab[MI(d-2,d-1)] = he = mfhei[fhe+3];
        mtab[MI(d-2,d-2)] = mhe[he].v0;
                  /* go north */
        for ( k = d+1; k < 2*d-1; k += 2 ) {
          he = mhe[mtab[MI(d-1,k-1)]].otherhalf;
          if ( mv[mtab[MI(d-2,k-1)]].degree != 4 || mv[mtab[MI(d,k-1)]].degree != 4 ||
               he < 0 )
            goto nextfacet;
          f = mhe[he].facetnum;
          if ( mfac[f].degree != 4 )
            goto nextfacet;
          mtab[MI(d-1,k)] = f;
          fhe = mfac[f].firsthalfedge;
          for ( l = 0; l < 4; l++ )
            if ( mfhei[fhe+l] == he )
              break;
          if ( l == 4 )
            goto failure;
          mtab[MI(d-2,k)] = mfhei[fhe+(l+1) % 4];
          mtab[MI(d-1,k+1)] = he = mfhei[fhe+(l+2) % 4];
          mtab[MI(d-2,k+1)] = mhe[he].v0;
          mtab[MI(d,k)] = he = mfhei[fhe+(l+3) %4];
          mtab[MI(d,k+1)] = mhe[he].v0;
        }
                  /* go south */
        for ( k = d-3; k > 0; k -= 2 ) {
          he = mhe[mtab[MI(d-1,k+1)]].otherhalf;
          if ( mv[mtab[MI(d-2,k+1)]].degree != 4 || mv[mtab[MI(d,k+1)]].degree != 4 ||
               he < 0 )
            goto nextfacet;
          f = mhe[he].facetnum;
          if ( mfac[f].degree != 4 )
            goto nextfacet;
          mtab[MI(d-1,k)] = f;
          fhe = mfac[f].firsthalfedge;
          for ( l = 0; l < 4; l++ )
            if ( mfhei[fhe+l] == he )
              break;
          if ( l == 4 )
            goto failure;
          mtab[MI(d,k)] = mfhei[fhe+(l+1) % 4];
          mtab[MI(d-1,k-1)] = he = mfhei[fhe+(l+2) % 4];
          mtab[MI(d,k-1)] = mhe[he].v0;
          mtab[MI(d-2,k)] = he = mfhei[fhe+(l+3) % 4];
          mtab[MI(d-2,k-1)] = mhe[he].v0;
        }
        for ( k = 1; k < 2*d-2; k += 2 ) {
                  /* go west */
          for ( j = d-3; j > 0; j -= 2 ) {
            he = mhe[mtab[MI(j+1,k)]].otherhalf;
            if ( he < 0 )
              goto nextfacet;
            f = mhe[he].facetnum;
            if ( mfac[f].degree != 4 )
              goto nextfacet;
            mtab[MI(j,k)] = f;
            fhe = mfac[f].firsthalfedge;
            for ( l = 0; l < 4; l++ )
              if ( mfhei[fhe+l] == he )
                break;
            if ( l == 4 )
              goto failure;
            mtab[MI(j,k-1)] = mfhei[fhe+(l+1) % 4];
            mtab[MI(j-1,k)] = he = mfhei[fhe+(l+2) % 4];
            mtab[MI(j-1,k-1)] = mhe[he].v0;
            mtab[MI(j,k+1)] = he = mfhei[fhe+(l+3) % 4];
            mtab[MI(j-1,k+1)] = mhe[he].v0;
          }
                  /* go east */
          for ( j = d+1; j < 2*d-1; j += 2 ) {
            he = mhe[mtab[MI(j-1,k)]].otherhalf;
            if ( he < 0 )
              goto nextfacet;
            f = mhe[he].facetnum;
            if ( mfac[f].degree != 4 )
              goto nextfacet;
            mtab[MI(j,k)] = f;
            fhe = mfac[f].firsthalfedge;
            for ( l = 0; l < 4; l++ )
              if ( mfhei[fhe+l] == he )
                break;
            if ( l == 4 )
              goto failure;
            mtab[MI(j,k+1)] = mfhei[fhe+(l+1) % 4];
            mtab[MI(j+1,k)] = he = mfhei[fhe+(l+2) % 4];
            mtab[MI(j+1,k+1)] = mhe[he].v0;
            mtab[MI(j,k-1)] = he = mfhei[fhe+(l+3) % 4];
            mtab[MI(j+1,k-1)] = mhe[he].v0;
          }
        }
                  /* verify the inner vertices */
        for ( j = 2; j < 2*d-3; j += 2 )
          for ( k = 2; k < 2*d-3; k += 2 ) {
            v = mtab[MI(j,k)];
            if ( mv[v].degree != 4 )
              goto nextfacet;
            fhe = mv[v].firsthalfedge;
            if ( mhe[mvhei[fhe+3]].otherhalf < 0 )  /* must be inner */
              goto nextfacet;
          }
                  /* output the subnet */
        for ( j = l = 0;  j < 2*d-1;  j += 2 )
          for ( k = 0;  k < 2*d-1;  k += 2, l ++ )
            vertnum[l] = mtab[MI(j,k)];
        output ( d, vertnum, mtab, usrptr );
      }
nextfacet:
      ;
    }
  }

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
#undef MI
} /*bsm_FindRegularSubnets*/ 

