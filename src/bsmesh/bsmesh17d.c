
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2010, 2012                            */
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
int bsm_ContractEdged ( int spdimen,
                        int inv, BSMvertex *imv, int *imvhei, double *iptc,
                        int inhe, BSMhalfedge *imhe,
                        int infac, BSMfacet *imfac, int *imfhei,
                        int nche,
                        int *onv, BSMvertex *omv, int *omvhei, double *optc,
                        int *onhe, BSMhalfedge *omhe,
                        int *onfac, BSMfacet *omfac, int *omfhei )
{
  void *sp;
  int  v0, v1, v2, v3, vn, oth, fac1, fac2;
  int  df1, df2, df, dv, le0, le1, le2, le3;
  int  i, j, k, l, m, n, o, fhe, fhe1, dv1;
  int  *newhenum;

  sp = pkv_GetScratchMemTop ();
  if ( nche < 0 || nche >= inhe )
    goto failure;
  le1 = v2 = 0;
  newhenum = pkv_GetScratchMemi ( inhe );
  if ( !newhenum )
    goto failure;
  memset ( newhenum, 0, inhe*sizeof(int) );
        /* copy the vertices data and compute the edge midpoint coordinates */
  *onv = inv-1;
  v0 = imhe[nche].v0;
  v1 = imhe[nche].v1;
  if ( v1 > 0 ) {
    memcpy ( omv, imv, v1*sizeof(BSMvertex) );
    memcpy ( optc, iptc, v1*spdimen*sizeof(double) );
  }
  if ( v1 < inv-1 ) {
    memcpy ( &omv[v1], &imv[v1+1], (inv-1-v1)*sizeof(BSMvertex) );
    memcpy ( &optc[v1*spdimen], &iptc[(v1+1)*spdimen],
             (inv-1-v1)*spdimen*sizeof(double) );
  }
  vn = v0 < v1 ? v0 : v0-1;
  pkn_AddMatrixd ( 1, spdimen, 0, &iptc[v0*spdimen], 0, &iptc[v1*spdimen],
                   0, &optc[vn*spdimen] );
  pkn_MultMatrixNumd ( 1, spdimen, 0, &optc[vn*spdimen], 0.5,
                       0, &optc[vn*spdimen] );
  
  if ( imhe[nche].otherhalf == -1 ) {  /* a boundary edge to be contracted */
    fac1 = imhe[nche].facetnum;
    df1 = imfac[fac1].degree;
    if ( df1 > 3 ) { /* all facets remain */
      *onhe = inhe-1;
      *onfac = infac;
                     /* setup new halfedge numbers */
      newhenum[nche] = -1;
      for ( i = j = 0; i < inhe; i++ )
        if ( newhenum[i] >= 0 )
          newhenum[i] = j++;
                     /* setup output facets */
      for ( i = k = 0; i < infac; i++ ) {
        df = imfac[i].degree;
        fhe = imfac[i].firsthalfedge;
        omfac[i].firsthalfedge = k;
        if ( i != fac1 ) {  /* just copy the data */
          omfac[i].degree = df;
          for ( j = 0; j < df; j++ )
            omfhei[k+j] = newhenum[imfhei[fhe+j]];
        }
        else {
          omfac[i].degree = df-1;
          for ( j = l = 0; j < df; j++ )
            if ( imfhei[fhe+j] != nche ) {
              omfhei[k+l] = newhenum[imfhei[fhe+j]];
              l ++;
            }
        }
        k += omfac[i].degree;
      }
            /* setup the output vertices */
      for ( i = j = k = 0; i < inv; i++ ) {
        if ( i == v0 ) {
          omv[j].firsthalfedge = k;
          omv[j].degree = imv[v0].degree + imv[v1].degree - 1;
          fhe = imv[v0].firsthalfedge;
          dv = imv[v0].degree-1;
          for ( l = 0; l < dv; l++ )
            omvhei[k+l] = newhenum[imvhei[fhe+l]];
          fhe = imv[v1].firsthalfedge;
          df = imv[v1].degree;
          for ( l = 0; l < df; l++ )
            omvhei[k+dv+l] = newhenum[imvhei[fhe+l]];
          k += omv[j].degree;
          j ++;
        }
        else if ( i != v1 ) {
          fhe = imv[i].firsthalfedge;
          omv[j].firsthalfedge = k;
          omv[j].degree = dv = imv[i].degree;
          for ( l = 0; l < dv; l++ )
            omvhei[k+l] = newhenum[imvhei[fhe+l]];
          k += dv;
          j ++;
        }
      }
            /* setup the output halfedges */
      for ( i = j = 0; i < inhe; i++ )
        if ( newhenum[i] >= 0 ) {
          omhe[j].facetnum = imhe[i].facetnum;
          omhe[j].otherhalf =
                 imhe[i].otherhalf < 0 ? -1 : newhenum[imhe[i].otherhalf];
          if ( imhe[i].v0 == v0 || imhe[i].v0 == v1 )
            omhe[j].v0 = v0 > v1 ? v0-1 : v0;
          else
            omhe[j].v0 = imhe[i].v0 > v1 ? imhe[i].v0-1 : imhe[i].v0;
          if ( imhe[i].v1 == v0 || imhe[i].v1 == v1 )
            omhe[j].v1 = v0 > v1 ? v0-1 : v0;
          else
            omhe[j].v1 = imhe[i].v1 > v1 ? imhe[i].v1-1 : imhe[i].v1;
          j ++;
        }
    }
    else {  /* the facet adjacent to the contracted edge is to be removed */
      *onhe = inhe-3;
      *onfac = infac-1;
              /* setup new halfedge numbers */
      fhe = imfac[fac1].firsthalfedge;
      for ( i = 0; i < df1; i++ )
        newhenum[imfhei[fhe+i]] = -1;
      for ( i = j = 0; i < inhe; i++ )
        if ( newhenum[i] >= 0 )
          newhenum[i] = j++;
             /* setup the output facets */
      for ( i = j = k = 0; i < infac; i++ )
        if ( i != fac1 ) {
          df = imfac[i].degree;
          fhe = imfac[i].firsthalfedge;
          omfac[j].degree = imfac[i].degree;
          omfac[j].firsthalfedge = k;
          for ( l = m = 0; l < df; l++ )
            if ( newhenum[imfhei[fhe+l]] >= 0 ) {
              omfhei[k+m] = newhenum[imfhei[fhe+l]];
              m ++;
            }
          j ++;
          k += m;
        }
              /* setup the output vertices */
      fhe = imfac[fac1].firsthalfedge;
      le0 = imfhei[fhe];
      le1 = imfhei[fhe+1];
      if ( le0 == nche ) le0 = imfhei[fhe+2];
      else if ( le1 == nche ) le1 = imfhei[fhe+2];
      if ( imhe[le0].v0 != v0 && imhe[le0].v0 != v1 )
        v2 = imhe[le0].v0;
      else
        v2 = imhe[le0].v1;
      for ( i = j = k = 0; i < inv; i++ ) {
        if ( i == v0 ) {
          omv[j].firsthalfedge = k;
          omv[j].degree = dv = imv[v0].degree + imv[v1].degree - 2;
          if ( dv > 127 )
            goto failure;
          fhe = imv[v0].firsthalfedge;
          dv = imv[v0].degree-1;
          for ( l = 0; l < dv; l++ )
            omvhei[k+l] = newhenum[imvhei[fhe+l]];
          fhe = imv[v1].firsthalfedge;
          df = imv[v1].degree-1;
          for ( l = 0; l < df; l++ )
            omvhei[k+dv+l] = newhenum[imvhei[fhe+l+1]];
          k += omv[j].degree;
          j ++;
        }
        else if ( i == v2 ) {
          fhe = imv[i].firsthalfedge;
          omv[j].degree = dv = imv[i].degree-1;
          if ( (dv < 3 && imhe[imvhei[fhe+dv]].otherhalf >= 0) || dv < 1 )
            goto failure;
          for ( l = m = 0; l <= dv; l++ )
            if ( newhenum[imvhei[fhe+l]] >= 0 ) {
              omvhei[k+m] = newhenum[imvhei[fhe+l]];
              m ++;
            }
          k += m;
          j ++;
        }
        else if ( i != v1 ) {
          fhe = imv[i].firsthalfedge;
          omv[j].firsthalfedge = k;
          omv[j].degree = dv = imv[i].degree;
          for ( l = 0; l < dv; l++ )
            omvhei[k+l] = newhenum[imvhei[fhe+l]];
          k += dv;
          j ++;
        }
      }
              /* setup the output halfedges */
      for ( i = j = 0; i < inhe; i++ )
        if ( newhenum[i] >= 0 ) {
          omhe[j].facetnum =
              imhe[i].facetnum < fac1 ? imhe[i].facetnum : imhe[i].facetnum-1;
          omhe[j].otherhalf =
              imhe[i].otherhalf < 0 ? -1 : newhenum[imhe[i].otherhalf];
          if ( imhe[i].v0 == v0 || imhe[i].v0 == v1 )
            omhe[j].v0 = v0 > v1 ? v0-1 : v0;
          else
            omhe[j].v0 = imhe[i].v0 > v1 ? imhe[i].v0-1 : imhe[i].v0;
          if ( imhe[i].v1 == v0 || imhe[i].v1 == v1 )
            omhe[j].v1 = v0 > v1 ? v0-1 : v0;
          else
            omhe[j].v1 = imhe[i].v1 > v1 ? imhe[i].v1-1 : imhe[i].v1;
          j ++;
        }
             /* link the pair of halfedges neighbouring the removed facet */
      le0 = imhe[le0].otherhalf;
      le1 = imhe[le1].otherhalf;
      if ( le0 >= 0 && le1 >= 0 ) {
        le0 = newhenum[le0];
        le1 = newhenum[le1];
        omhe[le0].otherhalf = le1;
        omhe[le1].otherhalf = le0;
      }
    }
  }
  else {  /* an inner edge to be contracted */
memset ( omv, 0, *onv*sizeof(BSMvertex) );
memset ( omvhei, 0, *onhe*sizeof(int) );
memset ( omhe, 0, *onhe*sizeof(BSMhalfedge) );
memset ( omfac, 0, *onfac*sizeof(BSMfacet) );
memset ( omfhei, 0, *onhe*sizeof(int) );

    oth = imhe[nche].otherhalf;
    fac1 = imhe[nche].facetnum;
    df1 = imfac[fac1].degree;
    fac2 = imhe[oth].facetnum;
    df2 = imfac[fac2].degree;
    if ( df1 < 3 || df2 < 3 )
      goto failure;
            /* cannot remove an inner edge, whose both vertices are */
            /* at the boundary */
    if ( imhe[imvhei[imv[v0].firsthalfedge+imv[v0].degree-1]].otherhalf == -1 &&
         imhe[imvhei[imv[v1].firsthalfedge+imv[v1].degree-1]].otherhalf == -1 )
      goto failure;
    if ( df1 > 3 && df2 > 3 ) {
              /* both adjacent facets remain */
                /* setup new halfedge numbers */
      newhenum[nche] = newhenum[oth] = -1;
      for ( i = j = 0; i < inhe; i++ )
        if ( newhenum[i] >= 0 )
          newhenum[i] = j++;
                /* setup output facets */
      for ( i = k = 0; i < infac; i++ ) {
        df = imfac[i].degree;
        fhe = imfac[i].firsthalfedge;
        omfac[i].firsthalfedge = k;
        if ( i == fac1 || i == fac2 ) {
          omfac[i].degree = df-1;
          for ( j = l = 0; j < df; j++ )
            if ( newhenum[imfhei[fhe+j]] >= 0 ) {
              omfhei[k+l] = newhenum[imfhei[fhe+j]];
              l ++;
            }
        }
        else {  /* just copy data */
          omfac[i].degree = df;
          for ( j = 0; j < df; j++ )
            omfhei[k+j] = newhenum[imfhei[fhe+j]];
        }
        k += omfac[i].degree;
      }
      *onfac = infac;
                /* setup the output vertices */
      for ( i = j = k = 0;  i < inv;  i++ )
        if ( i == v0 ) {
          omv[j].degree = dv = imv[v0].degree + imv[v1].degree - 2;
          if ( dv > 127 )
            goto failure;
          omv[j].firsthalfedge = k;
                    /* 3 cases possible: v0 at boundary, v1 at boundary, */
                    /* or none at boundary */
          fhe = imv[i].firsthalfedge;
          dv = imv[i].degree;
          fhe1 = imv[v1].firsthalfedge;
          dv1 = imv[v1].degree;
          if ( imhe[imvhei[fhe1+dv1-1]].otherhalf == -1 ) {
                      /* v1 at the boundary */
            for ( l = m = 0; l < dv1; l++ )
              if ( newhenum[imvhei[fhe1+l]] >= 0 ) {
                omvhei[k+m] = newhenum[imvhei[fhe1+l]];
                m ++;
              }
              else {
                for ( n = 0; n < dv; n++ )
                  if ( newhenum[imvhei[fhe+n]] == -1 )
                    break;
                if ( n >= dv )
                  goto failure;
                for ( o = 0, n = (n+1) % dv; o < dv-1;  o++, n = (n+1) % dv, m++ )
                  omvhei[k+m] = newhenum[imvhei[fhe+n]];
              }
          }
          else {
                      /* v1 inner, v0 inner or at the boundary  */
            for ( l = m = 0; l < dv; l++ )
              if ( newhenum[imvhei[fhe+l]] >= 0 ) {
                omvhei[k+m] = newhenum[imvhei[fhe+l]];
                m ++;
              }
              else {
                for ( n = 0; n < dv1; n++ )
                  if ( newhenum[imvhei[fhe1+n]] == -1 )
                    break;
                if ( n >= dv1 )
                  goto failure;
                for ( o = 0, n = (n+1) % dv1;  o < dv1-1;  o++, n = (n+1) % dv1, m++ )
                  omvhei[k+m] = newhenum[imvhei[fhe1+n]];
              }
          }
          k += omv[j].degree;
          j ++;
        }
        else if ( i != v1 ) {
          omv[j].degree = dv = imv[i].degree;
          omv[j].firsthalfedge = k;
          fhe = imv[i].firsthalfedge;
          for ( l = 0; l < dv; l++ )
            omvhei[k+l] = newhenum[imvhei[fhe+l]];
          k += dv;
          j ++;
        }
      *onv = inv-1;
                /* setup the output halfedges */
      for ( i = k = 0; i < inhe; i++ )
        if ( newhenum[i] >= 0 ) {
          omhe[k].v0 = imhe[i].v0 < v1 ? imhe[i].v0 :
                       (imhe[i].v0 == v1 ? vn : imhe[i].v0-1);
          omhe[k].v1 = imhe[i].v1 < v1 ? imhe[i].v1 :
                       (imhe[i].v1 == v1 ? vn : imhe[i].v1-1);
          omhe[k].facetnum = imhe[i].facetnum;
          omhe[k].otherhalf = imhe[i].otherhalf >= 0 ?
                                newhenum[imhe[i].otherhalf] : -1;
          k++;
        }
      *onhe = k;
    }
    else if ( df1 == 3 && df2 == 3 ) {
              /* both adjacent facets to be removed */
                /* setup new halfedge numbers */
      j = imfac[fac1].firsthalfedge;
      k = imfac[fac2].firsthalfedge;
      for ( i = 0; i < 3; i++ )
        newhenum[imfhei[j+i]] = newhenum[imfhei[k+i]] = -1;
      if ( imhe[imfhei[j]].v0 == v0 )      { le0 = imfhei[j+2];  le1 = imfhei[j+1]; }
      else if ( imhe[imfhei[j]].v0 == v1 ) { le0 = imfhei[j+1];  le1 = imfhei[j];   }
      else                                 { le0 = imfhei[j];    le1 = imfhei[j+2]; }
      v2 = imhe[le0].v0;
      if ( imhe[imfhei[k]].v0 == v0 )      { le2 = imfhei[k+1];  le3 = imfhei[k];   }
      else if ( imhe[imfhei[k]].v0 == v1 ) { le2 = imfhei[k+2];  le3 = imfhei[k+1]; }
      else                                 { le2 = imfhei[k];    le3 = imfhei[k+2]; }
      v3 = imhe[le2].v0;
      for ( i = j = 0; i < inhe; i++ )
        if ( newhenum[i] >= 0 )
          newhenum[i] = j++;
      *onhe = j;
                /* setup output facets */
      for ( i = j = k = 0; i < infac; i++ )
        if ( i != fac1 && i != fac2 ) {
          omfac[j].degree = df = imfac[i].degree;
          omfac[j].firsthalfedge = k;
          fhe = imfac[i].firsthalfedge;
          for ( l = 0; l < df; l++ )
            omfhei[k+l] = newhenum[imfhei[fhe+l]];
          k += df;
          j ++;
        }
      *onfac = j;
                /* setup the output vertices */
      for ( i = j = k = 0; i < inv; i++ ) {
        fhe = imv[i].firsthalfedge;
        dv = imv[i].degree;
        if ( i == v0 ) {
          omv[j].degree = imv[v0].degree + imv[v1].degree - 4;
          if ( dv > 127 )
            goto failure;
          omv[j].firsthalfedge = k;
          fhe1 = imv[v1].firsthalfedge;
          dv1 = imv[v1].degree;
          if ( imhe[imvhei[fhe1+dv1-1]].otherhalf == -1 ) {
                /* v1 at the boundary */
            for ( l = m = 0; l < dv1; l++ )
              if ( newhenum[imvhei[fhe1+l]] >= 0 ) {
                omvhei[k+m] = newhenum[imvhei[fhe1+l]];
                m ++;
              }
              else if ( imvhei[fhe1+l] == oth ) {
                for ( n = 0; n < dv; n++ )
                  if ( imvhei[fhe+n] == nche )
                    break;
                if ( n >= dv )
                  goto failure;
                for ( o = 0, n = (n+2) % dv;  o < dv-2;  o++, n = (n+1) % dv, m++ )
                  omvhei[k+m] = newhenum[imvhei[fhe+n]];
              }
          }
          else {
                /* v1 inner, v0 inner or at the boundary */
            for ( l = m = 0; l < dv; l++ )
              if ( newhenum[imvhei[fhe+l]] >= 0 ) {
                omvhei[k+m] = newhenum[imvhei[fhe+l]];
                m ++;
              }
              else if ( imvhei[fhe+l] == nche ) {
                for ( n = 0; n < dv1; n++ )
                  if ( imvhei[fhe1+n] == oth )
                    break;
                  if ( n >= dv1 )
                    goto failure;
                  for ( o = 0, n = (n+2) % dv1;  o < dv1-2;  o++, n = (n+1) % dv1, m++ )
                    omvhei[k+m] = newhenum[imvhei[fhe1+n]];
              }
          }
          k += omv[j].degree;
          j ++;
        }
        else if ( i == v2 || i == v3 ) {
          omv[j].degree = dv-1;
          omv[j].firsthalfedge = k;
          for ( l = m = 0; l < dv; l++ )
            if ( newhenum[imvhei[fhe+l]] >= 0 ) {
              omvhei[k+m] = newhenum[imvhei[fhe+l]];
              m ++;
            }
          k += dv-1;
          j ++;
        }
        else if ( i != v1 ) {
          omv[j].degree = dv;
          omv[j].firsthalfedge = k;
          for ( l = 0; l < dv; l++ )
            omvhei[k+l] = newhenum[imvhei[fhe+l]];
          k += dv;
          j ++;
        }
      }
      *onv = j;
                /* setup the output halfedges */
      for ( i = 0; i < inhe; i++ ) {
        j = newhenum[i];
        if ( j >= 0 ) {
          omhe[j] = imhe[i];
          if ( imhe[i].v0 == v1 ) omhe[j].v0 = vn;
          else if ( imhe[i].v0 > v1 ) omhe[j].v0 --;
          if ( imhe[i].v1 == v1 ) omhe[j].v1 = vn;
          else if ( imhe[i].v1 > v1 ) omhe[j].v1 --;
          if ( imhe[i].facetnum > fac1 ) omhe[j].facetnum --;
          if ( imhe[i].facetnum > fac2 ) omhe[j].facetnum --;
          omhe[j].otherhalf = imhe[i].otherhalf >= 0 ?
                              newhenum[imhe[i].otherhalf] : -1;
        }
      }
      le0 = imhe[le0].otherhalf;  if ( le0 >= 0 ) le0 = newhenum[le0];
      le1 = imhe[le1].otherhalf;  if ( le1 >= 0 ) le1 = newhenum[le1];
      le2 = imhe[le2].otherhalf;  if ( le2 >= 0 ) le2 = newhenum[le2];
      le3 = imhe[le3].otherhalf;  if ( le3 >= 0 ) le3 = newhenum[le3];
      if ( le0 >= 0 ) omhe[le0].otherhalf = le1;
      if ( le1 >= 0 ) omhe[le1].otherhalf = le0;
      if ( le2 >= 0 ) omhe[le2].otherhalf = le3;
      if ( le3 >= 0 ) omhe[le3].otherhalf = le2;
    }
    else {
              /* one facet adjacent to the edge has to be removed */
      if ( df1 == 3 ) {  /* exchange, so as to have the facet fac2 to remove */
        i = nche;  nche = oth;   oth = i;
        i = fac1;  fac1 = fac2;  fac2 = i;
        i = df1;   df1 = df2;    df2 = i;
      }
      fhe = imfac[fac2].firsthalfedge;
      if ( imfhei[fhe] == oth )        { le0 = imfhei[fhe+1];  le1 = imfhei[fhe+2]; }
      else if ( imfhei[fhe+1] == oth ) { le0 = imfhei[fhe+2];  le2 = imfhei[fhe];   }
      else                             { le0 = imfhei[fhe];    le1 = imfhei[fhe+1]; }
                /* setup new halfedge numbers */
      k = imfac[fac2].firsthalfedge;
      for ( i = 0; i < 3; i++ ) {
        j = imfhei[k+i];
        newhenum[j] = -1;
        if ( imhe[j].v0 != v0 && imhe[j].v0 != v1 )
          v2 = imhe[j].v0;
      }
      newhenum[nche] = newhenum[oth] = -1;
      for ( i = j = 0;i < inhe; i++ )
        if ( newhenum[i] >= 0 )
          newhenum[i] = j++;
                /* setup output facets */
      for ( i = j = k = 0; i < infac; i++ )
        if ( i == fac1 ) {
          omfac[j].degree = df = imfac[i].degree-1;
          omfac[j].firsthalfedge = k;
          fhe = imfac[i].firsthalfedge;
          for ( l = m = 0; l <= df; l++ )
            if ( newhenum[imfhei[fhe+l]] >= 0 ) {
              omfhei[k+m] = newhenum[imfhei[fhe+l]];
              m ++;
            }
          k += df;
          j ++;
        }
        else if ( i != fac2 ) {
          omfac[j].degree = df = imfac[i].degree;
          omfac[j].firsthalfedge = k;
          fhe = imfac[i].firsthalfedge;
          for ( l = 0; l < df; l++ )
            omfhei[k+l] = newhenum[imfhei[fhe+l]];
          k += df;
          j ++;
        }
      *onfac = infac-1;
                /* setup the output vertices */
      for ( i = j = k = 0; i < inv; i++ ) {
        dv = imv[i].degree;
        fhe = imv[i].firsthalfedge;
        if ( i == v0 ) {
          dv1 = imv[v1].degree;
          fhe1 = imv[v1].firsthalfedge;
          omv[j].degree = df = dv+dv1-3;
          if ( df > 127 )
            goto failure;
          omv[j].firsthalfedge = k;
          if ( imhe[imvhei[fhe1+dv1-1]].otherhalf == -1 ) {
                 /* v1 at the boundary */
            for ( l = m = 0; l < dv1; l++ )
              if ( newhenum[imvhei[fhe1+l]] >= 0 ) {
                omvhei[k+m] = newhenum[imvhei[fhe1+l]];
                m ++;
              }
              else if ( imvhei[fhe1+l] == nche || imvhei[fhe1+l] == oth ) {
                for ( n = 0; n < dv; n++ )
                  if ( imvhei[fhe+n] == oth || imvhei[fhe+n] == nche )
                    break;
                if ( n >= dv )
                  goto failure;
                for ( o = 0, n = (n+1) % dv;  o < dv-1;  o++, n = (n+1) % dv, m++ )
                  omvhei[k+m] = newhenum[imvhei[fhe+n]];
              }
          }
          else {
                 /* v1 inner */
            for ( l = m = 0; l < dv; l++ )
              if ( newhenum[imvhei[fhe+l]] >= 0 ) {
                omvhei[k+m] = newhenum[imvhei[fhe+l]];
                m ++;
              }
              else if ( imvhei[fhe+l] == nche || imvhei[fhe+l] == oth ) {
                for ( n = 0; n < dv1; n++ )
                  if ( imvhei[fhe1+n] == oth || imvhei[fhe1+n] == nche )
                    break;
                if ( n > dv1 )
                  goto failure;
                for ( o = 0, n = (n+1) % dv1;  o < dv1-1;  o++, n = (n+1) % dv1, m++ )
                  omvhei[k+m] = newhenum[imvhei[fhe1+n]];
              }
          }
          k += omv[j].degree;
          j ++;
        }
        else if ( i == v2 ) {
          omv[j].degree = dv-1;
          omv[j].firsthalfedge = k;
          for ( l = m = 0; l < dv; l++ )
            if ( newhenum[imvhei[fhe+l]] >= 0 ) {
              omvhei[k+m] = newhenum[imvhei[fhe+l]];
              m ++;
            }
          k += dv-1;
          j ++;
        }
        else if ( i != v1 ) {
          omv[j].degree = dv;
          omv[j].firsthalfedge = k;
          for ( l = 0; l < dv; l++ )
            omvhei[k+l] = newhenum[imvhei[fhe+l]];
          k += dv;
          j ++;
        }
      }
      *onv = inv-1;
                /* setup the output halfedges */
      for ( i = 0; i < inhe; i++ ) {
        j = newhenum[i];
        if ( j >= 0 ) {
          if ( imhe[i].v0 == v1 ) omhe[j].v0 = vn;
          else omhe[j].v0 = imhe[i].v0 > v1 ? imhe[i].v0-1 : imhe[i].v0;
          if ( imhe[i].v1 == v1 ) omhe[j].v1 = vn;
          else omhe[j].v1 = imhe[i].v1 > v1 ? imhe[i].v1-1 : imhe[i].v1;
          omhe[j].facetnum = imhe[i].facetnum > fac2 ?
                             imhe[i].facetnum-1 : imhe[i].facetnum;
          omhe[j].otherhalf = imhe[i].otherhalf >= 0 ?
                              newhenum[imhe[i].otherhalf] : -1;
        }
      }
      le0 = imhe[le0].otherhalf;  if ( le0 >= 0 ) le0 = newhenum[le0];
      le1 = imhe[le1].otherhalf;  if ( le1 >= 0 ) le1 = newhenum[le1];
      if ( le0 >= 0 ) omhe[le0].otherhalf = le1;
      if ( le1 >= 0 ) omhe[le1].otherhalf = le0;

      *onhe = inhe-4;
    }
  }

  pkv_SetScratchMemTop ( sp );
  return vn;

failure:
  *onv = *onhe = *onfac = -1;
  pkv_SetScratchMemTop ( sp );
  return -1;
} /*bsm_ContractEdged*/

