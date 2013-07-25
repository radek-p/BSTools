
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2010, 2011                            */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

#include <limits.h>
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


void (*g2mbl_outputnzdistr)( int nbl, int blnum, boolean final,
                             int nvcp, int n, byte *nzdistr ) = NULL;

/* ////////////////////////////////////////////////////////////////////////// */
boolean _g2mbl_OrderCPoints ( int nv, int nvcp, int fcp,
                              int nzcdsize, byte *nzcdistr,
                              int *nncpi, int *vncpi, vertex_desc *nvcpi,
                              int ndomelems, meshdom_elem *domelem, int *domelcpind,
                              int *domelind,
                              int blsize, int *hsize, int *hprof, int *vpermut )
{
  void *sp;
  int  nvneigh;
  int  i, j, k, l, vi, vb, ei, eb;
  int  *vneigh, *permut;

  sp = pkv_GetScratchMemTop ();
  nvcpi[0].fneigh = nvcpi[0].nncn = 0;
  nvneigh = nvcpi[0].nneigh;
  for ( i = 1; i < nvcp; i++ ) {
    nvcpi[i].fneigh = nvneigh;
    nvcpi[i].nncn = 0;
    nvneigh += nvcpi[i].nneigh;
  }
  vneigh = pkv_GetScratchMemi ( nvneigh );
  if ( !vneigh )
    goto failure;
  memset ( vneigh, 0xFF, nvneigh*sizeof(int) );
  memset ( nzcdistr, 0, nzcdsize );
  for ( i = 0; i < ndomelems; i++ ) {
    if ( domelind )
      j = domelind[i];
    else
      j = i;
    vi = domelem[j].ncp;
    vb = domelem[j].firstcpi;
    for ( j = 0; j < vi; j++ ) {
      ei = nncpi[domelcpind[vb+j]] - fcp;
      if ( ei >= 0 && ei < nvcp )
        for ( k = 0; k < vi; k++ ) {
          eb = nncpi[domelcpind[vb+k]] - fcp;
          if ( eb > ei && eb < nvcp )
            if ( !pkn_TMBTestAndSet ( nzcdistr, eb, ei ) ) {
              vneigh[nvcpi[eb].fneigh+nvcpi[eb].nncn] = ei;
              nvcpi[eb].nncn ++;
              vneigh[nvcpi[ei].fneigh+nvcpi[ei].nncn] = eb;
              nvcpi[ei].nncn ++;
            }
        }
    }
  }
  permut = pkv_GetScratchMemi ( nvcp );
  if ( !permut )
    goto failure;
  for ( i = 0; i < nvcp; i++ )
    permut[i] = i;
  j = 0;
  k = nvcpi[0].nncn;
  for ( i = 1; i < nvcp; i++ )
    if ( nvcpi[i].nncn < k ) { j = i;  k = nvcpi[i].nncn; }
  if ( j )
    { permut[0] = j;  permut[j] = 0; }
  nvcpi[j].firstsel = 0;
  for ( i = 0; i < nvcp-1; i++ ) {
    k = permut[i];
    ei = nvcpi[k].fneigh;
    eb = nvcpi[k].nneigh;
    for ( j = 0; j < eb; j++ ) {
      l = vneigh[ei+j];
      if ( pkn_TMBTestAndClear ( nzcdistr, k, l ) ) {
        nvcpi[l].firstsel = min ( nvcpi[l].firstsel, i );
        nvcpi[l].nncn --;
      }
    }
          /* selection */
    k = i+1;
    vi = permut[k];
    ei = nvcpi[vi].nncn;
    eb = nvcpi[vi].firstsel;
    for ( j = i+2; j < nvcp; j++ ) {
      l = permut[j];
      if ( nvcpi[l].firstsel < eb ||
           (nvcpi[l].firstsel == eb && nvcpi[l].nncn < ei ) ) {
        k = j;
        ei = nvcpi[l].nncn;
        eb = nvcpi[l].firstsel;
      }
    }
    if ( k != i+1 )
      { l = permut[i+1];  permut[i+1] = permut[k];  permut[k] = l; }
    if ( nvcpi[permut[i+1]].firstsel > i+1 )
      nvcpi[permut[i+1]].firstsel = i+1;
  }

        /* setup the Hessian matrix profile */
  if ( !vpermut ) {
    for ( i = k = 0;  i < nvcp;  i++, k += blsize ) {
      j = permut[i];
      vncpi[i] = nvcpi[j].vnum;
      nncpi[nvcpi[j].vnum] = i;
      hprof[k] = blsize*nvcpi[j].firstsel;
      for ( l = 1; l < blsize; l++ )
        hprof[k+l] = hprof[k];
    }
  }
  else {
    for ( i = k = 0;  i < nvcp;  i++, k += blsize ) {
      j = permut[i];
      vpermut[j] = i;
      vncpi[i] = nvcpi[j].vnum;
      hprof[k] = blsize*nvcpi[j].firstsel;
      for ( l = 1; l < blsize; l++ )
        hprof[k+l] = hprof[k];
    }
  }
  *hsize = pkn_NRBArraySize ( blsize*nvcp, hprof );
  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*_g2mbl_OrderCPoints*/

