
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2005, 2009                            */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include "pkvaria.h"
#include "pknum.h"
#include "pkgeom.h"
#include "multibs.h"

#include "eg1holef.h"
#include "eg1holed.h"
#include "eg2holef.h"
#include "eg2holed.h"

#include "eghprivated.h"

/* ///////////////////////////////////////////////////////////////////////// */
void gh_GetBspInd ( int hole_k, int i, int j, int *ind )
{
  int  l;
  int  ind1[16];

  switch ( j ) {
case 0:
    i = i % hole_k;
    ind[3] = 0;
    for ( j = 1; j <= 3; j++ ) {
      for ( l = 0; l <= 3; l++ )
        ind[(3-j)+4*l] = 12*i+3*l+j;
    }
    i = (i+1) % hole_k;
    for ( j = 1; j <= 3; j++ )
      ind[3+4*j] = 12*i+j;
    break;

case 1:
    gh_GetBspInd ( hole_k, i, 0, ind );
    gh_GetBspInd ( hole_k, i+1, 0, ind1 );
    for ( i = 0; i <= 3; i++ ) {
      memmove ( &ind[4*i], &ind[4*i+1], 3*sizeof(int) );
      ind[4*i+3] = ind1[7-i];
    }
    break;

case 2:
    gh_GetBspInd ( hole_k, i, 0, ind );
    gh_GetBspInd ( hole_k, i+1, 0, ind1 );
    for ( i = 0; i <= 3; i++ ) {
      memmove ( &ind[4*i], &ind[4*i+2], 2*sizeof(int) );
      ind[4*i+2] = ind1[7-i];
      ind[4*i+3] = ind1[11-i];
    }
    break;

default:
    exit ( 1 );
  }
} /*gh_GetBspInd*/

/* ////////////////////////////////////////////////////////////////////////// */
int gh_DrawBFcpn ( int hole_k, unsigned char *bfcpn )
{
  void *sp;
  int  *ind;
  int  i, j;

  sp = pkv_GetScratchMemTop ();
  ind = pkv_GetScratchMem ( 16*sizeof(int) );
  if ( !ind )
    goto failure;

  bfcpn[0] = 0;
  i = 1;
  for ( j = 0; j < hole_k; j++, i += 6 ) {
    gh_GetBspInd ( hole_k, j, 0, ind );
    bfcpn[i]   = (unsigned char)ind[2];
    bfcpn[i+1] = (unsigned char)ind[1];
    bfcpn[i+2] = (unsigned char)ind[6];
    bfcpn[i+3] = (unsigned char)ind[5];
    bfcpn[i+4] = (unsigned char)ind[10];
    bfcpn[i+5] = (unsigned char)ind[9];
  }
  pkv_SetScratchMemTop ( sp );
  return 6*hole_k+1;

failure:
  pkv_SetScratchMemTop ( sp );
  return 0;
} /*gh_DrawBFcpn*/

/* ///////////////////////////////////////////////////////////////////////// */
boolean _gh_FindBasisFuncbSupport ( int hole_k, int nfunc_b,
                                    unsigned short *support_b,
                                    unsigned char *bfcpn )
{
  void           *sp;
  int            i, j, k, l, n;
  unsigned short mask, m;
  int            *ind;

  sp = pkv_GetScratchMemTop ();
  ind = pkv_GetScratchMem ( 16*sizeof(int) );
  if ( !ind )
    goto failure;
  mask = m = 0x0001;
  for ( i = 0; i < hole_k; i++ ) {
    mask = (unsigned short)(mask | m);
    m = (unsigned short)(2*m);
  }
  for ( i = 0; i < nfunc_b; i++ )
    support_b[i] = mask;

  bfcpn[0] = 0;
  i = 1;
  for ( j = 0; j < hole_k; j++, i += 6 ) {
    k = (j + hole_k - 1) % hole_k;
    l = (j + hole_k - 2) % hole_k;
    n = (j + hole_k - 3) % hole_k;
    gh_GetBspInd ( hole_k, j, 0, ind );
    bfcpn[i]       = (unsigned char)ind[2];
    support_b[i]   = (unsigned short)((0x0001 << n) | (0x0001 << l) |
                                      (0x0001 << k) | (0x0001 << j));
    bfcpn[i+1]     = (unsigned char)ind[1];
    support_b[i+1] = (unsigned short)((0x0001 << l) | (0x0001 << k));
    bfcpn[i+2]     = (unsigned char)ind[6];
    support_b[i+2] = (unsigned short)((0x0001 << l) | (0x0001 << k) | (0x0001 << j));
    bfcpn[i+3]     = (unsigned char)ind[5];
    support_b[i+3] = (unsigned short)((0x0001 << l) | (0x0001 << k));
    bfcpn[i+4]     = (unsigned char)ind[10];
    support_b[i+4] = (unsigned short)((0x0001 << k) | (0x0001 << j));
    bfcpn[i+5]     = (unsigned char)ind[9];
    support_b[i+5] = (unsigned short)(0x0001 << k);
  }
  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*_gh_FindBasisFuncbSupport*/

