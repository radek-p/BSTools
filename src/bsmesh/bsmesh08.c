
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
boolean _bsm_RotateHalfedgei ( int degree, int *hei, int hn )
{
  void *sp;
  int  *whei, i;

  sp = pkv_GetScratchMemTop ();
  whei = pkv_GetScratchMemi ( degree );
  if ( !whei )
    goto failure;
  for ( i = 0; i < degree && hei[i] != hn; i++ )
    ;
  if ( i >= degree )
    goto failure;
  if ( i == degree-1 )
    goto success;
  memcpy ( whei, hei, degree*sizeof(int) );
  memcpy ( hei, &whei[i+1], (degree-i-1)*sizeof(int) );
  memcpy ( &hei[degree-i-1], whei, (i+1)*sizeof(int) );
success:
  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return true;
} /*_bsm_RotateHalfedgei*/

