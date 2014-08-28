
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2014                                  */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <malloc.h>
#include <string.h>
#include <ctype.h>

#include "pkvaria.h"
#include "pknum.h"
#include "pkgeom.h"
#include "camerad.h"
#include "multibs.h"
#include "bsmesh.h"
#include "bsfile.h"

#include "bsfprivate.h"

/* dependencies are specified by a keyword followed by a comma-separated   */
/* list of integer object identifiers in braces. The objects with these    */
/* identifiers ought to be present in the file, but the library procedures */
/* do not check it - it is left to be dealt with by the application.       */

boolean bsf_ReadDependencies ( bsf_UserReaders *readers )
{
  void *sp;
  int  depname;
  int  maxdep, ndep;
  int  *dep;

  sp = pkv_GetScratchMemTop ();
  maxdep = readers->maxdep;
  if ( maxdep > 0 )
    dep = pkv_GetScratchMemi ( maxdep );
  else
    dep = NULL;
  depname = bsf_nextsymbol;  /* it is supposed to be one of keywords */
  bsf_GetNextSymbol ();
  if ( bsf_nextsymbol != BSF_SYMB_LBRACE )
    goto failure;
  bsf_GetNextSymbol ();
                            /* now read the list of identifiers of */
                            /* objects with dependencies */
  ndep = 0;
  while ( bsf_nextsymbol == BSF_SYMB_INTEGER ) {
    if ( ndep >= maxdep )
      goto failure;
    if ( dep )
      dep[ndep] = bsf_nextint;
    ndep ++;
    bsf_GetNextSymbol ();
    switch ( bsf_nextsymbol ) {
  case BSF_SYMB_COMMA:
      bsf_GetNextSymbol ();
      break;
  case BSF_SYMB_RBRACE:
      bsf_GetNextSymbol ();
      goto finish;
  default:
      goto failure;
    }
  }
finish:
  if ( readers->DepReader )
    readers->DepReader ( readers->userData, depname, ndep, dep );

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*bsf_ReadDependencies*/

