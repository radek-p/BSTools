
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

/* The procedure bsf_ReadIdentifiers reads the entire .bs data file, scanning */
/* it to find all identifiers (ident n). It is supposed to be used by         */
/* applications before appending data to the .bs file in order to make all    */
/* identifiers unique through the .bs file. It is legal not to assign         */
/* identifiers to the objects (call the writing procedures with the parameter */
/* ident equal to -1) and ignore these identifiers. The identifiers are       */
/* intended to carry the information about dependencies among the objects.    */

/* the return value is true if the specified file is empty or it has been     */
/* successfully read. For each identifier the application procedure is called */
/* and it should make a list of the identifiers. Then, generating the         */
/* identifiers of objects to append the application ought to avoid the        */
/* numbers from the list.                                                     */
boolean bsf_ReadIdentifiers ( const char *filename, void *usrdata,  
                              boolean (*readident)( void *usrdata,
                                                    int objtype, int ident ) )
{
  void *sp;
  int  objtype;

  sp = pkv_GetScratchMemTop ();
  if ( !bsf_OpenInputFile ( filename ) )
    return false;
  objtype = -1;
  for (;;) {
    switch ( bsf_nextsymbol ) {
  case BSF_SYMB_BCURVE:
  case BSF_SYMB_BPATCH:
  case BSF_SYMB_BSCURVE:
  case BSF_SYMB_BSPATCH:
  case BSF_SYMB_BSMESH:
  case BSF_SYMB_BSHOLE:
      objtype = bsf_nextsymbol;
      break;
  case BSF_SYMB_IDENT:
      bsf_GetNextSymbol ();
      if ( bsf_nextsymbol != BSF_SYMB_INTEGER )
        goto failure;
      readident ( usrdata, objtype, bsf_nextint );
      break;
  case BSF_SYMB_ERROR:
      goto failure;
  case BSF_SYMB_EOF:
      goto finish;
  default:
      break;
    }
    bsf_GetNextSymbol ();
  }
finish:
  bsf_CloseInputFile ();
  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  bsf_CloseInputFile ();
  pkv_SetScratchMemTop ( sp );
  return false;
} /*bsf_ReadIdentifiers*/

