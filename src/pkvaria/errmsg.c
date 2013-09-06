
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2005, 2013                            */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#include "pkvaria.h"


static void (*ehnd) ( int module, const char *file, int line,
                      int errcode, const char *errstr ) = NULL;


void pkv_SignalError ( int module, const char *file, int line,
                       int errcode, const char *errstr )
{
  if ( ehnd )
    ehnd ( module, file, line, errcode, errstr );
  else {
    fprintf ( stderr, "Error in module %d, file %s, line %d: %s\n",
                      module, file, line, errstr );
    exit ( 1 );
  }
} /*pkv_SignalError*/

void pkv_SetErrorHandler (
        void (*ehandler)( int module, const char *file, int line,
                          int errcode, const char *errstr ) )
{
  ehnd = ehandler;
} /*pkv_SetErrorHandler*/

