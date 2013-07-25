
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2012                                  */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include "pkvaria.h"
#include "pkvthreads.h"
#include "msgpool.h"

/* this procedure tries to find the number of CPUs in the system, by reading */
/* the file /proc/cpuinfo and counting the number of appearances of the word */
/* "processor" in it; in case of failure (no such file exists or it does not */
/* contain this word) the return value is 1 */
int pkv_FindNCPU ( void )
{
  FILE *f;
  int  ncpu;
  char buf, theword[10] = "processor";
  int  i, j;

  f = fopen ( "/proc/cpuinfo", "r" );
  if ( f ) {
    ncpu = 0;
    for (;;) {
      for ( i = 0; theword[i]; i++ ) {
        j = fscanf ( f, "%c", &buf );
        if ( j <= 0 )
          goto endit;
        /* if a character in the file is ':', skip the rest of the line, */
        /* we are not interested in any appearances of the word */
        /* "processor" after the colon */
        if ( buf == ':' ) {
          do {
            j = fscanf ( f, "%c", &buf );
          } while ( j >= 0 && buf != '\n' );
          if ( j <= 0 )
            goto endit;
        }
        if ( buf != theword[i] )
          break;
      }
      if ( !theword[i] )
        ncpu ++;
    }
endit:
    fclose ( f );
    return max ( ncpu, 1 );
  }
  else
    return 1;
} /*pkv_FindNCPU*/

