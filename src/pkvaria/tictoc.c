
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2012                                  */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

#include <math.h>
#include <malloc.h>
#include <string.h>
#include <sys/times.h>
#include <unistd.h>

#include "pkvaria.h"

static clock_t tic;

void pkv_Tic ( clock_t *_tic )
{
  struct tms tt;

  if ( _tic )
    *_tic = times ( &tt );
  else
    tic = times ( &tt );
} /*pkv_Tic*/

int pkv_Toc ( clock_t *_tic )
{
  struct tms tt;
  clock_t    toc;

  toc = times ( &tt );
  if ( _tic )
    return toc - *_tic;
  else
    return toc-tic;
} /*pkv_Toc*/

float pkv_Seconds ( clock_t ticks )
{
  return (float)ticks/(float)sysconf( _SC_CLK_TCK );
} /*pkv_Seconds*/

