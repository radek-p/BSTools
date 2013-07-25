
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2005, 2012                            */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <memory.h>
#include <pthread.h>

#define CONST_

#include "pkvaria.h"
#include "pkvthreads.h"
#include "pkgeom.h"
#include "multibs.h"
#include "raybez.h"

pthread_mutex_t raybez_mutex;
boolean         raybez_use_mutex = false;

boolean raybez_InitMutex ( void )
{
  int rc;

  rc = pthread_mutex_init ( &raybez_mutex, NULL );
  return raybez_use_mutex = (rc == 0);
} /*raybez_InitMutex*/

void raybez_DestroyMutex ( void )
{
  pthread_mutex_destroy ( &raybez_mutex );
  raybez_use_mutex = false;
} /*raybez_DestroyMutex*/

