
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2012, 2013                            */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

/* Header file for the libraybez library of C procedures -               */
/* ray tracing for Bezier patches and curves (rendered as tubes)         */

#ifndef RAYBEZ_H
#define RAYBEZ_H

#ifndef _PTHREAD_H
#include <pthread.h>
#endif

#ifndef PKGEOM_H
#include "pkgeom.h"
#endif

#ifndef MULTIBS_H
#include "multibs.h"
#endif

#include "raybezf.h"
#include "raybezd.h"

extern pthread_mutex_t raybez_mutex;
extern boolean         raybez_use_mutex;

#ifdef __cplusplus   
extern "C" {
#endif

boolean raybez_InitMutex ( void );
void raybez_DestroyMutex ( void );

#ifdef __cplusplus
}
#endif

#endif

