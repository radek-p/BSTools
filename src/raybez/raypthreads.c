
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2015                                  */
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

#include "raybezprivate.h"

typedef struct {
    pthread_mutex_t mutex;
    pthread_cond_t  cond;
    void            *vertex;
    int             count;
    unsigned short  n_waiting_threads;
    unsigned short  errcode;
  } vertq;

static unsigned short  raybez_npthreads = 0;
static pthread_mutex_t raybez_mutex;
static vertq           *raybez_vert = NULL;
static unsigned short  *raybez_pt = NULL;
static unsigned short  raybez_cdiv = 0;

void raybez_DisablePThreads ( void )
{
  int i;

  if ( raybez_vert ) {
    for ( i = 0; i < raybez_npthreads; i++ ) {
      if ( !(raybez_vert[i].errcode & 0x01) )
        pthread_mutex_destroy ( &raybez_vert[i].mutex );
      if ( !(raybez_vert[i].errcode & 0x02 ) )
        pthread_cond_destroy ( &raybez_vert[i].cond );
    }
    PKV_FREE ( raybez_vert );
    pthread_mutex_destroy ( &raybez_mutex );
    raybez_npthreads = raybez_cdiv = 0;
    if ( raybez_pt ) PKV_FREE ( raybez_pt );
  }
} /*raybez_DisablePThreads*/

boolean raybez_EnablePThreads ( void )
{
  int i;

  raybez_DisablePThreads ();
  PKV_MALLOC ( raybez_vert, max_threads*sizeof(vertq) );
  PKV_MALLOC ( raybez_pt, max_threads*sizeof(unsigned short) );
  if ( raybez_vert && raybez_pt ) {
    memset ( raybez_vert, 0, max_threads*sizeof(vertq) );
    raybez_npthreads = max_threads;
    if ( pthread_mutex_init ( &raybez_mutex, NULL ) )
      goto failure;
    for ( i = 0; i < max_threads; i++ ) {
      raybez_vert[i].errcode = 0;
      if ( pthread_mutex_init ( &raybez_vert[i].mutex, NULL ) )
        raybez_vert[i].errcode = 1;
      if ( pthread_cond_init ( &raybez_vert[i].cond, NULL ) )
        raybez_vert[i].errcode += 2;
      if ( raybez_vert[i].errcode ) {
        raybez_npthreads = i;
        raybez_DisablePThreads ();
        return false;
      }
    }
    raybez_cdiv = 0;
    return true;
  }

failure:
  if ( raybez_pt )   PKV_FREE ( raybez_pt );
  if ( raybez_vert ) PKV_FREE ( raybez_vert );
  return false;
} /*raybez_EnablePThreads*/

void raybez_DivideTreeVertex ( void *tree, void *vertex, unsigned char *tag,
                               divide_vertex_proc DivideVertex )
{
  unsigned short my_id, its_id, vp;

  if ( raybez_npthreads == 0 ) {
    DivideVertex ( tree, vertex );
    *tag = 2;
  }
  else if ( *tag < 2 ) {
    my_id = pkv_PThreadMyPos ();
/*printf ( "%d (a)\n", my_id );*/
    pthread_mutex_lock ( &raybez_mutex );
    switch ( *tag ) {
  case 0:    /* the vertex is an intact leaf - divide it */
      *tag = 1;
      raybez_vert[my_id].count ++;
      raybez_vert[my_id].vertex = vertex;
      raybez_vert[my_id].n_waiting_threads = 0;
      raybez_pt[raybez_cdiv++] = my_id;
/*printf ( "%d (b)\n", my_id );*/
      pthread_mutex_unlock ( &raybez_mutex );

      DivideVertex ( tree, vertex );

/*printf ( "%d (c)\n", my_id );*/
      pthread_mutex_lock ( &raybez_mutex );
      *tag = 2;
      raybez_vert[my_id].vertex = NULL;
      raybez_cdiv --;
      for ( vp = 0; vp < raybez_cdiv; vp ++ )
        if ( raybez_pt[vp] == my_id ) {
          raybez_pt[vp] = raybez_pt[raybez_cdiv];
          break;
        }
      switch ( raybez_vert[my_id].n_waiting_threads ) {
    case 0:   /* none to wake up */
/*printf ( "%d (d)\n", my_id );*/
        pthread_mutex_unlock ( &raybez_mutex );
        break;
    case 1:   /* one to wake up */
/*printf ( "%d (e)\n", my_id );*/
        pthread_mutex_lock ( &raybez_vert[my_id].mutex );
        pthread_mutex_unlock ( &raybez_mutex );
/*printf ( "%d (f)\n", my_id );*/
        pthread_cond_signal ( &raybez_vert[my_id].cond );
/*printf ( "%d (g)\n", my_id );*/
        pthread_mutex_unlock ( &raybez_vert[my_id].mutex );
        break;
    default:  /* more than 1 to wake up */
/*printf ( "%d (h)\n", my_id );*/
        pthread_mutex_lock ( &raybez_vert[my_id].mutex );
        pthread_mutex_unlock ( &raybez_mutex );
/*printf ( "%d (i)\n", my_id );*/
        pthread_cond_broadcast ( &raybez_vert[my_id].cond );
/*printf ( "%d (j)\n", my_id );*/
        pthread_mutex_unlock ( &raybez_vert[my_id].mutex );
        break;
      }
/*printf ( "%d (k)\n", my_id );*/
      raybez_vert[my_id].n_waiting_threads = 0;
      break;

  case 1:    /* some thread is now dividing the leaf - wait for the result */
               /* linear search - perhaps to be changed later */
/*printf ( "%d (m)\n", my_id );*/
      for ( vp = 0; vp < raybez_cdiv; vp++ ) {
        its_id = raybez_pt[vp];
        if ( vertex == raybez_vert[its_id].vertex )
          goto hang;
      }
      printf ( "raybez_pthread error.\n" );
      exit ( 1 );
hang:
/*printf ( "%d (n)\n", my_id );*/
      pthread_mutex_lock ( &raybez_vert[its_id].mutex );
      raybez_vert[its_id].n_waiting_threads ++;
/*printf ( "%d (o)\n", my_id );*/
      pthread_mutex_unlock ( &raybez_mutex );
/*printf ( "%d (p)\n", my_id );*/
      pthread_cond_wait ( &raybez_vert[its_id].cond, &raybez_vert[its_id].mutex );
/*printf ( "%d (q)\n", my_id );*/
      pthread_mutex_unlock ( &raybez_vert[its_id].mutex );
/*printf ( "%d (r)\n", my_id );*/
      break;

  case 2:    /* meanwhile the leaf has been divided - enjoy */
/*printf ( "%d (s)\n", my_id );*/
      pthread_mutex_unlock ( &raybez_mutex );
      break;
    }
  }
/*printf ( "%d (t)\n", my_id );*/
} /*raybez_DivideTreeVertex*/

void raybez_GetPThreadCounts ( int *cnt )
{
  int i;

  if ( raybez_vert )
    for ( i = 0; i < raybez_npthreads; i++ )
      cnt[i] = raybez_vert[i].count;
} /*raybez_GetPThreadCounts*/

