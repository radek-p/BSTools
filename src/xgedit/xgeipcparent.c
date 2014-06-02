
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2009, 2014                            */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

#include <sys/types.h>
#include <unistd.h>
#include <fcntl.h>
#include <signal.h>
#include <sys/wait.h>

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <malloc.h>
#include <string.h>

#include <X11/Xlib.h>
#include <X11/Xutil.h>

#include "xgedit.h"

static boolean child_active = false;

/* ///////////////////////////////////////////////////////////////////////// */
static void SigChldHandler ( int sig )
{
  XEvent event;
  /*pid_t  pid;*/
  int    status;

  if ( sig == SIGCHLD ) {
    /*pid =*/ wait ( &status );
    child_active = false;
    close ( xge_pipe_in[1] );
    close ( xge_pipe_out[0] );
    event.type = ClientMessage;
    event.xclient.window = xgewindow;
    event.xclient.message_type = xgemsg_CHILD_FAILURE;
    event.xclient.format = 32;  
    event.xclient.data.l[0] = 0;
    event.xclient.data.l[1] = 0;
    event.xclient.data.l[2] = 0;
    XSendEvent ( xgedisplay, xgewindow, true, 0, &event );
  }
} /*SigChldHandler*/

boolean xge_ChildIsActive ( void )
{
  return child_active;
} /*xge_ChildIsActive*/

boolean xge_MakeTheChild ( const char *name, const char *suffix, int magic )
{
  pid_t pid;
  char procname[256], pin[16], pout[16], ppid[16], cpid[16], mag[16], win[16];
  int  l;
  int  flag;

  signal ( SIGCHLD, SigChldHandler );
  if ( pipe ( xge_pipe_in ) == -1 || pipe ( xge_pipe_out ) == -1 ) {
    printf ( "Error: cannot open pipes.\n" );
    return false;
  }
  pid = fork ();
  child_active = true;
  switch ( pid ) {
case -1: /* failure */
    child_active = false;
    return false;

case 0:  /* child process */
    close ( xge_pipe_in[1] );
    close ( xge_pipe_out[0] );
    xge_parent_pid = getppid ();
    xge_child_pid = getpid ();
    sprintf ( ppid, "%d", xge_parent_pid );
    sprintf ( cpid, "%d", xge_child_pid );
    sprintf ( pin, "%d", xge_pipe_in[0] );
    sprintf ( pout, "%d", xge_pipe_out[1] );
    sprintf ( win, "%d", (int)xgewindow );
    sprintf ( mag, "%d", magic );
/*
printf ( "0: parent window = %u\n", (int)xgewindow );
*/
    l = strlen ( name );
    if ( l < 255 - strlen(suffix) ) {
      strcpy ( procname, name );
      strcpy ( &procname[l], suffix  );
        /* the child process is executed with 6 parameters: */
        /* the parent and child pid, numbers of communication pipes */
        /* parent window identifier and the magic number, which */
        /* ought to be unique for each application */
      execl ( procname, ppid, cpid, pin, pout, win, mag, NULL );
      return false;
    }
    else
      return false;

default:  /* parent process */
    xgeparentwindow = xgewindow;
/*
printf ( "1: parent window = %d\n", (int)xgewindow );
*/
    close ( xge_pipe_in[0] );
    close ( xge_pipe_out[1] );
    xge_child_pid = pid;
    xge_parent_pid = getpid ();
          /* receive the number of the child window */
    flag = fcntl ( xge_pipe_out[0], F_GETFL, 0 );
    flag &= ~O_NONBLOCK;
    fcntl ( xge_pipe_out[0], F_SETFL, flag );
    if ( child_active )
      l = read ( xge_pipe_out[0], &xgechildwindow, sizeof(xgechildwindow) );
/*
printf ( "1: child window = %d, l = %d\n", (int)xgechildwindow, l );
*/
    return child_active;
  }
} /*xge_MakeTheChild*/

void xge_CallTheChild ( int cmd, int size )
{
  XEvent event;

  if ( child_active ) {
    event.type = ClientMessage;
    event.xclient.window = xgechildwindow;
    event.xclient.message_type = 1;
    event.xclient.format = 32;  
    event.xclient.data.l[0] = cmd;
    event.xclient.data.l[1] = size;
    event.xclient.data.l[2] = 0;
    XSendEvent ( xgedisplay, xgechildwindow, true, 0, &event );
  }
} /*xge_CallTheChild*/

int kill(pid_t pid, int sig);

void xge_SignalTheChild ( void )
{
  if ( child_active )
    kill ( xge_child_pid, SIGUSR1 );
} /*xge_SignalTheChild*/

void xge_ParentFlushPipe ( void )
{
  int flag, buf;

  if ( child_active ) {
        /* set the nonblocking mode for the input pipe */
    flag = fcntl ( xge_pipe_out[0], F_GETFL, 0 );
    flag |= O_NONBLOCK;
    fcntl ( xge_pipe_out[0], F_SETFL, flag );
        /* empty the pipe */
    while ( read ( xge_pipe_out[0], &buf, sizeof(buf) ) > 0 )
      ;
        /* restore the blocking mode */
    flag &= ~O_NONBLOCK;
    fcntl ( xge_pipe_out[0], F_SETFL, flag );
  }
} /*xge_ParentFlushPipe*/

