
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

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <malloc.h>
#include <string.h>

#include <X11/Xlib.h>
#include <X11/Xutil.h>

#include "xgedit.h"


void (*xge_childcallback) ( int msg, int size );

void xge_CallTheParent ( int cmd, int size )
{
  XEvent event;

  event.type = ClientMessage;
  event.xclient.window = xgeparentwindow;
  event.xclient.message_type = xgemsg_CHILD_MESSAGE;
  event.xclient.format = 32;
  event.xclient.data.l[0] = cmd;
  event.xclient.data.l[1] = size;
  event.xclient.data.l[2] = 0;
  XSendEvent ( xgedisplay, xgeparentwindow, true, 0, &event );
} /*xge_CallTheParent*/

void xge_ChildCallYourself ( int cmd )
{
  XEvent event;

  event.type = ClientMessage;
  event.xclient.window = xgechildwindow;
  event.xclient.message_type = 1;
  event.xclient.format = 32;
  event.xclient.data.l[0] = cmd;
  event.xclient.data.l[1] = 0;
  event.xclient.data.l[2] = 0;
  XSendEvent ( xgedisplay, xgechildwindow, true, 0, &event );
} /*xge_ChildCallYourself*/

void xge_ChildMessageLoop ( void )
{
  for (;;) {
    XNextEvent ( xgedisplay, &xgeevent );
    xge_childcallback ( xgeevent.xclient.data.l[0],
                        xgeevent.xclient.data.l[1] );
  }
} /*xge_ChildMessageLoop*/

void xge_ChildFlushPipe ( void )
{
  int flag, buf;

        /* set the nonblocking mode for the input pipe */
  flag = fcntl ( xge_pipe_in[0], F_GETFL, 0 );
  flag |= O_NONBLOCK;
  fcntl ( xge_pipe_in[0], F_SETFL, flag );
        /* empty the pipe */
  while ( read ( xge_pipe_in[0], &buf, sizeof(buf) ) > 0 )
    ;
        /* restore the blocking mode */
  flag &= ~O_NONBLOCK;
  fcntl ( xge_pipe_in[0], F_SETFL, flag );
} /*xge_ChildFlushPipe*/

boolean xge_ChildInit ( int argc, char **argv, int magic,
                        void (*callback) ( int msg, int size ) )
{
  int ppid, cpid, mag;
  int flag;

        /* the child process is executed with 6 parameters: */
        /* the parent and child pid, numbers of communication pipes */
        /* parent window identifier and the magic number, which */
        /* ought to be unique for each application */
  if ( argc != 6 )
    goto failure1;
        /* this program is supposed to work only if the */
        /* parent process passes the matching pids - */
        /* its own and the child's, and the magic number also matches */
  xge_parent_pid = getppid ();
  xge_child_pid = getpid ();
  ppid = atoi ( argv[0] );
  cpid = atoi ( argv[1] );
  mag = atoi ( argv[5] );
  if ( ppid != xge_parent_pid || cpid != xge_child_pid || mag != magic )
    goto failure1;
        /* prepare the communication pipes */
  xge_pipe_in[0] = atoi ( argv[2] );
  xge_pipe_out[1] = atoi ( argv[3] );
        /* the input pipe is read in the blocking mode */
  flag = fcntl ( xge_pipe_in[0], F_GETFL, 0 );
  flag &= ~O_NONBLOCK;
  fcntl ( xge_pipe_in[0], F_SETFL, flag );
        /* enable communication via XWindow */
  xgedisplay = XOpenDisplay ( "" );
  xgescreen = DefaultScreen ( xgedisplay );
  if ( !xgedisplay )
    goto failure2;
  xgeparentwindow = (Window)atol ( argv[4] );
/*
printf ( "2: parent window = %d\n", (int)xgeparentwindow );
*/
  xgechildwindow = XCreateWindow ( xgedisplay,
                      XRootWindow ( xgedisplay, xgescreen ), 0, 0, 1, 1, 0,
                      CopyFromParent, InputOnly, CopyFromParent, 0, NULL );
/*
printf ( "2: child window = %d\n", (int)xgechildwindow );
*/
        /* send the child window id as the very first thing */
  write ( xge_pipe_out[1], &xgechildwindow, sizeof(xgechildwindow) );

  xge_childcallback = callback;
  return true;

failure1:
  printf ( "Do not execute this program from a shell.\n" );
  return false;

failure2:
  printf ( "Error: Cannot open communication with XWindow\n" );
  return false;
} /*xge_ChildInit*/

