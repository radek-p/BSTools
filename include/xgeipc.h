
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2009                                  */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

#ifndef XGEIPC_H
#define XGEIPC_H

#ifndef _LIBC_LIMITS_H_
#include <limits.h>
#endif

#ifndef XGEDIT_H
#include "xgedit.h"
#endif

#ifdef __cplusplus
extern "C" {
#endif

extern pid_t xge_parent_pid, xge_child_pid;
extern int xge_pipe_in[2], xge_pipe_out[2];
extern Window xgeparentwindow, xgechildwindow;

extern void (*xge_childcallback) ( int msg, int size );

/* the procedures of the parent process */
boolean xge_MakeTheChild ( const char *name, const char *suffix, int magic );
boolean xge_ChildIsActive ( void );
void xge_CallTheChild ( int cmd, int size );
void xge_SignalTheChild ( void );
void xge_ParentFlushPipe ( void );

/* the procedures of the child process */
void xge_CallTheParent ( int cmd, int size );
void xge_ChildCallYourself ( int cmd );
void xge_ChildMessageLoop ( void );
void xge_ChildFlushPipe ( void );
boolean xge_ChildInit ( int argc, char **argv, int magic,
                        void (*callback) ( int msg, int size ) );

#ifdef __cplusplus
}
#endif

#endif

