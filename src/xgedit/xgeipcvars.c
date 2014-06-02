
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

#include <X11/Xlib.h>
#include <X11/Xutil.h>

#include "xgedit.h"


pid_t xge_parent_pid, xge_child_pid;
int xge_pipe_in[2], xge_pipe_out[2];
Window xgeparentwindow, xgechildwindow;

