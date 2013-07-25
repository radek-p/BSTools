
/* ///////////////////////////////////////////////////  */
/* This file is a part of the BSTools procedure package */
/* written by Przemyslaw Kiciak.                        */
/* ///////////////////////////////////////////////////  */

#include <sys/types.h>
#include <sys/times.h>
#include <signal.h>
#include <unistd.h>
#include <setjmp.h>

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <fpu_control.h>

#include <X11/Xlib.h>
#include <X11/Xutil.h>

#include "pkvaria.h"
#include "pknum.h"
#include "pkgeom.h"
#include "multibs.h"
#include "eg2holed.h"
#include "bsmesh.h"
#include "g1blendingd.h"
#include "g2blendingd.h"
#include "g2mblendingd.h"
#include "xgedit.h"
#include "xgeipc.h"

#define CHILD_SIDE
#include "pozwalajipc.h"
#undef CHILD_SIDE

#include "pozwalaj_proc.h"

/* /////////////////////////////////////////////////////////////////////////// */
void ResetIPCBuffer ( void )
{
  ipc_buf_count = ipc_data_size = 0;
} /*ResetIPCBuffer*/

boolean IPCAppendDataItem ( int desc, int size, void *ptr )
{
  if ( size >= 0 && ipc_buf_count < IPC_BUFFER_LENGTH ) {
    ipc_buffer[ipc_buf_count].desc = desc;
    ipc_buffer[ipc_buf_count].size = size;
    ipc_buffer[ipc_buf_count].ptr = ptr;
    ipc_data_size += size + 2*sizeof(int);
    ipc_buf_count ++;
    return true;
  }
  else
    return false;
} /*IPCAppendDataItem*/

void IPCSendData ( void )
{
  int i;

  for ( i = 0; i < ipc_buf_count; i++ ) {
    write ( xge_pipe_out[1], &ipc_buffer[i].desc, sizeof(int) );
    write ( xge_pipe_out[1], &ipc_buffer[i].size, sizeof(int) );
    if ( ipc_buffer[i].size > 0 )
    write ( xge_pipe_out[1], ipc_buffer[i].ptr, ipc_buffer[i].size );
  }
} /*IPCSendData*/

