
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
#include "mengerc.h"
#include "xgedit.h"

#define CHILD_SIDE
#include "pozwalajipc.h"

#include "pozwalaj_proc.h"

/* /////////////////////////////////////////////////////////////////////////// */
void GetData ( int size )
{
  int  cnt;
  int  idesc, isize, i;
  char b;

  cnt = 0;
  while ( cnt < size ) {
    read ( xge_pipe_in[0], &idesc, sizeof(int) );
    read ( xge_pipe_in[0], &isize, sizeof(int) );
/*printf ( "idesc = %d, isize = %d\n", idesc, isize ); */
    switch ( idesc ) {
  case ipcd_BLP_SIZE:
      ReadBLPSize ( isize );
      break;
  case ipcd_BLP_CPOINTS:
      ReadBLPCPoints ( isize );
      break;
  case ipcd_BLP_MKCP:
      ReadBLPMkcp ( isize );
      break;
  case ipcd_BLP_OPTIMIZE:
      ReadBLPOptimizeOptions ( isize );
      break;
  case ipcd_BSM_SIZE:
      ReadBSMSize ( isize );
      break;
  case ipcd_BSM_VERT:
      ReadBSMVert ( isize );
      break;
  case ipcd_BSM_VHE:
      ReadBSMVertHE ( isize );
      break;
  case ipcd_BSM_VERTC:
      ReadBSMVertPC ( isize );
      break;
  case ipcd_BSM_VERTMK:
      ReadBSMVertMK ( isize );
      break;
  case ipcd_BSM_HALFE:
      ReadBSMHalfedges ( isize );
      break;
  case ipcd_BSM_FAC:
      ReadBSMFacets ( isize );
      break;
  case ipcd_BSM_FHE:
      ReadBSMFacetHE ( isize );
      break;
  case ipcd_BSM_COARSE_SIZE:
      ReadBSMCSize ( isize );
      break;
  case ipcd_BSM_COARSE_VERT:
      ReadBSMCVert ( isize );
      break;
  case ipcd_BSM_COARSE_VHE:
      ReadBSMCVertHE ( isize );
      break;
  case ipcd_BSM_COARSE_HALFE:
      ReadBSMCHalfedges ( isize );
      break;
  case ipcd_BSM_COARSE_FAC:
      ReadBSMCFacets ( isize );
      break;
  case ipcd_BSM_COARSE_FHE:
      ReadBSMCFacetHE ( isize );
      break;
  case ipcd_BSM_OPTIMIZE:
      ReadBSMOptimizeOptions ( isize );
      break;
  case ipcd_BSC_SIZE:
      ReadBSCSize ( isize );
      break;
  case ipcd_BSC_KNOTS:
      ReadBSCKnots ( isize );
      break;
  case ipcd_BSC_CPOINTS:
      ReadBSCCPoints ( isize );
      break;
  case ipcd_BSC_MKCP:
      ReadBSCMkcp ( isize );
      break;
  case ipcd_BSC_OPTIMIZE:
      ReadBSCMCOptimizeOptions ( isize );
      break;
  default:  /* error - reject all pending data and signal error */
      for ( i = cnt; i < size; i++ )
        read ( xge_pipe_in[0], &b, 1 );
      cnt = size;
      xge_CallTheParent ( ipccmd_ERROR, 0 );
      return;
    }
    cnt += 2*sizeof(int) + isize;
/*printf ( "data item: %d size %d\n", idesc, isize ); */
  }
} /*GetData*/

boolean InvertPretransformation ( trans3d *pretrans )
{
  pretrans_inv = *pretrans;
  return InvertTrans3d ( &pretrans_inv );
} /*InvertPretransformation*/

void TransformCPoints ( trans3d *tr, int n, int m, int pitch, point3d *cp )
{
  int     i, j, k;
  point3d q;

  pitch /= 3;
  for ( i = k = 0;  i < n;  i++, k += pitch )
    for ( j = 0; j < m; j++ ) {
      TransPoint3d ( tr, &cp[k+j], &q );
      cp[k+j] = q;
    }
} /*TransformCPoints*/

