
/* //////////////////////////////////////////////////// */
/* This file is a part of the BSTools procedure package */
/* written by Przemyslaw Kiciak                         */
/* //////////////////////////////////////////////////// */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <malloc.h>
#include <string.h>

#include <X11/Xlib.h>
#include <X11/Xutil.h>

#include "pkvaria.h"
#include "pknum.h"
#include "pkgeom.h"
#include "camerad.h"
#include "multibs.h"
#include "eg1holed.h"
#include "eg2holed.h"
#include "xgedit.h"

#include "render.h"
#include "edpolicz.h"
#include "edwidgets.h"
#include "splhole.h"
#include "bsfile.h"


boolean FilenameCorrect ( const char *fn )
{
  return (boolean)(fn[0] != 0);
} /*FilenameCorrect*/

void BSHoleReader ( void *userData, const char *name, int holek,
                    const double *kn, const point2d *domaincp,
                    const point4d *holecp,
                    int spdimen, boolean rational )
{
  bsf_UserReaders *readers;

  if ( spdimen == 3 && !rational && holek >= 3 && holek <= GH_MAX_K ) {
    hole_k = holek;
    nctrlp = 12*hole_k+1;
    pkv_Selectd ( nctrlp, 3, 4, 3, holecp, hole_cp );
    memcpy ( domain_cp, domaincp, nctrlp*sizeof(point2d) );
    memset ( mkdcp, 0, nctrlp );
    memcpy ( knots, kn, 11*holek*sizeof(double) );
    readers = (bsf_UserReaders*)userData;    
    readers->done = true;
  }
} /*BSHoleReader*/

void OpenFile ( const char *fn )
{
  bsf_UserReaders readers;

  if ( RenderingIsOn )
    BreakRendering ( false );
  bsf_ClearReaders ( &readers );
  bsf_BSH4ReadFuncd ( &readers, BSHoleReader );
  readers.userData = &readers;
  readers.done = false;
  bsf_ReadBSFiled ( filename, &readers );
  if ( readers.done ) {
/*    memset ( mkhcp, 0, nctrlp );*/
    FindBoundingBox ( &swind.DefBBox );
    FindDomainBoundingBox ( &domwind.DefBBox );
    xge_3DwindSetupParProj ( &swind, &swind.DefBBox );
    swind.PerspBBox = swind.RefBBox;
    xge_3DwindSetupPerspProj ( &swind, false );
    ConfigureConstraintWidgets ( true );
    ResizeWinStatus ( win1 );
    view_surf_1 = view_surf_2 = swind_picture = false;
    options1.constr_matrix_valid = options2.constr_matrix_valid = false;
    InvalFinalSurfaces ();
    ProjectSurfaceNet ();
    ProjectDomainNet ();
    if ( !UpdateDomains () )
      goto failure;
    InitConstraintFrame ( 1 );
    InitConstraintFrame ( 2 );
    xge_RedrawAll ();
  }
  else {
failure:
    InitGHObject ( 5 );
    xge_DisplayErrorMessage ( ErrMsgFileReading, -1 );
  }
} /*OpenFile*/

boolean WriteHoleAttrib ( void *userData )
{
  void         *sp;
  int          ncp, i;
  unsigned int *mk;

  sp = pkv_GetScratchMemTop ();
  ncp = 12*hole_k+1;
  mk = (unsigned int*)pkv_GetScratchMemi ( ncp );
  if ( mk ) {
    for ( i = 0; i < ncp; i++ )
      mk[i] = mkhcp[i];
    bsf_WritePointsMK ( ncp, mk );
  }
  pkv_SetScratchMemTop ( sp );
  return true;
} /*WriteHoleAttrib*/

void SaveFile ( const char *fn )
{
  int lfn, lex;

  if ( !FilenameCorrect ( fn ) ) {
    xge_DisplayErrorMessage ( ErrMsgBadFilename, -1 );
    return;
  }
  lfn = strlen ( filename );
  lex = strlen ( file_ext );
  if ( strcmp ( file_ext, &filename[lfn-lex]) )
    strcpy ( &filename[lfn], file_ext );
  if ( bsf_OpenOutputFile ( filename, false ) ) {
    bsf_WriteBSplineHoled ( 3, 3, false, hole_k, knots, domain_cp, &hole_cp[0].x,
                            NULL, WriteHoleAttrib, NULL );
    bsf_CloseOutputFile ();
    return;
  }
  xge_DisplayErrorMessage ( ErrMsgFileWriting, -1 );
} /*SaveFile*/

