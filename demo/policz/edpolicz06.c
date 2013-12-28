
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

void OpenFile ( const char *fn )
{
  void    *sp;
  boolean result, rational;
  point4d *hcp;
  int     spdimen;

  sp = pkv_GetScratchMemTop ();
  hcp = pkv_GetScratchMem ( (12*GH_MAX_K+1)*sizeof(point4d) );
  if ( !hcp )
    goto failure;
  if ( RenderingIsOn )
    BreakRendering ( false );
  if ( bsf_OpenInputFile ( filename ) ) {
    result = bsf_ReadBSplineHoled ( GH_MAX_K, &hole_k, knots,
                                    domain_cp, hcp, &spdimen, &rational,
                                    mkhcp, NULL );
    bsf_CloseInputFile ();
  }
  else {
    bsf_PrintErrorLocation ();
    result = false;
  }
  if ( result && !rational ) {
    nctrlp = 12*hole_k+1;
    pkv_Selectd ( nctrlp, 3, 4, 3, hcp, hole_cp );
    memset ( mkdcp, 0, nctrlp );
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
  pkv_SetScratchMemTop ( sp );
} /*OpenFile*/

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
    bsf_WriteBSplineHoled ( hole_k, knots, domain_cp, hole_cp, mkhcp, NULL,
                            NULL, NULL );
    bsf_CloseOutputFile ();
    return;
  }
  xge_DisplayErrorMessage ( ErrMsgFileWriting, -1 );
} /*SaveFile*/

