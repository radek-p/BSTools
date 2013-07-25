
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2005, 2008                            */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

boolean g1h_ComputeBasisd ( GHoleDomaind *domain )
{
  void  *sp;

  sp = pkv_GetScratchMemTop ();

  if ( !gh_FindDomSurrndBezPatchesd ( domain ) )
    goto failure;
  if ( !_gh_FindDomSurrndBFuncPatchesd ( domain ) )
    goto failure;
  if ( !FindAuxDPatchesd ( domain ) )
    goto failure;
  if ( !FindJFunctionsd ( domain ) )
    goto failure;
  if ( !FindDiPatchesd ( domain ) )
    goto failure;
  if ( !FindBasisFunctionsd ( domain ) )
    goto failure;

  pkv_SetScratchMemTop ( sp );
  domain->basisG1 = true;
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  domain->basisG1 = false;
  return false;
} /*g1h_ComputeBasisd*/

