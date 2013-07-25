
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2005, 2008                            */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

boolean g1h_ComputeBasisf ( GHoleDomainf *domain )
{
  void  *sp;

  sp = pkv_GetScratchMemTop ();

  if ( !gh_FindDomSurrndBezPatchesf ( domain ) )
    goto failure;
  if ( !_gh_FindDomSurrndBFuncPatchesf ( domain ) )
    goto failure;
  if ( !FindAuxDPatchesf ( domain ) )
    goto failure;
  if ( !FindJFunctionsf ( domain ) )
    goto failure;
  if ( !FindDiPatchesf ( domain ) )
    goto failure;
  if ( !FindBasisFunctionsf ( domain ) )
    goto failure;

  pkv_SetScratchMemTop ( sp );
  domain->basisG1 = true;
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  domain->basisG1 = false;
  return false;
} /*g1h_ComputeBasisf*/

