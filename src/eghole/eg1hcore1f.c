
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2005, 2008                            */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

void g1h_SetOptionProcf ( GHoleDomainf *domain,  
    int (*OptionProc)( GHoleDomainf *domain, int query, int qn,
                       int *ndata, int **idata, float **fdata ) )
{
  G1HolePrivateRecf *privateG1;

  privateG1 = domain->privateG1;
  if ( !OptionProc )
    privateG1->GetOption = gh_GetDefaultOptionf;
  else
    privateG1->GetOption = OptionProc;
} /*g1h_SetOptionProcf*/

