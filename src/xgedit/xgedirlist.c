
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2007, 2011                            */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

#include <unistd.h>
#include <dirent.h>
#include <ctype.h>
#include <string.h>
#include <stdio.h>
#include <malloc.h>
#include <sys/stat.h>

#include <X11/Xlib.h>
#include <X11/Xutil.h>

#include "pkvaria.h"
#include "xgedit.h"

#include "xgeprivate.h"

/* ///////////////////////////////////////////////////////////////////////// */
static boolean xge_LessStr ( char *s1, char *s2 )
{
  char c1, c2, C1, C2;

        /* find the first position with different characters */
  while ( *s1 && *s2 && *s1 == *s2 ) { s1 ++;  s2 ++; }
  c1 = *s1;  c2 = *s2;
        /* letters follow everything else */
  if ( isalpha(c1) ) {
    if ( isalpha(c2) ) {
        /* both characters are letters; deal with the upper- and lowercase */
      C1 = (char)toupper(c1);  C2 = (char)toupper(c2);
      if ( C1 < C2 )
        return true;
      else if ( C1 > C2 )
        return false;
      else
        return (boolean)(c1 < c2);
    }
    else return false;
  }
  else {
    if ( isalpha(c2) )
      return true;
    else
      return (boolean)(c1 < c2);
  }
} /*xge_LessStr*/

static boolean xge_LessString ( int i, int j, void *usrptr )
{
  xge_listbox *lbox;

  lbox = (xge_listbox*)usrptr;
  return xge_LessStr ( &lbox->itemstr[lbox->itemind[i]],
                       &lbox->itemstr[lbox->itemind[j]] );
} /*xge_LessString*/

static void xge_SwapString ( int i, int j, void *usrptr )
{
  int         k;
  xge_listbox *lbox;

  lbox = (xge_listbox*)usrptr;
  k = lbox->itemind[i];  lbox->itemind[i] = lbox->itemind[j];
  lbox->itemind[j] = k;
} /*xge_SwapString*/

/* ///////////////////////////////////////////////////////////////////////// */
boolean xge_FilterMatches ( const char *name, const char *filter )
{
  int i, j, k, l;

  if ( !filter )
    return true;

  i = j = 0;
  while ( filter[j] ) {
    switch ( filter[j] ) {
  case '?':
      i ++;
      j ++;
      break;
  case '*':
      do j ++; while ( filter[j] == '*' );
      if ( !filter[j] )
        return true;
      for ( k = j+1; filter[k] && filter[k] != '*'; k++ )
        ;
      do {
        if ( filter[j] != '?' )
          for ( ; name[i] && name[i] != filter[j]; i++ )
            ;
        if ( !name[i] )
          return false;
        for ( l = 0; j+l < k; l++ ) {
          if ( filter[j+l] != '?' && name[i+l] != filter[j+l] ) {
            i ++;
            goto notfound;
          }
        }
        i += l;
        j += l;
        goto found;
notfound:
        if ( !name[i+l-1] )
          return false;
      } while ( 1 );  
found:
      if ( name[i] && !filter[j] )
        return false;
      break;
  default:
      if ( name[i] == filter[j] )
        { i ++;  j++; }
      else
        return false;
      break;
    }
  }
  return true;
} /*xge_FilterMatches*/

/* ///////////////////////////////////////////////////////////////////////// */
boolean xge_SetupFileList ( xge_listbox *lbox, const char *dir,
                            const char *filter )
{
  DIR           *dp;
  struct dirent *entry;
  struct stat   statbuf;
  int           lsum, nit, lgt;

  xge_ClearListBox ( lbox );
  if ( (dp = opendir ( dir ) ) == NULL ) {
    return false;
  }
  chdir ( dir );
  lsum = nit = 0;
  while ( (entry = readdir(dp)) != NULL ) {
    stat ( entry->d_name, &statbuf );
    if ( !(S_ISDIR(statbuf.st_mode)) ) {
      if ( xge_FilterMatches ( entry->d_name, filter ) ) {
        nit ++;
        lsum += strlen ( entry->d_name ) + 1;
      }
    }
  }
  closedir ( dp );
  lbox->itemind = malloc ( (nit+1)*sizeof(int) );
  lbox->itemstr = malloc ( lsum );
  if ( !lbox->itemind || !lbox->itemstr ) {
    xge_ClearListBox ( lbox );
    return false;
  }
  lbox->nitems = (short)nit;
  lsum = nit = 0;
  dp = opendir ( dir );
  while ( (entry = readdir(dp)) != NULL ) {
    stat ( entry->d_name, &statbuf );
    if ( !(S_ISDIR(statbuf.st_mode)) ) {
      if ( xge_FilterMatches ( entry->d_name, filter ) ) {
        lgt = strlen ( entry->d_name );
        memcpy ( &lbox->itemstr[lsum], entry->d_name, lgt+1 );
        lbox->itemind[nit] = lsum;
        lsum += lgt+1;
        nit++;
      }
    }
  }
  closedir ( dp );
  pkv_QuickSort ( nit, (void *)lbox, xge_LessString, xge_SwapString );
  return true;
} /*xge_SetupFileList*/

/* ///////////////////////////////////////////////////////////////////////// */
boolean xge_SetupDirList ( xge_listbox *lbox, const char *dir,
                           const char *filter, const char *prevdir )
{
  DIR           *dp;
  struct dirent *entry;
  struct stat   statbuf;
  int           lsum, nit, lgt, lgt1;
  char          *newdir;

  xge_ClearListBox ( lbox );
  if ( (dp = opendir ( dir ) ) == NULL ) {
    return false;
  }
  chdir ( dir );
  lsum = nit = 0;
  while ( (entry = readdir(dp)) != NULL ) {
    stat ( entry->d_name, &statbuf );
    if ( (S_ISDIR(statbuf.st_mode)) ) {
      if ( strcmp ( entry->d_name, ".") ) {
        if ( xge_FilterMatches ( entry->d_name, filter ) ) {
          nit ++;
          lsum += strlen ( entry->d_name ) + 1;
        }
      }
    }
  }
  closedir ( dp );
  lbox->itemind = malloc ( (nit+1)*sizeof(int) );
  lbox->itemstr = malloc ( lsum );
  if ( !lbox->itemind || !lbox->itemstr ) {
    xge_ClearListBox ( lbox );
    return false;
  }
  lbox->nitems = (short)nit;
  lsum = nit = 0;
  dp = opendir ( dir );
  while ( (entry = readdir(dp)) != NULL ) {
    stat ( entry->d_name, &statbuf );
    if ( (S_ISDIR(statbuf.st_mode)) ) {
      if ( strcmp ( entry->d_name, "." ) ) {
        if ( xge_FilterMatches ( entry->d_name, filter ) ) {
          lgt = strlen ( entry->d_name );
          memcpy ( &lbox->itemstr[lsum], entry->d_name, lgt+1 );
          lbox->itemind[nit] = lsum;
          lsum += lgt+1;
          nit++;
        }
      }
    }
  }
  closedir ( dp );
  pkv_QuickSort ( nit, (void*)lbox, xge_LessString, xge_SwapString );
        /* if possible, setup the current list position to the */
        /* directory just left */
  if ( prevdir ) {
    newdir = malloc ( 1024 );
    if ( newdir ) {
      getcwd ( newdir, 1023 );
      lgt = strlen ( newdir );
      lgt1 = strlen ( prevdir );
      if ( lgt1 > lgt && !strncmp ( prevdir, newdir, lgt ) ) {
        memcpy ( newdir, &prevdir[lgt+1], lgt1-lgt );
        for ( nit = lbox->nitems-1;  nit;  nit-- )
          if ( !strcmp ( newdir, &lbox->itemstr[lbox->itemind[nit]] ) )
            break;
        lbox->currentitem = nit;
        lbox->fditem = nit - lbox->dlistnpos/2;
        lbox->fditem = min ( lbox->fditem, lbox->nitems-lbox->dlistnpos );
        lbox->fditem = max ( lbox->fditem, 0 );
      }
      free ( newdir );
    }
  }
  return true;
} /*xge_SetupDirList*/

