
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
#include "xgedit.h"

#include "spltemplate.h"
#include "edtempwidgets.h"
#include "edtemplate.h"

int win0, win1;
xge_widget *menu0, *menu1, *menu2;
xge_widget *menu01alist;
xge_widget *popup00, *popup01, *popup02;
xge_widget *status0, *status0sw;

xge_widget *menu3, *menu4, *menu5;
xge_widget *menu14alist;
xge_widget *knwind;
xge_widget *status1, *status1sw;

xge_listbox   filelist, dirlist;
xge_string_ed filename_editor;
char current_directory[MAX_PATH_LGT+1] = "";
char filename[MAX_FILENAME_LGT+1] = "";
char file_filter[2] = "*";

boolean status_0_sw = false;
boolean status_1_sw = false;

char statustext0[STATUS_LINE_LENGTH+1] = "";
char statustext1[STATUS_LINE_LENGTH+1] = "";


int main ( int argc, char *argv[] )
{
  xge_Init ( argc, argv, CallBack, NULL );
  init_edwin ();
  xge_MessageLoop ();
  destroy_edwin ();
  xge_Cleanup ();
  exit ( 0 );
} /*main*/

