
/* //////////////////////////////////////////////////// */
/* This file is a part of the BSTools procedure package */
/* written by Przemyslaw Kiciak.                        */
/* //////////////////////////////////////////////////// */

#define IMWIDTH  1024
#define IMHEIGHT  768

#define MAX_PATH_LGT     1024
#define MAX_PATH_SHRT      64
#define MAX_FILENAME_LGT   64

#define POPUP01      100
#define btnP01LOAD   101
#define btnP01EXIT   102
#define lb01DIRLIST  103
#define lb01FILELIST 104

#define MAXSLIDES   1001
#define MAXCHAPTERS  101

#define ATIME 0.5  /* 0.5 s */
#define ANIMASTEP 16

void ReadSlOptFile ( void );
void ReadChapterFile ( void );
void MakeFileList ( void );

void ClearImage ( void );
boolean ReadImage ( char *fn );

void DrawImageWin ( xge_widget *er, boolean onscreen );

void InitAnimation ( int dir );
void ContinueAnimation ( void );
void ShowNextPicture ( boolean anim );
void ShowPreviousPicture ( boolean anim );
void ShowNextChapter ( void );
void ShowPreviousChapter ( void );

void DrawBackground ( void );

void ChangeDir ( void );

boolean ProcessKey ( int key );
boolean ImageWinMsg ( xge_widget *er, int msg, int key, short x, short y );
int CallBack ( xge_widget *er, int msg, int key, short x, short y );

void init_win ( void );
void destroy_win ( void );
int main ( int argc, char *argv[] );

