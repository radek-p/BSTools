
/* //////////////////////////////////////////////////// */
/* This file is a part of the BSTools procedure package */
/* written by Przemyslaw Kiciak.                        */
/* //////////////////////////////////////////////////// */

#define MAX_PATH_LGT     1024
#define MAX_PATH_SHRT      63
#define MAX_FILENAME_LGT   64

#define SIDE_MENU_WIDTH    80

#define WDT 10240/*8192*/
#define HGH 10240/*8192*/

#define NPARA 12

/* widget identifiers */
#define MENU0             0x100
#define M0_BTN_FILE       (MENU0+1)
#define M0_SW_SELECT      (MENU0+2)
#define M0_SW_TRANSLATE   (MENU0+3)
#define M0_SW_SCALE       (MENU0+4)
#define M0_SW_ROTATE      (MENU0+5)
#define M0_SW_SHEAR       (MENU0+6)
#define M0_SW_PANZOOM     (MENU0+7)
#define M0_SW_COORD       (MENU0+8)
#define M0_SW_ANTIALIAS   (MENU0+9)
#define M0_SWR0           (MENU0+10)

#define POPUP0            0x200
#define P0_BTN_OPEN       (POPUP0+1)
#define P0_BTN_SAVE       (POPUP0+2)
#define P0_BTN_EXIT       (POPUP0+3)

#define POPUP1            0x300
#define P1_BTN_OPEN       (POPUP1+1)
#define P1_BTN_CANCEL     (POPUP1+2)
#define P1_LB_DIRLIST     (POPUP1+3)
#define P1_LB_FILELIST    (POPUP1+4)

#define POPUP2            0x400
#define P2_BTN_SAVE       (POPUP2+1)
#define P2_BTN_CANCEL     (POPUP2+2)
#define P2_LB_DIRLIST     (POPUP2+3)
/*#define P2_LB_FILELIST    (POPUP2+4)*/
#define P2_TXT_SAVE_AS    (POPUP2+5)
#define P2_TXTED_FILENAME (POPUP2+6)

#define CWIND0            0x500


void SetupBitTranslation ( void );
boolean InitMyBitmap ( void );
void SetMyPixel ( short x, short y, char on );
char GetMyPixel ( short x, short y );
void SetMySubpixelAA ( short x, short y, char on );
char GetMySubpixelAA ( short x, short y );
char GetMyPixelAA ( short x, short y );

void MakeTr ( point2d *p, point2d *q, trans2d *tr,
              trans2d *trb, trans2d *traa );
void MakeAllTr ( void );
void IterIFS ( void );

void DrawParallelogram ( xge_2Dwind *_2Dwin, point2d *para,
                         xgecolour_int c, boolean *cpmk );
void DrawBitmap ( xge_widget *er );
void DrawBitmapAA ( xge_widget *er );
void RysujOkno ( xge_widget *er, boolean onscreen );

boolean FindNearestPoint ( xge_widget *er, short x, short y, short mindist );
void SetCPoint ( xge_widget *er, short x, short y );

void SelectPoints ( xge_widget *er, short x0, short x1, short y0, short y1,
                    boolean select );
void SavePoints ( void );
void TransformPoints ( xge_widget *er );

void SetCurrentDir ( void );
boolean OpenTheFile ( void );
boolean SaveTheFile ( void );
boolean ChangeDir ( xge_widget *popup,
                    xge_listbox *dirlist, xge_listbox *filelist );

int CallBack ( xge_widget *er, int msg, int key, short x, short y );

void init_edwin ( void );
void destroy_edwin ( void );

int main ( int argc, char *argv[] );

