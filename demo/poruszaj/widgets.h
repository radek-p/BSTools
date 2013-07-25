
/* //////////////////////////////////////////////////// */
/* This file is a part of the BSTools procedure package */
/* written by Przemyslaw Kiciak.                        */
/* //////////////////////////////////////////////////// */

#define MENUMASK 0xFFFF0000
#define MENUINC  0x00010000

/* top menu */
#define TOPMENU0              MENUINC
#define btnT0FILE        (TOPMENU0+1)
#define btnT0ARTICULATE  (TOPMENU0+2)
#define btnT0OPTIONS     (TOPMENU0+3)
#define btnT0LIGHT       (TOPMENU0+4)
#define btnT0ABOUT       (TOPMENU0+5)

/* side menu: a - articulate */
#define SIDEMENU0a          2*MENUINC
#define slS0aARTPARAM0 (SIDEMENU0a+1)
#define swS0aPANZOOM   (SIDEMENU0a+101)
#define swS0aCOORDS    (SIDEMENU0a+102)
/* b - options */ 
#define SIDEMENU0b          3*MENUINC
#define swS0bPALMLOD1  (SIDEMENU0b+1)
#define swS0bPALMLOD2  (SIDEMENU0b+2)
#define swS0bPALMLOD3  (SIDEMENU0b+3)
#define swS0bSPLHANDS  (SIDEMENU0b+4)
#define swS0bSPLNET    (SIDEMENU0b+5)
#define swS0bANTIALIAS (SIDEMENU0b+6)
/* c - light */
#define SIDEMENU0c          4*MENUINC
#define diS0cANG1      (SIDEMENU0c+1)
#define diS0cANG2      (SIDEMENU0c+2)
#define swS0cUSESPEC   (SIDEMENU0c+3)

/* geometry displaying windows */
#define GEOMMENU0           5*MENUINC

/* bottom menu */
#define BOTTOMMENU0         6*MENUINC
#define knwB0KNOTWIN  (BOTTOMMENU0+1)
#define btnB0PLAY     (BOTTOMMENU0+2)
#define btnB0PREVKEY  (BOTTOMMENU0+3)
#define btnB0SETPOSE  (BOTTOMMENU0+4)
#define btnB0NEXTKEY  (BOTTOMMENU0+5)
#define btnB0EDIT     (BOTTOMMENU0+6)
#define btnB0PLAYSTOP (BOTTOMMENU0+7)
#define swB0PERIODIC  (BOTTOMMENU0+8)

/* File popup */
#define POPUP0              7*MENUINC
#define btnP0OPEN          (POPUP0+1)
#define btnP0SAVE          (POPUP0+2)
#define btnP0SAVEAS        (POPUP0+3)
#define btnP0EXIT          (POPUP0+4)

/* Open file popup */
#define POPUP1              8*MENUINC
#define btnP1OPEN          (POPUP1+1)
#define btnP1CANCEL        (POPUP1+2)
#define lbP1DIRLIST        (POPUP1+3)
#define lbP1FILELIST       (POPUP1+4)

/* Save file as popup */
#define POPUP2              9*MENUINC
#define btnP2SAVE          (POPUP2+1)
#define btnP2CANCEL        (POPUP2+2)
#define lbP2DIRLIST        (POPUP2+3)
#define lbP2FILELIST       (POPUP2+4)
#define txtedP2FILENAME    (POPUP2+5)

/* Exit popup */
#define POPUP3             10*MENUINC
#define btnP3EXIT          (POPUP3+1)
#define btnP3SAVE          (POPUP3+2)
#define btnP3CANCEL        (POPUP3+3)


extern char txtFile[];
extern char txtArticulate[];
extern char txtOptions[];
extern char txtLight[];
extern char txtAbout[];
extern char txtOpen[];
extern char txtSave[];
extern char txtSaveAs[];
extern char txtCancel[];
extern char txtExit[];
extern char txtPanZoom[];
extern char txtCoordinates[];
extern char txtPalmLOD[];
extern char txtSplineHands[];
extern char txtDisplayNets[];
extern char txtEdit[];
extern char txtPlay[];
extern char txtSet[];
extern char txtPeriodic[];
extern char txtSpecular[];
extern char txtAntialias[];

extern char txtMsgReconsider[];

extern char ErrorMessageIncorrectFilename[];

extern char *InfoMsg[];

