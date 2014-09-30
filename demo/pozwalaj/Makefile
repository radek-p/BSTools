
# /////////////////////////////////////////////////////////////////////
# This file is a part of the BSTools package
# written by Przemyslaw Kiciak
# /////////////////////////////////////////////////////////////////////
# (C) Copyright by Przemyslaw Kiciak, 2010, 2014
# this package is distributed under the terms of the
# Lesser GNU Public License, see the file COPYING.LIB
# /////////////////////////////////////////////////////////////////////

HOME = /home/przemek
SRCDIR = .
BINDIR = ../bin
#BINDIR = .
INCDIR = ../../include
LIBDIR = -L/usr/X11R6/lib -L../../lib

include ../../options.mak
#include ../../gcc.debug.mak

CC = gcc
#CFLAGS = -g -ansi -pedantic -Wall
#LDFLAGS = -static
LDFLAGS =

LIBS = -lxgedit -lbsfile -lg1blending -lg2blending -leghole -lbsmesh -lmengerc \
  -lraybez -lmultibs -lcamera -lpkgeom -lpknum -lpkvaria \
  -lGLU -lGL -lX11 -lXext -lpthread -lm

OBJECTS = pozwalaj00.o pozwalaj01.o pozwalaj02.o statustext.o \
  top00menu.o side00menu.o side00bmenu.o side00emenu.o bottom00menu.o \
  geom00menu.o popup00menu.o popup01menu.o popup02menu.o popup03menu.o \
  top10menu.o side10menu.o side10bezc.o side10bezp.o \
  side10bsc.o side10bsp01.o side10bsp02.o side10bsp03.o \
  side10bsm01.o side10bsm02.o side10bsh.o bottom10menu.o geom10menu.o \
  popup10menu.o popup11menu.o popup12menu.o popup13menu.o popup14menu.o \
  popup15menu.o \
  arrowbtn.o setupobjmenu.o widgets.o calc.o \
  editor01.o editor02.o editor03.o editor04.o editor05.o editor06.o \
  editor07.o editor08.o editor09.o \
  edbezcurve.o edbscurve.o edbezpatch.o \
  edbspatch00.o edbspatch01.o edbspatch02.o edbspatch03.o edbspatch04.o \
  edbspatch05.o edbspatch06.o edbspatch07.o edbspatch08.o edbspatch09.o \
  edbspatch10.o \
  edbsmesh00.o edbsmesh01.o edbsmesh02.o edbsmesh03.o edbsmesh04.o \
  edbsmesh05.o edbsmesh06.o edbsmesh07.o edbsmesh08.o edbsmesh09.o \
  edbsmesh10.o edbsmesh11.o edbsmesh12.o \
  edbshole.o edrwfile.o \
  edlight.o render.o \
  pozwalajipc.o

INCL = $(SRCDIR)/pozwalaj.h $(SRCDIR)/pozwalaj0.h

EDINC = $(SRCDIR)/editor.h $(SRCDIR)/edcolours.h \
  $(SRCDIR)/editor_bezc.h $(SRCDIR)/editor_bsc.h $(SRCDIR)/editor_bezp.h \
  $(SRCDIR)/editor_bsp.h $(SRCDIR)/editor_bsm.h $(SRCDIR)/editor_bsh.h \
  $(SRCDIR)/arrowbtn.h

PROC_OBJECTS = pozwalaj_proc00.o pozwalaj_proc01.o pozwalaj_proc02.o \
   pozwalaj_proc03.o pozwalaj_proc04.o pozwalaj_proc05.o pozwalaj_proc06.o \
   pozwalaj_proc07.o pozwalaj_proc08.o


default: $(BINDIR)/pozwalaj $(BINDIR)/pozwalaj_proc

.c.o:
	$(CC) $(CFLAGS) -I$(INCDIR) -c $< -o $*.o

$(BINDIR)/pozwalaj: $(SRCDIR)/pozwalaj.c $(INCL) \
  $(EDINC) $(SRCDIR)/widgets.h $(OBJECTS)
	$(CC) $(LDFLAGS) $(CFLAGS) -I$(INCDIR) $(SRCDIR)/pozwalaj.c \
  $(OBJECTS) $(LIBDIR) $(LIBS) -o $(BINDIR)/pozwalaj
	rm -f $(SRCDIR)/pozwalaj
	ln -s $(BINDIR)/pozwalaj $(SRCDIR)/pozwalaj

pozwalaj00.o: $(SRCDIR)/pozwalaj00.c $(INCL) \
  $(EDINC) $(SRCDIR)/widgets.h

pozwalaj01.o: $(SRCDIR)/pozwalaj01.c $(INCL) $(SRCDIR)/pozwalajipc.h \
  $(SRCDIR)/render.h $(EDINC) $(SRCDIR)/widgets.h

pozwalaj02.o: $(SRCDIR)/pozwalaj02.c $(INCL) \
  $(EDINC) $(SRCDIR)/widgets.h

statustext.o: $(SRCDIR)/statustext.c $(INCL) \
  $(EDINC) $(SRCDIR)/widgets.h

arrowbtn.o: $(SRCDIR)/arrowbtn.c $(SRCDIR)/arrowbtn.h

calc.o: $(SRCDIR)/calc.c $(SRCDIR)/calc.h

top00menu.o: top00menu.c $(INCL) \
  $(EDINC) $(SRCDIR)/widgets.h

side00menu.o: side00menu.c $(INCL) \
  $(EDINC) $(SRCDIR)/widgets.h

side00bmenu.o: side00bmenu.c $(INCL) $(SRCDIR)/edlight.h \
  $(SRCDIR)/render.h $(EDINC) $(SRCDIR)/widgets.h

side00emenu.o: side00emenu.c $(INCL) \
  $(EDINC) $(SRCDIR)/widgets.h

bottom00menu.o: bottom00menu.c $(INCL) $(SRCDIR)/calc.h \
  $(EDINC) $(SRCDIR)/widgets.h

geom00menu.o: geom00menu.c $(INCL) \
  $(EDINC) $(SRCDIR)/widgets.h $(SRCDIR)/edlight.h $(SRCDIR)/render.h

popup00menu.o: popup00menu.c $(INCL) \
  $(EDINC) $(SRCDIR)/widgets.h

popup01menu.o: popup01menu.c $(INCL) \
  $(EDINC) $(SRCDIR)/widgets.h

popup02menu.o: popup02menu.c $(INCL) \
  $(EDINC) $(SRCDIR)/widgets.h

popup03menu.o: popup03menu.c $(INCL) \
  $(EDINC) $(SRCDIR)/widgets.h

top10menu.o: top10menu.c $(INCL) \
  $(EDINC) $(SRCDIR)/widgets.h

side10menu.o: side10menu.c $(INCL) \
  $(EDINC) $(SRCDIR)/widgets.h

side10bezc.o: side10bezc.c $(INCL) \
  $(SRCDIR)/editor.h $(SRCDIR)/widgets.h

side10bezp.o: side10bezp.c $(INCL) \
  $(EDINC) $(SRCDIR)/widgets.h

side10bsc.o: side10bsc.c $(INCL) \
  $(EDINC) $(SRCDIR)/widgets.h

side10bsp01.o: side10bsp01.c $(INCL) \
  $(EDINC) $(SRCDIR)/widgets.h

side10bsp02.o: side10bsp02.c $(INCL) \
  $(EDINC) $(SRCDIR)/widgets.h

side10bsp03.o: side10bsp03.c $(INCL) \
  $(EDINC) $(SRCDIR)/widgets.h

side10bsm01.o: side10bsm01.c $(INCL) $(SRCDIR)/pozwalajipc.h \
  $(EDINC) $(SRCDIR)/widgets.h

side10bsm02.o: side10bsm02.c $(INCL) $(SRCDIR)/pozwalajipc.h \
  $(EDINC) $(SRCDIR)/widgets.h

side10bsh.o: side10bsh.c $(INCL) \
  $(EDINC) $(SRCDIR)/widgets.h

bottom10menu.o: bottom10menu.c $(INCL) $(SRCDIR)/calc.h \
  $(EDINC) $(SRCDIR)/widgets.h

geom10menu.o: geom10menu.c $(INCL) \
  $(EDINC) $(SRCDIR)/widgets.h

popup10menu.o: popup10menu.c $(INCL) \
  $(EDINC) $(SRCDIR)/widgets.h

popup11menu.o: popup11menu.c $(INCL) \
  $(EDINC) $(SRCDIR)/widgets.h

popup12menu.o: popup12menu.c $(INCL) \
  $(EDINC) $(SRCDIR)/widgets.h $(SRCDIR)/render.h

popup13menu.o: popup13menu.c $(INCL) \
  $(EDINC) $(SRCDIR)/widgets.h $(SRCDIR)/render.h

popup14menu.o: popup14menu.c $(INCL) \
  $(EDINC) $(SRCDIR)/widgets.h $(SRCDIR)/render.h

popup15menu.o: popup15menu.c $(INCL) \
  $(EDINC) $(SRCDIR)/widgets.h

setupobjmenu.o: setupobjmenu.c $(INCL) \
  $(EDINC) $(SRCDIR)/widgets.h

widgets.o: $(SRCDIR)/widgets.c $(SRCDIR)/widgets.h

editor01.o: $(SRCDIR)/editor01.c $(EDINC)

editor02.o: $(SRCDIR)/editor02.c \
  $(SRCDIR)/editor.h $(SRCDIR)/editor_bsm.h

editor03.o: $(SRCDIR)/editor03.c $(EDINC)

editor04.o: $(SRCDIR)/editor04.c $(EDINC)

editor05.o: $(SRCDIR)/editor05.c $(EDINC)

editor06.o: $(SRCDIR)/editor06.c $(EDINC)

editor07.o: $(SRCDIR)/editor07.c $(EDINC)

editor08.o: $(SRCDIR)/editor08.c $(EDINC)

editor09.o: $(SRCDIR)/editor09.c $(EDINC)

edbezcurve.o: $(SRCDIR)/edbezcurve.c \
  $(SRCDIR)/editor.h $(SRCDIR)/editor_bezc.h $(SRCDIR)/render.h

edbscurve.o: $(SRCDIR)/edbscurve.c \
  $(SRCDIR)/editor.h $(SRCDIR)/editor_bsc.h $(SRCDIR)/render.h

edbezpatch.o: $(SRCDIR)/edbezpatch.c \
  $(SRCDIR)/editor.h $(SRCDIR)/editor_bezp.h $(SRCDIR)/render.h

edbspatch00.o: $(SRCDIR)/edbspatch00.c \
  $(SRCDIR)/editor.h $(SRCDIR)/editor_bsp.h

edbspatch01.o: $(SRCDIR)/edbspatch01.c \
  $(SRCDIR)/editor.h $(SRCDIR)/editor_bsp.h

edbspatch02.o: $(SRCDIR)/edbspatch02.c \
  $(SRCDIR)/editor.h $(SRCDIR)/editor_bsp.h

edbspatch03.o: $(SRCDIR)/edbspatch03.c \
  $(SRCDIR)/editor.h $(SRCDIR)/editor_bsp.h

edbspatch04.o: $(SRCDIR)/edbspatch04.c \
  $(SRCDIR)/editor.h $(SRCDIR)/editor_bsp.h

edbspatch05.o: $(SRCDIR)/edbspatch05.c \
  $(SRCDIR)/editor.h $(SRCDIR)/editor_bsp.h

edbspatch06.o: $(SRCDIR)/edbspatch06.c \
  $(SRCDIR)/editor.h $(SRCDIR)/editor_bsp.h $(SRCDIR)/render.h

edbspatch07.o: $(SRCDIR)/edbspatch07.c \
  $(SRCDIR)/editor.h $(SRCDIR)/editor_bsp.h $(SRCDIR)/render.h

edbspatch08.o: $(SRCDIR)/edbspatch08.c \
  $(SRCDIR)/editor.h $(SRCDIR)/editor_bsp.h $(SRCDIR)/render.h

edbspatch09.o: $(SRCDIR)/edbspatch09.c \
  $(SRCDIR)/editor.h $(SRCDIR)/editor_bsp.h $(SRCDIR)/render.h

edbspatch10.o: $(SRCDIR)/edbspatch10.c \
  $(SRCDIR)/editor.h $(SRCDIR)/editor_bsp.h $(SRCDIR)/render.h

edbsmesh00.o: $(SRCDIR)/edbsmesh00.c \
  $(SRCDIR)/editor.h $(SRCDIR)/editor_bsm.h

edbsmesh01.o: $(SRCDIR)/edbsmesh01.c \
  $(SRCDIR)/editor.h $(SRCDIR)/editor_bsm.h

edbsmesh02.o: $(SRCDIR)/edbsmesh02.c \
  $(SRCDIR)/editor.h $(SRCDIR)/editor_bsm.h

edbsmesh03.o: $(SRCDIR)/edbsmesh03.c \
  $(SRCDIR)/editor.h $(SRCDIR)/editor_bsm.h

edbsmesh04.o: $(SRCDIR)/edbsmesh04.c \
  $(SRCDIR)/editor.h $(SRCDIR)/editor_bsm.h

edbsmesh05.o: $(SRCDIR)/edbsmesh05.c \
  $(SRCDIR)/editor.h $(SRCDIR)/editor_bsm.h

edbsmesh06.o: $(SRCDIR)/edbsmesh06.c \
  $(SRCDIR)/editor.h $(SRCDIR)/editor_bsm.h

edbsmesh07.o: $(SRCDIR)/edbsmesh07.c \
  $(SRCDIR)/editor.h $(SRCDIR)/editor_bsm.h

edbsmesh08.o: $(SRCDIR)/edbsmesh08.c \
  $(SRCDIR)/editor.h $(SRCDIR)/editor_bsm.h $(SRCDIR)/render.h

edbsmesh09.o: $(SRCDIR)/edbsmesh09.c \
  $(SRCDIR)/editor.h $(SRCDIR)/editor_bsm.h

edbsmesh10.o: $(SRCDIR)/edbsmesh10.c \
  $(SRCDIR)/editor.h $(SRCDIR)/editor_bsm.h

edbsmesh11.o: $(SRCDIR)/edbsmesh11.c \
  $(SRCDIR)/editor.h $(SRCDIR)/editor_bsm.h

edbsmesh12.o: $(SRCDIR)/edbsmesh12.c \
  $(SRCDIR)/editor.h $(SRCDIR)/editor_bsm.h

edbshole.o: $(SRCDIR)/edbshole.c \
  $(SRCDIR)/editor.h $(SRCDIR)/editor_bsh.h

edrwfile.o: $(SRCDIR)/edrwfile.c $(EDINC)

edlight.o: $(SRCDIR)/edlight.c $(EDINC) $(SRCDIR)/widgets.h \
  $(SRCDIR)/edlight.h $(SRCDIR)/render.h $(INCL)

render.o: $(SRCDIR)/render.c $(SRCDIR)/render.h

pozwalaijpc.o: $(SRCDIR)/pozwalajipc.c $(SRCDIR)/widgets.h $(EDINC) \
  $(INCL) $(SRDIR)/pozwalajipc.h


$(BINDIR)/pozwalaj_proc: $(SRCDIR)/pozwalaj_proc.c $(SRCDIR)/pozwalajipc.h \
  $(SRCDIR)/pozwalaj_proc.h $(PROC_OBJECTS)
	$(CC) $(LDFLAGS) $(CFLAGS) -I$(INCDIR) $(SRCDIR)/pozwalaj_proc.c \
  $(PROC_OBJECTS) $(LIBDIR) $(LIBS) -o $(BINDIR)/pozwalaj_proc
	rm -f $(SRCDIR)/pozwalaj_proc
	ln -s $(BINDIR)/pozwalaj_proc $(SRCDIR)/pozwalaj_proc

pozwalaj_proc00.o: $(SRCDIR)/pozwalaj_proc00.c $(SRCDIR)/pozwalajipc.h \
  $(SRCDIR)/pozwalaj_proc.h

pozwalaj_proc01.o: $(SRCDIR)/pozwalaj_proc01.c $(SRCDIR)/pozwalajipc.h \
  $(SRCDIR)/pozwalaj_proc.h

pozwalaj_proc02.o: $(SRCDIR)/pozwalaj_proc02.c $(SRCDIR)/pozwalajipc.h \
  $(SRCDIR)/pozwalaj_proc.h

pozwalaj_proc03.o: $(SRCDIR)/pozwalaj_proc03.c $(SRCDIR)/pozwalajipc.h \
  $(SRCDIR)/pozwalaj_proc.h

pozwalaj_proc04.o: $(SRCDIR)/pozwalaj_proc04.c $(SRCDIR)/pozwalajipc.h \
  $(SRCDIR)/pozwalaj_proc.h

pozwalaj_proc05.o: $(SRCDIR)/pozwalaj_proc05.c $(SRCDIR)/pozwalajipc.h \
  $(SRCDIR)/pozwalaj_proc.h

pozwalaj_proc06.o: $(SRCDIR)/pozwalaj_proc06.c $(SRCDIR)/pozwalajipc.h \
  $(SRCDIR)/pozwalaj_proc.h

pozwalaj_proc07.o: $(SRCDIR)/pozwalaj_proc07.c $(SRCDIR)/pozwalajipc.h \
  $(SRCDIR)/pozwalaj_proc.h

pozwalaj_proc08.o: $(SRCDIR)/pozwalaj_proc08.c $(SRCDIR)/pozwalajipc.h \
  $(SRCDIR)/pozwalaj_proc.h


clean:
	rm -f *.o

mrproper:
	rm -f *.o
	rm -f $(BINDIR)/pozwalaj $(BINDIR)/pozwalaj_proc \
   $(SRCDIR)/pozwalaj $(SRCDIR)/pozwalaj_proc

