
# /////////////////////////////////////////////////////////////////////
# This file is a part of the BSTools package
# written by Przemyslaw Kiciak
# /////////////////////////////////////////////////////////////////////
# (C) Copyright by Przemyslaw Kiciak, 2007, 2014
# this package is distributed under the terms of the
# Lesser GNU Public License, see the file COPYING.LIB
# /////////////////////////////////////////////////////////////////////

HOME = /home/przemek
SRCDIR = .
INCDIR = ../../include
LIBDIR = -L/usr/X11R6/lib -L../../lib
BINDIR = ../bin

include ../../options.mak

#LDFLAGS = -static
LDFLAGS =

OBJECTS = $(SRCDIR)/spl3d00.o $(SRCDIR)/spl3d01.o $(SRCDIR)/spl3d02.o \
  $(SRCDIR)/spl3d03.o $(SRCDIR)/spl3d04.o $(SRCDIR)/spl3d05.o \
  $(SRCDIR)/spl3d06.o $(SRCDIR)/spl3d07.o $(SRCDIR)/spl3d08.o \
  $(SRCDIR)/spl3d09.o $(SRCDIR)/spl3d10.o $(SRCDIR)/spl3d11.o \
  $(SRCDIR)/spl3d12.o $(SRCDIR)/spl3d13.o $(SRCDIR)/spl3d14.o \
  $(SRCDIR)/spl3d15.o $(SRCDIR)/spl3d16.o $(SRCDIR)/spl3d17.o \
  $(SRCDIR)/spl3d18.o $(SRCDIR)/spl3d19.o \
  $(SRCDIR)/ed3ds01.o $(SRCDIR)/ed3ds02.o $(SRCDIR)/ed3ds03.o \
  $(SRCDIR)/ed3ds04.o $(SRCDIR)/ed3ds05.o $(SRCDIR)/ed3ds06.o \
  $(SRCDIR)/ed3ds07.o $(SRCDIR)/ed3ds08.o $(SRCDIR)/ed3ds09.o \
  $(SRCDIR)/ed3ds10.o $(SRCDIR)/ed3ds11.o $(SRCDIR)/ed3ds12.o \
  $(SRCDIR)/ed3ds13.o \
  $(SRCDIR)/ed3dswidgets.o $(SRCDIR)/render.o

EDINCL = $(SRCDIR)/spl3d.h $(SRCDIR)/ed3ds.h $(SRCDIR)/ed3dswidgets.h \
  $(SRCDIR)/render.h $(SRCDIR)/pomnijipc.h

OBJECTS_PROC = pomnij_proc01.o pomnij_proc02.o pomnij_proc03.o \
  pomnij_proc04.o pomnij_proc05.o proc_regmem.o

EDINCL_PROC = $(SRCDIR)/pomnijipc.h $(SRCDIR)/pomnij_proc.h \
  $(SRCDIR)/proc_regmem.h

LIBS = -lxgedit -lraybez -lg1blending -lg2blending -lbsfile -lmultibs \
  -lcamera -lpkgeom -lpknum -lpkvaria -lX11 -lXext -lpthread -lm

default: $(BINDIR)/pomnij $(BINDIR)/pomnij_proc

.c.o:
	$(CC) $(CFLAGS) -I$(INCDIR) -c $< -o $*.o

$(SRCDIR)/spl3d00.o: $(SRCDIR)/spl3d00.c $(SRCDIR)/spl3d.h

$(SRCDIR)/spl3d01.o: $(SRCDIR)/spl3d01.c $(SRCDIR)/spl3d.h

$(SRCDIR)/spl3d02.o: $(SRCDIR)/spl3d02.c $(SRCDIR)/spl3d.h

$(SRCDIR)/spl3d03.o: $(SRCDIR)/spl3d03.c $(EDINCL)

$(SRCDIR)/spl3d04.o: $(SRCDIR)/spl3d04.c $(SRCDIR)/spl3d.h

$(SRCDIR)/spl3d05.o: $(SRCDIR)/spl3d05.c $(SRCDIR)/spl3d.h

$(SRCDIR)/spl3d06.o: $(SRCDIR)/spl3d06.c $(SRCDIR)/spl3d.h

$(SRCDIR)/spl3d07.o: $(SRCDIR)/spl3d07.c $(SRCDIR)/spl3d.h

$(SRCDIR)/spl3d08.o: $(SRCDIR)/spl3d08.c $(SRCDIR)/spl3d.h

$(SRCDIR)/spl3d09.o: $(SRCDIR)/spl3d09.c $(SRCDIR)/spl3d.h

$(SRCDIR)/spl3d10.o: $(SRCDIR)/spl3d10.c $(SRCDIR)/spl3d.h

$(SRCDIR)/spl3d11.o: $(SRCDIR)/spl3d11.c $(SRCDIR)/spl3d.h

$(SRCDIR)/spl3d12.o: $(SRCDIR)/spl3d12.c $(EDINCL)

$(SRCDIR)/spl3d13.o: $(SRCDIR)/spl3d13.c $(SRCDIR)/spl3d.h

$(SRCDIR)/spl3d14.o: $(SRCDIR)/spl3d14.c $(SRCDIR)/spl3d.h

$(SRCDIR)/spl3d15.o: $(SRCDIR)/spl3d15.c $(SRCDIR)/spl3d.h

$(SRCDIR)/spl3d16.o: $(SRCDIR)/spl3d16.c $(SRCDIR)/spl3d.h

$(SRCDIR)/spl3d17.o: $(SRCDIR)/spl3d17.c $(SRCDIR)/spl3d.h $(SRCDIR)/pomnijipc.h

$(SRCDIR)/spl3d18.o: $(SRCDIR)/spl3d18.c $(SRCDIR)/spl3d.h

$(SRCDIR)/spl3d19.o: $(SRCDIR)/spl3d19.c $(SRCDIR)/spl3d.h

$(SRCDIR)/ed3ds01.o: $(SRCDIR)/ed3ds01.c $(EDINCL)

$(SRCDIR)/ed3ds02.o: $(SRCDIR)/ed3ds02.c $(EDINCL)

$(SRCDIR)/ed3ds03.o: $(SRCDIR)/ed3ds03.c $(EDINCL)

$(SRCDIR)/ed3ds04.o: $(SRCDIR)/ed3ds04.c $(EDINCL)

$(SRCDIR)/ed3ds05.o: $(SRCDIR)/ed3ds05.c $(EDINCL)

$(SRCDIR)/ed3ds06.o: $(SRCDIR)/ed3ds06.c $(EDINCL)

$(SRCDIR)/ed3ds07.o: $(SRCDIR)/ed3ds07.c $(EDINCL)

$(SRCDIR)/ed3ds08.o: $(SRCDIR)/ed3ds08.c $(EDINCL)

$(SRCDIR)/ed3ds09.o: $(SRCDIR)/ed3ds09.c $(EDINCL)

$(SRCDIR)/ed3ds10.o: $(SRCDIR)/ed3ds10.c $(EDINCL)

$(SRCDIR)/ed3ds11.o: $(SRCDIR)/ed3ds11.c $(EDINCL)

$(SRCDIR)/ed3ds12.o: $(SRCDIR)/ed3ds12.c $(EDINCL)

$(SRCDIR)/ed3dswidgets.o: $(SRCDIR)/ed3dswidgets.c $(EDINCL)

$(SRCDIR)/render.o: $(SRCDIR)/render.c $(SRCDIR)/render.h

$(BINDIR)/pomnij: $(SRCDIR)/ed3ds.c $(EDINCL) $(OBJECTS)
	$(CC) $(LDFLAGS) $(CFLAGS) -I$(INCDIR) $(SRCDIR)/ed3ds.c $(OBJECTS) \
  $(LIBDIR) $(LIBS) -o $(BINDIR)/pomnij
	rm -f $(SRCDIR)/pomnij
	ln -s $(BINDIR)/pomnij $(SRCDIR)/pomnij

$(SRCDIR)/proc_regmem.o: $(SRCDIR)/proc_regmem.c $(EDINCL_PROC)

pomnij_proc01.o: $(SRCDIR)/pomnij_proc01.c $(EDINCL_PROC)

pomnij_proc02.o: $(SRCDIR)/pomnij_proc02.c $(EDINCL_PROC)

pomnij_proc03.o: $(SRCDIR)/pomnij_proc03.c $(EDINCL_PROC)

pomnij_proc04.o: $(SRCDIR)/pomnij_proc04.c $(EDINCL_PROC)

pomnij_proc05.o: $(SRCDIR)/pomnij_proc05.c $(EDINCL_PROC)

$(BINDIR)/pomnij_proc: $(SRCDIR)/pomnij_proc.c $(EDINCL_PROC) $(OBJECTS_PROC)
	$(CC) $(LDFLAGS) $(CFLAGS) -I$(INCDIR) $(SRCDIR)/pomnij_proc.c \
  $(OBJECTS_PROC) $(LIBDIR) $(LIBS) -o $(BINDIR)/pomnij_proc
	rm -f $(SRCDIR)/pomnij_proc
	ln -s $(BINDIR)/pomnij_proc $(SRCDIR)/pomnij_proc

clean:
	rm -f *.o

mrproper:
	rm -f *.o
	rm -f $(BINDIR)/pomnij $(BINDIR)/pomnij_proc \
  $(SRCDIR)/pomnij $(SRCDIR)/pomnij_proc

