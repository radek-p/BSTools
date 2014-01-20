
# /////////////////////////////////////////////////////////////////////
# This file is a part of the BSTools package
# written by Przemyslaw Kiciak
# /////////////////////////////////////////////////////////////////////
# (C) Copyright by Przemyslaw Kiciak, 2008, 2014
# this package is distributed under the terms of the
# Lesser GNU Public License, see the file COPYING.LIB
# /////////////////////////////////////////////////////////////////////

HOME = /home/przemek
SRCDIR = .
INCDIR = ../../include
LIBDIR = -L/usr/X11R6/lib -L../../lib
BINDIR = ../bin

include ../../options.mak

OBJECTS = $(SRCDIR)/edwidgets.o $(SRCDIR)/edpolicz01.o $(SRCDIR)/edpolicz02.o \
  $(SRCDIR)/edpolicz03.o $(SRCDIR)/edpolicz04.o $(SRCDIR)/edpolicz05.o \
  $(SRCDIR)/edpolicz06.o $(SRCDIR)/edpolicz07.o $(SRCDIR)/edpolicz08.o \
  $(SRCDIR)/edpolicz09.o $(SRCDIR)/edpolicz10.o \
  $(SRCDIR)/splhole00.o $(SRCDIR)/splhole01.o $(SRCDIR)/splhole02.o \
  $(SRCDIR)/splhole03.o $(SRCDIR)/splhole04.o $(SRCDIR)/splhole05.o \
  $(SRCDIR)/splhole06.o $(SRCDIR)/splhole07.o \
  $(SRCDIR)/datagend.o  $(SRCDIR)/drawbezd.o $(SRCDIR)/render.o

LIBS = -lraybez -lxgedit -leghole -lbsfile -lmultibs -lcamera \
  -lpsout -lpkgeom -lpknum -lpkvaria -lX11 -lXext -lpthread -lm

EDINCL = $(SRCDIR)/splhole.h $(SRCDIR)/edpolicz.h $(SRCDIR)/edwidgets.h \
  $(SRCDIR)/datagend.h $(SRCDIR)/drawbezd.h $(SRCDIR)/render.h

default: $(BINDIR)/policz

.c.o:
	$(CC) $(CFLAGS) -I$(INCDIR) -c $< -o $*.o

$(SRCDIR)/edwidgets.o: $(SRCDIR)/edwidgets.c $(EDINCL)

$(SRCDIR)/edpolicz01.o: $(SRCDIR)/edpolicz01.c $(EDINCL)

$(SRCDIR)/edpolicz02.o: $(SRCDIR)/edpolicz02.c $(EDINCL)

$(SRCDIR)/edpolicz03.o: $(SRCDIR)/edpolicz03.c $(EDINCL)

$(SRCDIR)/edpolicz04.o: $(SRCDIR)/edpolicz04.c $(EDINCL)

$(SRCDIR)/edpolicz05.o: $(SRCDIR)/edpolicz05.c $(EDINCL)

$(SRCDIR)/edpolicz06.o: $(SRCDIR)/edpolicz06.c $(EDINCL)

$(SRCDIR)/edpolicz07.o: $(SRCDIR)/edpolicz07.c $(EDINCL)

$(SRCDIR)/edpolicz08.o: $(SRCDIR)/edpolicz08.c $(EDINCL)

$(SRCDIR)/edpolicz09.o: $(SRCDIR)/edpolicz09.c $(EDINCL)

$(SRCDIR)/edpolicz10.o: $(SRCDIR)/edpolicz10.c $(EDINCL)

$(SRCDIR)/splhole00.o: $(SRCDIR)/splhole00.c $(EDINCL)

$(SRCDIR)/splhole01.o: $(SRCDIR)/splhole01.c $(EDINCL)

$(SRCDIR)/splhole02.o: $(SRCDIR)/splhole02.c $(EDINCL)

$(SRCDIR)/splhole03.o: $(SRCDIR)/splhole03.c $(EDINCL)

$(SRCDIR)/splhole04.o: $(SRCDIR)/splhole04.c $(EDINCL)

$(SRCDIR)/splhole05.o: $(SRCDIR)/splhole05.c $(EDINCL)

$(SRCDIR)/splhole06.o: $(SRCDIR)/splhole06.c $(EDINCL)

$(SRCDIR)/splhole07.o: $(SRCDIR)/splhole07.c $(EDINCL)

$(SRCDIR)/datagend.o: $(SRCDIR)/datagend.c $(EDINCL)

$(SRCDIR)/drawbezd.o: $(SRCDIR)/drawbezd.c $(EDINCL)

$(SRCDIR)/render.o: $(SRCDIR)/render.c $(SRCDIR)/render.h

$(BINDIR)/policz: $(SRCDIR)/edpolicz.c $(EDINCL) $(OBJECTS)
	$(CC) $(LDFLAGS) $(CFLAGS) -I$(INCDIR) $(SRCDIR)/edpolicz.c \
  $(OBJECTS) $(LIBDIR) $(LIBS) -o $(BINDIR)/policz
	rm -f $(SRCDIR)/policz
	ln -s $(BINDIR)/policz $(SRCDIR)/policz

clean:
	rm -f $(SRCDIR)/*.o

mrproper:
	rm -f $(SRCDIR)/*.o
	rm -f $(BINDIR)/policz $(SRCDIR)/policz

