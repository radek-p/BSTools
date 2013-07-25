
# /////////////////////////////////////////////////////////////////////
# This file is a part of the BSTools package
# written by Przemyslaw Kiciak
# /////////////////////////////////////////////////////////////////////
# (C) Copyright by Przemyslaw Kiciak, 2005, 2009
# this package is distributed under the terms of the
# Lesser GNU Public License, see the file COPYING.LIB
# /////////////////////////////////////////////////////////////////////

HOME = /home/przemek
SRCDIR = .
INCDIR = ../../include
LIBDIR = ../../lib

include ../../options.mak

OBJECTS = ps_out1.o ps_out2.o ps_out3.o ps_out4.o ps_out5.o \
  ps_bmp1.o ps_bmp2.o ps_bmp3.o ps_bmp4.o \
  ps_lib0.o ps_lib1.o ps_lib2.o
TARGET = $(LIBDIR)/libpsout.a

default: $(TARGET)

$(TARGET): $(OBJECTS)
	$(AR) $(ARFLAGS) $(TARGET) $(OBJECTS)

.c.o:
	$(CC) $(CFLAGS) -I$(INCDIR) -c $< -o $*.o

ps_out1.o: ps_out1.c psprivate.h $(INCDIR)/psout.h $(INCDIR)/pkgeom.h
ps_out2.o: ps_out2.c psprivate.h $(INCDIR)/psout.h $(INCDIR)/pkgeom.h
ps_out3.o: ps_out3.c psprivate.h $(INCDIR)/psout.h $(INCDIR)/pkgeom.h
ps_out4.o: ps_out4.c psprivate.h $(INCDIR)/psout.h $(INCDIR)/pkgeom.h
ps_out5.o: ps_out5.c psprivate.h $(INCDIR)/psout.h $(INCDIR)/pkgeom.h
ps_bmp1.o: ps_bmp1.c psprivate.h $(INCDIR)/psout.h $(INCDIR)/pkgeom.h
ps_bmp2.o: ps_bmp2.c psprivate.h $(INCDIR)/psout.h $(INCDIR)/pkgeom.h
ps_bmp3.o: ps_bmp3.c psprivate.h $(INCDIR)/psout.h $(INCDIR)/pkgeom.h
ps_bmp4.o: ps_bmp4.c psprivate.h $(INCDIR)/psout.h $(INCDIR)/pkgeom.h
ps_lib0.o: ps_lib0.c psprivate.h $(INCDIR)/psout.h $(INCDIR)/pkgeom.h
ps_lib1.o: ps_lib1.c psprivate.h $(INCDIR)/psout.h $(INCDIR)/pkgeom.h
ps_lib2.o: ps_lib2.c psprivate.h $(INCDIR)/psout.h $(INCDIR)/pkgeom.h

clean:
	rm -f $(OBJECTS)

mrproper:
	rm -f $(OBJECTS)
	rm -f $(TARGET)
