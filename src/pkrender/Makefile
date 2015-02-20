
# /////////////////////////////////////////////////////////////////////
# This file is a part of the BSTools package
# written by Przemyslaw Kiciak
# /////////////////////////////////////////////////////////////////////
# (C) Copyright by Przemyslaw Kiciak, 2015
# this package is distributed under the terms of the
# Lesser GNU Public License, see the file COPYING.LIB
# /////////////////////////////////////////////////////////////////////

HOME = /home/przemek
SRCDIR = .
INCDIR = ../../include
LIBDIR = ../../lib

include ../../options.mak

OBJECTS = $(SRCDIR)/pkrender01.o $(SRCDIR)/pkrender02.o $(SRCDIR)/pkrender03.o \
  $(SRCDIR)/pkrender04.o $(SRCDIR)/pkrender05.o $(SRCDIR)/pkrender06.o \
  $(SRCDIR)/pkrender07.o $(SRCDIR)/pkrender08.o $(SRCDIR)/pkrender09.o \
  $(SRCDIR)/pkrender10.o $(SRCDIR)/pkrender11.o $(SRCDIR)/pkrender12.o

INCLUDES = $(INCDIR)/pkrender.h $(SRCDIR)/pkrenderprivate.h    

TARGET = $(LIBDIR)/libpkrender.a

default: $(TARGET)

$(TARGET): $(OBJECTS)
	$(AR) $(ARFLAGS) $(TARGET) $(OBJECTS)

.c.o:
	$(CC) $(CFLAGS) -I$(INCDIR) -c $< -o $*.o

$(SRCDIR)/pkrender01.o: $(SRCDIR)/pkrender01.c $(INCLUDES)

$(SRCDIR)/pkrender02.o: $(SRCDIR)/pkrender02.c $(INCLUDES)

$(SRCDIR)/pkrender03.o: $(SRCDIR)/pkrender03.c $(INCLUDES)

$(SRCDIR)/pkrender04.o: $(SRCDIR)/pkrender04.c $(INCLUDES)

$(SRCDIR)/pkrender05.o: $(SRCDIR)/pkrender05.c $(INCLUDES)

$(SRCDIR)/pkrender06.o: $(SRCDIR)/pkrender06.c $(INCLUDES)

$(SRCDIR)/pkrender07.o: $(SRCDIR)/pkrender07.c $(INCLUDES)

$(SRCDIR)/pkrender08.o: $(SRCDIR)/pkrender08.c $(INCLUDES)

$(SRCDIR)/pkrender09.o: $(SRCDIR)/pkrender09.c $(INCLUDES)

$(SRCDIR)/pkrender10.o: $(SRCDIR)/pkrender10.c $(INCLUDES)

$(SRCDIR)/pkrender11.o: $(SRCDIR)/pkrender11.c $(INCLUDES)

$(SRCDIR)/pkrender12.o: $(SRCDIR)/pkrender12.c $(INCLUDES)

clean:
	rm -f $(OBJECTS)

mrproper:
	rm -f $(OBJECTS)
	rm -f $(TARGET)

