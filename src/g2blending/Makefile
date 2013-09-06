
# /////////////////////////////////////////////////////////////////////
# This file is a part of the BSTools package
# written by Przemyslaw Kiciak
# /////////////////////////////////////////////////////////////////////
# (C) Copyright by Przemyslaw Kiciak, 2009, 2013
# this package is distributed under the terms of the
# Lesser GNU Public License, see the file COPYING.LIB
# /////////////////////////////////////////////////////////////////////

HOME = /home/przemek
SRCDIR = .
INCDIR = ../../include
LIBDIR = ../../lib

include ../../options.mak

OBJECTS = $(SRCDIR)/g2blsize.o $(SRCDIR)/g2blhprofile.o \
  $(SRCDIR)/g2blhcprofile.o $(SRCDIR)/g2mblsize.o $(SRCDIR)/g2mblsdist.o

OBJECTSF = $(SRCDIR)/g2bltriharmf.o $(SRCDIR)/g2blctriharmf.o \
  $(SRCDIR)/g2bluconstrf.o

OBJECTSD = $(SRCDIR)/g2bltriharmd.o $(SRCDIR)/g2blctriharmd.o \
  $(SRCDIR)/g2bluconstrd.o $(SRCDIR)/g2blcuconstrd.o $(SRCDIR)/g2blutabfd.o \
  $(SRCDIR)/g2blufgh01d.o $(SRCDIR)/g2blufgh02d.o $(SRCDIR)/g2blufgh03d.o \
  $(SRCDIR)/g2blopt01d.o $(SRCDIR)/g2blopt02d.o $(SRCDIR)/g2blopt03d.o \
  $(SRCDIR)/g2blopt04d.o $(SRCDIR)/g2blopt05d.o $(SRCDIR)/g2blopt06d.o \
  $(SRCDIR)/g2blopt07d.o $(SRCDIR)/g2blopt08d.o $(SRCDIR)/g2bloptuconstrd.o \
  $(SRCDIR)/g2bloptuclconstrd.o $(SRCDIR)/g2bltsqfd.o $(SRCDIR)/g2mblholedomd.o \
  $(SRCDIR)/g2mblufgh01d.o $(SRCDIR)/g2mblufgh02d.o $(SRCDIR)/g2mblufgh03d.o \
  $(SRCDIR)/g2mblufgh04d.o $(SRCDIR)/g2mblufgh05d.o \
  $(SRCDIR)/g2mblopt01d.o $(SRCDIR)/g2mblopt02d.o $(SRCDIR)/g2mblopt03d.o \
  $(SRCDIR)/g2mblopt04d.o $(SRCDIR)/g2mblopt05d.o $(SRCDIR)/g2mblopt06d.o \
  $(SRCDIR)/g2mblopt07d.o $(SRCDIR)/g2mblopt08d.o $(SRCDIR)/g2mblopt09d.o \
  $(SRCDIR)/g2mblopt10d.o $(SRCDIR)/g2mblcpopt01d.o \
  $(SRCDIR)/g2mblnvcp.o \
  $(SRCDIR)/g2mbltestfghd.o \
  $(SRCDIR)/g2mblmlopt01d.o $(SRCDIR)/g2mblmlopt02d.o $(SRCDIR)/g2mblmlopt03d.o \
  $(SRCDIR)/g2mblmlopt04d.o $(SRCDIR)/g2mblmlopt05d.o $(SRCDIR)/g2mblmlopt06d.o \
  $(SRCDIR)/g2mblmlopt07d.o $(SRCDIR)/g2mblmlopt08d.o $(SRCDIR)/g2mblmlopt09d.o \
  $(SRCDIR)/g2mblmlopt10d.o $(SRCDIR)/g2mblmlopt11d.o $(SRCDIR)/g2mblmlopt12d.o \
  $(SRCDIR)/g2mblmlopt13d.o $(SRCDIR)/g2mblmlopt14d.o $(SRCDIR)/g2mblmlopt15d.o \
  $(SRCDIR)/g2mblmlopt16d.o $(SRCDIR)/g2mblmlopt17d.o $(SRCDIR)/g2mblmlopt18d.o \
  $(SRCDIR)/g2mblmlopt19d.o \
  $(SRCDIR)/g2mblmlauxd.o $(SRCDIR)/g2mblmlfghd.o \
  $(SRCDIR)/g2mblmlnlevels.o $(SRCDIR)/g2mblmlsnlevels.o \
  $(SRCDIR)/g2mblmlcpnlevels.o $(SRCDIR)/g2mblmlcpsnlevels.o \
  $(SRCDIR)/g2mblsfgh01d.o $(SRCDIR)/g2mblsfgh02d.o \
  $(SRCDIR)/g2mblmlsopt01d.o $(SRCDIR)/g2mblmlsopt02d.o $(SRCDIR)/g2mblmlsopt03d.o \
  $(SRCDIR)/g2mblmlsopt04d.o $(SRCDIR)/g2mblmlsopt05d.o $(SRCDIR)/g2mblmlsopt06d.o \
  $(SRCDIR)/g2mblmlsopt07d.o  $(SRCDIR)/g2mblmlsopt08d.o \
  $(SRCDIR)/g2mblmlcpopt01d.o $(SRCDIR)/g2mblmlcpopt02d.o \
  $(SRCDIR)/g2mblmlcpopt03d.o $(SRCDIR)/g2mblmlcpopt04d.o \
  $(SRCDIR)/g2mblmlcpsopt01d.o $(SRCDIR)/g2mblmlcpsopt02d.o \
  $(SRCDIR)/g2mblmltime.o

TARGET = $(LIBDIR)/libg2blending.a

default: $(TARGET)

$(TARGET): $(OBJECTS) $(OBJECTSF) $(OBJECTSD)
	$(AR) $(ARFLAGS) $(TARGET) $(OBJECTS) $(OBJECTSF) $(OBJECTSD)

.c.o:
	$(CC) $(CFLAGS) -I$(INCDIR) -c $< -o $*.o

$(SRCDIR)/g2bltriharmf.o: $(SRCDIR)/g2bltriharmf.c $(INCDIR)/g2blendingf.h \
  $(INCDIR)/msgpool.h

$(SRCDIR)/g2blsize.o: $(SRCDIR)/g2blsize.c \
  $(INCDIR)/g2blendingf.h $(INCDIR)/g2blendingd.h $(SRCDIR)/g2blprivated.h \
  $(INCDIR)/msgpool.h

$(SRCDIR)/g2blhprofile.o: $(SRCDIR)/g2blhprofile.c \
  $(INCDIR)/g2blendingf.h $(INCDIR)/g2blendingd.h $(SRCDIR)/g2blprivated.h \
  $(INCDIR)/msgpool.h

$(SRCDIR)/g2bluconstrf.o: $(SRCDIR)/g2bluconstrf.c $(INCDIR)/g2blendingf.h \
  $(INCDIR)/msgpool.h

$(SRCDIR)/g2bloptuconstrf.o: $(SRCDIR)/g2bloptuconstrf.c $(INCDIR)/g2blendingf.h \
  $(SRCDIR)/g2blprivatef.h $(INCDIR)/msgpool.h

$(SRCDIR)/g2bltriharmd.o: $(SRCDIR)/g2bltriharmd.c $(INCDIR)/g2blendingd.h \
  $(SRCDIR)/g2blprivated.h $(INCDIR)/msgpool.h

$(SRCDIR)/g2blutabfd.o: $(SRCDIR)/g2blutabfd.c $(INCDIR)/g2blendingd.h \
  $(SRCDIR)/g2blprivated.h $(INCDIR)/msgpool.h

$(SRCDIR)/g2blufgh01d.o: $(SRCDIR)/g2blufgh01d.c $(INCDIR)/g2blendingd.h \
  $(SRCDIR)/g2blprivated.h $(INCDIR)/msgpool.h

$(SRCDIR)/g2blufgh02d.o: $(SRCDIR)/g2blufgh02d.c $(INCDIR)/g2blendingd.h \
  $(SRCDIR)/g2blprivated.h $(INCDIR)/msgpool.h

$(SRCDIR)/g2blufgh03d.o: $(SRCDIR)/g2blufgh03d.c $(INCDIR)/g2blendingd.h \
  $(SRCDIR)/g2blprivated.h $(INCDIR)/msgpool.h

$(SRCDIR)/g2blopt01d.o: $(SRCDIR)/g2blopt01d.c $(INCDIR)/g2blendingd.h \
  $(SRCDIR)/g2blprivated.h $(INCDIR)/msgpool.h

$(SRCDIR)/g2blopt02d.o: $(SRCDIR)/g2blopt02d.c $(INCDIR)/g2blendingd.h \
  $(SRCDIR)/g2blprivated.h $(INCDIR)/msgpool.h

$(SRCDIR)/g2blopt03d.o: $(SRCDIR)/g2blopt03d.c $(INCDIR)/g2blendingd.h \
  $(SRCDIR)/g2blprivated.h $(INCDIR)/msgpool.h

$(SRCDIR)/g2blopt04d.o: $(SRCDIR)/g2blopt04d.c $(INCDIR)/g2blendingd.h \
  $(SRCDIR)/g2blprivated.h $(INCDIR)/msgpool.h

$(SRCDIR)/g2blopt05d.o: $(SRCDIR)/g2blopt05d.c $(INCDIR)/g2blendingd.h \
  $(SRCDIR)/g2blprivated.h $(INCDIR)/msgpool.h

$(SRCDIR)/g2blopt06d.o: $(SRCDIR)/g2blopt06d.c $(INCDIR)/g2blendingd.h \
  $(SRCDIR)/g2blprivated.h $(INCDIR)/msgpool.h

$(SRCDIR)/g2blopt07d.o: $(SRCDIR)/g2blopt07d.c $(INCDIR)/g2blendingd.h \
  $(SRCDIR)/g2blprivated.h $(INCDIR)/msgpool.h

$(SRCDIR)/g2blopt08d.o: $(SRCDIR)/g2blopt08d.c $(INCDIR)/g2blendingd.h \
  $(SRCDIR)/g2blprivated.h $(INCDIR)/msgpool.h

$(SRCDIR)/g2bluconstrd.o: $(SRCDIR)/g2bluconstrd.c $(INCDIR)/g2blendingd.h \
  $(SRCDIR)/g2blprivated.h $(INCDIR)/msgpool.h

$(SRCDIR)/g2blcuconstrd.o: $(SRCDIR)/g2blcuconstrd.c $(INCDIR)/g2blendingd.h \
  $(SRCDIR)/g2blprivated.h $(INCDIR)/msgpool.h

$(SRCDIR)/g2bloptuconstrd.o: $(SRCDIR)/g2bloptuconstrd.c $(INCDIR)/g2blendingd.h \
  $(SRCDIR)/g2blprivated.h $(INCDIR)/msgpool.h

$(SRCDIR)/g2bltsqfd.o: $(SRCDIR)/g2bltsqfd.c $(INCDIR)/g2blendingd.h \
  $(SRCDIR)/g2blprivated.h $(INCDIR)/msgpool.h


$(SRCDIR)/g2mblsize.o: $(SRCDIR)/g2mblsize.c $(INCDIR)/g2mblendingd.h \
  $(SRCDIR)/g2mblprivated.h $(SRCDIR)/g2blprivated.h $(INCDIR)/msgpool.h

$(SRCDIR)/g2mblholedomd.o: $(SRCDIR)/g2mblholedomd.c $(INCDIR)/g2mblendingd.h \
  $(SRCDIR)/g2mblprivated.h $(SRCDIR)/g2blprivated.h $(INCDIR)/msgpool.h

$(SRCDIR)/g2mblufgh01d.o: $(SRCDIR)/g2mblufgh01d.c $(INCDIR)/g2mblendingd.h \
  $(SRCDIR)/g2mblprivated.h $(SRCDIR)/g2blprivated.h $(INCDIR)/msgpool.h

$(SRCDIR)/g2mblufgh02d.o: $(SRCDIR)/g2mblufgh02d.c $(INCDIR)/g2mblendingd.h \
  $(SRCDIR)/g2mblprivated.h $(SRCDIR)/g2blprivated.h $(INCDIR)/msgpool.h

$(SRCDIR)/g2mblufgh03d.o: $(SRCDIR)/g2mblufgh03d.c $(INCDIR)/g2mblendingd.h \
  $(SRCDIR)/g2mblprivated.h $(SRCDIR)/g2blprivated.h $(INCDIR)/msgpool.h

$(SRCDIR)/g2mblufgh04d.o: $(SRCDIR)/g2mblufgh04d.c $(INCDIR)/g2mblendingd.h \
  $(SRCDIR)/g2mblprivated.h $(SRCDIR)/g2blprivated.h $(INCDIR)/msgpool.h

$(SRCDIR)/g2mblufgh05d.o: $(SRCDIR)/g2mblufgh05d.c $(INCDIR)/g2mblendingd.h \
  $(SRCDIR)/g2mblprivated.h $(SRCDIR)/g2blprivated.h $(INCDIR)/msgpool.h

$(SRCDIR)/g2mblopt01d.o: $(SRCDIR)/g2mblopt01d.c $(INCDIR)/g2mblendingd.h \
  $(SRCDIR)/g2mblprivated.h $(SRCDIR)/g2blprivated.h $(INCDIR)/msgpool.h

$(SRCDIR)/g2mblopt02d.o: $(SRCDIR)/g2mblopt02d.c $(INCDIR)/g2mblendingd.h \
  $(SRCDIR)/g2mblprivated.h $(SRCDIR)/g2blprivated.h $(INCDIR)/msgpool.h

$(SRCDIR)/g2mblopt03d.o: $(SRCDIR)/g2mblopt03d.c $(INCDIR)/g2mblendingd.h \
  $(SRCDIR)/g2mblprivated.h $(SRCDIR)/g2blprivated.h $(INCDIR)/msgpool.h

$(SRCDIR)/g2mblopt04d.o: $(SRCDIR)/g2mblopt04d.c $(INCDIR)/g2mblendingd.h \
  $(SRCDIR)/g2mblprivated.h $(SRCDIR)/g2blprivated.h $(INCDIR)/msgpool.h

$(SRCDIR)/g2mblopt05d.o: $(SRCDIR)/g2mblopt05d.c $(INCDIR)/g2mblendingd.h \
  $(SRCDIR)/g2mblprivated.h $(SRCDIR)/g2blprivated.h $(INCDIR)/msgpool.h

$(SRCDIR)/g2mblopt06d.o: $(SRCDIR)/g2mblopt06d.c $(INCDIR)/g2mblendingd.h \
  $(SRCDIR)/g2mblprivated.h $(SRCDIR)/g2blprivated.h $(INCDIR)/msgpool.h

$(SRCDIR)/g2mblopt07d.o: $(SRCDIR)/g2mblopt07d.c $(INCDIR)/g2mblendingd.h \
  $(SRCDIR)/g2mblprivated.h $(SRCDIR)/g2blprivated.h $(INCDIR)/msgpool.h

$(SRCDIR)/g2mblopt08d.o: $(SRCDIR)/g2mblopt08d.c $(INCDIR)/g2mblendingd.h \
  $(SRCDIR)/g2mblprivated.h $(SRCDIR)/g2blprivated.h $(INCDIR)/msgpool.h

$(SRCDIR)/g2mblopt09d.o: $(SRCDIR)/g2mblopt09d.c $(INCDIR)/g2mblendingd.h \
  $(SRCDIR)/g2mblprivated.h $(SRCDIR)/g2blprivated.h $(INCDIR)/msgpool.h

$(SRCDIR)/g2mblopt10d.o: $(SRCDIR)/g2mblopt10d.c $(INCDIR)/g2mblendingd.h \
  $(SRCDIR)/g2mblprivated.h $(SRCDIR)/g2blprivated.h $(INCDIR)/msgpool.h

$(SRCDIR)/g2mblcpopt01d.o: $(SRCDIR)/g2mblcpopt01d.c $(INCDIR)/g2mblendingd.h \
  $(SRCDIR)/g2mblprivated.h $(SRCDIR)/g2blprivated.h $(INCDIR)/msgpool.h

$(SRCDIR)/g2mblsdist.o: $(SRCDIR)/g2mblsdist.c $(INCDIR)/g2mblendingd.h \
  $(SRCDIR)/g2mblprivated.h $(SRCDIR)/g2blprivated.h $(INCDIR)/msgpool.h 

$(SRCDIR)/g2mbltestfghd.o: $(SRCDIR)/g2mbltestfghd.c $(INCDIR)/g2mblendingd.h \
  $(SRCDIR)/g2mblprivated.h $(SRCDIR)/g2blprivated.h $(INCDIR)/msgpool.h

$(SRCDIR)/g2mblnvcp.o: $(SRCDIR)/g2mblnvcp.c $(INCDIR)/g2mblendingd.h \
  $(SRCDIR)/g2mblprivated.h $(SRCDIR)/g2blprivated.h $(INCDIR)/msgpool.h


$(SRCDIR)/g2mblmlopt01d.o: $(SRCDIR)/g2mblmlopt01d.c $(INCDIR)/g2mblendingd.h \
  $(SRCDIR)/g2mblprivated.h $(SRCDIR)/g2mblmlprivated.h \
  $(SRCDIR)/g2blprivated.h $(INCDIR)/msgpool.h

$(SRCDIR)/g2mblmlopt02d.o: $(SRCDIR)/g2mblmlopt02d.c $(INCDIR)/g2mblendingd.h \
  $(SRCDIR)/g2mblprivated.h $(SRCDIR)/g2mblmlprivated.h \
  $(SRCDIR)/g2blprivated.h $(INCDIR)/msgpool.h

$(SRCDIR)/g2mblmlopt03d.o: $(SRCDIR)/g2mblmlopt03d.c $(INCDIR)/g2mblendingd.h \
  $(SRCDIR)/g2mblprivated.h $(SRCDIR)/g2mblmlprivated.h \
  $(SRCDIR)/g2blprivated.h $(INCDIR)/msgpool.h

$(SRCDIR)/g2mblmlopt04d.o: $(SRCDIR)/g2mblmlopt04d.c $(INCDIR)/g2mblendingd.h \
  $(SRCDIR)/g2mblprivated.h $(SRCDIR)/g2mblmlprivated.h \
  $(SRCDIR)/g2blprivated.h $(INCDIR)/msgpool.h

$(SRCDIR)/g2mblmlopt05d.o: $(SRCDIR)/g2mblmlopt05d.c $(INCDIR)/g2mblendingd.h \
  $(SRCDIR)/g2mblprivated.h $(SRCDIR)/g2mblmlprivated.h \
  $(SRCDIR)/g2blprivated.h $(INCDIR)/msgpool.h

$(SRCDIR)/g2mblmlopt06d.o: $(SRCDIR)/g2mblmlopt06d.c $(INCDIR)/g2mblendingd.h \
  $(SRCDIR)/g2mblprivated.h $(SRCDIR)/g2mblmlprivated.h \
  $(SRCDIR)/g2blprivated.h $(INCDIR)/msgpool.h

$(SRCDIR)/g2mblmlopt07d.o: $(SRCDIR)/g2mblmlopt07d.c $(INCDIR)/g2mblendingd.h \
  $(SRCDIR)/g2mblprivated.h $(SRCDIR)/g2mblmlprivated.h \
  $(SRCDIR)/g2blprivated.h $(INCDIR)/msgpool.h

$(SRCDIR)/g2mblmlopt08d.o: $(SRCDIR)/g2mblmlopt08d.c $(INCDIR)/g2mblendingd.h \
  $(SRCDIR)/g2mblprivated.h $(SRCDIR)/g2mblmlprivated.h \
  $(SRCDIR)/g2blprivated.h $(INCDIR)/msgpool.h

$(SRCDIR)/g2mblmlopt09d.o: $(SRCDIR)/g2mblmlopt09d.c $(INCDIR)/g2mblendingd.h \
  $(SRCDIR)/g2mblprivated.h $(SRCDIR)/g2mblmlprivated.h \
  $(SRCDIR)/g2blprivated.h $(INCDIR)/msgpool.h

$(SRCDIR)/g2mblmlopt10d.o: $(SRCDIR)/g2mblmlopt10d.c $(INCDIR)/g2mblendingd.h \
  $(SRCDIR)/g2mblprivated.h $(SRCDIR)/g2mblmlprivated.h \
  $(SRCDIR)/g2blprivated.h $(INCDIR)/msgpool.h

$(SRCDIR)/g2mblmlopt11d.o: $(SRCDIR)/g2mblmlopt11d.c $(INCDIR)/g2mblendingd.h \
  $(SRCDIR)/g2mblprivated.h $(SRCDIR)/g2mblmlprivated.h \
  $(SRCDIR)/g2blprivated.h $(INCDIR)/msgpool.h

$(SRCDIR)/g2mblmlopt12d.o: $(SRCDIR)/g2mblmlopt12d.c $(INCDIR)/g2mblendingd.h \
  $(SRCDIR)/g2mblprivated.h $(SRCDIR)/g2mblmlprivated.h \
  $(SRCDIR)/g2blprivated.h $(INCDIR)/msgpool.h

$(SRCDIR)/g2mblmlopt13d.o: $(SRCDIR)/g2mblmlopt13d.c $(INCDIR)/g2mblendingd.h \
  $(SRCDIR)/g2mblprivated.h $(SRCDIR)/g2mblmlprivated.h \
  $(SRCDIR)/g2blprivated.h $(INCDIR)/msgpool.h

$(SRCDIR)/g2mblmlopt14d.o: $(SRCDIR)/g2mblmlopt14d.c $(INCDIR)/g2mblendingd.h \
  $(SRCDIR)/g2mblprivated.h $(SRCDIR)/g2mblmlprivated.h \
  $(SRCDIR)/g2blprivated.h $(INCDIR)/msgpool.h

$(SRCDIR)/g2mblmlopt15d.o: $(SRCDIR)/g2mblmlopt15d.c $(INCDIR)/g2mblendingd.h \
  $(SRCDIR)/g2mblprivated.h $(SRCDIR)/g2mblmlprivated.h \
  $(SRCDIR)/g2blprivated.h $(INCDIR)/msgpool.h

$(SRCDIR)/g2mblmlopt16d.o: $(SRCDIR)/g2mblmlopt16d.c $(INCDIR)/g2mblendingd.h \
  $(SRCDIR)/g2mblprivated.h $(SRCDIR)/g2mblmlprivated.h \
  $(SRCDIR)/g2blprivated.h $(INCDIR)/msgpool.h

$(SRCDIR)/g2mblmlopt17d.o: $(SRCDIR)/g2mblmlopt17d.c $(INCDIR)/g2mblendingd.h \
  $(SRCDIR)/g2mblprivated.h $(SRCDIR)/g2mblmlprivated.h \
  $(SRCDIR)/g2blprivated.h $(INCDIR)/msgpool.h

$(SRCDIR)/g2mblmlopt18d.o: $(SRCDIR)/g2mblmlopt18d.c $(INCDIR)/g2mblendingd.h \
  $(SRCDIR)/g2mblprivated.h $(SRCDIR)/g2mblmlprivated.h \
  $(SRCDIR)/g2blprivated.h $(INCDIR)/msgpool.h

$(SRCDIR)/g2mblmlfghd.o: $(SRCDIR)/g2mblmlfghd.c \
  $(INCDIR)/g2mblendingd.h \
  $(SRCDIR)/g2mblprivated.h $(SRCDIR)/g2mblmlprivated.h \
  $(SRCDIR)/g2blprivated.h $(INCDIR)/msgpool.h

$(SRCDIR)/g2mblmlauxd.o: $(SRCDIR)/g2mblmlauxd.c $(INCDIR)/g2mblendingd.h \
  $(SRCDIR)/g2mblprivated.h $(SRCDIR)/g2mblmlprivated.h \
  $(SRCDIR)/g2blprivated.h $(INCDIR)/msgpool.h

$(SRCDIR)/g2mblmlnlevels.o: $(SRCDIR)/g2mblmlnlevels.c $(INCDIR)/g2mblendingd.h \
  $(SRCDIR)/g2mblprivated.h $(SRCDIR)/g2mblmlprivated.h \
  $(SRCDIR)/g2blprivated.h $(INCDIR)/msgpool.h


$(SRCDIR)/g2mblsfgh01d.o: $(SRCDIR)/g2mblsfgh01d.c $(INCDIR)/g2mblendingd.h \
  $(SRCDIR)/g2mblprivated.h $(SRCDIR)/g2mblmlprivated.h \
  $(SRCDIR)/g2blprivated.h $(INCDIR)/msgpool.h

$(SRCDIR)/g2mblsfgh02d.o: $(SRCDIR)/g2mblsfgh02d.c $(INCDIR)/g2mblendingd.h \
  $(SRCDIR)/g2mblprivated.h $(SRCDIR)/g2mblmlprivated.h \
  $(SRCDIR)/g2blprivated.h $(INCDIR)/msgpool.h

$(SRCDIR)/g2mblmlsnlevels.o: $(SRCDIR)/g2mblmlsnlevels.c $(INCDIR)/g2mblendingd.h \
  $(SRCDIR)/g2mblprivated.h $(SRCDIR)/g2mblmlprivated.h \
  $(SRCDIR)/g2blprivated.h $(INCDIR)/msgpool.h

$(SRCDIR)/g2mblmlsopt01d.o: $(SRCDIR)/g2mblmlsopt01d.c $(INCDIR)/g2mblendingd.h \
  $(SRCDIR)/g2mblprivated.h $(SRCDIR)/g2mblmlprivated.h \
  $(SRCDIR)/g2blprivated.h $(INCDIR)/msgpool.h

$(SRCDIR)/g2mblmlsopt02d.o: $(SRCDIR)/g2mblmlsopt02d.c $(INCDIR)/g2mblendingd.h \
  $(SRCDIR)/g2mblprivated.h $(SRCDIR)/g2mblmlprivated.h \
  $(SRCDIR)/g2blprivated.h $(INCDIR)/msgpool.h

$(SRCDIR)/g2mblmlsopt03d.o: $(SRCDIR)/g2mblmlsopt03d.c $(INCDIR)/g2mblendingd.h \
  $(SRCDIR)/g2mblprivated.h $(SRCDIR)/g2mblmlprivated.h \
  $(SRCDIR)/g2blprivated.h $(INCDIR)/msgpool.h

$(SRCDIR)/g2mblmlsopt04d.o: $(SRCDIR)/g2mblmlsopt04d.c $(INCDIR)/g2mblendingd.h \
  $(SRCDIR)/g2mblprivated.h $(SRCDIR)/g2mblmlprivated.h \
  $(SRCDIR)/g2blprivated.h $(INCDIR)/msgpool.h

$(SRCDIR)/g2mblmlsopt05d.o: $(SRCDIR)/g2mblmlsopt05d.c $(INCDIR)/g2mblendingd.h \
  $(SRCDIR)/g2mblprivated.h $(SRCDIR)/g2mblmlprivated.h \
  $(SRCDIR)/g2blprivated.h $(INCDIR)/msgpool.h

$(SRCDIR)/g2mblmlsopt06d.o: $(SRCDIR)/g2mblmlsopt06d.c $(INCDIR)/g2mblendingd.h \
  $(SRCDIR)/g2mblprivated.h $(SRCDIR)/g2mblmlprivated.h \
  $(SRCDIR)/g2blprivated.h $(INCDIR)/msgpool.h

$(SRCDIR)/g2mblmlsopt07d.o: $(SRCDIR)/g2mblmlsopt07d.c $(INCDIR)/g2mblendingd.h \
  $(SRCDIR)/g2mblprivated.h $(SRCDIR)/g2mblmlprivated.h \
  $(SRCDIR)/g2blprivated.h $(INCDIR)/msgpool.h

$(SRCDIR)/g2mblmlsopt08d.o: $(SRCDIR)/g2mblmlsopt08d.c $(INCDIR)/g2mblendingd.h \
  $(SRCDIR)/g2mblprivated.h $(SRCDIR)/g2mblmlprivated.h \
  $(SRCDIR)/g2blprivated.h $(INCDIR)/msgpool.h


$(SRCDIR)/g2mblmlcpopt01d.o: $(SRCDIR)/g2mblmlcpopt01d.c $(INCDIR)/g2mblendingd.h \
  $(SRCDIR)/g2mblprivated.h $(SRCDIR)/g2mblmlprivated.h \
  $(SRCDIR)/g2blprivated.h $(INCDIR)/msgpool.h

$(SRCDIR)/g2mblmlcpopt02d.o: $(SRCDIR)/g2mblmlcpopt02d.c $(INCDIR)/g2mblendingd.h \
  $(SRCDIR)/g2mblprivated.h $(SRCDIR)/g2mblmlprivated.h \
  $(SRCDIR)/g2blprivated.h $(INCDIR)/msgpool.h

$(SRCDIR)/g2mblmlcpopt03d.o: $(SRCDIR)/g2mblmlcpopt03d.c $(INCDIR)/g2mblendingd.h \
  $(SRCDIR)/g2mblprivated.h $(SRCDIR)/g2mblmlprivated.h \
  $(SRCDIR)/g2blprivated.h $(INCDIR)/msgpool.h

$(SRCDIR)/g2mblmlcpopt04d.o: $(SRCDIR)/g2mblmlcpopt04d.c $(INCDIR)/g2mblendingd.h \
  $(SRCDIR)/g2mblprivated.h $(SRCDIR)/g2mblmlprivated.h \
  $(SRCDIR)/g2blprivated.h $(INCDIR)/msgpool.h

$(SRCDIR)/g2mblmlcpsopt01d.o: $(SRCDIR)/g2mblmlcpsopt01d.c $(INCDIR)/g2mblendingd.h \
  $(SRCDIR)/g2mblprivated.h $(SRCDIR)/g2mblmlprivated.h \
  $(SRCDIR)/g2blprivated.h $(INCDIR)/msgpool.h

$(SRCDIR)/g2mblmlcpsopt02d.o: $(SRCDIR)/g2mblmlcpsopt02d.c $(INCDIR)/g2mblendingd.h \
  $(SRCDIR)/g2mblprivated.h $(SRCDIR)/g2mblmlprivated.h \
  $(SRCDIR)/g2blprivated.h $(INCDIR)/msgpool.h

$(SRCDIR)/g2mblmlcpnlevels.o: $(SRCDIR)/g2mblmlcpnlevels.c $(INCDIR)/g2mblendingd.h \
  $(SRCDIR)/g2mblprivated.h $(SRCDIR)/g2mblmlprivated.h \
  $(SRCDIR)/g2blprivated.h $(INCDIR)/msgpool.h

$(SRCDIR)/g2mblmlcpsnlevels.o: $(SRCDIR)/g2mblmlcpsnlevels.c $(INCDIR)/g2mblendingd.h \
  $(SRCDIR)/g2mblprivated.h $(SRCDIR)/g2mblmlprivated.h \
  $(SRCDIR)/g2blprivated.h $(INCDIR)/msgpool.h


$(SRCDIR)/g2mblmltime.o: $(SRCDIR)/g2mblmltime.c $(INCDIR)/g2mblendingd.h \
  $(SRCDIR)/g2mblprivated.h $(SRCDIR)/g2mblmlprivated.h \
  $(SRCDIR)/g2blprivated.h


clean:
	rm -f $(OBJECTS) $(OBJECTSF) $(OBJECTSD)

mrproper:
	rm -f $(OBJECTS) $(OBJECTSF) $(OBJECTSD)
	rm -f $(TARGET)

