
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

OBJECTS = camera01f.o camera02f.o camera03f.o camera04f.o camera05f.o \
  camera06f.o camera07f.o camera08f.o camera09f.o camera10f.o \
  stereo01f.o stereo02f.o \
  camera01d.o camera02d.o camera03d.o camera04d.o camera05d.o \
  camera06d.o camera07d.o camera08d.o camera09d.o camera10d.o \
  stereo01d.o stereo02d.o

TARGET = $(LIBDIR)/libcamera.a

default: $(TARGET)

$(TARGET): $(OBJECTS)
	$(AR) $(ARFLAGS) $(TARGET) $(OBJECTS)

.c.o:
	$(CC) $(CFLAGS) -I$(INCDIR) -c $< -o $*.o

camera01f.o: camera01f.c cprivate.h msgpool.h \
  $(INCDIR)/camera.h $(INCDIR)/cameraf.h \
  $(INCDIR)/camerad.h
camera02f.o: camera02f.c cprivate.h msgpool.h \
  $(INCDIR)/camera.h $(INCDIR)/cameraf.h \
  $(INCDIR)/camerad.h
camera03f.o: camera03f.c cprivate.h msgpool.h \
  $(INCDIR)/camera.h $(INCDIR)/cameraf.h \
  $(INCDIR)/camerad.h
camera04f.o: camera04f.c cprivate.h msgpool.h \
  $(INCDIR)/camera.h $(INCDIR)/cameraf.h \
  $(INCDIR)/camerad.h
camera05f.o: camera05f.c cprivate.h msgpool.h \
  $(INCDIR)/camera.h $(INCDIR)/cameraf.h \
  $(INCDIR)/camerad.h
camera06f.o: camera06f.c cprivate.h msgpool.h \
  $(INCDIR)/camera.h $(INCDIR)/cameraf.h \
  $(INCDIR)/camerad.h
camera07f.o: camera07f.c cprivate.h msgpool.h \
  $(INCDIR)/camera.h $(INCDIR)/cameraf.h \
  $(INCDIR)/camerad.h
camera08f.o: camera08f.c cprivate.h msgpool.h \
  $(INCDIR)/camera.h $(INCDIR)/cameraf.h \
  $(INCDIR)/camerad.h
camera09f.o: camera09f.c cprivate.h msgpool.h \
  $(INCDIR)/camera.h $(INCDIR)/cameraf.h \
  $(INCDIR)/camerad.h
camera10f.o: camera10f.c cprivate.h msgpool.h \
  $(INCDIR)/camera.h $(INCDIR)/cameraf.h \
  $(INCDIR)/camerad.h

camera01d.o: camera01d.c cprivate.h msgpool.h \
  $(INCDIR)/camera.h $(INCDIR)/cameraf.h \
  $(INCDIR)/camerad.h
camera02d.o: camera02d.c cprivate.h msgpool.h \
  $(INCDIR)/camera.h $(INCDIR)/cameraf.h \
  $(INCDIR)/camerad.h
camera03d.o: camera03d.c cprivate.h msgpool.h \
  $(INCDIR)/camera.h $(INCDIR)/cameraf.h \
  $(INCDIR)/camerad.h
camera04d.o: camera04d.c cprivate.h msgpool.h \
  $(INCDIR)/camera.h $(INCDIR)/cameraf.h \
  $(INCDIR)/camerad.h
camera05d.o: camera05d.c cprivate.h msgpool.h \
  $(INCDIR)/camera.h $(INCDIR)/cameraf.h \
  $(INCDIR)/camerad.h
camera06d.o: camera06d.c cprivate.h msgpool.h \
  $(INCDIR)/camera.h $(INCDIR)/cameraf.h \
  $(INCDIR)/camerad.h
camera07d.o: camera07d.c cprivate.h msgpool.h \
  $(INCDIR)/camera.h $(INCDIR)/cameraf.h \
  $(INCDIR)/camerad.h
camera08d.o: camera08d.c cprivate.h msgpool.h \
  $(INCDIR)/camera.h $(INCDIR)/cameraf.h \
  $(INCDIR)/camerad.h
camera09d.o: camera09d.c cprivate.h msgpool.h \
  $(INCDIR)/camera.h $(INCDIR)/cameraf.h \
  $(INCDIR)/camerad.h
camera10d.o: camera10d.c cprivate.h msgpool.h \
  $(INCDIR)/camera.h $(INCDIR)/cameraf.h \
  $(INCDIR)/camerad.h

stereo01f.o: stereo01f.c cprivate.h msgpool.h \
  $(INCDIR)/stereo.h $(INCDIR)/stereof.h \
  $(INCDIR)/stereod.h $(INCDIR)/camera.h \
  $(INCDIR)/cameraf.h $(INCDIR)/camerad.h \
 
stereo02f.o: stereo02f.c cprivate.h msgpool.h \
  $(INCDIR)/stereo.h $(INCDIR)/stereof.h \
  $(INCDIR)/stereod.h $(INCDIR)/camera.h \
  $(INCDIR)/cameraf.h $(INCDIR)/camerad.h \
 

stereo01d.o: stereo01d.c cprivate.h msgpool.h \
  $(INCDIR)/stereo.h $(INCDIR)/stereof.h \
  $(INCDIR)/stereod.h $(INCDIR)/camera.h \
  $(INCDIR)/cameraf.h $(INCDIR)/camerad.h \
 
stereo02d.o: stereo02d.c cprivate.h msgpool.h \
  $(INCDIR)/stereo.h $(INCDIR)/stereof.h \
  $(INCDIR)/stereod.h $(INCDIR)/camera.h \
  $(INCDIR)/cameraf.h $(INCDIR)/camerad.h \
 

clean:
	rm -f $(OBJECTS)

mrproper:
	rm -f $(OBJECTS)
	rm -f $(TARGET)

