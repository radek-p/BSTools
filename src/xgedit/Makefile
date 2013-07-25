
# /////////////////////////////////////////////////////////////////////
# This file is a part of the BSTools package
# written by Przemyslaw Kiciak
# /////////////////////////////////////////////////////////////////////
# (C) Copyright by Przemyslaw Kiciak, 2005, 2012
# this package is distributed under the terms of the
# Lesser GNU Public License, see the file COPYING.LIB
# /////////////////////////////////////////////////////////////////////

HOME = /home/przemek
SRCDIR = .
INCDIR = ../../include
LIBDIR = ../../lib

include ../../options.mak

INCLF = $(INCDIR)/xgedit.h $(INCDIR)/xgergb.h $(SRCDIR)/xgeprivate.h
INCGLF = $(INCDIR)/xgedit.h $(INCDIR)/xgledit.h $(INCDIR)/xgergb.h \
  $(SRCDIR)/xgeprivate.h

OBJECTS = xgevars.o xgeinit.o xgekernel.o xgedmsg.o xgergb.o xgecnames.o \
  xgepixel.o xgefocus.o \
  xgedrect.o xgerepos.o xgeempty.o xgemenu.o xgebutton.o xgeswitch.o \
  xgeslidebarf.o xgeslidebard.o xgevslidebarf.o xgevslidebard.o \
  xgeslidebar2f.o xgeslidebar2d.o xgevslidebar2f.o xgevslidebar2d.o \
  xgeslidebarrgbf.o xgeslidebarrgbd.o xgedialf.o xgediald.o \
  xgetextout.o xgefour.o xgestringed.o xgepopup.o xgeintwid.o \
  xgecursor.o xgelistbox.o xgedirlist.o xge2dwinf.o xge2dwind.o xge3dwinf.o \
  xge3dwind.o xgeknotwinf.o xgeknotwind.o xget2knotwinf.o xget2knotwind.o \
  xgescrollwidg.o xgegeomwinback.o xgeoutpix.o xgedrawbezf.o xgedrawbezd.o \
  xgergbsamplef.o xgergbsampled.o xgeshrect.o \
  xgequatrotball.o xgequatrotballf.o xgequatrotballd.o xgeflush.o \
  xgeipcchild.o xgeipcparent.o xgeipcvars.o \
  xglevars.o xgleinit.o xglergb.o xgleidentmap.o xglecameraf.o xglecamerad.o \
  xgleacameraf.o xgleacamerad.o \
  xglegeomwinback.o xglestring.o xgledraw.o xgle3dwinf.o xgle3dwind.o \
  xgle2dwinf.o xgle2dwind.o xglet2knotwinf.o xglet2knotwind.o \
  xgletrans3d.o xgletrans3f.o

TARGET = $(LIBDIR)/libxgedit.a

default: $(TARGET)

$(TARGET): $(OBJECTS)
	$(AR) $(ARFLAGS) $(TARGET) $(OBJECTS)

.c.o:
	$(CC) $(CFLAGS) -I$(INCDIR) -c $< -o $*.o

xgevars.o: $(SRCDIR)/xgevars.c $(INCLF)
xgeinit.o: $(SRCDIR)/xgeinit.c $(INCLF)
xgekernel.o: $(SRCDIR)/xgekernel.c $(INCLF)
xgedmsg.o: $(SRCDIR)/xgedmsg.c $(INCLF)
xgergb.o: $(SRCDIR)/xgergb.c $(SRCDIR)/xge_rgb.c $(INCLF)
xgecnames.o: $(SRCDIR)/xgecnames.c $(INCLF)
xgepixel.o: $(SRCDIR)/xgepixel.c $(INCLF)
xgefocus.o: $(SRCDIR)/xgefocus.c $(INCLF)
xgedrect.o: $(SRCDIR)/xgedrect.c $(INCLF)
xgerepos.o: $(SRCDIR)/xgerepos.c $(INCLF)
xgeempty.o: $(SRCDIR)/xgeempty.c $(INCLF)
xgemenu.o: $(SRCDIR)/xgemenu.c $(INCLF)
xgebutton.o: $(SRCDIR)/xgebutton.c $(INCLF)
xgeswitch.o: $(SRCDIR)/xgeswitch.c $(INCLF)
xgeslidebarf.o: $(SRCDIR)/xgeslidebarf.c $(INCLF)
xgeslidebard.o: $(SRCDIR)/xgeslidebard.c $(INCLF)
xgeslidebarrgbf.o: $(SRCDIR)/xgeslidebarrgbf.c $(INCLF)
xgeslidebarrgbd.o: $(SRCDIR)/xgeslidebarrgbd.c $(INCLF)
xgevslidebarf.o: $(SRCDIR)/xgevslidebarf.c $(INCLF)
xgevslidebard.o: $(SRCDIR)/xgevslidebard.c $(INCLF)
xgeslidebar2f.o: $(SRCDIR)/xgeslidebar2f.c $(INCLF)
xgeslidebar2d.o: $(SRCDIR)/xgeslidebar2d.c $(INCLF)
xgevslidebar2f.o: $(SRCDIR)/xgevslidebar2f.c $(INCLF)
xgevslidebar2d.o: $(SRCDIR)/xgevslidebar2d.c $(INCLF)
xgedialf.o: $(SRCDIR)/xgedialf.c $(INCLF)
xgediald.o: $(SRCDIR)/xgediald.c $(INCLF)
xgetextout.o: $(SRCDIR)/xgetextout.c $(INCLF)
xgefour.o: $(SRCDIR)/xgefour.c $(INCLF)
xgestringed.o: $(SRCDIR)/xgestringed.c $(INCLF)
xgepopup.o: $(SRCDIR)/xgepopup.c $(INCLF)
xgeintwid.o: $(SRCDIR)/xgeintwid.c $(INCLF)
xgelistbox.o: $(SRCDIR)/xgelistbox.c $(INCLF)
xgedirlist.o: $(SRCDIR)/xgedirlist.c $(INCLF)
xge2dwinf.o: $(SRCDIR)/xge2dwinf.c $(INCLF)
xge2dwind.o: $(SRCDIR)/xge2dwind.c $(INCLF)
xge3dwinf.o: $(SRCDIR)/xge3dwinf.c $(INCLF)
xge3dwind.o: $(SRCDIR)/xge3dwind.c $(INCLF)
xgeknotwinf.o: $(SRCDIR)/xgeknotwinf.c $(INCLF)
xgeknotwind.o: $(SRCDIR)/xgeknotwind.c $(INCLF)
xget2knotwinf.o: $(SRCDIR)/xget2knotwinf.c $(INCLF)
xget2knotwind.o: $(SRCDIR)/xget2knotwind.c $(INCLF)
xgegeomwinback.o: $(SRCDIR)/xgegeomwinback.c $(INCLF)
xgeoutpix.o: $(SRCDIR)/xgeoutpix.c $(INCLF)
xgedrawbezf.o: $(SRCDIR)/xgedrawbezf.c $(INCLF)
xgedrawbezd.o: $(SRCDIR)/xgedrawbezd.c $(INCLF)
xgescrollwidg.o: $(SRCDIR)/xgescrollwidg.c $(INCLF)
xgergbsamplef.o: $(SRCDIR)/xgergbsamplef.c $(INCLF)
xgergbsampled.o: $(SRCDIR)/xgergbsampled.c $(INCLF)
xgeshrect.o: $(SRCDIR)/xgeshrect.c $(INCLF)
xgequatrotball.o: $(SRCDIR)/xgequatrotball.c $(INCLF)
xgequatrotballf.o: $(SRCDIR)/xgequatrotballf.c $(INCLF)
xgequatrotballd.o: $(SRCDIR)/xgequatrotballd.c $(INCLF)
xgeflush.o: $(SRCDIR)/xgeflush.c $(INCLF)

xgeipcchild.o: $(SRCDIR)/xgeipcchild.c $(INCLF) $(INCDIR)/xgeipc.h
xgeipcparent.o: $(SRCDIR)/xgeipcparent.c $(INCLF) $(INCDIR)/xgeipc.h
xgeipcvars.o: $(SRCDIR)/xgeipcvars.c $(INCLF) $(INCDIR)/xgeipc.h

xglevars.o: $(SRCDIR)/xglevars.c $(INCGLF)
xgleinit.o: $(SRCDIR)/xgleinit.c $(INCGLF)
xglergb.o: $(SRCDIR)/xglergb.c $(SRCDIR)/xge_rgb.c $(INCGLF)
xgleidentmap.o: $(SRCDIR)/xgleidentmap.c $(INCGLF)
xglegeomwinback.o: $(SRCDIR)/xglegeomwinback.c $(INCGLF)
xglestring.o: $(SRCDIR)/xglestring.c $(INCGLF)
xgledraw.o: $(SRCDIR)/xgledraw.c $(INCGLF)
xgle3fwind.o: $(SRCDIR)/xgle3fwind.c $(INCGLF)
xgle3dwind.o: $(SRCDIR)/xgle3dwind.c $(INCGLF)
xgle2dwinf.o: $(SRCDIR)/xgle2dwinf.c $(INCGLF)
xgle2dwind.o: $(SRCDIR)/xgle2dwind.c $(INCGLF)
xglecameraf.o: $(SRCDIR)/xglecameraf.c $(INCGLF)
xglecamerad.o: $(SRCDIR)/xglecamerad.c $(INCGLF)
xglet2knotwinf.o: $(SRCDIR)/xglet2knotwinf.c $(INCGLF)
xglet2knotwind.o: $(SRCDIR)/xglet2knotwind.c $(INCGLF)
xgletrans3f.o: $(SRCDIR)/xgletrans3f.c $(INCGLF)
xgletrans3d.o: $(SRCDIR)/xgletrans3d.c $(INCGLF)

clean:
	rm -f $(OBJECTS)

mrproper:
	rm -f $(OBJECTS)
	rm -f $(TARGET)

