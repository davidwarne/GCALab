#!/bin/make

MAKE=make
LIBMESHDIR=libMesh
LIBGCADIR=libGCA
LIBBITMAPDIR=libBitMap
GCALABDIR=GCALab
DOXYGENCONF=GCALab_doxygen.conf

LIBS= $(LIBMESHDIR)/libmesh.so $(LIBGCADIR)/libGCA.so $(LIBBITMAPDIR)/libbitmap.so

all:
	$(MAKE) $(LIBMESHDIR)/libmesh.so
	$(MAKE) $(LIBGCADIR)/libGCA.so
	$(MAKE) $(LIBBITMAPDIR)/libbitmap.so
	$(MAKE) $(GCALABDIR)/GCALab

$(LIBMESHDIR)/libmesh.so: 
	$(MAKE) -C $(LIBMESHDIR)

$(LIBGCADIR)/libGCA.so: $(LIBMESHDIR)/libmesh.so
	$(MAKE) -C $(LIBGCADIR)

$(LIBBITMAPDIR)/libbitmap.so:
	$(MAKE) -C $(LIBBITMAPDIR)

$(GCALABDIR)/GCALab: $(LIBS)
	$(MAKE) -C $(GCALABDIR)

docs:
	doxygen $(DOXYGENCONF)

clean:
	rm */*.o */*.a */*.so
