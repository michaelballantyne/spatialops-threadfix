## Makefile.am -- edit me!
# Process this file with 'automake' to produce Makefile.in
# Then ./configure to produce Makefile

AM_CPPFLAGS	= -I$(top_srcdir)/include -I$(top_srcdir)/include/trilinos -I$(top_srcdir)/include/boost-1_33_1 $(SOPS_CPPFLAGS) #preprocessor flags go here
EPETRA_LIBS 	= -lepetra -lepetraext -lblas -llapack
SOPS_LIBDIR	= -L$(top_srcdir)/lib

SUBDIRS	= test
lib_LIBRARIES	= libspatialops.a
libspatialops_a_SOURCES	=  LinAlgTrilinos.cpp FVStaggeredSpatialOps.cpp LinearSystem.cpp

 #any other object files to pack in the lib archive

#bunch of junk for experimenting; ignore this:
#
#bin_PROGRAMS 	= LinAlgTrilinos FVStaggeredSpatialOps LinearSystem
#
#LinAlgTrilinos_SOURCES 	= LinAlgTrilinos.cpp ../include/LinAlgTrilinos.h
#FVStaggeredSpatialOps_SOURCES = FVStaggeredSpatialOps.cpp ../include/SpatialField.h ../include/SpatialOperator.h
#LinearSystem_SOURCES 		= LinearSystem.cpp ../include/SpatialOperator.h ../include/SpatialField.h ../include/LinearSystem.h 
#
#link with these
#LinAlgTrilinos_LDADD		= -L../lib $(EPETRA_LIBS) $(SOPS_LIBS) 
#FVStaggeredSpatialOps_LDADD	= -L../lib -lepetra -lblas -lmpi $(SOPS_LIBS)
#LinearSystem_LDADD		= -L../lib $(SOPS_LIBS)
#
#
#
#false.o: true.c
#	$(COMPILE)
