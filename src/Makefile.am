## Makefile.am -- edit me!
# Process this file with 'automake' to produce Makefile.in
# Then ./configure to produce Makefile

AM_CPPFLAGS	= -I$(top_srcdir)/include -I$(top_srcdir)/include/daixtrose -I$(top_srcdir)/include/trilinos -I$(top_srcdir)/include/boost $(SOPS_CPPFLAGS) #preprocessor flags go here
libspatialops_a_CXXFLAGS	= -O3 -Wall -fexpensive-optimizations -funroll-loops
EPETRA_LIBS 	= -lepetra -lepetraext -lblas -llapack
SOPS_LIBDIR	= -L$(top_srcdir)/lib

SUBDIRS	= test
lib_LIBRARIES	= libspatialops.a
libspatialops_a_SOURCES	=  LinAlgTrilinos.cpp FVStaggeredSpatialOps.cpp LinearSystem.cpp
