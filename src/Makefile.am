## Makefile.am -- edit me!
# Process this file with 'automake' to produce Makefile.in
# Then ./configure to produce Makefile

#main include directory
SOPS_INCDIR	= $(top_srcdir)/include

#arch-dependent include directory (i.e. Daixtrose, since it's prebuilt at config time)
DEPS_INCDIR	= $(top_builddir)/include

#preprocessor flags go here
AM_CPPFLAGS = -O4 -Wall -fexpensive-optimizations -funroll-loops \
	-DBOOST_UBLAS_SHALLOW_ARRAY_ADAPTOR -DNDEBUG -DUINTAH_FIELD_TYPES \
	-I$(SOPS_INCDIR) -I$(DEPS_INCDIR)/daixtrose

lib_LIBRARIES	= libspatialops.a

#Optimize flags go here.
libspatialops_a_CXXFLAGS	= -O3 -Wall -fexpensive-optimizations -funroll-loops

#Compile these files and pack them into a static library.
#Sources not listed here might not get packed into the source distribution.
libspatialops_a_SOURCES	=  \
	FVOneDimensional.cpp \
	LagrangePoly.cpp \
	LinAlgTrilinos.cpp \
	LinAlgUBlas.cpp \
	LinearSystem.cpp

#libspatialops_a_HEADERS = #headers to be installed

#libspatialops_a_noinst_HEADERS = #headers just used for compiling

#where is an extra_INSTALL target to specify libteuchos.a, etc. ?
#How about an extra_UNINSTALL target to clean up subpackages' headers?
