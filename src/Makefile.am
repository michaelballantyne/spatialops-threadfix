## Makefile.am -- edit me!
# Process this file with 'automake' to produce Makefile.in
# Then ./configure to produce Makefile

#SpatialOps sources are not always near the build location.
#For example, we might build each arch in its own directory.
#We need to refer to headers realtive to the sources' path,
#not the build path. ($srcdir instead of $top_srcdir)
SOPS_INCDIR	= $(srcdir)/../include

+#XXX change project tree structure?
+#(first) move all headers from /include to /src
+#(second) rename /src to /spatialops
+#(third)  change all includes to #include <spatialops/foo.h>
+#(fourth) make buildsys install headers to /usr/include/spatialops
+# Now we can use the same header paths to build the library,
+# as we can to link with it later on.

#Dependency package headers _are_ near our build location
# (since they're installed on a per-arch basis.)  Therefore
# they're relative to the build path, not the source path.
DEPS_INCDIR	= ../include

#preprocessor flags go here
AM_CPPFLAGS	= \
	-I$(SOPS_INCDIR) \
	-I$(DEPS_INCDIR) \
	-I$(DEPS_INCDIR)/daixtrose \
	-I$(DEPS_INCDIR)/trilinos \
	-I$(DEPS_INCDIR)/boost \
	$(SOPS_CPPFLAGS)

EPETRA_LIBS 	= -lepetra -lepetraext -lblas -llapack
SOPS_LIBDIR	= -L$(top_srcdir)/lib

SUBDIRS	= test
lib_LIBRARIES	= libspatialops.a

#Optimize flags go here.
libspatialops_a_CXXFLAGS	= -O3 -Wall -fexpensive-optimizations -funroll-loops

#Compile these files and pack them into a static library.
#Sources not listed here might not get packed into the source distribution.
libspatialops_a_SOURCES	=  \
	LinAlgTrilinos.cpp \
	FVStaggeredSpatialOps.cpp \
	LinearSystem.cpp \
	$(SOPS_INCDIR)/SpatialField.h \
	$(SOPS_INCDIR)/SpatialOpsDefs.h \
	$(SOPS_INCDIR)/SpatialFieldStore.h \
	$(SOPS_INCDIR)/SpatialOperator.h \
	$(SOPS_INCDIR)/FVStaggeredBCTools.h \
	$(SOPS_INCDIR)/FVStaggeredSpatialOps.h \
	$(SOPS_INCDIR)/LinAlgTrilinos.h \
	$(SOPS_INCDIR)/FV2ndOrderTypes.h \
	$(SOPS_INCDIR)/LinearSystem.h

#libspatialops_a_HEADERS = #headers to be installed
#libspatialops_a_noinst_HEADERS = #headers just used for compiling

#where is an extra_INSTALL target to specify libteuchos.a, etc. ?
#How about an extra_UNINSTALL target to clean up subpackages' headers?
