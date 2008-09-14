#Makefile.am - dinky makefile input for test suite
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

#note: linker is very picky about library precedence order, so leave these as shown
SOPS_LIBS	= -lspatialops -laztecoo -lteuchos -lepetra -lepetraext -lblas -llapack
SOPS_LIBDIR 	= -L../ -L$(top_builddir)/src 

#targets to build when 'make check' is called
check_PROGRAMS	= simplebc testnew poisson lagrange testOneDim 
#targets to execute when 'make check' is called
#'testnew' requires interactive input, so it's disabled
TESTS 		= simplebc poisson lagrange testOneDim
#testsuites that are expected to fail when run
FAIL_TESTS	= #none

simplebc_SOURCES  =bcTest.cpp
simplebc_LDADD    =$(SOPS_LIBS) $(SOPS_LIBDIR)

testnew_SOURCES   =testNew.cpp buildOps.cpp Grid.h Functions.h buildOps.h
testnew_LDADD     =$(SOPS_LIBS) $(SOPS_LIBDIR)

poisson_SOURCES   =testPoisson.cpp buildOps.cpp
poisson_LDADD     =$(SOPS_LIBS) $(SOPS_LIBDIR)

lagrange_SOURCES  =testLagrange.cpp
lagrange_LDADD    =$(SOPS_LIBS) $(SOPS_LIBDIR)

testOneDim_SOURCES=testOneDim.cpp
testOneDim_LDADD  =$(SOPS_LIBS) $(SOPS_LIBDIR)

#disabled; not up-to-date and doesn't compile
#smalltest_SOURCES=test.cpp
#smalltest_LDADD	= $(SOPS_LIBS) $(SOPS_LIBDIR)

