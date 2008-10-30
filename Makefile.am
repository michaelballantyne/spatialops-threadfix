#Makefile.am - edit me! 
# Process this file with 'automake' to produce Makefile.in
# Then ./configure to produce Makefile

# ==== Miscellaneous AutoMake Config

#main include directory
SOPS_INCDIR	= -I$(top_srcdir)/include -I$(top_srcdir)/test
#arch-dependent include directory (i.e. Daixtrose, since it's prebuilt at config time)
DEPS_INCDIR	= -I$(top_builddir)/include

EXTRA_DIST 	=	\
	debian/daixtrose-0.0.3.tar.gz \
	debian/changelog\
	debian/compat	\
	debian/control	\
	debian/copyright\
	debian/rules

EXTRA_HEADERS	=
library_include_HEADERS =@SOPS_HEADERS@
library_includedir	= $(includedir)/@PACKAGE@

#distribute the library pkg-config file.
pkgconfigdir 	= $(libdir)/pkgconfig
pkgconfig_DATA 	= spatialops-1.0.pc

nobase_include_HEADERS = \
	daixtrose/BinOps.h \
	daixtrose/ChangeDisambiguation.h \
	daixtrose/CompileTimeChecks.h \
	daixtrose/ConstRef.h \
	daixtrose/CountLeaves.h \
	daixtrose/CountOccurence.h \
	daixtrose/Daixt.h \
	daixtrose/DefaultOps.h \
	daixtrose/Differentiation.h \
	daixtrose/Disambiguation.h \
	daixtrose/Expr.h \
	daixtrose/ExtractOp.h \
	daixtrose/MatrixSelect.h \
	daixtrose/NeutralElements.h \
	daixtrose/PrettyPrinter.h \
	daixtrose/Scalar.h \
	daixtrose/Simplify.h \
	daixtrose/UnOps.h

#Distclean-local: a hook to make daixtrose clean up after itself
distclean-local:
	cd daixtrose-0.0.3; make uninstall; make distclean
	rm -f test/converge.out
	rm -f converge.out
	rm -rf include/daixtrose
	rm -f daixtrose
	rm -r daixtrose-0.0.3
	rm -rf include/tiny include/home include/error_handling include/linalg include/demos
	if test `ls include | wc -l` -eq 0; then rm -rf include; fi

# ==== Library Config ====

lib_LIBRARIES	= libspatialops.a

#preprocessor flags go here
AM_CPPFLAGS = -DBOOST_UBLAS_SHALLOW_ARRAY_ADAPTOR -DNDEBUG -DUINTAH_FIELD_TYPES \
	$(SOPS_INCDIR) $(DEPS_INCDIR) @BOOST_CPPFLAGS@
AM_LDFLAGS  = @BOOST_LDFLAGS@

#linker flags go here
#note: linker is very picky about library precedence order, so leave these as shown
SOPS_LIBS	= -lspatialops
SOPS_LIBDIR 	= -L$(top_builddir) 

#Optimize flags go here.
libspatialops_a_CXXFLAGS= -O4 -Wall -fexpensive-optimizations -funroll-loops

#Compile these files and pack them into a static library.
#Sources not listed here might not get packed into the source distribution.
libspatialops_a_SOURCES	=  \
	src/FVOneDimensional.cpp \
	src/LagrangePoly.cpp \
	src/LinAlgTrilinos.cpp \
	src/LinAlgUBlas.cpp \
	src/LinearSystem.cpp

#libspatialops_a_HEADERS = #headers to be installed
#libspatialops_a_noinst_HEADERS = #headers just used for compiling


# ==== Test Suite Config ====

test:	check
#targets to build when 'make check' is called
check_PROGRAMS	= simplebc testnew poisson lagrange testOneDim 
#targets to execute when 'make check' is called
#'testnew' requires interactive input, so it's disabled
TESTS 		= simplebc poisson lagrange testOneDim
#testsuites that are eXpected to fail when run
XFAIL_TESTS	= testOneDim

simplebc_SOURCES  =test/bcTest.cpp
simplebc_LDADD    =$(SOPS_LIBS) $(SOPS_LIBDIR)

testnew_SOURCES   =test/testNew.cpp test/buildOps.cpp test/Grid.h test/Functions.h test/buildOps.h
testnew_LDADD     =$(SOPS_LIBS) $(SOPS_LIBDIR)

poisson_SOURCES   =test/testPoisson.cpp test/buildOps.cpp
poisson_LDADD     =$(SOPS_LIBS) $(SOPS_LIBDIR)

lagrange_SOURCES  =test/testLagrange.cpp
lagrange_LDADD    =$(SOPS_LIBS) $(SOPS_LIBDIR)

testOneDim_SOURCES=test/testOneDim.cpp
testOneDim_LDADD  =$(SOPS_LIBS) $(SOPS_LIBDIR)

#disabled; not up-to-date and doesn't compile
#smalltest_SOURCES=test/test.cpp
#smalltest_LDAD   =$(SOPS_LIBS) $(SOPS_LIBDIR)

