#Makefile.am - dinky makefile input for test suite
## Makefile.am -- edit me!
# Process this file with 'automake' to produce Makefile.in
# Then ./configure to produce Makefile

#SpatialOps sources are not always near the build location.
#For example, we might build each arch in its own directory.
#We need to refer to headers realtive to the sources' path,
#not the build path. ($srcdir instead of $top_srcdir)
SOPS_INCDIR	= $(srcdir)/../../include

#Dependency package headers _are_ near our build location
# (since they're installed on a per-arch basis.)  Therefore
# they're relative to the build path, not the source path.
DEPS_INCDIR	= ../../include

#preprocessor flags go here
AM_CPPFLAGS	= \
	-I$(SOPS_INCDIR) \
	-I$(DEPS_INCDIR) \
	-I$(DEPS_INCDIR)/daixtrose \
	-I$(DEPS_INCDIR)/trilinos \
	-I$(DEPS_INCDIR)/boost \
	$(SOPS_CPPFLAGS)

TESTS 		= test smalltest testnew #programs to automatically run when 'make check'
check_PROGRAMS	= test smalltest testnew #targets to build when 'make check' called

EXTRA_DIST	= test.input #non-source file that needs to be packed up too

#Right now, Alberto's test looks for its input in the source tree
#instead of the arch-build tree.  So we'll just ignore it for now...
XFAIL_TESTS	= test #tests to ignore

smalltest_SOURCES=test.cpp
smalltest_LDADD	= $(SOPS_LIBS) $(SOPS_LIBDIR)

test_SOURCES 	= test_alberto.cpp test_alberto/analytical_class.h test_alberto/grid.h test_alberto/test_functions.h test_alberto/utilities.h
test_LDADD	= $(SOPS_LIBS) $(SOPS_LIBDIR)
test_CPPFLAGS	= -I$(srcdir)/test_alberto $(AM_CPPFLAGS)

testnew_SOURCES	= testNew.cpp Functions.h Grid.h
testnew_LDADD	= $(SOPS_LIBS) $(SOPS_LIBDIR)

SOPS_LIBDIR 	= -L../ -L$(DEPS_INCDIR)/../lib 
SOPS_LIBS	= -lspatialops -laztecoo -lteuchos -lepetra -lepetraext -lblas -llapack
#note: linker is very picky about library precedence order, so leave these as shown
