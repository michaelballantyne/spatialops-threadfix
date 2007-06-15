#Makefile.am - dinky makefile input for test suite

AM_CPPFLAGS     = -I$(top_srcdir)/include -I$(top_srcdir)/include/trilinos -I$(top_srcdir)/include/daixtrose #preprocessor flags go here

TESTS 		= test smalltest #programs to automatically run when 'make check'
check_PROGRAMS	= test smalltest #targets to build when 'make check' called

#Right now, Alberto's test looks for its input in the source tree
#instead of the arch-build tree.  So we'll just ignore it for now...
XFAIL_TESTS	= test #tests to ignore

smalltest_SOURCES=test.cpp
smalltest_LDADD	= $(SOPS_LIBS) $(SOPS_LIBDIR)

test_SOURCES 	= test_alberto.cpp
test_LDADD	= $(SOPS_LIBS) $(SOPS_LIBDIR)
test_CPPFLAGS	= -I./test_alberto $(AM_CPPFLAGS)

SOPS_LIBDIR 	= -L../ -L$(top_srcdir)/lib 
SOPS_LIBS	= -lspatialops -laztecoo -lteuchos -lepetra -lepetraext -lblas -llapack
#note: linker is very picky about library precedence order, so leave these as shown
