#Makefile.am - dinky makefile input for test suite
check_PROGRAMS 	=test
test_SOURCES 	=test.cpp 
test_LDADD 	= $(SOPS_LIBDIR) $(EPETRA_LIBS) $(LIBOBJS)  #define required libraries in configure.ac ??
test_CPPFLAGS 	= -I$(top_srcdir)/include $(SOPS_CPPFLAGS)

EPETRA_LIBS 	= -lepetra -lepetraext -lblas -llapack -lspatialops -laztecoo -lteuchos
SOPS_LIBDIR	= -L$(top_srcdir)/lib -L../

