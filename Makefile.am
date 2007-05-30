#bin_PROGRAMS = test
#test_SOURCES = i
#test_CFLAGS = 
SUBDIRS = 	src \
		src/test \
		# trilinos-7.0.6
#test_LDFLAGS = @LDFLAGS@ @GTK_LIBS@ -lgthread
#INCLUDES = @GTK_CFLAGS@
BOOST_VERSION=1_33_1
boost: jam
	cd boost_${BOOST_VERSION}; ./bjam -sTOOLS=gcc --prefix=`pwd`/../ install

boost-clean:;

jam:
	cd boost_${BOOST_VERSION}/tools/build/jam_src; sh build.sh gcc; cp `find -name bjam` ../../../

SHELL=/bin/bash
TRILINOS_VERSION=7.0.6
trilinos:
#wget http://software.sandia.gov/trilinos/downloads/files/trilinos-$(TRILINOS_VERSION).tar.gz
#trilinos-$(TRILINOS_VERSION)/configure
	cd trilinos-$(TRILINOS_VERSION) #; $(MAKE)
	if test ! -e config.status; \
	then \
	  ./configure; \
	fi
	$(MAKE)

trilinos-clean:
	$(MAKE) -C trilinos/trilinos-$(TRILINOS_VERSION) clean


