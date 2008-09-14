SUBDIRS = src test 

#tarfiles.
EXTRA_DIST = tarfiles/daixtrose-0.0.3.tar.gz

SOPS_INCDIR=$(top_srcdir)/include
library_includedir=$(includedir)/spatialops
#maybe @AC_PACKAGE_STRING@ ? 

EXTRA_HEADERS	=
library_include_HEADERS =@SOPS_HEADERS@

#distribute the library pkg-config file.
pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = spatialops-1.0.pc

#Clean-local: a hook to make daixtrose clean up after itself
clean-local:
	cd daixtrose-0.0.3; make uninstall; make distclean
	rm -f test/converge.out
#	cd trilinos-7.0.6; make uninstall; make distclean

distclean-local:
	rm -r daixtrose-0.0.3

