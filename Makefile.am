SUBDIRS = src 

#tarfiles.
EXTRA_DIST = tarfiles/trilinos-7.0.6.tar.gz \
	tarfiles/boost_1_33_1.tar.gz  \
	tarfiles/daixtrose-0.0.3.tar.gz


SOPS_INCDIR=$(top_srcdir)/include
library_includedir=$(includedir)/spatialops
#maybe @AC_PACKAGE_STRING@ ? 
EXTRA_HEADERS	=
library_include_HEADERS =@SOPS_HEADERS@
old_stuffders = \
        $(SOPS_INCDIR)/SpatialField.h \
        $(SOPS_INCDIR)/SpatialOpsDefs.h \
        $(SOPS_INCDIR)/SpatialFieldStore.h \
        $(SOPS_INCDIR)/SpatialOperator.h \
        $(SOPS_INCDIR)/FVStaggeredBCTools.h \
        $(SOPS_INCDIR)/FVStaggeredSpatialOps.h \
        $(SOPS_INCDIR)/LinAlgTrilinos.h \
        $(SOPS_INCDIR)/FV2ndOrderTypes.h \
        $(SOPS_INCDIR)/LinearSystem.h

#distribute the library pkg-config file.
pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = spatialops-1.0.pc

#Special code to run during distcheck.
#This script executed one level up from the source tree
#that distcheck creates for itself.
disabled-distcheck-hook:
	echo "faking distcheck dependencies."
	chmod a+w spatialops-1.0
	(find trilinos-7.0.6 -name config.status |xargs rm) ||true
	(find daixtrose-0.0.3 -name config.status |xargs rm)||true
	echo "Reusing sources directories:"
	ln -s ../../trilinos-7.0.6 spatialops-1.0/
	ln -s ../../boost_1_33_1 spatialops-1.0/
	ln -s ../../daixtrose-0.0.3 spatialops-1.0/
	echo "Reusing build directories:"
	ln -s ../../trilinos-7.0.6 spatialops-1.0/_build
	ln -s ../../boost-1_33_1 spatialops-1.0/_build
	ln -s ../../daixtrose-0.0.3 spatialops-1.0/_build

#Clean-local: a hook to make trilinos, boost, etc clean up after themselves
#FIXME: daixtrose barfs on this currently
clean-local:
	cd trilinos-7.0.6; make uninstall; make distclean
	cd daixtrose-0.0.3; make uninstall; make distclean

#Extra install targets.
#Could be trilinos libraries, etc.
#XXX need to fix distclean so it recurses too.
EXTRA_INSTALLLIBS	= `cd lib; ls *`
EXTRA_INSTALLHEADERS	= `cd include; ls *config.h daixtrose/daixtrose/* trilinos/*`	#just install anything that winds up in our built lib directory
EXTRA_INSTALL		= $(EXTRA_INSTALLHEADERS) $(EXTRA_INSTALLLIBS)
libLIBRARIES_INSTALL 	= $(INSTALL_DATA)


install-exec-local:
	$(mkinstalldirs) $(DESTDIR)
	for file in $(EXTRA_INSTALLLIBS); do \
	  echo "$(libLIBRARIES_INSTALL) lib/$$file $(DESTDIR)$(libdir)/$$file"; \
	  $(libLIBRARIES_INSTALL) lib/$$file $(DESTDIR)$(libdir)/$$file; \
	done
	for file in $(EXTRA_INSTALLHEADERS); do \
	  echo "(cd include; $(INSTALL_DATA) -D $$file $(DESTDIR)$(pkgincludedir)/$$file)"; \
	  (cd include; $(INSTALL_DATA) -D $$file $(DESTDIR)$(pkgincludedir)/$$file); \
	done

uninstall-local:
	for file in $(EXTRA_INSTALLLIBS); do \
	  echo "rm -f $(DESTDIR)$(libdir)/$$file"; \
	  rm -f $(DESTDIR)$(libdir)/$$file; \
	done
	for file in $(EXTRA_INSTALLHEADERS); do \
	  echo "rm -f $(DESTDIR)$(pkgincludedir)/$$file"; \
	  rm -f $(DESTDIR)$(pkgincludedir)/$$file; \
	done
