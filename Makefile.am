SUBDIRS = src 

SOPS_INCDIR=$(top_srcdir)/include
library_includedir=$(includedir)/spatialops
#maybe @AC_PACKAGE_STRING@ ? 
library_include_HEADERS = \
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
distcheck-hook:
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


#Extra install targets.
#Could be trilinos libraries, etc.
#XXX need to fix distclean so it recurses too.
EXTRA_INSTALLLIBS	= `ls lib/*`
EXTRA_INSTALLHEADERS	= `cd include; ls *config.h daixtrose/daixtrose/* trilinos/*`	#just install anything that winds up in our built lib directory
EXTRA_INSTALL		= $(EXTRA_INSTALLHEADERS) $(EXTRA_INSTALLLIBS)
libLIBRARIES_INSTALL 	= $(INSTALL_DATA)


install-exec-local:
	$(mkinstalldirs) $(DESTDIR)
	for file in $(EXTRA_INSTALLLIBS); do \
	  $(libLIBRARIES_INSTALL) $$file $(DESTDIR)$(libdir); \
	done
	for file in $(EXTRA_INSTALLHEADERS); do \
	  echo "$(INSTALL_DATA) -D $$file $(DESTDIR)$(pkgincludedir)/$$file"; \
	  (cd include; $(INSTALL_DATA) -D $$file $(DESTDIR)$(pkgincludedir)/$$file); \
	done

uninstall-local:
	for file in $(EXTRA_INSTALLLIBS); do \
	  rm -f $(DESTDIR)$(libdir)/$$file; \
	done
	for file in $(EXTRA_INSTALLHEADERS); do \
	  rm -f $(DESTDIR)$(pkgincludedir)/$$file; \
	done
