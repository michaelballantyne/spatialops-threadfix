SUBDIRS = src 

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
EXTRA_INSTALL 		= `ls lib`	#just install anything that winds up in our built lib directory
libLIBRARIES_INSTALL 	= $(INSTALL_DATA)


install-exec-local:
	$(mkinstalldirs) $(DESTDIR)$(libdir)
	for file in $(EXTRA_INSTALL); do \
	  echo " $(libLIBRARIES_INSTALL) lib/$$file \\"; \
	  echo "        $(DESTDIR)$(libdir)"; \
	  $(libLIBRARIES_INSTALL) lib/$$file $(DESTDIR)$(libdir); \
	done

uninstall-local:
	for file in $(EXTRA_INSTALL); do \
	  rm -f $(DESTDIR)$(libdir)/$$file; \
	done
