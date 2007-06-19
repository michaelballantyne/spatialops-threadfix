SUBDIRS = src src/test 

distcheck-hook:
	echo "faking distcheck dependencies."
	chmod a+w spatialops-1.0
	(find trilinos-7.0.6 -name config.status |xargs rm) ||true
	(find daixtrose-0.0.3 -name config.status |xargs rm)||true
	ln -s ../trilinos-7.0.6 spatialops-1.0/
	ln -s ../boost_1_33_1 spatialops-1.0/
	ln -s ../daixtrose-0.0.3 spatialops-1.0/
