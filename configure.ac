AC_INIT([SpatialOps], [1.0], [james.sutherland@utah.edu])
AM_INIT_AUTOMAKE([-Wall -Werror foreign])
#check if the compiler works.
#Use Fink defaults if available
if test `uname` = "Darwin" -a -e "/sw/bin/gcc-4"; then
	echo "Using Fink compilers."
	export CC="gcc-4"
	export CXX="g++-4"
	export F77="gfortran"
	export with_blas=/System/Library/Frameworks/vecLib.framework/vecLib
fi

AC_PROG_CXX
AC_PROG_CC
AC_PROG_F77

AC_PROG_RANLIB

#check for libraries.
#AC_CHECK_LIB (libblas, function)
#AC_CHECK_LIB (liblapack, function)
#AC_CHECK_LIB (library, function)


#check that the compiler isn't insane.
AC_HEADER_STDBOOL
AC_C_CONST

#check system library functions.
AC_CHECK_FUNCS([sqrt])

#Generate these makefiles when all done.
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_FILES([
  Makefile
  src/Makefile
  src/test/Makefile
])

AC_OUTPUT

TRILINOS_VERSION=7.0.6
TRILINOS_FLAGS=" --prefix=`pwd`\
  --includedir=`pwd`/include/trilinos \
  --disable-dependency-tracking \
  --disable-default-packages \
  --enable-aztecoo \
  --enable-epetra \
  --enable-epetraext \
  --disable-examples \
  --disable-tests \
"
#use this to invoke configure on subpackages (trilinos boost etc.)
if test ! -e "trilinos-7.0.6/config.status"; then
  BUILD_CONFIG_SUBDIR(trilinos-7.0.6,[ ${TRILINOS_FLAGS} ])
  cd trilinos-7.0.6
  make install
  cd ..
else
  echo "Trilinos seems to be configured already.  Skipping."
fi

#   AC_CONFIG_SUBDIRS(["trilinos-7.0.6"])
#fi
#echo "Configuring subpackage trilinos..."
#  cd trilinos-7.0.6
#  ./configure --prefix=`pwd`/../ \
#  make
#  cd ..


BOOST_VERSION=1_33_1
echo "Configuring subpackage boost."
if test ! -d include/boost; then
  cp -r boost_1_33_1/boost include/.
fi
#BUILD_CONFIG_SUBDIR(boost_1_33_1, [ --prefix=`pwd`/..])
#cd boost_1_33_1
#make install
#cd ..

#cd boost_1_33_1
#if test ! -e "bjam"; then
#   cd tools/build/jam_src
#   ./build.sh
#   cp `find . |grep bjam` ../../../
#   cd ../../../
#fi
#   ./bjam "-sTOOLS=gcc" --without-python --prefix=`pwd`/../ install
#   cd ../
#fi
echo "All done configuring boost."

DAIXTROSE_VERSION=0.0.3
DAIXTROSE_FLAGS="    \
  --prefix=`pwd` \
  --srcdir='' \
  --includedir=`pwd`/include/daixtrose
  --with-boost=`pwd`/include   \
"

if test ! -e "daixtrose-0.0.3/config.status"; then
  echo "Configuring subpackage daixtrose."
  BUILD_CONFIG_SUBDIR("daixtrose-0.0.3",[ ${DAIXTROSE_FLAGS} ])
  cd daixtrose-0.0.3
  make install
  cd ..
  echo "All done configuring daixtrose."
else
  echo "Daixtrose seems to be configured already.  Skipping."
fi

CFLAGS= #nothing; cleanup stray -O2 from children above?
CXXFLAGS= #nothing
#what are these for?
#PKG_CHECK_MODULES(SOPS)
#AC_SUBST(SOPS_CPPFLAGS)
#AC_SUBST(SOPS_LIBS)

echo "All done.  Type 'make <enter>' to build SpatialOps."
