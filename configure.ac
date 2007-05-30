echo "Alakazaam!  Configuring the Spiffy Spatial Ops Build System 0.1"

AC_INIT([SpatialOps], [1.0], [james.sutherland@utah.edu])
AM_INIT_AUTOMAKE([-Wall -Werror foreign])

#check if the compiler works.
AC_PROG_CXX
AC_PROG_CC
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

AC_CONFIG_HEADERS([config.h])
AC_CONFIG_FILES([
 Makefile
 src/Makefile
 src/test/Makefile
])
#trilinos-7.0.6/Makefile
#])

TRILINOS_VERSION=7.0.6
#use this to invoke configure on subpackages (trilinos boost etc.)
if test ! -e "trilinos-7.0.6/config.status"; then
#   AC_CONFIG_SUBDIRS(["trilinos-7.0.6"])
#fi
echo "Configuring subpackage trilinos..."
  cd trilinos-7.0.6
  ./configure --prefix=`pwd`/../
  make
  cd ..
else
echo "Trilinos seems to be configured already.  Skipping."
fi


BOOST_VERSION=1_33_1

echo "Configuring subpackage boost."
cd boost_1_33_1
if test ! -e "bjam"; then
   cd tools/build/jam_src
   ./build.sh
   cp `find -name bjam` ../../../
   cd ../../../
fi
   ./bjam "-sTOOLS=gcc" --prefix=`pwd`/../ install
   cd ../
#fi
echo "All done configuring boost."

#PKG_CHECK_MODULES(SOPS)
AC_SUBST(SOPS_CPPFLAGS)
AC_SUBST(SOPS_LIBS)
AC_OUTPUT

echo "All done.  Type 'make <enter>' to build SpatialOps."
