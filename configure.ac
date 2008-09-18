AC_INIT([SpatialOps], [0.3a], [james.sutherland@utah.edu])

#show where to get all our pieces
#All the buildsystem pieces are tucked away neatly
# out of sight in the "buildsys" directory
AC_CONFIG_AUX_DIR(build-aux)

AM_INIT_AUTOMAKE([-Wall -Werror foreign])


if test -x "./configure"; then
	AC_MSG_ERROR([You cannot configure in the source tree.  Make an empty build directory and configure there instead.])
fi

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
AC_CHECK_LIB(blas,  dgemm_ )
AC_CHECK_LIB(lapack,dgecon_)


AX_BOOST_BASE(1.33.1)

#AX_TRILINOS_BASE(7.0.6, [] , [AC_MSG_ERROR([Could not find Trilinos.])])  
#Requisite Trilinos Libraries >=7.0.6 

#set our buildsystem to use C++ tests.
#does this screw up C++ header filenames when it comes time to compile?
# e.g. .h vs .hpp
AC_LANG([C++])

AC_ARG_WITH([trilinos], [AS_HELP_STRING([--with-trilinos], [Trilinos package, needed for linear algebra routines @<:@default=check@:>@])],
  [TRILINOS_INCLUDE=$withval/include
   TRILINOS_LIBS=$withval/lib
   CPPFLAGS="$CPPFLAGS -I$TRILINOS_INCLUDE"
   LDFLAGS="$LDFLAGS -L$TRILINOS_LIBS"]
)

trilinos_fail=no
AC_CHECK_LIB(epetra,	exit,		[], [trilinos_fail=yes])
AC_CHECK_LIB(epetraext,	main,		[], [trilinos_fail=yes])
AC_CHECK_LIB(teuchos,	Teuchos_exit_helper,[], [trilinos_fail=yes])
AC_CHECK_LIB(triutils,	main,		[], [trilinos_fail=yes])
AC_CHECK_LIB(aztecoo,	AZ_broadcast,	[], [trilinos_fail=yes])

AS_IF([test "x$trilinos_fail" != "xno"],
  [AC_MSG_FAILURE([Trilinos seems to be improperly installed.])
])

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
  spatialops-1.0.pc
])

#List all the headers in /include
#so they can get packed into the archive.
SOPS_HEADERS="`cd $srcdir/include; ls -D *.h|sed -e 's/^/$(srcdir)\/include\//g'| tr '\n' ' '`"
AC_SUBST(SOPS_HEADERS)

AC_OUTPUT

DAIXTROSE_VERSION=0.0.3
DAIXTROSE_FLAGS="    \
  --prefix=`pwd` \
  --srcdir='' \
  --includedir=`pwd`/include
  --with-boost=`pwd`/include   \
"

if test ! -e "daixtrose-0.0.3/config.status"; then
  echo "Configuring subpackage daixtrose."
  pushd $srcdir; chmod u+w .; tar xzf tarfiles/daixtrose-0.0.3.tar.gz; popd
  mkdir_p='mkdir -p'
  BUILD_CONFIG_SUBDIR("daixtrose-0.0.3",[ ${DAIXTROSE_FLAGS} ])
  cd daixtrose-0.0.3
  make install
  cd ..
  ln -s include/daixtrose daixtrose
  echo "All done configuring daixtrose."
else
  echo "Daixtrose seems to be configured already.  Skipping."
fi

CFLAGS= #nothing; cleanup stray -O2 from children above?
CXXFLAGS= #nothing

echo "All done.  Type 'make <enter>' to build SpatialOps."
