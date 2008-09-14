dnl AX_TRILINOS_BASE
dnl http://autoconf-archive.cryp.to/ax_trilinos_base.html
dnl checks version of Trilinos libraries

AC_DEFUN([AX_TRILINOS_BASE],
[
AC_REQUIRE([AC_LANG_CPLUSPLUS])
AC_ARG_VAR(TRILINOS_HOME,[Location of Trilinos installation])

AC_ARG_WITH(trilinos, [AS_HELP_STRING([--with-trilinos[=DIR]],[root directory of Trilinos installation])],[
with_trilinos=$withval
if test "${with_trilinos}" != yes; then
    TRILINOS_HOME=$withval
        trilinos_include="$withval/include"
        trilinos_libdir="$withval/lib"
fi
],[
with_trilinos=$withval
if test "x${TRILINOS_HOME}" != "x"; then
        trilinos_include="${TRILINOS_HOME}/include"
        trilinos_libdir="${TRILINOS_HOME}/lib"
fi
])

AC_ARG_WITH(trilinos-include,
[AS_HELP_STRING([--with-trilinos-include=DIR],[specify exact directory for Trilinos headers])],[
if test -d $withval; then
    trilinos_include="$withval"
else
    AC_MSG_ERROR([--with-trilinos-include expected directory name])
fi
])

AC_ARG_WITH(trilinos-libdir, [AS_HELP_STRING([--with-trilinos-libdir=DIR],[specify exact directory for Trilinos libraries])],[
if test -d $withval; then
    trilinos_libdir="$withval"
else
    AC_MSG_ERROR([--with-trilinos-libdir expected directory name])
fi
])

if test "${with_trilinos}" != no ; then
        OLD_LIBS=$LIBS
        OLD_LDFLAGS=$LDFLAGS
        OLD_CFLAGS=$CFLAGS
        OLD_CPPFLAGS=$CPPFLAGS

        if test "${trilinos_libdir}" ; then
                LDFLAGS="$LDFLAGS -L${trilinos_libdir}"
        fi
        if test "${trilinos_include}" ; then
                CPPFLAGS="$CPPFLAGS -I${trilinos_include}"
                CFLAGS="$CFLAGS -I${trilinos_include}"
        fi

    succeeded=no
        AC_CHECK_HEADER([Trilinos_version.h],found_header=yes,found_header=no)
        if test "$found_header" = yes; then
        dnl Patterned after AX_BOOST_BASE
        trilinos_lib_version_req=ifelse([$1],,8.0.0,$1)
        trilinos_lib_version_req_shorten=`expr $trilinos_lib_version_req : '\([[0-9]]*\.[[0-9]]*\)'`
        trilinos_lib_version_req_major=`expr $trilinos_lib_version_req : '\([[0-9]]*\)'`
        trilinos_lib_version_req_minor=`expr $trilinos_lib_version_req : '[[0-9]]*\.\([[0-9]]*\)'`
        trilinos_lib_version_req_sub_minor=`expr $trilinos_lib_version_req : '[[0-9]]*\.[[0-9]]*\.\([[0-9]]*\)'`
        if test "x$trilinos_lib_version_req_sub_minor" = "x" ; then
                trilinos_lib_version_req_sub_minor="0"
        fi
        WANT_TRILINOS_VERSION=`expr $trilinos_lib_version_req_major \* 10000 \+  $trilinos_lib_version_req_minor \* 100 \+ $trilinos_lib_version_req_sub_minor`
        AC_MSG_CHECKING(for Trilinos release >= $trilinos_lib_version_req)
        AC_LANG_PUSH(C++)
        AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[
@%:@include <Trilinos_version.h>
]], [[
#if (TRILINOS_MAJOR_VERSION == 0)
/* Development branch has zero major version. Everything is okay. */
#elif (TRILINOS_MAJOR_MINOR_VERSION >= $WANT_TRILINOS_VERSION)
/* Stable release of appropriate version. Everything is okay. */
#else
#  error Trilinos version is too old
#endif
        ]])],[
            AC_MSG_RESULT(yes)
            succeeded=yes
            ],[
            AC_MSG_RESULT(no)
            ])
        AC_LANG_POP([C++])
    fi

    if test "$succeeded" = no; then
                LIBS=$OLD_LIBS
                LDFLAGS=$OLD_LDFLAGS
                CPPFLAGS=$OLD_CPPFLAGS
                CFLAGS=$OLD_CFLAGS

                ifelse([$3], , , [$3])
        else
                AC_DEFINE(HAVE_TRILINOS,1,[Define if Trilinos is available])
                ifelse([$2], , , [$2])
        fi

fi

])


dnl AX_BOOST_BASE
dnl http://autoconf-archive.cryp.to/ax_boost_base.html
dnl checks version of Boost C++ headers
AC_DEFUN([AX_BOOST_BASE],
[
AC_ARG_WITH([boost],
        AS_HELP_STRING([--with-boost@<:@=DIR@:>@], [use boost (default is yes) - it is possible to specify the root directory for boost (optional)]),
        [
    if test "$withval" = "no"; then
                want_boost="no"
    elif test "$withval" = "yes"; then
        want_boost="yes"
        ac_boost_path=""
    else
            want_boost="yes"
        ac_boost_path="$withval"
        fi
    ],
    [want_boost="yes"])


AC_ARG_WITH([boost-libdir],
        AS_HELP_STRING([--with-boost-libdir=LIB_DIR],
        [Force given directory for boost libraries. Note that this will overwrite library path detection, so use this parameter only if default library detection fails and you know exactly where your boost libraries are located.]),
        [
        if test -d $withval
        then
                ac_boost_lib_path="$withval"
        else
                AC_MSG_ERROR(--with-boost-libdir expected directory name)
        fi
        ],
        [ac_boost_lib_path=""]
)

if test "x$want_boost" = "xyes"; then
        boost_lib_version_req=ifelse([$1], ,1.20.0,$1)
        boost_lib_version_req_shorten=`expr $boost_lib_version_req : '\([[0-9]]*\.[[0-9]]*\)'`
        boost_lib_version_req_major=`expr $boost_lib_version_req : '\([[0-9]]*\)'`
        boost_lib_version_req_minor=`expr $boost_lib_version_req : '[[0-9]]*\.\([[0-9]]*\)'`
        boost_lib_version_req_sub_minor=`expr $boost_lib_version_req : '[[0-9]]*\.[[0-9]]*\.\([[0-9]]*\)'`
        if test "x$boost_lib_version_req_sub_minor" = "x" ; then
                boost_lib_version_req_sub_minor="0"
        fi
        WANT_BOOST_VERSION=`expr $boost_lib_version_req_major \* 100000 \+  $boost_lib_version_req_minor \* 100 \+ $boost_lib_version_req_sub_minor`
        AC_MSG_CHECKING(for boostlib >= $boost_lib_version_req)
        succeeded=no

        dnl first we check the system location for boost libraries
        dnl this location ist chosen if boost libraries are installed with the --layout=system option
        dnl or if you install boost with RPM
        if test "$ac_boost_path" != ""; then
                BOOST_LDFLAGS="-L$ac_boost_path/lib"
                BOOST_CPPFLAGS="-I$ac_boost_path/include"
        else
                for ac_boost_path_tmp in /usr /usr/local /opt /opt/local ; do
                        if test -d "$ac_boost_path_tmp/include/boost" && test -r "$ac_boost_path_tmp/include/boost"; then
                                BOOST_LDFLAGS="-L$ac_boost_path_tmp/lib"
                                BOOST_CPPFLAGS="-I$ac_boost_path_tmp/include"
                                break;
                        fi
                done
        fi

    dnl overwrite ld flags if we have required special directory with
    dnl --with-boost-libdir parameter
    if test "$ac_boost_lib_path" != ""; then
       BOOST_LDFLAGS="-L$ac_boost_lib_path"
    fi

        CPPFLAGS_SAVED="$CPPFLAGS"
        CPPFLAGS="$CPPFLAGS $BOOST_CPPFLAGS"
        export CPPFLAGS

        LDFLAGS_SAVED="$LDFLAGS"
        LDFLAGS="$LDFLAGS $BOOST_LDFLAGS"
        export LDFLAGS

        AC_LANG_PUSH(C++)
        AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[
        @%:@include <boost/version.hpp>
        ]], [[
        #if BOOST_VERSION >= $WANT_BOOST_VERSION
        // Everything is okay
        #else
        #  error Boost version is too old
        #endif
        ]])],[
        AC_MSG_RESULT(yes)
        succeeded=yes
        found_system=yes
        ],[
        ])
        AC_LANG_POP([C++])



        dnl if we found no boost with system layout we search for boost libraries
        dnl built and installed without the --layout=system option or for a staged(not installed) version
        if test "x$succeeded" != "xyes"; then
                _version=0
                if test "$ac_boost_path" != ""; then
                        if test -d "$ac_boost_path" && test -r "$ac_boost_path"; then
                                for i in `ls -d $ac_boost_path/include/boost-* 2>/dev/null`; do
                                        _version_tmp=`echo $i | sed "s#$ac_boost_path##" | sed 's/\/include\/boost-//' | sed 's/_/./'`
                                        V_CHECK=`expr $_version_tmp \> $_version`
                                        if test "$V_CHECK" = "1" ; then
                                                _version=$_version_tmp
                                        fi
                                        VERSION_UNDERSCORE=`echo $_version | sed 's/\./_/'`
                                        BOOST_CPPFLAGS="-I$ac_boost_path/include/boost-$VERSION_UNDERSCORE"
                                done
                        fi
                else
                        for ac_boost_path in /usr /usr/local /opt /opt/local ; do
                                if test -d "$ac_boost_path" && test -r "$ac_boost_path"; then
                                        for i in `ls -d $ac_boost_path/include/boost-* 2>/dev/null`; do
                                                _version_tmp=`echo $i | sed "s#$ac_boost_path##" | sed 's/\/include\/boost-//' | sed 's/_/./'`
                                                V_CHECK=`expr $_version_tmp \> $_version`
                                                if test "$V_CHECK" = "1" ; then
                                                        _version=$_version_tmp
                                                        best_path=$ac_boost_path
                                                fi
                                        done
                                fi
                        done

                        VERSION_UNDERSCORE=`echo $_version | sed 's/\./_/'`
                        BOOST_CPPFLAGS="-I$best_path/include/boost-$VERSION_UNDERSCORE"
            if test "$ac_boost_lib_path" = ""
            then
               BOOST_LDFLAGS="-L$best_path/lib"
            fi

                        if test "x$BOOST_ROOT" != "x"; then
                                if test -d "$BOOST_ROOT" && test -r "$BOOST_ROOT" && test -d "$BOOST_ROOT/stage/lib" && test -r "$BOOST_ROOT/stage/lib"; then
                                        version_dir=`expr //$BOOST_ROOT : '.*/\(.*\)'`
                                        stage_version=`echo $version_dir | sed 's/boost_//' | sed 's/_/./g'`
                                        stage_version_shorten=`expr $stage_version : '\([[0-9]]*\.[[0-9]]*\)'`
                                        V_CHECK=`expr $stage_version_shorten \>\= $_version`
                    if test "$V_CHECK" = "1" -a "$ac_boost_lib_path" = "" ; then
                                                AC_MSG_NOTICE(We will use a staged boost library from $BOOST_ROOT)
                                                BOOST_CPPFLAGS="-I$BOOST_ROOT"
                                                BOOST_LDFLAGS="-L$BOOST_ROOT/stage/lib"
                                        fi
                                fi
                        fi
                fi

                CPPFLAGS="$CPPFLAGS $BOOST_CPPFLAGS"
                export CPPFLAGS
                LDFLAGS="$LDFLAGS $BOOST_LDFLAGS"
                export LDFLAGS

                AC_LANG_PUSH(C++)
                AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[
                @%:@include <boost/version.hpp>
                ]], [[
                #if BOOST_VERSION >= $WANT_BOOST_VERSION
                // Everything is okay
                #else
                #  error Boost version is too old
                #endif
                ]])],[
                AC_MSG_RESULT(yes)
                succeeded=yes
                found_system=yes
                ],[
                ])
                AC_LANG_POP([C++])
        fi

        if test "$succeeded" != "yes" ; then
                if test "$_version" = "0" ; then
                        AC_MSG_ERROR([[We could not detect the boost libraries (version $boost_lib_version_req_shorten or higher). If you have a staged boost library (still not installed) please specify \$BOOST_ROOT in your environment and do not give a PATH to --with-boost option.  If you are sure you have boost installed, then check your version number looking in <boost/version.hpp>. See http://randspringer.de/boost for more documentation.]])
                else
                        AC_MSG_NOTICE([Your boost libraries seems to old (version $_version).])
                fi
        else
                AC_SUBST(BOOST_CPPFLAGS)
                AC_SUBST(BOOST_LDFLAGS)
                AC_DEFINE(HAVE_BOOST,,[define if the Boost library is available])
        fi

        CPPFLAGS="$CPPFLAGS_SAVED"
        LDFLAGS="$LDFLAGS_SAVED"
fi

])




dnl
dnl BUILD_CONFIG_SUBDIR(dir [, sub-package-cmdline-args, args-to-drop])
dnl
dnl This was stolen from Apache2.2's APR_SUBDIR_CONFIG in file
dnl    httpd-2.2.4/build/apr_common.m4
dnl which is in turn a hacked macro to replace AutoConf's AC_CONFIG_SUBDIRS.
dnl "Your mileage may vary."
dnl
dnl dir: directory to find configure in
dnl sub-package-cmdline-args: arguments to add to the invocation (optional)
dnl args-to-drop: arguments to drop from the invocation (optional)
dnl
dnl Note: This macro relies on ac_configure_args being set properly.
dnl
dnl The args-to-drop argument is shoved into a case statement, so
dnl multiple arguments can be separated with a |.
dnl
dnl Note: Older versions of autoconf do not single-quote args, while 2.54+
dnl places quotes around every argument.  So, if you want to drop the
dnl argument called --enable-layout, you must pass the third argument as:
dnl [--enable-layout=*|\'--enable-layout=*]
dnl
dnl Trying to optimize this is left as an exercise to the reader who wants
dnl to put up with more autoconf craziness.  I give up.
dnl
AC_DEFUN([BUILD_CONFIG_SUBDIR], [
  # save our work to this point; this allows the sub-package to use it
  AC_CACHE_SAVE

  echo "configuring package in $1 now"
  ac_popdir=`pwd`
  apr_config_subdirs="$1"
  test -d $1 || $mkdir_p $1
  ac_abs_srcdir=`(cd $srcdir/$1 && pwd)`
  cd $1

changequote(, )dnl
      # A "../" for each directory in /$config_subdirs.
      ac_dots=`echo $apr_config_subdirs|sed -e 's%^\./%%' -e 's%[^/]$%&/%' -e 's%[^/]*/%../%g'`
changequote([, ])dnl

  # Make the cache file pathname absolute for the subdirs
  # required to correctly handle subdirs that might actually
  # be symlinks
  case "$cache_file" in
  /*) # already absolute
    ac_sub_cache_file=$cache_file ;;
  *)  # Was relative path.
    ac_sub_cache_file="$ac_popdir/$cache_file" ;;
  esac

  ifelse($3, [], [apr_configure_args=$ac_configure_args],[
  apr_configure_args=
  apr_sep=
  for apr_configure_arg in $ac_configure_args
  do
    case "$apr_configure_arg" in
      $3)
        continue ;;
    esac
    apr_configure_args="$apr_configure_args$apr_sep'$apr_configure_arg'"
    apr_sep=" "
  done
  ])

  # autoconf doesn't add --silent to ac_configure_args; explicitly pass it
  test "x$silent" = "xyes" && apr_configure_args="$apr_configure_args --silent"

  dnl The eval makes quoting arguments work - specifically $2 where the
  dnl quoting mechanisms used is "" rather than [].
  dnl
  dnl We need to execute another shell because some autoconf/shell combinations
  dnl will choke after doing repeated BUILD_CONFIG_SUBDIR()s.  (Namely Solaris
  dnl and autoconf-2.54+)
  echo executing $SHELL $ac_abs_srcdir/configure $apr_configure_args --cache-file=$ac_sub_cache_file --srcdir=$ac_abs_srcdir $2
  if eval $SHELL $ac_abs_srcdir/configure $apr_configure_args --cache-file=$ac_sub_cache_file --srcdir=$ac_abs_srcdir $2
  then :
    echo "$1 configured properly"
  else
    echo "configure failed for $1"
    exit 1
  fi

  cd $ac_popdir

  # grab any updates from the sub-package
  AC_CACHE_LOAD
])dnl


