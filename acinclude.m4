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


