#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.
AC_PREREQ(2.59)
AC_INIT(sga, 0.10.10, js18@sanger.ac.uk)
# move aux files to subdir
AC_CONFIG_AUX_DIR([auxfiles])
AM_INIT_AUTOMAKE(foreign)
AC_CONFIG_SRCDIR([SGA/sga.cpp])
AC_CONFIG_HEADER([config.h])
#m4_include([m4/m4_ax_openmp.m4])
AC_CONFIG_MACRO_DIR([m4])

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC
#AC_PROG_RANLIB

# don't build shared libs for now; must appear before AC_PROG_LIBTOOL to take effect
AC_DISABLE_SHARED
# check for libtool
AC_PROG_LIBTOOL

# Check for headers
AC_LANG([C++])
AC_CHECK_HEADER([zlib.h])

# Check for libraries
AC_SEARCH_LIBS([gzopen],[z],,[AC_MSG_ERROR([libz not found, please install zlib (http://www.zlib.net/)])])
AC_SEARCH_LIBS([clock_gettime], [rt], [AC_DEFINE([HAVE_CLOCK_GETTIME], [1], [clock_getttime found])], )

# Check for openmp
#AX_OPENMP([openmp_cppflags="-fopenmp" AC_DEFINE(HAVE_OPENMP,1,[Define if OpenMP is enabled])])
AC_OPENMP([AC_DEFINE(HAVE_OPENMP, 1, [OpenMP is enabled])])
AC_SUBST(OPENMP_CXXFLAGS)

# Check for the bamtools library path
# Bamtools has two different path formats
# If it is built in place the files will be in {path}/lib/ and {path}/include/
# If it is installed (with make install) it will be in {path}/include/bamtools
# and {path}/lib/bamtools. Both variants are handled here
AC_ARG_WITH(bamtools, AS_HELP_STRING([--with-bamtools=PATH],
	[specify directory containing the bamtools library (http://github.com/pezmaster31/bamtools)]))

if test "$with_bamtools" -a -d "$with_bamtools"; then
    #bamtools_ldflags="-Wl,-rpath,$with_bamtools/lib -L$with_bamtools/lib -Wl,-rpath,$with_bamtools/lib/bamtools -L$with_bamtools/lib/bamtools"
    #bamtools_include="-I$with_bamtools/include -I$with_bamtools/include/bamtools"
    AC_SUBST(BAMTOOLS_CPPFLAGS, ["-I$with_bamtools/include -I$with_bamtools/include/bamtools]")
    AC_SUBST(BAMTOOLS_LDFLAGS, ["-Wl,-rpath,$with_bamtools/lib -L$with_bamtools/lib -Wl,-rpath,$with_bamtools/lib/bamtools -L$with_bamtools/lib/bamtools"])
    AC_SUBST(BAMTOOLS_LIBS, ["-lbamtools"])
fi

# Support for tcmalloc/jemalloc/hoard

# Check for the jemalloc memory allocator
AC_ARG_WITH(jemalloc, AS_HELP_STRING([--with-jemalloc=PATH],
	[specify directory containing the jemalloc library]))

# Check for the tcmalloc
AC_ARG_WITH(tcmalloc, AS_HELP_STRING([--with-tcmalloc=PATH],
	[specify directory containing the tcmalloc library]))

# Check for the hoard memory allocator
AC_ARG_WITH(hoard, AS_HELP_STRING([--with-hoard=PATH],
	[specify directory containing the hoard memory allocator library]))

# Set library path to user-selected allocator
if test "$with_hoard" -a -d "$with_hoard"; then
#	external_malloc_ldflags="-Wl,-rpath,$with_hoard -L$with_hoard"
    AC_SUBST(EXTERNAL_MALLOC_LDFLAGS, ["-Wl,-rpath,$with_hoard -L$with_hoard"])
    enable_hoard=1
fi

if test "$with_tcmalloc" -a -d "$with_tcmalloc"; then
#	external_malloc_ldflags="-Wl,-rpath,$with_tcmalloc -L$with_tcmalloc"
    AC_SUBST(EXTERNAL_MALLOC_LDFLAGS, ["-Wl,-rpath,$with_tcmalloc -L$with_tcmalloc"])
    enable_tcmalloc=1
fi

if test "$with_jemalloc" -a -d "$with_jemalloc"; then
#	external_malloc_ldflags="-Wl,-rpath,$with_jemalloc -L$with_jemalloc"
    AC_SUBST(EXTERNAL_MALLOC_LDFLAGS, ["-Wl,-rpath,$with_jemalloc -L$with_jemalloc"])
    enable_jemalloc=1
fi

# Check for the google sparse hash
AC_ARG_WITH(sparsehash, AS_HELP_STRING([--with-sparsehash=PATH],
	[specify directory containing the google sparsehash headers http://code.google.com/p/google-sparsehash/)]))

if test "$with_sparsehash" -a -d "$with_sparsehash"; then
#    sparsehash_include="-I$with_sparsehash/include"
    AC_SUBST(SPARSEHASH_CPPFLAGS, ["-I$with_sparsehash/include"])
fi

# Warn that multithreading is not available on macosx, since it does not implement unnamed semaphores
AC_MSG_CHECKING(for host type)
host="`uname -a | awk '{print $1}'`";
if test "$host" = Darwin;then
  AC_MSG_RESULT(warning: multi-threading is not available since OSX does not support un-named pthread semaphores.)
else
  AC_MSG_RESULT(you are not using osx so multi-threading should work.);
fi

# Check for boost libs.
BOOST_REQUIRE([1.55])

# The BOOST_IOSTREAMS macro is broken: it wrongfully picks up the first version
# it finds, even if it is an old version from a standard location, and the user
# specifically requested a new one.
#
# If boost path is specifed (either as BOOST_ROOT or --with-boost), we define
# BOOST_IOSTREAMS_LDFLAGS to include that path, and check that the
# boost_iostreams library can indeed be found at that specific location.
AC_SUBST(BOOST_IOSTREAMS_LIBS, ["-lboost_iostreams"])
AC_SUBST(BOOST_PROGRAM_OPTIONS_LIBS, ["-lboost_program_options"])
if test "$with_boost"; then
    AC_SUBST(BOOST_LDFLAGS, ["-L$with_boost/lib -Wl,--rpath -Wl,$with_boost/lib"])
    boost_iostreams_libs_test="$with_boost/lib/libboost_iostreams.so"
    boost_program_options_libs_test="$with_boost/lib/libboost_program_options.so"
else
    boost_iostreams_libs_test="-lboost_iostreams"
    boost_program_options_libs_test="-lboost_program_options"
fi
# check for BOOST_IOSTREAMS
AC_MSG_CHECKING([for boost_iostreams lib])
cat >conftest.$ac_ext <<_ACEOF
#include <boost/iostreams/device/file_descriptor.hpp>
int main ()
{
    boost::iostreams::file_descriptor fd;
    fd.close();
    return 0;
}
_ACEOF
save_CPPFLAGS=$CPPFLAGS
save_LDFLAGS=$LDFLAGS
save_LIBS=$LIBS
CPPFLAGS="$CPPFLAGS $BOOST_CPPFLAGS"
LDFLAGS="$LDFLAGS $BOOST_LDFLAGS"
LIBS="$LIBS $boost_iostreams_libs_test"
if ac_fn_cxx_try_link "$LINENO"; then
   AC_MSG_RESULT([$boost_iostreams_libs_test])
else
   AC_MSG_ERROR([boost_iostreams not found])
fi
CPPFLAGS=$save_CPPFLAGS
LDFLAGS=$save_LDFLAGS
LIBS=$save_LIBS
rm -f conftest.$ac_ext
# check for BOOST_PROGRAM_OPTIONS
AC_MSG_CHECKING([for boost_program_options lib])
cat >conftest.$ac_ext <<_ACEOF
#include <boost/program_options.hpp>
int main ()
{
    boost::program_options::options_description d("test");
    return 0;
}
_ACEOF
save_CPPFLAGS=$CPPFLAGS
save_LDFLAGS=$LDFLAGS
save_LIBS=$LIBS
CPPFLAGS="$CPPFLAGS $BOOST_CPPFLAGS"
LDFLAGS="$LDFLAGS $BOOST_LDFLAGS"
LIBS="$LIBS $boost_program_options_libs_test"
if ac_fn_cxx_try_link "$LINENO"; then
   AC_MSG_RESULT([$boost_program_options_libs_test])
else
   AC_MSG_ERROR([boost_iostreams not found])
fi
CPPFLAGS=$save_CPPFLAGS
LDFLAGS=$save_LDFLAGS
LIBS=$save_LIBS
rm -f conftest.$ac_ext

# To run tests, we need to add things to CPPFLAGS and LDFLAGS; first we save their original values.
save_CPPFLAGS=$CPPFLAGS
CPPFLAGS="$CPPFLAGS $BAMTOOLS_CPPFLAGS $SPARSEHASH_CPPFLAGS"
save_LDFLAGS=$LDFLAGS
LDFLAGS="$LDFLAGS $EXTERNAL_MALLOC_LDFLAGS $BAMTOOLS_LDFLAGS"

#
if test "$enable_hoard"; then
    AC_CHECK_LIB([hoard], [malloc])
fi

if test "$enable_tcmalloc"; then
    AC_CHECK_LIB([tcmalloc], [malloc])
fi

if test "$enable_jemalloc"; then
    AC_CHECK_LIB([jemalloc], [malloc])
fi


# Ensure the sparse hash is available
AC_CHECK_HEADERS([google/sparse_hash_set],,[AC_MSG_ERROR([google sparse hash library is required: http://code.google.com/p/google-sparsehash/])])

# Check for hash table headers
AC_CHECK_HEADERS([ \
	google/sparse_hash_set google/sparse_hash_map \
	unordered_map tr1/unordered_map ext/hash_map \
	unordered_set tr1/unordered_set ext/hash_set \
])

# Make sure the bamtools headers can be found
AC_CHECK_HEADERS([api/BamReader.h],,[AC_MSG_ERROR([The bamtools library must be installed (http://github.com/pezmaster31/bamtools). You can specify its path with the --with-bamtools=PATH option])])

# restore CPPFLAGS and LDFLAGS
CPPFLAGS=$save_CPPFLAGS
LDFLAGS=$save_LDFLAGS

# Set compiler flags.
AC_SUBST(AM_CFLAGS, ["-Wall -Wextra -pedantic $OPENMP_CXXFLAGS"])
AC_SUBST(AM_CXXFLAGS, ["-Wall -Wextra -pedantic $OPENMP_CXXFLAGS"])
AC_SUBST(AM_LDFLAGS, ["-Wl,--rpath=\${LD_RUN_PATH} $EXTERNAL_MALLOC_LDFLAGS"])
AC_SUBST(AM_CPPFLAGS, ['$(BOOST_CPPFLAGS) $(SPARSEHASH_CPPFLAGS) $(BAMTOOLS_CPPFLAGS) \
    -I$(top_srcdir)/Thirdparty \
    -I$(top_srcdir)/Util \
    -I$(top_srcdir)/SQG \
    -I$(top_srcdir)/Bigraph \
    -I$(top_srcdir)/Algorithm \
    -I$(top_srcdir)/StringGraph \
    -I$(top_srcdir)/Concurrency \
    -I$(top_srcdir)/SuffixTools \
    -I$(top_srcdir)/Scaffold \
    -I$(top_srcdir)/GraphDiff'])

AC_CONFIG_FILES([Makefile
    bin/Makefile
    Thirdparty/Makefile
    Util/Makefile
    SQG/Makefile
    Bigraph/Makefile
    Algorithm/Makefile
    StringGraph/Makefile
    Concurrency/Makefile
    SuffixTools/Makefile
    GraphDiff/Makefile
    Scaffold/Makefile
    SGA/Makefile
    Dev/Makefile])

AC_OUTPUT
