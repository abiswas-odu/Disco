#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT([DISCO], [1.0.2], [biswasa@ornl.gov])
#AC_CONFIG_MACRO_DIR([m4])
m4_include([m4/ax_prog_cxx_mpi.m4])
m4_include([m4/ax_check_zlib.m4])
m4_include([m4/ax_mpi.m4])

systemuname=`uname`

DISCOVERSION=1.0.2

# Checks for programs.
AC_PROG_AWK
AC_PROG_CXX
AC_PROG_CC

# check for gcc version if using gcc
if test x$GCC = xyes; then
  gccversion=`$CC -dumpversion`
  gxxversion=`$CXX -dumpversion`
  if test "x${gccversion}" != "x${gxxversion}"; then
    echo "gcc is version ${gccversion} while g++ is version ${gxxversion}"
    echo "This is not a sane environment, both versions should be the same (>=4.6.1). Stopping here."
    exit 10
  fi
  GCC_VERSION_MAJOR=$(echo $gccversion | cut -d'.' -f1)
  GCC_VERSION_MINOR=$(echo $gccversion | cut -d'.' -f2)
  gccgoodversion="no"
  if test ${GCC_VERSION_MAJOR} -ge 4; then
    if test ${GCC_VERSION_MAJOR} -ge 5; then
      gccgoodversion="yes"
    elif test ${GCC_VERSION_MINOR} -ge 8; then
      gccgoodversion="yes"
    fi
  fi
  if test "${gccgoodversion}" = "no"; then
    echo
    echo "When using gcc, you need to have gcc >= 4.8.x. You have only" $gccversion", aborting."
    exit 10
  fi
fi

AC_PROG_MAKE_SET

AC_LANG_PUSH([C++])


# Checks for libraries.

#AX_PROG_CXX_MPI([],has_mpi="yes",has_mpi="no")

AX_MPI(has_mpi="yes",has_mpi="no")

AX_CHECK_ZLIB(has_zlib="yes",has_zlib="no")

#if test "${HAVE_MPI}" = "1"; then
#   AC_MSG_WARN([No MPI compiler found, won't use MPI.])
#   has_mpi="yes"
#   mpi_compiler="${MPICXX}"
#fi

#AC_ARG_WITH(mpi, [AS_HELP_STRING([--with-mpi],
#      [compile with MPI (parallelization) support. If none is found,
#      MPI is not used. Default: auto])
#  ],,[with_mpi=auto])

#AX_PROG_CXX_MPI([test x"$with_mpi" != xno],[use_mpi=yes],[
#    use_mpi=no
#    if test x"$with_mpi" = xyes; then
#      AC_MSG_FAILURE([MPI compiler requested, but couldn't use MPI.])
#    else
#      AC_MSG_WARN([No MPI compiler found, won't use MPI.])
#    fi
#  ])




# Checks for header files.
AC_CHECK_HEADERS([inttypes.h limits.h memory.h stdint.h stdlib.h string.h sys/time.h unistd.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_CHECK_HEADER_STDBOOL
AC_C_INLINE
AC_TYPE_PID_T
AC_TYPE_SIZE_T
AC_TYPE_UINT32_T
AC_TYPE_UINT64_T
AC_TYPE_UINT8_T
AC_CHECK_TYPE(int64_t,has_int64_t="yes",has_int64_t="no")

# Checks for library functions.
AC_FUNC_MALLOC
AC_FUNC_REALLOC
AC_CHECK_FUNCS([gettimeofday localtime_r memset pow sqrt strchr])


echo "Building version ................................. ${DISCOVERSION}"

echo
echo "We are building on ............................... ${systemuname}"
echo "CPU supports 64 bit? ............................. ${has_int64_t}"
echo "Supports MPI library? ............................ ${has_mpi}"
echo "Supports ZLIB library? ........................... ${has_zlib}"
echo "Compiler ......................................... ${CC}"

if test ${has_int64_t} = "no"; then
	echo
	AC_MSG_FAILURE([ERROR! 64bit Integer support required for proper installation.])
fi

if test x$GCC != xyes; then
  possibleproblems="yes"
  echo
  echo "Compiler is not GCC. This *MIGHT* be a problem! (more compilation tests"
  echo "needed, please report success or failure)"
  echo
fi

if test "${has_mpi}" = "yes"; then
   echo "MPI compiler? .................................... ${MPICXX}"
   echo
   echo "Can Disco single computer mode be built?.......... yes"
   echo "Can Disco distributed computation mode?........... yes"
   echo "Can Disco distributed memory mode?................ yes"
   echo
   if test "${has_zlib}" = "yes"; then
	   echo "Type 'make all READGZ=1' to make the all the binaries with gzip compressed FASTA/FASTQ (*.gz) file reading capability."
   else
       echo "Type 'make all' to make the all the binaries."
	   echo "WARNING! Input of compressed FASTA/FASTQ (*.gz) is not possible as ZLIB library is missing."
   fi
   echo
else
   echo
   echo "Can Disco single computer mode be built?.......... yes"
   echo "Can Disco distributed computation mode?........... no"
   echo "Can Disco distributed memory mode?................ no"
   echo
   if test "${has_zlib}" = "yes"; then
       echo "Type 'make openmp READGZ=1' to make the all the binaries with gzip compressed FASTA/FASTQ (*.gz) file reading capability."
   else
       echo "Type 'make openmp' to make the all the binaries." 
       echo "WARNING! Input of compressed FASTA/FASTQ (*.gz) is not possible as ZLIB library is missing."
   fi
fi
AC_OUTPUT
