# ----------------------------------------------------------------------
# MNDO installation and configuration script.
#
# First version by Ernst-Udo Wallenborn in 1999 at Zurich University.
# Final version by Axel Koslowski in 2003-2013 at MPI Muelheim.
#
# Process this file with autoconf to produce a configure script.
# The first version of this script has been rewritten because it
# did not work with autoconf 2.53 any more. The current script
# requires a sufficiently recent version of autoconf.
# ----------------------------------------------------------------------

# Autoconf initialization.
AC_PREREQ(2.57)
AC_INIT(MNDO2020, 8.0)

# Get build, host, and target information.
AC_CANONICAL_TARGET

# Automake initialization.
AM_INIT_AUTOMAKE(no-dependencies no-exeext)

# Check for ranlib.
AC_PROG_RANLIB

# ----------------------------------------------------------------------
# The following section contains the list of the optional subdirectories
# and the corresponding option strings.
# ----------------------------------------------------------------------

# List of optional subdirectories.
SUBDIRS='PS GUGACI PSEX HDLCOPT CGDMS SPARSE SCOSMO DFTB VB ES'

# List of corresponding options.
OPTIONS='analytic guga surface hdlc lindir sparse scosmo dftb vb es'

# List of options that are enabled by default.
DEFAULT_OPTIONS='analytic,guga,surface,hdlc,lindir,sparse,vb,es'

# ----------------------------------------------------------------------
# End of the optional subdirectory section.
# ----------------------------------------------------------------------

# Initialize variables OPTION_PS, OPTION_GUGACI etc. to obtain the
# corresponding option string from the subdirectory name.

ni=0
for i in $SUBDIRS; do
  let ni=ni+1
  nj=0
  for j in $OPTIONS; do
    let nj=nj+1
    if test "$ni" = "$nj"; then
      eval OPTION_$i=$j
    fi
  done
done

# Print subdirectory names and corresponding option strings.
if test -n "$SCRIPTDEBUG"; then
  for i in $SUBDIRS; do
    eval echo $i \$OPTION_$i
  done
fi

# Default libraries. The least dependent ones go to the right.
# More libraries will be appended as needed by the options.
MNDOLIBS='-lmndo'

# The executable depends on the following archives.
# More dependencies will be appended as needed by the options.
ARCHIVES='libmndo.a'

# Flags concerning conditional compilation.
for opt in $SUBDIRS; do
  eval NEED_$opt='no'
  eval NEED_NO$opt='no'
done
NEED_BLAS='no'
NEED_BLAS3='no'
NEED_LAPACK='no'
NEED_MACHINE='no'
NEED_DEFMEM='no'
NEED_DMR='no'
NEED_MEMX='no'
NEED_GEMMBL='no'
NEED_NOPGPLOT='no'
HAVE_CXX='no'
HAVE_CUDA='no'
HAVE_CSRMM='no'

# Non-standard files to delete on 'make clean'.
CLEANFILES='*.mod *.d *.stb *.a'


###############################################################################
# Check library, parallelizing, debugging, and profiling options.
# See './configure --help' for full description.
# Overview:
#  --enable-subdirs=ARG     Enable code in selected optional subdirectories.
#  --disable-subdirs        Disable code in all optional subdirectories.
#  --disable-analytic       Disable analytic derivatives and NMR chemical shifts.
#  --disable-guga           Disable GUGA-CI code.
#  --disable-surface        Disable surface crossing/hopping code.
#  --disable-hdlc           Disable the HDLC optimizer.
#  --disable-lindir         Disable full-matrix density matrix search.
#  --disable-sparse         Disable sparse density matrix search.
#  --disable-scosmo         Disable Smooth COSMO solvation model.
#  --disable-dftb           Disable the DFTB module.
#  --disable-vb             Disable semiempirical valence-bond module.
#  --enable-blas=system     Use system blas library.
#  --enable-math=system     Use system blas and lapack libraries.
#  --enable-blas=mkl        Use blas routines from Intel's Math Kernel Library.
#  --enable-math=mkl        Use Intel's MKL blas and lapack routines.
#  --enable-blas=acml       Use blas routines from AMD Core Math Library.
#  --enable-math=acml       Use ACML blas and lapack routines.
#  --enable-blas=atlas      Use ATLAS blas library.
#  --enable-math=atlas      Use ATLAS blas and lapack libraries.
#  --enable-blas=generic    Use generic blas library (-lblas).
#  --enable-math=generic    Use generic math libraries (-lblas and -llapack).
#  --enable-lapack=private  Use private lapack library.
#  --enable-math=private    Use private blas and lapack libraries.
#  --enable-math=gemmbl     Include blas level 3 routines modified by Per Ling.
#  --enable-blas=ARG        Use blas library -lARG.
#  --enable-lapack=ARG      Use lapack library -lARG.
#  --enable-math=ARG        Use blas and lapack library -lARG.
#  --enable-pgplot=ARG      Use pgplot library -lARG.
#  --disable-pgplot         Do not use pgplot library.
#  --enable-openmp[=ARG]    Parallelize using OpenMP.
#  --enable-debug[=ARG]     Set debugging flag in selected subdirectories.
#  --enable-bound[=ARG]     Enable boundary checking in selected subdirs.
#  --enable-diag[=ARG]      Enable full runtime diagnosis in selected subdirs.
#  --enable-prof[=ARG]      Turn on profiling.
#
# MNDO default settings:
#  --enable-subdirs=analytic,guga,surface,hdlc,lindir,sparse,dftb,vb,es
#  --enable-math=system
#  --enable-pgplot=pgplot
#  --disable-openmp
#  --disable-debug
#  --disable-bound
#  --disable-diag
#  --disable-prof
###############################################################################

# --enable-subdirs
AC_ARG_ENABLE(subdirs,dnl
[[  --enable-subdirs=ARG    Enable code in selected optional subdirectories.
                          ARG is a comma-separated list of identifiers for the
                          subdirectories, namely analytic, guga, surface, hdlc,
                          lindir, sparse, scosmo, dftb, vb and es. By default
                          all optional subdirectories which are present in the
                          distribution except SCOSMO and DFTB are enabled.
                          '--disable-subdirs' disables all subdirectories.
                          'make dist' includes only enabled subdirectories
                          in the distribution.]],[dnl
if test "$enable_subdirs" = "yes"; then
  # --enable-subdirs has been specified without argument.
  AC_MSG_RESULT([
ERROR: --enable-subdirs requires an argument.
])
  exit
elif test "$enable_subdirs" != "no"; then
  TMP='no'
  for i in $SUBDIRS; do
    if test -d "$i"  &&  echo "$enable_subdirs" | eval grep -i "\$OPTION_$i" > /dev/null; then
      if test "$TMP" = "no"; then
        TMP=$(eval echo \$OPTION_$i)
      else
        TMP=$(eval echo \$TMP','\$OPTION_$i)
      fi
    fi
  done
  enable_subdirs=$TMP
fi],[dnl
  # Default if --enable-subdirs has not been specified.
  enable_subdirs='no'
  for i in $SUBDIRS; do
    if test -d "$i"  &&  echo "$DEFAULT_OPTIONS" | eval grep -i "\$OPTION_$i" > /dev/null; then
      if test "$enable_subdirs" = "no"; then
        enable_subdirs=$(eval echo \$OPTION_$i)
      else
        enable_subdirs=$(eval echo \$enable_subdirs','\$OPTION_$i)
      fi
    fi
  done
])

# Define help text macro.
AC_DEFUN([AC_HELPTEXT],[dnl
AC_ARG_ENABLE($1,$2,[dnl
if test "$enable_$1" = "yes"; then
  # --enable-$1 has been specified without argument.
  # Add keyword $1 to the string of subdirs if not yet present.
  if echo "$enable_subdirs" | grep "$1" > /dev/null; then
    # Keyword $1 already present in the string. Do nothing.
    :
  elif test "$enable_subdirs" = "no"; then
    # No keywords present in the string so far.
    # Replace the current value.
    enable_subdirs='$1'
  else
    # Append keyword $1 to the current string.
    enable_subdirs=$enable_subdirs',$1'
  fi
elif test "$enable_$1" = "no"; then
  # Option --disable-$1 has been specified.
  # Rebuild option string.
  TMP='no'
  for i in $OPTIONS; do
    if test "$i" != "$1"  &&  echo "$enable_subdirs" | grep "$i" > /dev/null; then
      if test "$TMP" = "no"; then
        TMP=$i
      else
        TMP=$TMP','$i
      fi
    fi
  done
  enable_subdirs=$TMP
else
  # Option --enable-$1 has been specified with an argument.
  AC_MSG_RESULT([
ERROR: --enable-$1 cannot have an argument.
])
  exit
fi])])

# ----------------------------------------------------------------------
# The following section contains the help texts for ./configure --help.
# The correct number of spaces at the beginning of the lines is crucial
# for a clean alignment of the output.
# ----------------------------------------------------------------------

# --enable-analytic
AC_HELPTEXT(analytic,dnl
[[  --disable-analytic      Disable analytic derivatives and NMR chemical
                          shifts, and exclude the PS directory from the
                          distribution when invoking 'make dist'.]])

# --enable-guga
AC_HELPTEXT(guga,dnl
[[  --disable-guga          Disable GUGA-CI code and exclude the GUGACI and
                          GTO directories from the distribution when invoking
                          'make dist'.]])

# --enable-es
AC_HELPTEXT(es,dnl
[[  --disable-es            Disable excited-state module and exclude
                          the ES directory from the distribution when
                          invoking 'make dist'.]])

# --enable-surface
AC_HELPTEXT(surface,dnl
[[  --disable-surface       Disable surface crossing and surface hopping module
                          and exclude the PSEX directory from the distribution
                          when invoking 'make dist'.]])

# --enable-hdlc
AC_HELPTEXT(hdlc,dnl
[[  --disable-hdlc          Disable the HDLC optimizer and exclude the HDLCOPT
                          directory from the distribution when invoking
                          'make dist'.]])

# --enable-lindir
AC_HELPTEXT(lindir,dnl
[[  --disable-lindir        Disable full-matrix density matrix search and
                          exclude the subdirectory CGDMS from the distribution
                          when invoking 'make dist'.]])

# --enable-sparse
AC_HELPTEXT(sparse,dnl
[[  --disable-sparse        Disable sparse density matrix search and exclude
                          the subdirectory SPARSE from the distribution when
                          invoking 'make dist'.]])

# --enable-scosmo
AC_HELPTEXT(scosmo,dnl
[[  --enable-scosmo         Enable Smooth COSMO solvation model and include the
                          subdirectory SCOSMO in the distribution when invoking
                          'make dist'. *** INCOMPLETE IMPLEMENTATION ***]])

# --enable-dftb
AC_HELPTEXT(dftb,dnl
[[  --enable-dftb           Enable DFTB module and include the subdirectory
                          DFTB in the distribution when invoking 'make dist'.
                          *** NOT PRESENT IN THE OPEN SOURCE DISTRIBUTION ***]])

# --enable-vb
AC_HELPTEXT(vb,dnl
[[  --disable-vb            Disable semiempirical valence-bond module and
                          exclude the VB directory from the distribution
                          when invoking 'make dist'.]])

# ----------------------------------------------------------------------
# End of help text section.
# ----------------------------------------------------------------------

# Print initial list of options.
if test -n "$SCRIPTDEBUG"; then
  echo $enable_subdirs
fi

# Print list of enabled options.
AC_MSG_CHECKING([which optional subdirs to enable])
AC_MSG_RESULT([$enable_subdirs])

# Assemble list of enabled subdirectory names.
ENABLED_SUBDIRS='none'
for i in $SUBDIRS; do
  if echo "$enable_subdirs" | eval grep "\$OPTION_$i" > /dev/null; then
    if test "$ENABLED_SUBDIRS" = "none"; then
      ENABLED_SUBDIRS=$i
    else
      ENABLED_SUBDIRS=$ENABLED_SUBDIRS' '$i
    fi
  fi
done

# Set flags according to list of enabled subdirs.
# Please note that the order of the libraries in MNDOLIBS is crucial.
# The least dependent libraries go to the right.
# There is no interdependence between the optional subdirectories
# so that their relative order does not matter.
for i in $SUBDIRS; do
  if echo "$enable_subdirs" | eval grep "\$OPTION_$i" > /dev/null; then
    eval "NEED_$i='yes'"
    MNDOLIBS=$(eval echo \$MNDOLIBS -l\$OPTION_$i)
    ARCHIVES=$(eval echo \$ARCHIVES lib\$OPTION_$i.a)
  else
    eval "NEED_NO$i='yes'"
    NEED_MACHINE='yes'
  fi
done

# --enable-math
AC_ARG_ENABLE(math,dnl
[[  --enable-math=ARG       Specify which math library (with blas and lapack
                          routines) shall be used. If the argument equals
                          'system' the configure script will look for a
                          vendor-specific library. This is the default in
                          MNDO. On an Intel-compatible Linux system, the
                          following libraries are searched for by default:
                          Intel's Math Kernel Library (MKL), the AMD Core
                          Math Library (ACML), the ATLAS blas library (from
                          http://math-atlas.sourceforge.net) together with a
                          patched lapack library named -llamndo), and generic
                          blas and lapack libraries. Use of a specific library
                          may be requested by one of the arguments mkl, acml,
                          atlas, or generic. If no optimized library is
                          found, a subset of the blas and lapack routines
                          included in the MNDO distribution is used which can
                          be obtained from http://www.netlib.org/lapack/double.
                          These routines may be requested directly by the
                          argument 'private'. If the argument is 'gemmbl' the
                          original netlib blas level 3 routines are replaced by
                          versions modified by Per Ling. For other values of ARG
                          the executable will be linked against -lARG. If some
                          of the required library files are located in non-
                          standard directories the paths must be specified on
                          the command line using LDFLAGS='-L<dir1> -L<dir2> ...'.
                          To specify individual names for the blas and lapack
                          libraries use the --enable-blas and --enable-lapack
                          options. It is not possible to leave out the argument
                          or to disable the math library.]],[dnl
if test "$enable_math" = "yes"; then
  AC_MSG_RESULT([
ERROR: --enable-math requires an argument.
])
  exit
elif test "$enable_math" = "no"; then
  AC_MSG_RESULT([
ERROR: The math library cannot be disabled.
])
  exit
elif test "$enable_math" = "system"; then
  BLAS=''
  LAPACK=''
elif test "$enable_math" = "mkl"; then
  BLAS='mkl_intel'
  LAPACK='mkl_intel'
  F2C='-lmkl_sequential -lmkl_core'
elif test "$enable_math" = "acml"; then
  BLAS='acml'
  LAPACK='acml'
  F2C='-lacml_mv'
elif test "$enable_math" = "atlas"; then
  BLAS='f77blas'
  LAPACK='lamndo'
  if test -z "$F2C"; then
    F2C='-lcblas -latlas'
  else
    F2C='-lcblas -latlas '$F2C
  fi
elif test "$enable_math" = "generic"; then
  BLAS='blas'
  LAPACK='lapack'
elif test "$enable_math" = "private"; then
  BLAS='-lbl2011'
  LAPACK='-lla350'
  NEED_BLAS='yes'
  NEED_BLAS3='yes'
  NEED_LAPACK='yes'
elif test "$enable_math" = "gemmbl"; then
  BLAS='-lbl2011'
  LAPACK='-lla350'
  NEED_BLAS='yes'
  NEED_GEMMBL='yes'
  NEED_LAPACK='yes'
  NEED_MACHINE='yes'
else
  BLAS=${enable_math#-l}
  LAPACK=${enable_math#-l}
fi],[dnl
  # Default if --enable-math is not specified.
  BLAS=''
  LAPACK=''
  enable_math='system'
])

# --enable-blas
AC_ARG_ENABLE(blas,dnl
[[  --enable-blas=ARG       Link the executable against the blas library -lARG.
                          Possible values for ARG with a special meaning are
                          system, mkl, acml, gsl, atlas, generic, private, and
                          gemmbl as explained above. By default the lapack
                          library from the MNDO distribution is used if a
                          specific blas library but no lapack library has been
                          requested by the user.]],[dnl
if test "$enable_blas" = "yes"; then
  AC_MSG_RESULT([
ERROR: --enable-blas requires an argument.
  exit
])
elif test "$enable_blas" = "no"; then
  AC_MSG_RESULT([
ERROR: The blas library cannot be disabled.
])
  exit
fi
if test "$enable_math" != "system"; then
  AC_MSG_RESULT([
ERROR: --enable-math and --enable-blas are mutually exclusive.
])
  exit
fi
if test "$enable_blas" = "system"; then
  BLAS=''
elif test "$enable_blas" = "mkl"; then
  BLAS='mkl_intel'
  LAPACK='-lla350'
  NEED_LAPACK='yes'
elif test "$enable_blas" = "acml"; then
  BLAS='acml'
  LAPACK='-lla350'
  F2C='-lacml_mv'
  NEED_LAPACK='yes'
elif test "$enable_blas" = "atlas"; then
  BLAS='f77blas'
  LAPACK='-lla350'
  F2C='-lcblas -latlas'
  NEED_LAPACK='yes'
elif test "$enable_blas" = "generic"; then
  BLAS='blas'
  LAPACK='-lla350'
  NEED_LAPACK='yes'
elif test "$enable_blas" = "private"; then
  BLAS='-lbl2011'
  LAPACK='-lla350'
  NEED_BLAS='yes'
  NEED_BLAS3='yes'
  NEED_LAPACK='yes'
elif test "$enable_blas" = "gemmbl"; then
  BLAS='-lbl2011'
  LAPACK='-lla350'
  NEED_BLAS='yes'
  NEED_GEMMBL='yes'
  NEED_LAPACK='yes'
  NEED_MACHINE='yes'
else
  BLAS=${enable_blas#-l}
  LAPACK='-lla350'
  NEED_LAPACK='yes'
fi],[dnl
  # Default if --enable-blas is not specified.
  enable_blas=${enable_math#-l}
])

# --enable-lapack
AC_ARG_ENABLE(lapack,dnl
[[  --enable-lapack=ARG     Link the executable against the lapack library
                          -lARG. Possible values for ARG with a special
                          meaning are 'system' and 'private' as explained for
                          --enable-math. A system blas library may be combined
                          with the lapack library from the MNDO distribution
                          by specifying --enable-lapack=private.]],[dnl
if test "$enable_lapack" = "yes"; then
  AC_MSG_RESULT([
ERROR: --enable-lapack requires an argument.
])
  exit
elif test "$enable_lapack" = "no"; then
  AC_MSG_RESULT([
ERROR: The lapack library cannot be disabled.
])
  exit
fi
if test "$enable_math" != "system"; then
  AC_MSG_RESULT([
ERROR: --enable-math and --enable-lapack are mutually exclusive.
])
  exit
fi
if test "$enable_lapack" = "system"; then
  LAPACK=''
  NEED_LAPACK='no'
elif test "$enable_lapack" = "private"; then
  LAPACK='-lla350'
  NEED_LAPACK='yes'
else
  LAPACK=${enable_lapack#-l}
  NEED_LAPACK='no'
fi],[dnl
  # This branch is executed if --enable-lapack is not specified.
  enable_lapack=${enable_math#-l}
])

AC_MSG_CHECKING([which blas library to use])
if test -z "$BLAS"; then
  AC_MSG_RESULT([$enable_blas])
elif echo "$BLAS" | grep -- "-l" > /dev/null; then
  AC_MSG_RESULT([$BLAS])
else
  AC_MSG_RESULT([-l$BLAS])
fi

AC_MSG_CHECKING([which lapack library to use])
if test -z "$LAPACK"; then
  AC_MSG_RESULT([$enable_lapack])
elif echo "$LAPACK" | grep -- "-l" > /dev/null; then
  AC_MSG_RESULT([$LAPACK])
else
  AC_MSG_RESULT([-l$LAPACK])
fi

# --disable-cuda
AC_ARG_ENABLE(cuda,dnl
[[  --disable-cuda          Do not compile any CUDA files (*.cu) and do not
                          search for the NVidia CUDA compiler nvcc. Instead,
                          compile C++ files (*.cxx) using the GNU compiler.]],[dnl
if test "$enable_cuda" = "yes"; then
  CXXCOMP=$CXX' nvcc g++'
elif test "$enable_cuda" = "no"; then
  CXXCOMP=$CXX' g++'
else
  CXXCOMP="${enable_cuda} $CXX"
fi],[dnl
  # Default if --enable-cuda is not specified.
  CXXCOMP=$CXX' nvcc g++'
])

# --disable-cxx
AC_ARG_ENABLE(cxx,dnl
[[  --disable-cxx           Do not compile any C++ or CUDA files (*.cxx or *.cu)
                          and do not search for a corresponding compiler.]],[dnl
if test "$enable_cxx" = "yes"; then
  CXXCOMP=$CXX' g++'
elif test "$enable_cxx" = "no"; then
  CXXCOMP=''
  CXX=''
else
  CXXCOMP="${enable_cxx} $CXX"
fi])

# --disable-pgplot
AC_ARG_ENABLE(pgplot,dnl
[[  --disable-pgplot        Do not search for the pgplot library. By default,
                          the pgplot library is included as -lpgplot if
                          available. If your library has a different name,
                          use --enable-pgplot=mypgplot (without the -l flag).
                          If the library file is not in a standard directory,
                          you must specify its path on the command line using
                          LDFLAGS='-L<path>'. The pgplot library is used in the
                          GUGA-CI code to generate postscript files with plots of
                          the Shavitt graph and of loops. This feature is present
                          mainly for debugging purposes and the corresponding
                          code will not be called unless ciplot.gt.0 is specified
                          in the input.]],[dnl
if test "$enable_pgplot" = "yes"; then
  PGPLOT='pgplot'
elif test "$enable_pgplot" = "no"; then
  PGPLOT=''
else
  PGPLOT=${enable_pgplot#-l}
fi],[dnl
  # Default if --enable-pgplot is not specified.
  PGPLOT='pgplot'
])

AC_MSG_CHECKING([which pgplot library to use])
if test -z "$PGPLOT"; then
  AC_MSG_RESULT([none])
else
  AC_MSG_RESULT([-l$PGPLOT])
fi

# --enable-openmp
AC_ARG_ENABLE(openmp,dnl
[[  --enable-openmp[=ARG]   Parallelize using OpenMP. The use of a parallel
                          system math library by default is implied.]],[dnl
if test "$enable_openmp" != "yes" && test "$enable_openmp" != "no"; then
  OPENMP=$enable_openmp
  enable_openmp='yes'
fi],[dnl
  # Default if --enable-openmp is not specified.
  enable_openmp='no'
])

AC_MSG_CHECKING([whether to parallelize using OpenMP])
AC_MSG_RESULT([$enable_openmp])

# --enable-bound
AC_ARG_ENABLE(bound,dnl
[[  --enable-bound[=ARG]    Enable boundary checking for selected subdirecto-
                          ries. ARG is a comma-separated list of identifiers
                          for the subdirectories, namely mndo (for the main
                          program and the 2020 subdirectory), analytic, guga,
                          surface, hdlc, lindir, sparse, scosmo, dftb, vb, es,
                          blas, lapack, or machine. If the option is specified
                          without argument the list is mndo,analytic,guga,
                          surface,hdlc,lindir,sparse,scosmo,dftb,vb es. If ARG
                          equals 'all', boundary checking is turned on
                          throughout the program. If --enable-debug is not
                          given the debugging flag is implied. The optimization
                          level is reduced if necessary.]],[dnl
if test "$enable_bound" = "yes"; then
  enable_bound='mndo'
  for i in $OPTIONS; do
    if echo "$enable_subdirs" | grep "$i" > /dev/null; then
      enable_bound=$enable_bound','$i
    fi
  done
fi],[dnl
  # Default if --enable-bound is not specified.
  enable_bound='no'
])

AC_MSG_CHECKING([whether to use boundary checking])
AC_MSG_RESULT([$enable_bound])

# --enable-diag
AC_ARG_ENABLE(diag,dnl
[[  --enable-diag[=ARG]     Enable full run-time diagnostics for selected
                          subdirectories (see --enable-bound) if supported
                          by the compiler. Otherwise boundary checking is
                          activated. The debugging flag is implied by
                          default. The optimization level is reduced if
                          necessary. Using this option is not recommended
                          for subdirectories other than 'guga'.]],[dnl
if test "$enable_diag" = "yes"; then
  enable_diag='mndo'
  for i in $OPTIONS; do
    if echo "$enable_subdirs" | grep "$i" > /dev/null; then
      enable_diag=$enable_diag','$i
    fi
  done
fi],[dnl
  # Default if --enable-diag is not specified.
  enable_diag='no'
])

AC_MSG_CHECKING([whether to use full run-time diagnostics])
AC_MSG_RESULT([$enable_diag])

# --enable-debug
AC_ARG_ENABLE(debug,dnl
[[  --enable-debug[=ARG]    Include debugging symbols in the executable for
                          selected subdirectories (see --enable-bound).
                          The optimization level is reduced if necessary.]],[dnl
if test "$enable_debug" = "yes"; then
  enable_debug='mndo'
  for i in $OPTIONS; do
    if echo "$enable_subdirs" | grep "$i" > /dev/null; then
      enable_debug=$enable_debug','$i
    fi
  done
  enable_debug=$enable_debug',blas,lapack,machine'
fi],[dnl
# Default settings if --enable-debug is not specified.
# Join the lists of --enable-bound and --enable-diag.
enable_debug=''
for i in mndo $OPTIONS blas lapack machine; do
  if echo "$enable_bound" | grep -- "$i" > /dev/null  ||
     echo "$enable_diag"  | grep -- "$i" > /dev/null; then
    if test -n "$enable_debug"; then
      enable_debug=$enable_debug','$i
    else
      enable_debug=$i
    fi
  fi
done
if test -z "$enable_debug"; then
  enable_debug='no'
fi])

AC_MSG_CHECKING([whether to use the debugging flag])
AC_MSG_RESULT([$enable_debug])

# --enable-prof
AC_ARG_ENABLE(prof,dnl
[[  --enable-prof[=ARG]     Turn on profiling. A specific flag may be given
                          as ARG. Different flags for the compiler and the
                          linker may be specified by setting the variables
                          FPROF and LDPROF on the command line. If debug
                          symbols are needed they must be enabled explicitly
                          by the user. Please check the documentation of your
                          compiler.]],[dnl
if test "$enable_prof" != "yes" && test "$enable_prof" != "no"; then
  FPROF=$enable_prof
  LDPROF=$enable_prof
  enable_prof='yes'
fi],[dnl
  # Default settings if --enable-prof is not specified.
  enable_prof='no'
])

AC_MSG_CHECKING([whether to enable profiling])
AC_MSG_RESULT([$enable_prof])

AC_MSG_RESULT([
We will now search for a Fortran 90 and a C compiler. If several different
compilers are installed on your machine, we will select those which we
believe are best suited for the MNDO program. If you wish to force usage
of a specific compiler, run './configure F90=myf90 CC=mycc'.
])

# Initialize FIXED, FREE, FFLAGS1, LDFLAGS, and LIBS1.
FIXED=''
FREE=''
FFLAGS1=''
if test -z "$LDFLAGS"; then
  LDFLAGS='-L.'
else
  LDFLAGS=$LDFLAGS' -L.'
fi
LIBS1=''

###########################################################################
### System-specific compiler search lists and general compiler options. ###
###########################################################################

if echo "$target_os" | grep osf > /dev/null; then
  # OSF, i.e. Tru64
  echo 'Configuring a Hewlett-Packard Tru64 system...'
  FCOMP=$F90' f95'
  CCOMP=$CC' cc gcc'
  FFLAGS1='-fpe1'
elif echo "$target_os" | grep hpux > /dev/null; then
  # HP-UX/Itanium-2
  echo 'Configuring a Hewlett-Packard HP-UX/Itanium-2 system...'
  FCOMP=$F90' f95'
  CCOMP=$CC' cc gcc'
  FFLAGS1='-fpe1'
elif echo "$target_os" | grep irix > /dev/null; then
  # SGI IRIX
  echo 'Configuring an SGI IRIX system...'
  FCOMP=$F90' f90'
  CCOMP=$CC' cc gcc'
  FFLAGS1='-n32 -mips4 -r10000'
  LDFLAGS='-n32 -mips4 -r10000 '$LDFLAGS
elif echo "$target_os" | grep aix > /dev/null; then
  # IBM AIX
  echo 'Configuring an IBM AIX system...'
  echo 'WARNING: AIX support still experimental!'
  FCOMP=$F90' xlf95_r'
  CCOMP=$CC' xlc gcc'
  FIXED='-qfixed '
  FREE='-qsuffix=f=f90 '
elif echo "$target" | grep i686 > /dev/null; then
  # 32-bit Intel-compatible system
  echo 'Configuring a 32-bit Intel-compatible system...'
  AC_CHECK_PROGS([FCOMP],[$F90 ifort pathf95 pgf90])
  if test "$FCOMP" = "pathf95"; then
    # General options for the PathScale Fortran compiler (pathf95).
    LDFLAGS='-static '$LDFLAGS
  elif test "$FCOMP" = "ifort"; then
    # General options for the Intel Fortran compiler 9.1 (ifort).
    # LDFLAGS='-static '$LDFLAGS
    # LIBS1='-liomp5 -lpthread'
    LIBS1=''
  elif test "$FCOMP" = "pgf90"; then
    # General options for the Portland Group Fortran compiler 5.2 (pgf90).
    LDFLAGS='-Bstatic '$LDFLAGS
  fi;
  CCOMP=$CC' gcc icc pgcc'
  # TO DO: Select library path depending on the compiler.
  # WARNING: Support of 32-bit ACML library still experimental.
  LDFLAGS=$LDFLAGS' -L/opt/acml4.4.0/ifort32/lib'
elif echo "$target" | grep x86_64 > /dev/null; then
  # AMD64 or EM64T system.
  echo 'Configuring an x86-64 system...'
  # Try to read vendor id from /proc/cpuinfo
  if test -f /proc/cpuinfo; then
    CPUVENDOR=$(grep vendor_id /proc/cpuinfo | head -1)
    CPUVENDOR=${CPUVENDOR/vendor_id/}
    CPUVENDOR=${CPUVENDOR/:/}
    CPUVENDOR=$(echo $CPUVENDOR)
    echo "CPU vendor id: $CPUVENDOR"
  fi
  # Set MATHLIB.
  MATHLIB=''
  if   test "$enable_blas" = "mkl"; then
    MATHLIB='MKL'
  elif test "$enable_blas" = "acml"; then
    MATHLIB='ACML'
  elif test "$enable_blas" = "system"; then
    if   test "$CPUVENDOR" = "GenuineIntel"; then
      MATHLIB='MKL'
    elif test "$CPUVENDOR" = "AuthenticAMD"; then
      MATHLIB='ACML'
    else
      echo "No optimized math library for this CPU."
      exit
    fi
  fi
  # Find compiler.
  AC_CHECK_PROGS([FCOMP],[$F90 ifort pgf95 pathf95])
  if test "$FCOMP" = "ifort"; then
    # General options for the Intel Fortran compiler 12.1 (ifort).
    # LDFLAGS='-static '$LDFLAGS
    # LIBS1='-liomp5 -lpthread'
    # LIBS1='-lifcore'
    # ifort finds the paths to MKL using enviroment variables,
    # so there is no need to include MKL in MATHPATHS.
    MATHPATHS=''
    if test "$MATHLIB" = "ACML" -o -z "$MATHLIB"; then
      MATHPATHS=$MATHPATHS'/opt/acml5.0.0/ifort64/lib '
      MATHPATHS=$MATHPATHS'/opt/acml4.4.0/ifort64/lib '
    fi
  elif test "${FCOMP:0:3}" = "pgf"; then
    # General options for the Portland Group Fortran compiler 10.9 (pgf95).
    # LDFLAGS='-Bstatic '$LDFLAGS
    MATHPATHS=''
    if test "$MATHLIB" = "MKL" -o -z "$MATHLIB"; then
      MATHPATHS=$MATHPATHS'/opt/intel/mkl/lib/intel64 '
      MATHPATHS=$MATHPATHS'/opt/intel/mkl/10.2.7.041/lib/em64t '
    fi
    if test "$MATHLIB" = "ACML" -o -z "$MATHLIB"; then
      MATHPATHS=$MATHPATHS'/opt/acml5.0.0/pgi64/lib '
      MATHPATHS=$MATHPATHS'/opt/acml4.4.0/pgi64/lib '
    fi
  elif test "$FCOMP" = "pathf95"; then
    # General options for the PathScale Fortran compiler (pathf95).
    LDFLAGS='-static '$LDFLAGS
    echo Pathscale compiler currently not supported.
    exit
  fi
  CCOMP=$CC' gcc icc pgcc'
  # Find installed math library.
  MATHPATH=''
  for p in $MATHPATHS; do
    if test -d "$p"; then
      if test -z "$MATHPATH"; then
        MATHPATH="$p"
      fi
    fi
  done
  if echo "$MATHPATH" | grep acml5 > /dev/null; then
    F2C=''
  fi
  if test -n "$MATHPATH"; then
    LDFLAGS=$LDFLAGS" -L$MATHPATH"
  fi
else
  # Generic system.
  echo 'Configuring a generic system...'
  FCOMP=$F90' f95 f90'
  CCOMP=$CC' gcc cc'
fi

# Concat FFLAGS1 and FFLAGS.
if test -n "$FFLAGS1"  &&  test -n "$FFLAGS"; then
  FFLAGS=$FFLAGS1' '$FFLAGS
else
  FFLAGS=$FFLAGS1$FFLAGS
fi

# Concat LIBS1 and LIBS.
if   test -n "$LIBS1"  &&  test -n "$LIBS"; then
  LIBS=' '$LIBS1' '$LIBS
elif test -n "$LIBS1"  ||  test -n "$LIBS"; then
  LIBS=' '$LIBS1$LIBS
fi

# AC_PROG_F77 and AC_PROG_CC produce different results if they appear in
# the configure script more than once. Thus we use the variables FCOMP
# and CCOMP to determine which compilers will be searched for so we get
# along with a single instance of these macros.
FFLAGS_BAK=$FFLAGS
AC_PROG_FC($FCOMP)
FFLAGS=$FFLAGS_BAK

CFLAGS_BAK=$CFLAGS
AC_PROG_CC($CCOMP)
CFLAGS=$CFLAGS_BAK

if test -n "$CXXCOMP"; then
CXXFLAGS_BAK=$CXXFLAGS
AC_PROG_CXX($CXXCOMP)
CXXFLAGS=$CXXFLAGS_BAK
fi

if test "$CXX" = "nvcc"; then
  LDFLAGS=$LDFLAGS' -L/usr/local/cuda/lib64'
  CXXLIBS=' -lcusparse -lcublas -lcudart'
  HAVE_CUDA='yes'
fi

if test -n "$CXX"; then
  CXXLIBS=$CXXLIBS' -lrt -lstdc++'
  HAVE_CXX='yes'
fi

# The F77 compiler is in fact an F90!
F77=$FC
F90=$FC

#######################################################################
### System- and compiler-specific optimization flags and libraries. ###
#######################################################################

AC_LANG(Fortran 77)
FLIBS=''
if echo "$target_os" | grep osf > /dev/null; then
  ############################
  ## DEC OSF, i.e. HP Tru64 ##
  ############################
  if test "$enable_debug" != "no"  &&  test -z "$DEBUG"; then
    # Check if -g3 is supported.
    FFLAGS_BAK=$FFLAGS
    FFLAGS='-g3'
    AC_MSG_CHECKING([whether $F90 accepts -g3])
    AC_COMPILE_IFELSE([[      PROGRAM TESTG3
      END]],[DEBUG='-g3'],[DEBUG='-g'])
    AC_MSG_RESULT([$DEBUG])
    FFLAGS=$FFLAGS_BAK
  fi
  if test -z "$CDEBUG";  then  CDEBUG='-g';   fi
  if test -z "$BOUND";   then   BOUND='-C';   fi
  if test -z "$DIAG";    then    DIAG='-C';   fi
  if test -z "$OFLAGS1"; then OFLAGS1='-O4';  fi
  if test -z "$OFLAGS2"; then
    OFLAGS2='-O5 -align records -align dcommons -align sequence -assume noaccuracy_sensitive -assume bigarrays -assume nozsize -math_library fast -tune host'
  fi
  if test -z "$OFLAGSC"; then OFLAGSC='-O2';  fi
  if test -z "$LDPROF";  then  LDPROF='-p';   fi
  if test -z "$OPENMP";  then  OPENMP='-omp'; fi
  if test "$DEBUG" = "-g3"; then
    OPT_WITH_DEBUG='yes'
  else
    OPT_WITH_DEBUG='no'
  fi
  OPT_WITH_BOUND='yes'
  OPT_WITH_DIAG='yes'
  #
  # Check user-specified blas and lapack libraries. If the blas library
  # is not found use the default libraries. If the lapack library is not
  # found use the lapack library from the MNDO distribution.
  #
  if test -n "$F2C"; then F2C=' '$F2C; fi
  if test -n "$BLAS"  &&  echo "$BLAS" | grep -v -- "-l" > /dev/null; then
    AC_CHECK_LIB([$BLAS],[dgemm],[BLAS='-l'$BLAS],[BLAS=''
      LAPACK=''; F2C=''],[$F2C])
  fi
  if test -n "$LAPACK"  &&  echo "$LAPACK" | grep -v -- "-l" > /dev/null; then
    AC_CHECK_LIB([$LAPACK],[dsyevx],[LAPACK='-l'$LAPACK],[LAPACK='-lla350'
      NEED_LAPACK='yes'],[$BLAS$F2C])
  fi
  if test "$enable_openmp" = "no"; then
    # Search for serial math library versions.
    # Search for cxml library if unsuccessful so far.
    if test -z "$BLAS"; then
      AC_CHECK_LIB([cxml],[dgemm],[BLAS='-lcxml'],[BLAS=''])
    fi
    if test -z "$LAPACK"; then
      AC_CHECK_LIB([cxml],[dsyevx],[LAPACK='-lcxml'],[LAPACK=''])
    fi
    # Search for dxml library if unsuccessful.
    if test -z "$BLAS"; then
      AC_CHECK_LIB([dxml],[dgemm],[BLAS='-ldxml'],[BLAS=''])
    fi
    if test -z "$LAPACK"; then
      AC_CHECK_LIB([dxml],[dsyevx],[LAPACK='-ldxml'],[LAPACK=''])
    fi
  else
    # Search for parallel math library versions.
    # Search for cxmlp library if unsuccessful so far.
    if test -z "$BLAS"; then
      AC_CHECK_LIB([cxmlp],[dgemm],[BLAS='-lcxmlp'],[BLAS=''])
    fi
    if test -z "$LAPACK"; then
      AC_CHECK_LIB([cxmlp],[dsyevx],[LAPACK='-lcxmlp'],[LAPACK=''])
    fi
    # Search for dxmlp library if unsuccessful.
    if test -z "$BLAS"; then
      AC_CHECK_LIB([dxmlp],[dgemm],[BLAS='-ldxmlp'],[BLAS=''])
    fi
    if test -z "$LAPACK"; then
      AC_CHECK_LIB([dxmlp],[dsyevx],[LAPACK='-ldxmlp'],[LAPACK=''])
    fi
  fi
elif echo "$target_os" | grep hpux > /dev/null; then
  #####################
  ## HP-UX/Itanium-2 ##
  #####################
  echo HP-UX not yet implemented.
  exit
elif echo "$target_os" | grep irix > /dev/null; then
  ##############
  ## SGI IRIX ##
  ##############
  if test -z "$DEBUG";   then   DEBUG='-g'; fi
  if test -z "$CDEBUG";  then  CDEBUG='-g'; fi
  if test -z "$BOUND";   then   BOUND='-C'; fi
  if test -z "$DIAG";    then    DIAG='-C'; fi
  if test -z "$OFLAGS1"; then OFLAGS1='-O2 -OPT:Olimit=0'; fi
  if test -z "$OFLAGS2"; then OFLAGS2='-O2 -OPT:Olimit=0'; fi
  if test -z "$OFLAGSC"; then OFLAGSC='-O'; fi
  OPT_WITH_DEBUG='no'
  OPT_WITH_BOUND='yes'
  OPT_WITH_DIAG='yes'
  # Check user-specified blas and lapack libraries. If the blas library
  # is not found use the default libraries. If the lapack library is not
  # found use the lapack library from the MNDO distribution.
  if test -n "$F2C"; then F2C=' '$F2C; fi
  if test -n "$BLAS"  &&  echo "$BLAS" | grep -v -- "-l" > /dev/null; then
    AC_CHECK_LIB([$BLAS],[dgemm],[BLAS='-l'$BLAS],[BLAS=''
      LAPACK=''; F2C=''],[$F2C])
  fi
  if test -n "$LAPACK"  &&  echo "$LAPACK" | grep -v -- "-l" > /dev/null; then
    AC_CHECK_LIB([$LAPACK],[dsyevx],[LAPACK='-l'$LAPACK],[LAPACK='-lla350'
      NEED_LAPACK='yes'],[$BLAS$F2C])
  fi
  # Search for complib.sgimath library if unsuccessful.
  if test -z "$BLAS"; then
    AC_CHECK_LIB([complib.sgimath],[dgemm],dnl
      [BLAS='-lcomplib.sgimath'],[BLAS=''])
  fi
  if test -z "$LAPACK"; then
    AC_CHECK_LIB([complib.sgimath],[dsyevx],dnl
      [LAPACK='-lcomplib.sgimath'],[LAPACK=''])
  fi
elif echo "$target_os" | grep aix > /dev/null; then
  #############
  ## IBM AIX ##
  #############
  if test -z "$DEBUG";   then   DEBUG='-g'; fi
  if test -z "$CDEBUG";  then  CDEBUG='-g'; fi
  if test -z "$BOUND";   then   BOUND='-C'; fi
  if test -z "$DIAG";    then    DIAG='-C'; fi
  if test -z "$OFLAGS1"; then OFLAGS1='-O'; fi
  if test -z "$OFLAGS2"; then OFLAGS2='-O'; fi
  if test -z "$OFLAGSC"; then OFLAGSC='-O'; fi
  OPT_WITH_DEBUG='yes'
  OPT_WITH_BOUND='yes'
  OPT_WITH_DIAG='yes'
  # Check user-specified blas and lapack libraries. If the blas library
  # is not found use the default libraries. If the lapack library is not
  # found use the lapack library from the MNDO distribution.
  if test -n "$F2C"; then F2C=' '$F2C; fi
  if test -n "$BLAS"  &&  echo "$BLAS" | grep -v -- "-l" > /dev/null; then
    AC_CHECK_LIB([$BLAS],[dgemm],[BLAS='-l'$BLAS],[BLAS=''
      LAPACK=''; F2C=''],[$F2C])
  fi
  if test -n "$LAPACK"  &&  echo "$LAPACK" | grep -v -- "-l" > /dev/null; then
    AC_CHECK_LIB([$LAPACK],[dsyevx],[LAPACK='-l'$LAPACK],[LAPACK='-lla350'
      NEED_LAPACK='yes'],[$BLAS$F2C])
  fi
  # Search for essl library if unsuccessful.
  if test -z "$BLAS"; then
    AC_CHECK_LIB([essl],[dgemm],[BLAS='-lessl'],[BLAS=''])
  fi
  if test -z "$LAPACK"; then
    AC_CHECK_LIB([lapack-essl],[dsyevx],dnl
      [LAPACK='-llapack-essl'],[LAPACK=''],[-lessl])
  fi
elif echo "$target" | grep i686 > /dev/null; then
  ###########################################################
  ## 32-bit Intel-compatible Linux system, compiler flags. ##
  ###########################################################
  if test "$F90" = "pathf95"; then
    ###############################################
    ## PathScale Fortran compiler 2.0 for IA-32. ##
    ###############################################
    if test -z "$DEBUG";   then   DEBUG='-g';      fi
    if test -z "$CDEBUG";  then  CDEBUG='-g';      fi
    if test -z "$BOUND";   then   BOUND='-C';      fi
    if test -z "$DIAG";    then    DIAG='-C';      fi
    if test -z "$OFLAGS1"; then OFLAGS1='-O0';     fi
    if test -z "$OFLAGS2"; then OFLAGS2='-O0';     fi
    if test -z "$OFLAGSC"; then OFLAGSC='-O';      fi
    if test -z "$FPROF";   then   FPROF='-p';      fi
    if test -z "$LDPROF";  then  LDPROF='-p';      fi
    if test -z "$OPENMP";  then  OPENMP='';        fi
    OPT_WITH_DEBUG='yes'
    OPT_WITH_BOUND='yes'
    OPT_WITH_DIAG='yes'
  elif test "$F90" = "ifort"; then
    ###########################################
    ## Intel Fortran compiler 8.1 for IA-32. ##
    ###########################################
    if test -z "$DEBUG";   then   DEBUG='-g';      fi
    if test -z "$CDEBUG";  then  CDEBUG='-g';      fi
    if test -z "$BOUND";   then   BOUND='-C';      fi
    if test -z "$DIAG";    then    DIAG='-C';      fi
    if test -z "$OFLAGS1"; then OFLAGS1='-O1 -fp-model strict'; fi
    if test -z "$OFLAGS2"; then OFLAGS2='-O3';     fi
    if test -z "$OFLAGSC"; then OFLAGSC='-O';      fi
    if test -z "$FPROF";   then   FPROF='-p';      fi
    if test -z "$LDPROF";  then  LDPROF='-p';      fi
    if test -z "$OPENMP";  then  OPENMP='-openmp -fpp2'; fi
    OPT_WITH_DEBUG='yes'
    OPT_WITH_BOUND='no'
    OPT_WITH_DIAG='no'
  elif test "$F90" = "pgf90"; then
    ##############################################
    ## Portland Fortran compiler 5.2 for IA-32. ##
    ##############################################
    if test -z "$DEBUG";   then   DEBUG='-g';     fi
    if test -z "$CDEBUG";  then  CDEBUG='-g';     fi
    if test -z "$BOUND";   then   BOUND='-C';     fi
    if test -z "$DIAG";    then    DIAG='-C';     fi
    if test -z "$OFLAGS1"; then OFLAGS1='-O0';    fi
    if test -z "$OFLAGS2"; then OFLAGS2='-O0';    fi
    if test -z "$OFLAGSC"; then OFLAGSC='-O';     fi
    if test -z "$FPROF";   then   FPROF='-Mprof'; fi
    if test -z "$LDPROF";  then  LDPROF='-Mprof'; fi
    if test -z "$OPENMP";  then  OPENMP='-mp';    fi
    OPT_WITH_DEBUG='yes'
    OPT_WITH_BOUND='yes'
    OPT_WITH_DIAG='yes'
  fi
  ###########################################################
  ## 32-bit Intel-compatible Linux system, math libraries. ##
  ###########################################################
  # Check user-specified blas library.
  # If not found use the default search order.
  if test -n "$F2C"; then F2C=' '$F2C; fi
  if test -n "$BLAS"  &&  echo "$BLAS" | grep -v -- "-l" > /dev/null; then
    AC_CHECK_LIB([$BLAS],[dgemm],[BLAS='-l'$BLAS],[BLAS=''
      LAPACK=''; F2C=''],[$F2C])
  fi
  # Search for Intel's Math Kernel Library (MKL) if unsuccessful.
  if test -z "$BLAS"; then
    AC_CHECK_LIB([mkl_intel],[dgemm],[BLAS='-lmkl_intel'; LAPACK='mkl_intel'
      F2C=' -lmkl_sequential -lmkl_core'],[BLAS=''; LAPACK=''
      F2C=''],[-lmkl_sequential -lmkl_core])
  fi
  # Search for AMD Core Math Library (ACML) if unsuccessful.
  if test -z "$BLAS"; then
    AC_CHECK_LIB([acml],[dgemm],[BLAS='-lacml'; LAPACK='acml'
      F2C=' -lacml_mv'],[BLAS=''; LAPACK=''; F2C=''],[-lacml_mv])
  fi
  # Search for ATLAS library if unsuccessful.
  if test -z "$BLAS"; then
    AC_CHECK_LIB([atlas],[dgemm],[BLAS='-latlas'],[BLAS=''])
  fi
  # Search for generic blas library if unsuccessful.
  if test -z "$BLAS"; then
    AC_CHECK_LIB([blas],[dgemm],[BLAS='-lblas'],[BLAS=''])
  fi
  # Test if BLAS library contains MKL sparse BLAS.
  if test -n "$BLAS"; then
    BLASNAME=${BLAS#-l}
    AC_CHECK_LIB([$BLASNAME],[mkl_dcsrmm],dnl
      [HAVE_CSRMM='yes'],[HAVE_CSRMM='no'],[$F2C])
  fi
  # Check lapack library. If not found use MNDO's private routines.
  if test -n "$LAPACK"  &&  echo "$LAPACK" | grep -v -- "-l" > /dev/null; then
    AC_CHECK_LIB([$LAPACK],[dsyevx],[LAPACK='-l'$LAPACK],[LAPACK='-lla350'
      NEED_LAPACK='yes'],[$BLAS$F2C])
  fi
elif echo "$target" | grep x86_64 > /dev/null; then
  ##################################################
  ## AMD64 or EM64T Linux system, compiler flags. ##
  ##################################################
  if test "$F90" = "ifort"; then
    ########################################
    ## Intel Fortran compiler for x86-64. ##
    ########################################
    if test -z "$DEBUG";   then   DEBUG='-g';      fi
    if test -z "$CDEBUG";  then  CDEBUG='-g';      fi
    if test -z "$BOUND";   then   BOUND='-CB';     fi
    if test -z "$DIAG";    then    DIAG='-C';      fi
    if test -z "$OFLAGS1"; then OFLAGS1='-diag-disable 8290,8291 -O1'; fi
    if test -z "$OFLAGS2"; then OFLAGS2='-diag-disable 8290,8291 -O3'; fi
    if test -z "$OFLAGSC"; then OFLAGSC='-O';      fi
    if test -z "$FPROF";   then   FPROF='-p';      fi
    if test -z "$LDPROF";  then  LDPROF='-p';      fi
    if test -z "$OPENMP";  then  OPENMP='-openmp -fpp2'; fi
    OPT_WITH_DEBUG='yes'
    OPT_WITH_BOUND='no'
    OPT_WITH_DIAG='no'
  elif test "${F90:0:3}" = "pgf"; then
    ######################################################
    ## Portland Group Fortran compiler 10.9 for x86-64. ##
    ######################################################
    if test -z "$DEBUG";   then   DEBUG='-g';     fi
    if test -z "$CDEBUG";  then  CDEBUG='-g';     fi
    if test -z "$BOUND";   then   BOUND='-C';     fi
    if test -z "$DIAG";    then    DIAG='-C';     fi
    if test -z "$OFLAGS1"; then OFLAGS1='-O1';    fi
    if test -z "$OFLAGS2"; then OFLAGS2='-O4';    fi
    if test -z "$OFLAGSC"; then OFLAGSC='-O';     fi
    if test -z "$FPROF";   then   FPROF='-Mprof'; fi
    if test -z "$LDPROF";  then  LDPROF='-Mprof'; fi
    if test -z "$OPENMP";  then  OPENMP='-mp';    fi
    OPT_WITH_DEBUG='yes'
    OPT_WITH_BOUND='yes'
    OPT_WITH_DIAG='yes'
  elif test "$F90" = "pathf95"; then
    ################################################
    ## PathScale Fortran compiler 2.0 for x86-64. ##
    ################################################
    if test -z "$DEBUG";   then   DEBUG='-g';      fi
    if test -z "$CDEBUG";  then  CDEBUG='-g';      fi
    if test -z "$BOUND";   then   BOUND='-C';      fi
    if test -z "$DIAG";    then    DIAG='-C';      fi
    if test -z "$OFLAGS1"; then OFLAGS1='-O2';     fi
    if test -z "$OFLAGS2"; then OFLAGS2='-O2';     fi
    if test -z "$OFLAGSC"; then OFLAGSC='-O';      fi
    if test -z "$FPROF";   then   FPROF='-p';      fi
    if test -z "$LDPROF";  then  LDPROF='-p';      fi
    OPT_WITH_DEBUG='yes'
    OPT_WITH_BOUND='no'
    OPT_WITH_DIAG='no'
  fi
  ##################################################
  ## AMD64 or EM64T Linux system, math libraries. ##
  ##################################################
  # Check user-specified blas library.
  # If not found use the default search order.
  if test -n "$F2C"; then F2C=' '$F2C; fi
  if test -n "$BLAS"  &&  echo "$BLAS" | grep -v -- "-l" > /dev/null; then
    AC_CHECK_LIB([$BLAS],[dgemm],[BLAS='-l'$BLAS],[BLAS=''
      LAPACK=''; F2C=''],[$F2C])
  fi
  # Search for Intel's Math Kernel Library (MKL) if unsuccessful.
  if test -z "$BLAS" -a \( "$MATHLIB" = "MKL" -o -z "$MATHLIB" \); then
    AC_CHECK_LIB([mkl_intel_lp64],[dgemm],[HAVE_MKL='yes'
      BLAS='-lmkl_intel_lp64'; LAPACK='mkl_intel_lp64'
      F2C=' -lmkl_sequential -lmkl_core'],[BLAS=''
      LAPACK=''; F2C=''],[-lmkl_sequential -lmkl_core])
  fi
  # Search for AMD Core Math Library (ACML) if unsuccessful.
  # Test without -lacml_mv (acml5.0.0) first.
  if test -z "$BLAS" -a \( "$MATHLIB" = "ACML" -o -z "$MATHLIB" \); then
    AC_CHECK_LIB([acml],[dgemm],[BLAS='-lacml'; LAPACK='acml'
      F2C=''],[BLAS=''; LAPACK=''; F2C=''])
  fi
  # Test with -lacml_mv (acml4.4.0) if unsuccessful.
  if test -z "$BLAS" -a \( "$MATHLIB" = "ACML" -o -z "$MATHLIB" \); then
    AC_CHECK_LIB([acml],[dgemm],[BLAS='-lacml'; LAPACK='acml'
      F2C=' -lacml_mv'],[BLAS=''; LAPACK=''; F2C=''],[-lacml_mv])
  fi
  # Search for ATLAS library if unsuccessful.
  if test -z "$BLAS" -a -z "$MATHLIB"; then
    AC_CHECK_LIB([atlas],[dgemm],[BLAS='-latlas'],[BLAS=''])
  fi
  # Search for generic blas library if unsuccessful.
  if test -z "$BLAS" -a -z "$MATHLIB"; then
    AC_CHECK_LIB([blas],[dgemm],[BLAS='-lblas'],[BLAS=''])
  fi
  # Test if BLAS library contains MKL sparse BLAS.
  if test -n "$BLAS"; then
    BLASNAME=${BLAS#-l}
    AC_CHECK_LIB([$BLASNAME],[mkl_dcsrmm],dnl
      [HAVE_CSRMM='yes'],[HAVE_CSRMM='no'],[$F2C])
  fi
  # Check lapack library. If not found use MNDO's private routines.
  if test -n "$LAPACK"  &&  echo "$LAPACK" | grep -v -- "-l" > /dev/null; then
    AC_CHECK_LIB([$LAPACK],[dsyevx],[LAPACK='-l'$LAPACK],[LAPACK='-lla350'
      NEED_LAPACK='yes'],[$BLAS$F2C])
  fi
else
  #####################
  ## Generic system. ##
  #####################
  if test -z "$DEBUG";   then   DEBUG='-g'; fi
  if test -z "$CDEBUG";  then  CDEBUG='-g'; fi
  if test -z "$BOUND";   then   BOUND='-C'; fi
  if test -z "$DIAG";    then    DIAG='-C'; fi
  if test -z "$OFLAGS2"; then OFLAGS2='-O'; fi
  if test -z "$OFLAGSC"; then OFLAGSC='-O'; fi
  OPT_WITH_DEBUG='no'
  OPT_WITH_BOUND='no'
  OPT_WITH_DIAG='no'
fi

######################################################################
## By default, assume g++ and nvcc for C++ and CUDA on all systems. ##
######################################################################
if test -z "$CXXDEBUG";  then  CXXDEBUG='-g'; fi
if test -z "$NVDEBUG";   then   NVDEBUG='-g'; fi
if test -z "$OFLAGSCXX"; then OFLAGSCXX='-O2'; fi
if test -z "$OFLAGSNV";  then  OFLAGSNV='-arch=sm_20'; fi

# Use math libraries from MNDO distribution if none have been found so far.
if test -z "$BLAS"; then
  BLAS='-lbl2011'
  LAPACK='-lla350'
  NEED_BLAS='yes'
  NEED_BLAS3='yes'
  NEED_LAPACK='yes'
  # 
fi
if test -z "$LAPACK"; then
  LAPACK='-lla350'
  NEED_LAPACK='yes'
fi

# cpudavliu.cxx in the GUGACI module requires mkl_dcsrmm,
# so if mkl_dcsrmm is not available, cpudavliu.cxx must
# not be compiled (prevented by disabling C++).
if test "$HAVE_CSRMM" = 'no'; then
  HAVE_CXX='no'
fi

# Check PGPLOT library.
if test -n "$PGPLOT"; then
  AC_CHECK_LIB([$PGPLOT],[pgbeg],[PGPLOT='-l'$PGPLOT],[PGPLOT=''])
fi
if test -z "$PGPLOT"; then
  NEED_NOPGPLOT='yes'
  NEED_MACHINE='yes'
fi

# Check if the optional content of libblas.a matches the current needs and
# delete this library file otherwise. The tests on libmachine.a have been
# abandoned because the AC_CHECK_LIB calls cannot be arranged in a loop.
# Always delete libmachine.a and libguga.a.
#
if test -f "BLAS/libbl2011.a"; then
  AC_CHECK_LIB([bl2011],[dgemm],[HAVE_BLAS='yes'],[HAVE_BLAS='no'])
  if test "$HAVE_BLAS" = "no"; then rm -f BLAS/libbl2011.a libbl2011.a; fi
  AC_CHECK_LIB([bl2011],[dsymm],[HAVE_BLAS3='yes'],[HAVE_BLAS3='no'])
else
  HAVE_BLAS='no'
  HAVE_BLAS3='no'
fi
if test "$NEED_BLAS"   = "yes"  &&  \
   test "$HAVE_BLAS"   = "yes"  &&  \
   test "$NEED_BLAS3" != "$HAVE_BLAS3"; then
  rm -f BLAS/libbl2011.a
fi
rm -f MACHINE/libmachine.a
rm -f GUGACI/libguga.a

# Check required C functions.
AC_LANG(C)
AC_CHECK_FUNCS([gettimeofday])

# Define macro to derive FFLAGSn_ARG from $OFLAGSn and other variables.
# AC_SET_FFLAGS(FFLAGSn_ARG, arg, OFLAGSn)
AC_DEFUN([AC_SET_FFLAGS],dnl
[{
  # Clear optimization flags if they don't work with -g or -C.
  eval $1="\$$3"
  if (echo "$enable_bound" | eval grep "$2" > /dev/null)  &&  \
      test -n "$BOUND"  &&  test "$OPT_WITH_BOUND" = "no"; then
    eval $1="''"
  fi
  if (echo "$enable_debug" | eval grep "$2" > /dev/null)  &&  \
      test -n "$DEBUG"  &&  test "$OPT_WITH_DEBUG" = "no"; then
    eval $1="''"
  fi
  if (echo "$enable_diag" | eval grep "$2" > /dev/null)  &&  \
      test -n "$DIAG"  &&  test "$OPT_WITH_DIAG" = "no"; then
    eval $1="''"
  fi
  # Add FFLAGS if length is non-zero.
  if test -n "$FFLAGS"  &&  test -n "$(eval echo \$$1)"; then
    eval $1="\$FFLAGS' '\$$1"
  else
    eval $1="\$FFLAGS\$$1"
  fi
  # Add debug flags.
  if (echo "$enable_debug" | eval grep "$2" > /dev/null)  &&  \
      test -n "$DEBUG"; then
        if test -n "$(eval echo \$$1)"; then
          eval $1="\$DEBUG' '\$$1"
        else
          eval $1="\$DEBUG"
        fi
  fi
  # Add diagnostics flags.
  if (echo "$enable_diag" | eval grep "$2" > /dev/null)  &&  \
      test -n "$DIAG"; then
        if test -n "$(eval echo \$$1)"; then
          eval $1="\$DIAG' '\$$1"
        else
          eval $1="\$DIAG"
        fi
  elif (echo "$enable_bound" | eval grep "$2" > /dev/null)  &&  \
        test -n "$BOUND"; then
          if test -n "$(eval echo \$$1)"; then
            eval $1="\$BOUND' '\$$1"
          else
            eval $1="\$BOUND"
          fi
  fi
  # Add profiling flags.
  if test "$enable_prof" != "no"; then
    if test "$enable_prof" = "yes"; then
      if test -n "$FPROF"  &&  test -n "$(eval echo  \$$1)"; then
        eval $1="\$FPROF' '\$$1"
      else
        eval $1="\$FPROF\$$1"
      fi
    else
      eval $1="\$enable_prof' '\$$1"
    fi
  fi
  # Add OpenMP flags.
  if test "$enable_openmp" != "no"; then
    if test "$enable_openmp" = "yes"; then
      if test -n "$OPENMP"  && test -n "$(eval echo \$$1)"; then
        eval $1="\$OPENMP' '\$$1"
      else
        eval $1="\$OPENMP\$$1"
      fi
    else
      eval $1="\$enable_openmp' '\$$1"
    fi
  fi
}])

# Define macro to derive CFLAGS_ARG from $OFLAGSC and other variables.
# AC_SET_CFLAGS(CFLAGS_ARG, arg)
AC_DEFUN([AC_SET_CFLAGS],
[{
  $1=$OFLAGSC
  if (echo "$enable_debug" | grep "$2" > /dev/null)  &&  \
      test -n "$CDEBUG"; then
    $1=$CDEBUG
  fi
  if test -n "$CFLAGS"  &&  test -n "$$1"; then
    $1=$CFLAGS$' '$$1
  else
    $1=$CFLAGS$$1
  fi
}])

# Define macro to derive CXXFLAGS_ARG from $OFLAGSCXX and other variables.
# AC_SET_CXXFLAGS(CXXFLAGS_ARG, arg)
AC_DEFUN([AC_SET_CXXFLAGS],
[{
  $1=$OFLAGSCXX
  if (echo "$enable_debug" | grep "$2" > /dev/null)  &&  \
      test -n "$CXXDEBUG"; then
    $1=$CXXDEBUG
  fi
  if test -n "$CXXFLAGS"  &&  test -n "$$1"; then
    $1=$CXXFLAGS$' '$$1
  else
    $1=$CXXFLAGS$$1
  fi
}])

# Define macro to derive NVFLAGS_ARG from $OFLAGSNV and other variables.
# AC_SET_NVFLAGS(NVFLAGS_ARG, arg)
AC_DEFUN([AC_SET_NVFLAGS],
[{
  $1=$OFLAGSNV
  if (echo "$enable_debug" | grep "$2" > /dev/null)  &&  \
      test -n "$NVDEBUG"; then
    $1=$NVDEBUG
  fi
  if test -n "$NVFLAGS"  &&  test -n "$$1"; then
    $1=$NVFLAGS$' '$$1
  else
    $1=$NVFLAGS$$1
  fi
}])

# Derive all FFLAGS and CFLAGS variables.
AC_SET_FFLAGS([FFLAGS1_MNDO],      [mndo],        [OFLAGS1])
AC_SET_FFLAGS([FFLAGS2_MNDO],      [mndo],        [OFLAGS2])
for i in $SUBDIRS; do
AC_SET_FFLAGS([FFLAGS1_$i],        [\$OPTION_$i], [OFLAGS1])
AC_SET_FFLAGS([FFLAGS2_$i],        [\$OPTION_$i], [OFLAGS2])
done
AC_SET_FFLAGS([FFLAGS1_BLAS],      [blas],        [OFLAGS1])
AC_SET_FFLAGS([FFLAGS2_BLAS],      [blas],        [OFLAGS2])
AC_SET_FFLAGS([FFLAGS1_LAPACK],    [lapack],      [OFLAGS1])
AC_SET_FFLAGS([FFLAGS2_LAPACK],    [lapack],      [OFLAGS2])
AC_SET_FFLAGS([FFLAGS1_MACHINE],   [machine],     [OFLAGS1])
AC_SET_FFLAGS([FFLAGS2_MACHINE],   [machine],     [OFLAGS2])
AC_SET_CFLAGS([CFLAGS_GUGACI],     [guga])
AC_SET_CFLAGS([CFLAGS_MACHINE],    [machine])
AC_SET_CXXFLAGS([CXXFLAGS_GUGACI], [guga])
AC_SET_NVFLAGS([NVFLAGS_GUGACI],   [guga])
AC_SET_NVFLAGS([NVFLAGS_MACHINE],  [machine])

# Assemble ARCHIVES, LDFLAGS, and LIBS which are used in the Makefile.
if test "$BLAS" = "-lbl2011"  &&  test "$LAPACK" = "-lla350"; then F2C=''; fi
if test "$NEED_BLAS" = "yes"; then
  ARCHIVES=$ARCHIVES' libbl2011.a'
fi
if test "$NEED_LAPACK" = "yes"; then
  ARCHIVES=$ARCHIVES' libla350.a'
fi
if test "$NEED_MACHINE" = "yes"; then
  MNDOLIBS=$MNDOLIBS' -lmachine'
  ARCHIVES=$ARCHIVES' libmachine.a'
fi
if test "$enable_prof" != "no"; then
  if test "$enable_prof" = "yes"; then
    if test -n "$LDPROF"; then
      LDFLAGS=$LDPROF' '$LDFLAGS
    fi
  else
    LDFLAGS=$enable_prof' '$LDFLAGS
  fi
fi
if test "$enable_openmp" != "no"; then
  if test "$enable_openmp" = "yes"; then
    if test -n "$OPENMP"; then
      LDFLAGS=$OPENMP' '$LDFLAGS
    fi
  else
    LDFLAGS=$enable_openmp' '$LDFLAGS
  fi
fi
if test -n "$PGPLOT";  then MNDOLIBS=$MNDOLIBS' '$PGPLOT; fi
if test -n "$LAPACK";  then MNDOLIBS=$MNDOLIBS' '$LAPACK; fi
if test -n "$BLAS" && test "$BLAS" != "$LAPACK"
                       then MNDOLIBS=$MNDOLIBS' '$BLAS;   fi
LIBS=$MNDOLIBS$FLIBS$F2C$CXXLIBS$LIBS

# ----------------------------------------------------------------------
# Automake conditionals for various subdirectories.
# This section cannot be arranged in a loop because the
# AM_CONDITIONAL macros are interpreted by automake.
# ----------------------------------------------------------------------

AM_CONDITIONAL([COND_PS],       [test "$NEED_PS"        = "yes"])
AM_CONDITIONAL([COND_GUGACI],   [test "$NEED_GUGACI"    = "yes"])
AM_CONDITIONAL([COND_HDLCOPT],  [test "$NEED_HDLCOPT"   = "yes"])
AM_CONDITIONAL([COND_PSEX],     [test "$NEED_PSEX"      = "yes"])
AM_CONDITIONAL([COND_CGDMS],    [test "$NEED_CGDMS"     = "yes"])
AM_CONDITIONAL([COND_SPARSE],   [test "$NEED_SPARSE"    = "yes"])
AM_CONDITIONAL([COND_SCOSMO],   [test "$NEED_SCOSMO"    = "yes"])
AM_CONDITIONAL([COND_DFTB],     [test "$NEED_DFTB"      = "yes"])
AM_CONDITIONAL([COND_VB],       [test "$NEED_VB"        = "yes"])
AM_CONDITIONAL([COND_ES],       [test "$NEED_ES"        = "yes"])
AM_CONDITIONAL([COND_NOPS],     [test "$NEED_NOPS"      = "yes"])
AM_CONDITIONAL([COND_NOGUGACI], [test "$NEED_NOGUGACI"  = "yes"])
AM_CONDITIONAL([COND_NOPSEX],   [test "$NEED_NOPSEX"    = "yes"])
AM_CONDITIONAL([COND_NOHDLCOPT],[test "$NEED_NOHDLCOPT" = "yes"])
AM_CONDITIONAL([COND_NOCGDMS],  [test "$NEED_NOCGDMS"   = "yes"])
AM_CONDITIONAL([COND_NOSPARSE], [test "$NEED_NOSPARSE"  = "yes"])
AM_CONDITIONAL([COND_NOSCOSMO], [test "$NEED_NOSCOSMO"  = "yes"])
AM_CONDITIONAL([COND_NODFTB],   [test "$NEED_NODFTB"    = "yes"])
AM_CONDITIONAL([COND_NOVB],     [test "$NEED_NOVB"      = "yes"])
AM_CONDITIONAL([COND_NOES],     [test "$NEED_NOES"      = "yes"])
AM_CONDITIONAL([COND_NOPGPLOT], [test "$NEED_NOPGPLOT"  = "yes"])
AM_CONDITIONAL([COND_BLAS],     [test "$NEED_BLAS"      = "yes"])
AM_CONDITIONAL([COND_BLAS3],    [test "$NEED_BLAS3"     = "yes"])
AM_CONDITIONAL([COND_LAPACK],   [test "$NEED_LAPACK"    = "yes"])
AM_CONDITIONAL([COND_MACHINE],  [test "$NEED_MACHINE"   = "yes"])
AM_CONDITIONAL([COND_DEFMEM],   [test "$NEED_DEFMEM"    = "yes"])
AM_CONDITIONAL([COND_DMR],      [test "$NEED_DMR"       = "yes"])
AM_CONDITIONAL([COND_MEMX],     [test "$NEED_MEMX"      = "yes"])
AM_CONDITIONAL([COND_GEMMBL],   [test "$NEED_GEMMBL"    = "yes"])
AM_CONDITIONAL([COND_CXX],      [test "$HAVE_CXX"       = "yes"])
AM_CONDITIONAL([COND_CUDA],     [test "$HAVE_CUDA"      = "yes"])
AM_CONDITIONAL([COND_CSRMM],    [test "$HAVE_CSRMM"     = "yes"])

# ----------------------------------------------------------------------
# Variables not automatically substituted in the Makefiles.
# This section cannot be arranged in a loop because the
# AC_SUBST macros are interpreted by automake.
# ----------------------------------------------------------------------

AC_SUBST(F90)
AC_SUBST(FFLAGS1_MNDO)
AC_SUBST(FFLAGS2_MNDO)
AC_SUBST(FFLAGS1_PS)
AC_SUBST(FFLAGS2_PS)
AC_SUBST(FFLAGS1_GUGACI)
AC_SUBST(FFLAGS2_GUGACI)
AC_SUBST(FFLAGS1_PSEX)
AC_SUBST(FFLAGS2_PSEX)
AC_SUBST(FFLAGS1_HDLCOPT)
AC_SUBST(FFLAGS2_HDLCOPT)
AC_SUBST(FFLAGS1_CGDMS)
AC_SUBST(FFLAGS2_CGDMS)
AC_SUBST(FFLAGS1_SPARSE)
AC_SUBST(FFLAGS2_SPARSE)
AC_SUBST(FFLAGS1_SCOSMO)
AC_SUBST(FFLAGS2_SCOSMO)
AC_SUBST(FFLAGS1_DFTB)
AC_SUBST(FFLAGS2_DFTB)
AC_SUBST(FFLAGS1_VB)
AC_SUBST(FFLAGS2_VB)
AC_SUBST(FFLAGS1_ES)
AC_SUBST(FFLAGS2_ES)
AC_SUBST(FFLAGS1_BLAS)
AC_SUBST(FFLAGS2_BLAS)
AC_SUBST(FFLAGS1_LAPACK)
AC_SUBST(FFLAGS2_LAPACK)
AC_SUBST(FFLAGS1_MACHINE)
AC_SUBST(FFLAGS2_MACHINE)
AC_SUBST(CFLAGS_GUGACI)
AC_SUBST(CFLAGS_MACHINE)
AC_SUBST(CXXFLAGS_GUGACI)
AC_SUBST(NVFLAGS_GUGACI)
AC_SUBST(NVFLAGS_MACHINE)
AC_SUBST(ARCHIVES)
AC_SUBST(CLEANFILES)
AC_SUBST(OPTDIR_PS)
AC_SUBST(OPTDIR_GUGACI)
AC_SUBST(OPTDIR_PSEX)
AC_SUBST(OPTDIR_HDLCOPT)
AC_SUBST(OPTDIR_CGDMS)
AC_SUBST(OPTDIR_SPARSE)
AC_SUBST(OPTDIR_SCOSMO)
AC_SUBST(OPTDIR_DFTB)
AC_SUBST(OPTDIR_VB)
AC_SUBST(OPTDIR_ES)

# ----------------------------------------------------------------------
# Generate Makefiles for all existing directories.
# This section cannot be arranged in a loop because
# the AC_CONFIG_FILES macros are interpreted by automake.
# ----------------------------------------------------------------------

AC_CONFIG_FILES(Makefile)
AC_CONFIG_FILES(2020/Makefile)
AC_CONFIG_FILES(BLAS/Makefile)
AC_CONFIG_FILES(LAPACK/Makefile)
AC_CONFIG_FILES(MACHINE/Makefile)

if test -d "PS"; then
  OPTDIR_PS='PS'
  AC_CONFIG_FILES(PS/Makefile)
else
  OPTDIR_PS=''
fi

if test -d "GUGACI"; then
  OPTDIR_GUGACI='GUGACI'
  AC_CONFIG_FILES(GUGACI/Makefile)
else
  OPTDIR_GUGACI=''
fi

if test -d "PSEX"; then
  OPTDIR_PSEX='PSEX'
  AC_CONFIG_FILES(PSEX/Makefile)
else
  OPTDIR_PSEX=''
fi

if test -d "HDLCOPT"; then
  OPTDIR_HDLCOPT='HDLCOPT'
  AC_CONFIG_FILES(HDLCOPT/Makefile)
else
  OPTDIR_HDLCOPT=''
fi

if test -d "CGDMS"; then
  OPTDIR_CGDMS='CGDMS'
  AC_CONFIG_FILES(CGDMS/Makefile)
else
  OPTDIR_CGDMS=''
fi

if test -d "SPARSE"; then
  OPTDIR_SPARSE='SPARSE'
  AC_CONFIG_FILES(SPARSE/Makefile)
else
  OPTDIR_SPARSE=''
fi

if test -d "SCOSMO"; then
  OPTDIR_SCOSMO='SCOSMO'
  AC_CONFIG_FILES(SCOSMO/Makefile)
else
  OPTDIR_SCOSMO=''
fi

if test -d "DFTB"; then
  OPTDIR_DFTB='DFTB'
  AC_CONFIG_FILES(DFTB/Makefile)
else
  OPTDIR_DFTB=''
fi

if test -d "VB"; then
  OPTDIR_VB='VB'
  AC_CONFIG_FILES(VB/Makefile)
else
  OPTDIR_VB=''
fi

if test -d "ES"; then
  OPTDIR_ES='ES'
  AC_CONFIG_FILES(ES/Makefile)
else
  OPTDIR_ES=''
fi

# ----------------------------------------------------------------------
# End of subdirectory-dependent sections.
# ----------------------------------------------------------------------

# Generate all files!
AC_OUTPUT

# We're all set! Print a summary of the options used.
AC_MSG_RESULT([
MNDO has been configured successfully.])

AC_MSG_RESULT([Enabled subdirectories: $ENABLED_SUBDIRS])

if test "$BLAS" = "$LAPACK"; then
  AC_MSG_RESULT([The math library $BLAS will be used.])
else
  if test "$BLAS" = "-lbl2011"; then
    if test "$NEED_GEMMBL" = "yes"; then
      AC_MSG_RESULT([MNDO's private blas library with gemmbl.f will be used.])
    else
      AC_MSG_RESULT([MNDO's private blas library will be used.])
    fi
  else
    AC_MSG_RESULT([The blas library $BLAS will be used.])
  fi
  if test "$LAPACK" = "-lla350"; then
    AC_MSG_RESULT([MNDO's private lapack library will be used.])
  else
    AC_MSG_RESULT([The lapack library $LAPACK will be used.])
  fi
fi

if test -n "$PGPLOT"; then
  AC_MSG_RESULT([The pgplot library $PGPLOT will be used.])
else
  AC_MSG_RESULT([The pgplot library will not be used.])
fi

if test "$enable_openmp" != "no"; then
  AC_MSG_RESULT([The code will be parallelized using OpenMP.])
fi

if test "$enable_bound" != "no"; then
  AC_MSG_RESULT([Boundary checking: $enable_bound.])
fi

if test "$enable_diag" != "no"; then
  AC_MSG_RESULT([Full runtime diagnosis: $enable_diag.])
fi

if test "$enable_debug" != "no"; then
  AC_MSG_RESULT([Debugging symbols: $enable_debug.])
fi

if test "$enable_prof" != "no"; then
  AC_MSG_RESULT([Profiling: $enable_prof])
fi

AC_MSG_RESULT([
See './configure --help' for an explanation of all
available options. The next step is to run 'make'.])
