# VolEsti (volume computation and sampling library)
# Copyright (c) 2012-2020 Vissarion Fisikopoulos
# Copyright (c) 2018-2020 Apostolos Chalkis
# Contributed and/or modified by Vaibhav Thakkar
# Licensed under GNU LGPL.3, see LICENCE file

project( VolEsti )

enable_testing()

CMAKE_MINIMUM_REQUIRED(VERSION 3.11)

set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)

# Locate Intel MKL root (in case it is enabled)

if (APPLE)
	set(MKLROOT /opt/intel/oneapi/mkl/latest)
elseif(UNIX)
	set(MKLROOT $ENV{HOME}/intel/mkl)
endif()

option(DISABLE_NLP_ORACLES "Disable non-linear oracles (used in collocation)" ON)
option(BUILTIN_EIGEN "Use eigen from ../external" OFF)
option(USE_MKL "Use MKL library to build eigen" OFF)

if(DISABLE_NLP_ORACLES)
  add_definitions(-DDISABLE_NLP_ORACLES)
else()
  find_library(IFOPT NAMES libifopt_core.so PATHS /usr/local/lib)
  find_library(IFOPT_IPOPT NAMES libifopt_ipopt.so PATHS /usr/local/lib)
  find_library(GMP NAMES libgmp.so PATHS /usr/lib/x86_64-linux-gnu /usr/lib/i386-linux-gnu)
  find_library(MPSOLVE NAMES libmps.so PATHS /usr/local/lib)
  find_library(PTHREAD NAMES libpthread.so PATHS /usr/lib/x86_64-linux-gnu /usr/lib/i386-linux-gnu)
  find_library(FFTW3 NAMES libfftw3.so.3 PATHS /usr/lib/x86_64-linux-gnu /usr/lib/i386-linux-gnu)

  if (NOT IFOPT)

    message(FATAL_ERROR "This program requires the ifopt library, and will not be compiled.")

  elseif (NOT GMP)

    message(FATAL_ERROR "This program requires the gmp library, and will not be compiled.")

  elseif (NOT MPSOLVE)

    message(FATAL_ERROR "This program requires the mpsolve library, and will not be compiled.")

  elseif (NOT FFTW3)

    message(FATAL_ERROR "This program requires the fftw3 library, and will not be compiled.")

  else()
    message(STATUS "Library ifopt found: ${IFOPT}")
    message(STATUS "Library gmp found: ${GMP}")
    message(STATUS "Library mpsolve found: ${MPSOLVE}")
    message(STATUS "Library fftw3 found:" ${FFTW3})

  endif(NOT IFOPT)

endif(DISABLE_NLP_ORACLES)

option(DISABLE_NLP_ORACLES "Disable non-linear oracles (used in collocation)" ON)

if(DISABLE_NLP_ORACLES)
  add_definitions(-DDISABLE_NLP_ORACLES)
else()
  find_library(IFOPT NAMES libifopt_core.so PATHS /usr/local/lib)
  find_library(IFOPT_IPOPT NAMES libifopt_ipopt.so PATHS /usr/local/lib)
  find_library(GMP NAMES libgmp.so PATHS /usr/lib/x86_64-linux-gnu /usr/lib/i386-linux-gnu)
  find_library(MPSOLVE NAMES libmps.so PATHS /usr/local/lib)
  find_library(PTHREAD NAMES libpthread.so PATHS /usr/lib/x86_64-linux-gnu /usr/lib/i386-linux-gnu)
  find_library(FFTW3 NAMES libfftw3.so.3 PATHS /usr/lib/x86_64-linux-gnu /usr/lib/i386-linux-gnu)

  if (NOT IFOPT)

    message(FATAL_ERROR "This program requires the ifopt library, and will not be compiled.")

  elseif (NOT GMP)

    message(FATAL_ERROR "This program requires the gmp library, and will not be compiled.")

  elseif (NOT MPSOLVE)

    message(FATAL_ERROR "This program requires the mpsolve library, and will not be compiled.")

  elseif (NOT FFTW3)

    message(FATAL_ERROR "This program requires the fftw3 library, and will not be compiled.")

  else()
    message(STATUS "Library ifopt found: ${IFOPT}")
    message(STATUS "Library gmp found: ${GMP}")
    message(STATUS "Library mpsolve found: ${MPSOLVE}")
    message(STATUS "Library fftw3 found:" ${FFTW3})

  endif(NOT IFOPT)

endif(DISABLE_NLP_ORACLES)

if(COMMAND cmake_policy)
       cmake_policy(SET CMP0003 NEW)
endif(COMMAND cmake_policy)

include("lpSolver.cmake")
GetLPsolver()

include("Eigen.cmake")
GetEigen()

#if(NOT DEFINED BOOST)
#  message(FATAL_ERROR "This program requires the boost library, and will not be compiled.")
#else()
#  message(STATUS "BOOST is DEFINED")
#endif()

include("Boost.cmake")
GetBoost()

find_library(LP_SOLVE NAMES liblpsolve55.so liblpsolve55.dylib PATHS /usr/lib/lp_solve /opt/local/lib)

if (NOT LP_SOLVE)

  message(FATAL_ERROR "This program requires the lp_solve library, and will not be compiled.")

else ()

  message(STATUS "Library lp_solve found: ${LP_SOLVE}")
  set(CMAKE_EXPORT_COMPILE_COMMANDS "ON")

  if (BUILTIN_EIGEN)
      include_directories (BEFORE ../external/Eigen)
  else ()
      include_directories(BEFORE /usr/include/eigen3)
  endif(BUILTIN_EIGEN)

  if (USE_MKL)
    find_library(BLAS NAMES libblas.so libblas.dylib PATHS /usr/lib/x86_64-linux-gnu /usr/lib/i386-linux-gnu /usr/lib)
    find_library(GFORTRAN NAME libgfortran.dylib PATHS /usr/local/Cellar/gcc/10.2.0_4/lib/gcc/10)
    find_library(LAPACK NAME liblapack.dylib PATHS /usr/lib)
    include_directories (BEFORE ${MKLROOT}/include)
    set(PROJECT_LIBS ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} ${GFORTRAN_LIBRARIES})
    add_definitions(-DEIGEN_USE_MKL_ALL)
  endif(USE_MKL)

  include_directories (BEFORE ../external)
  include_directories (BEFORE ../external/minimum_ellipsoid)
  include_directories (BEFORE ../include/generators)
  include_directories (BEFORE ../include/volume)
  include_directories (BEFORE ../include)
  include_directories (BEFORE ../include/convex_bodies)
  include_directories (BEFORE ../include/convex_bodies/spectrahedra)
  include_directories (BEFORE ../include/annealing)
  include_directories (BEFORE ../include/ode_solvers)
  include_directories (BEFORE ../include/root_finders)
  include_directories (BEFORE ../include/samplers)
  include_directories (BEFORE ../include/lp_oracles)
  include_directories (BEFORE ../include/nlp_oracles)
  include_directories (BEFORE ../include/misc)

  #for Eigen
  if (${CMAKE_VERSION} VERSION_LESS "3.12.0")
    add_compile_options(-D "EIGEN_NO_DEBUG")
  else ()
    add_compile_definitions("EIGEN_NO_DEBUG")
  endif ()

  add_definitions(${CMAKE_CXX_FLAGS} "-std=c++11")  # enable C++11 standard
  add_definitions(${CMAKE_CXX_FLAGS} "-g")  # enable debuger
  #add_definitions(${CMAKE_CXX_FLAGS} "-Wint-in-bool-context")
  #add_definitions(${CMAKE_CXX_FLAGS} "-Wall")

  add_definitions(${CMAKE_CXX_FLAGS} "-O3")  # optimization of the compiler
  #add_definitions(${CXX_COVERAGE_COMPILE_FLAGS} "-lgsl")
  add_definitions(${CXX_COVERAGE_COMPILE_FLAGS} "-lm")
  add_definitions(${CXX_COVERAGE_COMPILE_FLAGS} "-ldl")
  add_definitions(${CXX_COVERAGE_COMPILE_FLAGS} "-DBOOST_NO_AUTO_PTR")
  #add_definitions(${CXX_COVERAGE_COMPILE_FLAGS} "-lgslcblas")
  #add_definitions( "-O3 -lgsl -lm -ldl -lgslcblas" )

  add_executable (new_volume_example new_volume_example.cpp)
  add_executable (benchmarks_sob benchmarks_sob.cpp)
  add_executable (benchmarks_cg benchmarks_cg.cpp)
  add_executable (benchmarks_cb benchmarks_cb.cpp)

  add_library(test_main OBJECT test_main.cpp)

  add_executable (mcmc_diagnostics_test mcmc_diagnostics_test.cpp $<TARGET_OBJECTS:test_main>)
  add_test(NAME test_psrf COMMAND mcmc_diagnostics_test -tc=psrf)
  add_test(NAME test_geweke COMMAND mcmc_diagnostics_test -tc=geweke)
  add_test(NAME test_raftery COMMAND mcmc_diagnostics_test -tc=raftery)

  add_executable (ode_solvers_test ode_solvers_test.cpp $<TARGET_OBJECTS:test_main>)
  add_test(NAME ode_solvers_test_first_order
          COMMAND ode_solvers_test -tc=first_order)
  add_test(NAME ode_solvers_test_second_order
          COMMAND ode_solvers_test -tc=second_order)

  add_executable (root_finders_test root_finders_test.cpp $<TARGET_OBJECTS:test_main>)
  add_test(NAME root_finders_test_root_finders
          COMMAND root_finders_test -tc=root_finders)

  add_executable (boundary_oracles_test boundary_oracles_test.cpp $<TARGET_OBJECTS:test_main>)
  add_test(NAME boundary_oracles_test_h_poly_oracles
    COMMAND boundary_oracles_test -tc=h_poly_oracles)

  add_executable(test_sdpa_format test_sdpa_format.cpp $<TARGET_OBJECTS:test_main>)
  add_test(NAME test_sdpa_format COMMAND test_sdpa_format -tc=sdpa_format_parser)

  add_executable (volume_sob_hpolytope volume_sob_hpolytope.cpp $<TARGET_OBJECTS:test_main>)
  add_test(NAME volume_sob_hpolytope_cube COMMAND volume_sob_hpolytope -tc=cube)
  add_test(NAME volume_sob_hpolytope_cross COMMAND volume_sob_hpolytope -tc=cross)
  add_test(NAME volume_sob_hpolytope_birkhoff COMMAND volume_sob_hpolytope -tc=birk)
  add_test(NAME volume_sob_hpolytope_prod_simplex COMMAND volume_sob_hpolytope -tc=prod_simplex)
  add_test(NAME volume_sob_hpolytope_simplex COMMAND volume_sob_hpolytope -tc=simplex)
  add_test(NAME volume_sob_hpolytope_skinny_cube COMMAND volume_sob_hpolytope -tc=skinny_cube)

  add_executable (volume_sob_vpolytope volume_sob_vpolytope.cpp $<TARGET_OBJECTS:test_main>)
  add_test(NAME volume_sob_vpolytope_cube COMMAND volume_sob_vpolytope -tc=cube)
  add_test(NAME volume_sob_vpolytope_cross COMMAND volume_sob_vpolytope -tc=cross)
  add_test(NAME volume_sob_vpolytope_simplex COMMAND volume_sob_vpolytope -tc=simplex)

  add_executable (volume_cg_hpolytope volume_cg_hpolytope.cpp $<TARGET_OBJECTS:test_main>)
  add_test(NAME volume_cg_hpolytope_cube COMMAND volume_cg_hpolytope -tc=cube)
  add_test(NAME volume_cg_hpolytope_cross COMMAND volume_cg_hpolytope -tc=cross)
  add_test(NAME volume_cg_hpolytope_birkhoff COMMAND volume_cg_hpolytope -tc=birk)
  add_test(NAME volume_cg_hpolytope_prod_simplex COMMAND volume_cg_hpolytope -tc=prod_simplex)
  add_test(NAME volume_cg_hpolytope_simplex COMMAND volume_cg_hpolytope -tc=simplex)
  add_test(NAME volume_cg_hpolytope_skinny_cube COMMAND volume_cg_hpolytope -tc=skinny_cube)

  add_executable (volume_cg_vpolytope volume_cg_vpolytope.cpp $<TARGET_OBJECTS:test_main>)
  add_test(NAME volume_cg_vpolytope_cube COMMAND volume_cg_vpolytope -tc=cube)
  add_test(NAME volume_cg_vpolytope_cross COMMAND volume_cg_vpolytope -tc=cross)
  add_test(NAME volume_cg_vpolytope_simplex COMMAND volume_cg_vpolytope -tc=simplex)

  add_executable (volume_cb_hpolytope volume_cb_hpolytope.cpp $<TARGET_OBJECTS:test_main>)
  add_test(NAME volume_cb_hpolytope_cube COMMAND volume_cb_hpolytope -tc=cube)
  add_test(NAME volume_cb_hpolytope_cross COMMAND volume_cb_hpolytope -tc=cross)
  add_test(NAME volume_cb_hpolytope_birkhoff COMMAND volume_cb_hpolytope -tc=birk)
  add_test(NAME volume_cb_hpolytope_prod_simplex COMMAND volume_cb_hpolytope -tc=prod_simplex)
  add_test(NAME volume_cb_hpolytope_simplex COMMAND volume_cb_hpolytope -tc=simplex)
  add_test(NAME volume_cb_hpolytope_skinny_cube COMMAND volume_cb_hpolytope -tc=skinny_cube)

  add_executable (volume_cb_vpolytope volume_cb_vpolytope.cpp $<TARGET_OBJECTS:test_main>)
  add_test(NAME volume_cb_vpolytope_cube COMMAND volume_cb_vpolytope -tc=cube)
  add_test(NAME volume_cb_vpolytope_cross COMMAND volume_cb_vpolytope -tc=cross)
  add_test(NAME volume_cb_vpolytope_simplex COMMAND volume_cb_vpolytope -tc=simplex)

  add_executable (volume_cb_zonotopes volume_cb_zonotopes.cpp $<TARGET_OBJECTS:test_main>)
  add_test(NAME volume_cb_zonotopes_uniform_zonotopes
           COMMAND volume_cb_zonotopes -tc=uniform_zonotopes)

  add_executable (volume_cb_vpoly_intersection_vpoly volume_cb_vpoly_intersection_vpoly.cpp $<TARGET_OBJECTS:test_main>)
  add_test(NAME volume_cb_vpoly_intersection_vpoly_random_vpoly_sphere
           COMMAND volume_cb_vpoly_intersection_vpoly -tc=random_vpoly_sphere)

  add_executable (new_rounding_test new_rounding_test.cpp $<TARGET_OBJECTS:test_main>)
  add_test(NAME new_rounding_test_round_skinny_cube
           COMMAND new_rounding_test -tc=round_skinny_cube)

  add_executable (logconcave_sampling_test logconcave_sampling_test.cpp $<TARGET_OBJECTS:test_main>)
  add_test(NAME logconcave_sampling_test_hmc
          COMMAND logconcave_sampling_test -tc=hmc)
  add_test(NAME logconcave_sampling_test_uld
          COMMAND logconcave_sampling_test -tc=uld)

  TARGET_LINK_LIBRARIES(new_volume_example ${LP_SOLVE})
  TARGET_LINK_LIBRARIES(volume_sob_hpolytope ${LP_SOLVE})
  TARGET_LINK_LIBRARIES(volume_sob_vpolytope ${LP_SOLVE})
  TARGET_LINK_LIBRARIES(volume_cg_hpolytope ${LP_SOLVE})
  TARGET_LINK_LIBRARIES(volume_cg_vpolytope ${LP_SOLVE})
  TARGET_LINK_LIBRARIES(volume_cb_hpolytope ${LP_SOLVE})
  TARGET_LINK_LIBRARIES(volume_cb_vpolytope ${LP_SOLVE})
  TARGET_LINK_LIBRARIES(volume_cb_zonotopes ${LP_SOLVE})
  TARGET_LINK_LIBRARIES(volume_cb_vpoly_intersection_vpoly ${LP_SOLVE})
  TARGET_LINK_LIBRARIES(new_rounding_test ${LP_SOLVE})
  TARGET_LINK_LIBRARIES(mcmc_diagnostics_test ${LP_SOLVE})
  TARGET_LINK_LIBRARIES(benchmarks_sob ${LP_SOLVE})
  TARGET_LINK_LIBRARIES(benchmarks_cg ${LP_SOLVE})
  TARGET_LINK_LIBRARIES(benchmarks_cb ${LP_SOLVE})
  TARGET_LINK_LIBRARIES(ode_solvers_test ${LP_SOLVE} ${IFOPT} ${IFOPT_IPOPT} ${PTHREAD} ${GMP} ${MPSOLVE} ${FFTW3})
  TARGET_LINK_LIBRARIES(boundary_oracles_test ${LP_SOLVE} ${IFOPT} ${IFOPT_IPOPT} ${PTHREAD} ${GMP} ${MPSOLVE} ${FFTW3})
  TARGET_LINK_LIBRARIES(root_finders_test ${PTHREAD} ${GMP} ${MPSOLVE} ${FFTW3})

  if (USE_MKL)
      TARGET_LINK_LIBRARIES(logconcave_sampling_test ${LP_SOLVE} ${IFOPT} ${IFOPT_IPOPT} ${PTHREAD} ${GMP} ${MPSOLVE} ${FFTW3} ${BLAS} "-L${MKLROOT}/lib/intel64 -Wl,--no-as-needed -lmkl_intel_ilp64 -lmkl_gnu_thread -lmkl_core -lgomp -lpthread -lm -ldl")
  else()
      TARGET_LINK_LIBRARIES(logconcave_sampling_test ${LP_SOLVE} ${IFOPT} ${IFOPT_IPOPT} ${PTHREAD} ${GMP} ${MPSOLVE} ${FFTW3} ${BLAS})
  endif()

endif()
