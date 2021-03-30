# VolEsti (volume computation and sampling library)
# Copyright (c) 2012-2018 Vissarion Fisikopoulos
# Copyright (c) 2018 Apostolos Chalkis
# Contributed and/or modified by Marios Papachristou, as part of Google Summer of Code 2020 program.
# Licensed under GNU LGPL.3, see LICENCE file


project( VolEsti )


CMAKE_MINIMUM_REQUIRED(VERSION 2.4.5)

set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)

if(COMMAND cmake_policy)
       cmake_policy(SET CMP0003 NEW)
endif(COMMAND cmake_policy)


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

include("Eigen.cmake")
GetEigen()

include("Boost.cmake")
GetBoost()

# Find lpsolve library
find_library(LP_SOLVE NAMES liblpsolve55.so PATHS /usr/lib/lp_solve)

if (NOT LP_SOLVE)
  message(FATAL_ERROR "This program requires the lp_solve library, and will not be compiled.")
else ()
  message(STATUS "Library lp_solve found: ${LP_SOLVE}")

  set(CMAKE_EXPORT_COMPILE_COMMANDS "ON")

  include_directories (BEFORE ../../external/Eigen)
  include_directories (BEFORE ../../external)
  include_directories (BEFORE ../../external/minimum_ellipsoid)
  #include_directories (BEFORE ../../include/cartesian_geom)
  #include_directories (BEFORE ../../include/convex_bodies)
  include_directories (BEFORE ../../external/LPsolve_src/run_headers)
  include_directories (BEFORE ../../external/boost)
  #include_directories (BEFORE BOOST)
  include_directories (BEFORE ../../include/generators)
  include_directories (BEFORE ../../include/volume)
  include_directories (BEFORE ../../include)
  include_directories (BEFORE ../../include/lp_oracles)
  include_directories (BEFORE ../../include/nlp_oracles)
  include_directories (BEFORE ../../include/convex_bodies)
  include_directories (BEFORE ../../include/random_walks)
  include_directories (BEFORE ../../include/annealing)
  include_directories (BEFORE ../../include/ode_solvers)
  include_directories (BEFORE ../../include/root_finders)
  include_directories (BEFORE ../../include/samplers)
  include_directories (BEFORE ../../include/misc)
  include_directories (BEFORE ../../include/optimization)

  # for Eigen
  if (${CMAKE_VERSION} VERSION_LESS "3.12.0")
    add_compile_options(-D "EIGEN_NO_DEBUG")
  else ()
    add_compile_definitions("EIGEN_NO_DEBUG")
  endif ()


  add_definitions(${CMAKE_CXX_FLAGS} "-std=c++11")  # enable C++11 standard
  add_definitions(${CMAKE_CXX_FLAGS} "-O3")  # optimization of the compiler
  #add_definitions(${CXX_COVERAGE_COMPILE_FLAGS} "-lgsl")
  add_definitions(${CXX_COVERAGE_COMPILE_FLAGS} "-lm")
  add_definitions(${CXX_COVERAGE_COMPILE_FLAGS} "-ldl")
  add_definitions(${CXX_COVERAGE_COMPILE_FLAGS} "-DBOOST_NO_AUTO_PTR")
  #add_definitions(${CXX_COVERAGE_COMPILE_FLAGS} "-lgslcblas")
  #add_definitions( "-O3 -lgsl -lm -ldl -lgslcblas" )
  
  find_package(OpenMP REQUIRED)

  add_executable (skinny_cube_10_dim skinny_cube_10_dim.cpp)
  add_executable (random_hpoly_50_dim random_hpoly_50_dim.cpp)
  add_executable (random_hpoly_50_dim_parallel random_hpoly_50_dim_parallel.cpp)

  TARGET_LINK_LIBRARIES(skinny_cube_10_dim ${LP_SOLVE})
  TARGET_LINK_LIBRARIES(random_hpoly_50_dim ${LP_SOLVE})
  TARGET_LINK_LIBRARIES(random_hpoly_50_dim_parallel ${LP_SOLVE} OpenMP::OpenMP_CXX)


endif()
