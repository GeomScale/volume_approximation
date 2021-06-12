function(GetLPsolver)
    find_path(LPSOLVE_DIR NAMES LPsolve_src PATHS ../../external)

  if (NOT LPSOLVE_DIR) 
    include(FetchContent)
    FetchContent_Declare(
        lpsolve
        #GIT_REPOSITORY https://github.com/gaborcsardi/lpSolve.git 
        #GIT_TAG 5.6.13.3
        URL https://sourceforge.net/projects/lpsolve/files/lpsolve/5.5.2.11/lp_solve_5.5.2.11_source.tar.gz/download 
        URL_HASH MD5=a829a8d9c60ff81dc72ff52363703886
    )

    FetchContent_GetProperties(lpSolve)
    
    if(NOT lpsolve_POPULATED)
      message(STATUS "lpSolve not found locally, downloading it.")
      FetchContent_Populate(lpsolve)
    endif()

    set(LPSOLVE_DIR "${lpsolve_SOURCE_DIR}")
    message(STATUS "Using downloaded lpSolve library at: ${LPSOLVE_DIR}")

  else ()

    message(STATUS "Eigen Library found: ${LPSOLVE_DIR}")

endif()

  include_directories(${LPSOLVE_DIR})

endfunction()
