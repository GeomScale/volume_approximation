function(GetEigen)
  find_path(EIGEN_DIR NAMES Eigen PATHS ../../external)

  if (NOT EIGEN_DIR) 
    include(FetchContent)
    FetchContent_Declare(
      eigen
      GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
      GIT_TAG 3.3.9
    )

    FetchContent_GetProperties(eigen)
    
    if(NOT eigen_POPULATED)
      message(STATUS "Eigen library not found locally, downloading it.")
      FetchContent_Populate(eigen)
    endif()

    set(EIGEN_DIR ${eigen_SOURCE_DIR})
    message(STATUS "Using downloaded Eigen library at: ${EIGEN_DIR}")

  else ()

    message(STATUS "Eigen Library found: ${EIGEN_DIR}")

  endif()

  include_directories(${EIGEN_DIR})

endfunction()