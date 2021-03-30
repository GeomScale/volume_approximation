function(GetBoost)

	find_path(BOOST_DIR NAMES boost PATHS ../external/)

	if (NOT BOOST_DIR) 
  
    set(BOOST_URL "https://dl.bintray.com/boostorg/release/1.75.0/source/boost_1_75_0.tar.bz2" CACHE STRING "Boost download URL")
    set(BOOST_URL_SHA256 "953db31e016db7bb207f11432bef7df100516eeb746843fa0486a222e3fd49cb" CACHE STRING "Boost download URL SHA256 checksum")

    include(FetchContent)
    FetchContent_Declare(
      Boost
      URL ${BOOST_URL}
      URL_HASH SHA256=${BOOST_URL_SHA256}
    )
    FetchContent_GetProperties(Boost)

    if(NOT Boost_POPULATED)
      message(STATUS "Fetching Boost")
      FetchContent_Populate(Boost)
      message(STATUS "Fetching Boost - done")
      set(BOOST_DIR ${boost_SOURCE_DIR})
    endif()

    message(STATUS "Using downloaded Boost library at ${BOOST_DIR}")

  else ()
    message(STATUS "Boost Library found: ${BOOST_DIR}")

  endif()

  include_directories(${BOOST_DIR})
  include_directories(${BOOST_DIR}/boost)

endfunction()
