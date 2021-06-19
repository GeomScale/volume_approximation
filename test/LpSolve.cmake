function(GetLpSolve)
    find_path(LPSOLVE_DIR NAMES lp_solve_5.5 PATHS ../../external)

    if(NOT LPSOLVE_DIR)

        include(FetchContent)
            FetchContent_Declare(
                lpSolve
                URL https://sourceforge.net/projects/lpsolve/files/lpsolve/5.5.2.0/lp_solve_5.5.2.0_source.tar.gz/download
            )

            FetchContent_GetProperties(lpSolve)

            if(NOT lpSolve_POPULATE)
                message(STATUS "Downloading lp_solve.")
                FetchContent_Populate(lpSolve)
                message(STATUS "Download complete")

                set(LPSOLVE_DIR ${lpSolve_SOURCE_DIR})
                message(STATUS "Using downloaded lp_solve at: ${LPSOLVE_DIR}")

            else()
                message(STATUS "lp_solve found:  ${LPSOLVE_DIR}")

            endif()

            set(LPSOLVE_DIR "_deps/lpsolve-src")
            message(STATUS "Using downloaded lp_solve at: ${LPSOLVE_DIR}")

            # target bar is only build when `make bar` is issued
            add_custom_target(liblpsolve55.dylib
                COMMAND sh ccc.osx
                WORKING_DIRECTORY ${LP_SOLVE_DIR}/lpsolve55
            )


    endif()

        include_directories(${LPSOLVE_DIR})

endfunction()
