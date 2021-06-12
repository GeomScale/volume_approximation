function(Get_LP_Solver)
    find_path(LPSOLVE_DIR NAMES LPsolve_src PATHS ../../external)
    
    if(NOT LPSOLVE_DIR)
    
        include(FetchContent)
            FetchContent_Declare(
                lpSolve 
                URL https://sourceforge.net/projects/lpsolve/files/lpsolve/5.5.2.11/lp_solve_5.5.2.11_source.tar.gz/download 
                URL_HASH MD5=a829a8d9c60ff81dc72ff52363703886
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

    endif()

        include_directories(${LPSOLVE_DIR})
        
endfunction()
