function(GetLPSolve)
    find_path(LP_SOLVE_DIR NAMES LPsolve_src PATHS ../../external)

    if (NOT LP_SOLVE_DIR)    
        include(FetchContent)
        FetchContent_Declare(
            lpsolve 
            URL https://cran.r-project.org/src/contrib/lpSolve_5.6.15.tar.gz
            URL_HASH MD5=3b8d780f703e0da2e4863939add164bc
        )
        
        FetchContent_GetProperties(lpsolve)
            
        if(NOT lpsolve_POPULATED)
            message(STATUS "lp_solve library not found locally, downloading it.")
            FetchContent_Populate(lpsolve)
        endif()
        
        set(LP_SOLVE_DIR "${lpsolve_SOURCE_DIR}/src")
        message(STATUS "Using downloaded lp_solve at: ${LP_SOLVE_DIR}")

    else()

        message(STATUS "lp_solve library found: ${LP_SOLVE_DIR}")
    
    endif()

    include_directories(${LP_SOLVE_DIR})
        
endfunction()
