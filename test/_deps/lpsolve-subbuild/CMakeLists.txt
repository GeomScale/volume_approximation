# Distributed under the OSI-approved BSD 3-Clause License.  See accompanying
# file Copyright.txt or https://cmake.org/licensing for details.

cmake_minimum_required(VERSION 3.16.3)

# We name the project and the target for the ExternalProject_Add() call
# to something that will highlight to the user what we are working on if
# something goes wrong and an error message is produced.

project(lpsolve-populate NONE)

include(ExternalProject)
ExternalProject_Add(lpsolve-populate
                     "UPDATE_DISCONNECTED" "False" "URL" "https://sourceforge.net/projects/lpsolve/files/lpsolve/5.5.2.11/lp_solve_5.5.2.11_source.tar.gz/download" "URL_HASH" "MD5=a829a8d9c60ff81dc72ff52363703886"
                    SOURCE_DIR          "/home/manav/volesti_fork/test/_deps/lpsolve-src"
                    BINARY_DIR          "/home/manav/volesti_fork/test/_deps/lpsolve-build"
                    CONFIGURE_COMMAND   ""
                    BUILD_COMMAND       ""
                    INSTALL_COMMAND     ""
                    TEST_COMMAND        ""
                    USES_TERMINAL_DOWNLOAD  YES
                    USES_TERMINAL_UPDATE    YES
)
