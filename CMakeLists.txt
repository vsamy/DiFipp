# Version minimum
cmake_minimum_required(VERSION 3.0.2)

include(cmake/base.cmake)
include(cmake/eigen.cmake)
include(cmake/msvc-specific.cmake)

set(PROJECT_NAME fratio)
set(PROJECT_DESCRIPTION "Filter using rational transfer function")
set(PROJECT_URL "...")
SET(PROJECT_DEBUG_POSTFIX "_d")

#SET(CXX_DISABLE_WERROR True)
set(DOXYGEN_USE_MATHJAX "NO")

project(${PROJECT_NAME} CXX)
set(CMAKE_CXX_STANDARD 14)

setup_project()

option(DISABLE_TESTS "Disable unit tests." OFF)

# for MSVC
if(MSVC)
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /W1 /MP")
    set(CMAKE_MSVCIDE_RUN_PATH "\$(SolutionDir)/src/\$(Configuration)")
endif()

add_compile_options("-D_USE_MATH_DEFINES")

# Eigen
set(Eigen_REQUIRED "eigen3 >= 3.3")
search_for_eigen()

add_subdirectory(include)

if(NOT ${DISABLE_TESTS})
    add_subdirectory(tests)
endif()

pkg_config_append_libs(${PROJECT_NAME})
setup_project_finalize()