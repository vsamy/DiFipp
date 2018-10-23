# Version minimum
cmake_minimum_required(VERSION 3)

include(cmake/base.cmake)
include(cmake/eigen.cmake)

set(PROJECT_NAME fratio)
set(PROJECT_DESCRIPTION "Filter using rational transfer function")
set(PROJECT_URL "...")

#SET(CXX_DISABLE_WERROR True)
set(DOXYGEN_USE_MATHJAX "NO")

project(${PROJECT_NAME} CXX)
set(CMAKE_CXX_STANDARD 14)

setup_project()

# for MSVC
if(MSVC)
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /W1 /MP")
endif()

# Eigen
set(Eigen_REQUIRED "eigen3 >= 3.3")
search_for_eigen()

add_subdirectory(src)

setup_project_finalize()