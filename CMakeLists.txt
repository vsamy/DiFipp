# Copyright (c) 2019, Vincent SAMY
# All rights reserved.

# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#     * Redistributions of source code must retain the above copyright
#       notice, this list of conditions and the following disclaimer.
#     * Redistributions in binary form must reproduce the above copyright
#       notice, this list of conditions and the following disclaimer in the
#       documentation and/or other materials provided with the distribution.
#     * Neither the name of the AIST nor the
#       names of its contributors may be used to endorse or promote products
#       derived from this software without specific prior written permission.

# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
# ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL <COPYRIGHT HOLDER> BE LIABLE FOR ANY
# DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

# Version minimum
cmake_minimum_required(VERSION 3.0.2)

include(cmake/base.cmake)
include(cmake/eigen.cmake)
include(cmake/msvc-specific.cmake)

set(PROJECT_NAME DiFi++)
set(PROJECT_DESCRIPTION "Filter using rational transfer function")
set(PROJECT_URL "...")
SET(PROJECT_DEBUG_POSTFIX "_d")

#SET(CXX_DISABLE_WERROR True)
set(DOXYGEN_USE_MATHJAX "NO")

project(${PROJECT_NAME} CXX)
set(CMAKE_CXX_STANDARD 14)

setup_project()

option(DISABLE_TESTS "Disable unit tests." OFF)
option(THROW_ON_CONTRACT_VIOLATION "Throw an error when program fails." ON)
option(TERMINATE_ON_CONTRACT_VIOLATION "Terminate program when an error occurs. (Default)" OFF)
option(UNENFORCED_ON_CONTRACT_VIOLATION "Do not perform any check." OFF)

# Handle contracts specifications
if(${THROW_ON_CONTRACT_VIOLATION})
    add_compile_options(-DGSL_THROW_ON_CONTRACT_VIOLATION)
elseif(${TERMINATE_ON_CONTRACT_VIOLATION})
    add_compile_options(-DGSL_TERMINATE_ON_CONTRACT_VIOLATION)
elseif()
    add_compile_options(-DGSL_UNENFORCED_ON_CONTRACT_VIOLATION)
endif()

# for MSVC
if(MSVC)
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /W1 /MP")
    set(CMAKE_MSVCIDE_RUN_PATH "\$(SolutionDir)/src/\$(Configuration)")
endif()

# Eigen
set(Eigen_REQUIRED "eigen3 >= 3.3")
search_for_eigen()

add_subdirectory(include)

if(NOT ${DISABLE_TESTS})
    add_subdirectory(tests)
endif()

pkg_config_append_libs(${PROJECT_NAME})
setup_project_finalize()