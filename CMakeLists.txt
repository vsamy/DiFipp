# Copyright (c) 2019, Vincent SAMY
# All rights reserved.

# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met: 

# 1. Redistributions of source code must retain the above copyright notice,
#    this list of conditions and the following disclaimer. 
# 2. Redistributions in binary form must reproduce the above copyright notice,
#    this list of conditions and the following disclaimer in the documentation
#    and/or other materials provided with the distribution. 

# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
# ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
# ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

# The views and conclusions contained in the software and documentation are those
# of the authors and should not be interpreted as representing official policies, 
# either expressed or implied, of the FreeBSD Project.

# Version minimum
cmake_minimum_required(VERSION 3.8.2)

set(PROJECT_NAME difi++)
set(PROJECT_DESCRIPTION "Filter using rational transfer function")
set(PROJECT_URL "https://github.com/vsamy/DiFipp")
SET(PROJECT_DEBUG_POSTFIX "_d")
set(INSTALL_GENERATED_HEADERS OFF)

set(DOXYGEN_USE_MATHJAX "NO")
set(CMAKE_CXX_STANDARD 17)

option(BUILD_TESTING "Disable unit tests." ON)
option(THROW_ON_CONTRACT_VIOLATION "Throw an error when program fails." ON)
option(TERMINATE_ON_CONTRACT_VIOLATION "Terminate program when an error occurs. (Default)" OFF)
option(UNENFORCED_ON_CONTRACT_VIOLATION "Do not perform any check." OFF)

include(cmake/base.cmake)
include(cmake/msvc-specific.cmake)

project(${PROJECT_NAME} CXX)

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
    set(CMAKE_MSVCIDE_RUN_PATH "\$(SolutionDir)/src/\$(Configuration)")
endif()

# Eigen
set(Eigen_REQUIRED "eigen3 >= 3.3")
add_project_dependency(Eigen3 REQUIRED)

add_subdirectory(include)

if(${BUILD_TESTING})
    add_subdirectory(tests)
endif()
