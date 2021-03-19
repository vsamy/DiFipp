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

enable_testing()

macro(addTest testName)
    add_executable(${testName} ${testName}.cpp)
    if (MSVC)
        target_compile_definitions(${testName} PUBLIC _SILENCE_ALL_CXX17_DEPRECATION_WARNINGS)
    endif()
    target_compile_definitions(${testName} PUBLIC DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN)
    target_link_libraries(${testName} PUBLIC Eigen3::Eigen)
    # Adding a project configuration file (for MSVC only)
    generate_msvc_dot_user_file(${testName})

    add_test(${testName}Unit ${testName})
endmacro(addTest)

addTest(GenericFilterTests)
addTest(polynome_functions_tests)
addTest(DigitalFilterTests)
addTest(MovingAverageFilterTests)
addTest(ButterworthFilterTests)

# Differentiators
addTest(differentiator_tests)