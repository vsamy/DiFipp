include(${CMAKE_SOURCE_DIR}/cmake/boost.cmake)

enable_testing()

set(BOOST_COMPONENTS unit_test_framework)
search_for_boost()

add_definitions(-DBOOST_TEST_DYN_LINK)
include_directories(${Boost_INCLUDE_DIRS})
include_directories(${CMAKE_SOURCE_DIR}/include)
link_directories(${Boost_LIBRARY_DIR} ${Boost_LIBRARY_DIRS})

macro(addTest testName)
    add_executable(${testName} ${testName}.cpp)
    if(NOT MSVC)
        target_link_libraries(${testName} PUBLIC ${Boost_LIBRARIES})
    endif()
    add_test(${testName}Unit ${testName})
    # Adding a project configuration file (for MSVC only)
    GENERATE_MSVC_DOT_USER_FILE(${testName})
endmacro(addTest)

addTest(ButterWorthFilterTests)
addTest(GenericFilterTests)
addTest(DigitalFilterTests)
addTest(MovingAverageFilterTests)
addTest(polynome_functions_tests)