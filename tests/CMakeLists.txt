enable_testing()

set(BOOST_COMPONENTS system unit_test_framework)
search_for_boost()

add_definitions(-DBOOST_TEST_DYN_LINK)
include_directories(${Boost_INCLUDE_DIRS})

macro(addTest testName)
    add_executable(${testName} ${testName}.cpp systems.h tools.h)
    if(NOT MSVC)
        target_link_libraries(${name} PUBLIC ${Boost_LIBRARIES})
    endif(MSVC)
    add_test(${testName}Unit ${testName})
    # Adding a project configuration file (for MSVC only)
    GENERATE_MSVC_DOT_USER_FILE(${name})
endmacro(addTest)

addTest(FilterTest)