set(HEADERS
    Butterworth.h
    DigitalFilter.h
    GenericFilter.h
    MovingAverage.h
    polynome_functions.h
)

add_library(${PROJECT_NAME} INTERFACE)
target_include_directories(${PROJECT_NAME} INTERFACE ${HEADERS})
install(TARGETS ${PROJECT_NAME}
    RUNTIME DESTINATION bin
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib)