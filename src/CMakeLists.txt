set(HEADERS
    ../include/GenericFilter.h
)

set(SRC
    GenericFilter.cpp
)

add_library(${PROJECT_NAME} ${HEADERS} ${SRC})
install(TARGETS ${PROJECT_NAME}
    RUNTIME DESTINATION bin
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib)