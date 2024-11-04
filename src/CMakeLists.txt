file(COPY ../vendor/opal DESTINATION .)
file(COPY ./opal DESTINATION .)

add_subdirectory(${CMAKE_CURRENT_BINARY_DIR}/opal ${CMAKE_CURRENT_BINARY_DIR}/opal.build)
add_subdirectory(pyopal)