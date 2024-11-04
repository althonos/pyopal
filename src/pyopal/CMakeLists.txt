add_subdirectory(platform)

cython_extension(lib SOURCES lib.pyx LINKS opal)
install(FILES ${PROJECT_SOURCE_DIR}/vendor/opal/src/opal.h DESTINATION ${CMAKE_INSTALL_PREFIX}/pyopal)
