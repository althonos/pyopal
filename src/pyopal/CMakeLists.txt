add_subdirectory(platform)

cython_extension(lib SOURCES lib.pyx LINKS opal)