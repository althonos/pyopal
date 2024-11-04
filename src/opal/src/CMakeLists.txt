add_library(opal STATIC opal.h ScoreMatrix.cpp)
target_include_directories(opal PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
