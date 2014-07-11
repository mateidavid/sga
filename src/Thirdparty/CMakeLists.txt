message(STATUS "Processing: ${CMAKE_CURRENT_SOURCE_DIR}")

# disable errors from warnings for external sources
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wno-error")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-error")

add_library(thirdparty STATIC
    gzstream.C
    stdaln.c
    overlapper.cpp
    multiple_alignment.cpp
    bcr.c
    MurmurHash3.cpp
    count_min_sketch.cpp
)
