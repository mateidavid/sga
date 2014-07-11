message(STATUS "Processing: ${CMAKE_CURRENT_SOURCE_DIR}")

add_library(concurrency STATIC
    OverlapProcess.cpp
    RmdupProcess.cpp
)
