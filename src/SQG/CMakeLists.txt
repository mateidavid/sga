message(STATUS "Processing: ${CMAKE_CURRENT_SOURCE_DIR}")

add_library(sqg STATIC
    SQG.cpp
    ASQG.cpp
)
