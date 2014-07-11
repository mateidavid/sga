message(STATUS "Processing: ${CMAKE_CURRENT_SOURCE_DIR}")

add_library(bigraph STATIC
    Bigraph.cpp
    Vertex.cpp
    Edge.cpp
    EdgeDesc.cpp
)
