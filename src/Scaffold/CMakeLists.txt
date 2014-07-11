message(STATUS "Processing: ${CMAKE_CURRENT_SOURCE_DIR}")

add_library(scaffold STATIC
    ScaffoldVertex.cpp
    ScaffoldEdge.cpp
    ScaffoldGraph.cpp
    ScaffoldVisitors.cpp
    ScaffoldLink.cpp
    ScaffoldRecord.cpp
    ScaffoldGroup.cpp
    ScaffoldSearch.cpp
    ScaffoldWalk.cpp
    ScaffoldAlgorithms.cpp
    ScaffoldSequenceCollection.cpp
)
