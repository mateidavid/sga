message(STATUS "Processing: ${CMAKE_CURRENT_SOURCE_DIR}")

add_library(stringgraph STATIC
    SGUtil.cpp
    SGAlgorithms.cpp
    SGVisitors.cpp
    CompleteOverlapSet.cpp
    RemovalAlgorithm.cpp
    SGSearch.cpp
    SGWalk.cpp
)
