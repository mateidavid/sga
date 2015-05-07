message(STATUS "Processing: ${CMAKE_CURRENT_SOURCE_DIR}")

# adjust compile flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

add_executable(BFTrim
    BFTrim.cpp
    )
target_link_libraries(BFTrim
    util
    thirdparty
    ${ZLIB_LIBRARIES}
    )
