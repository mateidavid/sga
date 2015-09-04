message(STATUS "Processing: ${CMAKE_CURRENT_SOURCE_DIR}")

# adjust compile flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

add_executable(BFTrim BFTrim.cpp)
target_link_libraries(BFTrim
    util
    thirdparty
    version
    ${ZLIB_LIBRARIES}
    )

if(NOT ${CMAKE_BUILD_TYPE} STREQUAL "Release")
    message(STATUS "Including debug and test files and targets; CMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}.")
    add_executable(test-kmer_gen test-kmer_gen.cpp)
endif()
