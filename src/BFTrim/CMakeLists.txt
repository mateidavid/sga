message(STATUS "Processing: ${CMAKE_CURRENT_SOURCE_DIR}")

# adjust compile flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

add_executable(BFTrim
    BFTrim.cpp
    )
target_link_libraries(BFTrim
    util
    thirdparty
    ${Boost_IOSTREAMS_LIBRARY}
    ${Boost_PROGRAM_OPTIONS_LIBRARY}
    ${ZLIB_LIBRARIES}
    )

if(OPENMP_FOUND)
    add_executable(BFTrim-omp
        BFTrim-omp.cpp
        )
    target_link_libraries(BFTrim-omp
        util
        thirdparty
        ${Boost_IOSTREAMS_LIBRARY}
        ${Boost_PROGRAM_OPTIONS_LIBRARY}
        ${ZLIB_LIBRARIES}
        )
else()
    message(STATUS "Cannot build BFTrim-omp, need OpenMP")
endif()
