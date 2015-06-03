message(STATUS "Processing: ${CMAKE_CURRENT_SOURCE_DIR}")

# adjust compile flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
# order matters: boost-intrusive-augmented must appear before general boost
include_directories(BEFORE ${CMAKE_SOURCE_DIR}/interval-tree/include)
include_directories(BEFORE ${CMAKE_SOURCE_DIR}/boost-intrusive-augmented/include) # <-first

add_library(mac STATIC
    Mutation.cpp
    Mutation_Cont.cpp
    Read_Chunk_Pos.cpp
    Read_Chunk.cpp
    Read_Chunk_Ptr_Cont.cpp
    Read_Chunk_Cont.cpp
    Read_Entry.cpp
    Contig_Entry.cpp
    Graph.cpp
    Unmap_Mut_Clusters.cpp
    Validate_Variations.cpp
    Hap_Map.cpp
    CLI.cpp
    )

# workaround library&executable having the same name
add_executable(mac-bin mac.cpp)
set_target_properties(mac-bin PROPERTIES OUTPUT_NAME mac)
add_dependencies(mac-bin version)
target_link_libraries(mac-bin
    mac
    suffixtools
    util
    thirdparty
    ${ZLIB_LIBRARIES}
    )

add_executable(test-itree-and-factory test-itree-and-factory.cpp)

add_executable(test-Read_Chunk-traversal test-Read_Chunk-traversal.cpp)
target_link_libraries(test-Read_Chunk-traversal
    mac
    thirdparty
    ${ZLIB_LIBRARIES}
    )

add_executable(test-affine-overlapper test-affine-overlapper.cpp)
target_link_libraries(test-affine-overlapper
    mac
    thirdparty
    ${ZLIB_LIBRARIES}
    )

add_executable(test-RC_Sequence test-RC_Sequence.cpp)
target_link_libraries(test-RC_Sequence
    util
    thirdparty
    ${ZLIB_LIBRARIES}
    )

if(NOT ${CMAKE_BUILD_TYPE} STREQUAL "Release")
    message(STATUS "Including debug and test files and targets; CMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}.")
    file(COPY
        ${CMAKE_SOURCE_DIR}/interval-tree/include/boost/intrusive/itree.gdb
        factory.gdb
        mac-gdb.gdb
        test-itree-and-factory-gdb.gdb
        DESTINATION .)
    add_executable(test-Cigar_Op test-Cigar_Op.cpp)
    add_executable(test-Cigar test-Cigar.cpp)
    add_test(test-Cigar test-Cigar)
    add_test(test-itree-and-factory test-itree-and-factory)
    add_test(test-Read_Chunk-traversal test-Read_Chunk-traversal)
    add_test(test-RC_Sequence test-RC_Sequence)
    add_test(test-affine-overlapper test-affine-overlapper)
    file(GLOB test_graphs test_graphs/test.*)
    file(GLOB test_graphs_real test_graphs/test-real.*)
    foreach(f ${test_graphs})
        get_filename_component(f_name ${f} NAME)
        add_test(mac-${f_name} mac -i ${f} -s)
    endforeach()
    foreach(f ${test_graphs_real})
        get_filename_component(f_name ${f} NAME)
        add_test(mac-${f_name} mac -i ${f} -s --unmap-single-chunks --unmap-mut-clusters)
    endforeach()
    add_executable(sample-logger sample-logger.cpp)
    add_executable(sample-pfor sample-pfor.cpp)
endif()
