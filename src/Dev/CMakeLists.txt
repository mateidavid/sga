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
    Unmapper.cpp
    Unmap_Mut_Clusters.cpp
    Validate_Variations.cpp
    Hap_Map.cpp
    Read_Merger.cpp
    Hap_Tree.cpp
    Anchor_Support.cpp
    CLI.cpp
    MAC_explicit_instantiations.cpp
    )

# workaround library&executable having the same name
add_executable(mac-bin mac.cpp)
set_target_properties(mac-bin PROPERTIES OUTPUT_NAME mac)
if(NOT ${CMAKE_BUILD_TYPE} STREQUAL "Release")
    target_link_libraries(mac-bin -Wl,--whole-archive)
endif()
target_link_libraries(mac-bin
    mac
    suffixtools
    util
    thirdparty
    version
    )
if(NOT ${CMAKE_BUILD_TYPE} STREQUAL "Release")
    target_link_libraries(mac-bin -Wl,--no-whole-archive)
endif()
target_link_libraries(mac-bin
    ${ZLIB_LIBRARIES}
    rt
    )

if(NOT ${CMAKE_BUILD_TYPE} STREQUAL "Release")
    # extract active_ptr addresses for reliable gdb operation
    add_custom_command(TARGET mac-bin POST_BUILD
        COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/scripts/get_active_ptr_addresses mac >mac-addresses.gdb
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        )
endif()

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
    add_executable(test_algorithm_extra test_algorithm_extra.cpp)
    add_test(test_algorithm_extra test_algorithm_extra)
    add_executable(sample-consensus-alignment sample-consensus-alignment.cpp)
    target_link_libraries(sample-consensus-alignment rt)
endif()
