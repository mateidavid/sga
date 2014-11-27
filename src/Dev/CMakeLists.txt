message(STATUS "Processing: ${CMAKE_CURRENT_SOURCE_DIR}")

# adjust compile flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
# order matters: local boost-intrusive must appear before general boost
include_directories(BEFORE ${CMAKE_SOURCE_DIR}/interval-tree/include)
include_directories(BEFORE ${CMAKE_SOURCE_DIR}/boost-intrusive/include) # <-first

add_library(mac STATIC
    global_assert.cpp
    Mutation.cpp
    Mutation_Cont.cpp
    Read_Chunk_Pos.cpp
    Read_Chunk.cpp
    Read_Chunk_Cont.cpp
    Read_Entry.cpp
    Contig_Entry.cpp
    Graph.cpp
    Unmap_Mut_Clusters.cpp
)
#set_target_properties(mac PROPERTIES LINK_FLAGS "${BOOST_LDFLAGS}")

# workaround library&executable having the same name
add_executable(mac-bin mac.cpp)
set_target_properties(mac-bin PROPERTIES OUTPUT_NAME mac)
add_dependencies(mac-bin version)
#set_target_properties(mac-bin PROPERTIES LINK_FLAGS "${BOOST_LDFLAGS} -lboost_iostreams -lboost_program_options")
#target_link_libraries(mac-bin mac ${CMAKE_BINARY_DIR}/Util/libutil.a ${CMAKE_BINARY_DIR}/Thirdparty/libthirdparty.a)
target_link_libraries(mac-bin
    mac
    thirdparty
    ${Boost_IOSTREAMS_LIBRARY}
    ${Boost_PROGRAM_OPTIONS_LIBRARY}
    ${ZLIB_LIBRARIES}
)

add_executable(test-itree-and-factory
    test-itree-and-factory.cpp
    global_assert.cpp
)
target_link_libraries(test-itree-and-factory
    ${Boost_PROGRAM_OPTIONS_LIBRARY}
)

add_executable(test-Read_Chunk-traversal
    test-Read_Chunk-traversal.cpp
    globals.cpp
    global_assert.cpp
)
target_link_libraries(test-Read_Chunk-traversal
    mac
    thirdparty
    ${Boost_PROGRAM_OPTIONS_LIBRARY}
    ${ZLIB_LIBRARIES}
)

add_executable(check-assert check-assert.cpp global_assert.cpp)

add_executable(test-affine-overlapper test-affine-overlapper.cpp)
target_link_libraries(test-affine-overlapper
    mac
    thirdparty
    ${Boost_PROGRAM_OPTIONS_LIBRARY}
    ${ZLIB_LIBRARIES}
)

add_executable(test-DNA_Sequence test-DNA_Sequence.cpp)
target_link_libraries(test-DNA_Sequence
    util
    thirdparty
    ${Boost_PROGRAM_OPTIONS_LIBRARY}
    ${ZLIB_LIBRARIES}
)

if(NOT ${CMAKE_BUILD_TYPE} STREQUAL "Release")
  message(STATUS "Including debug and test files and targets; CMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}.")
  file(COPY factory.gdb mac-gdb.gdb DESTINATION .)
  add_test(test-itree-and-factory test-itree-and-factory)
  add_test(test-Read_Chunk-traversal test-Read_Chunk-traversal)
  add_test(test-DNA_Sequence test-DNA_Sequence)
  add_test(test-affine-overlapper test-affine-overlapper)
  file(GLOB test_graphs test_graphs/test.*)
  file(GLOB test_graphs_real test_graphs/test-real.*)
  foreach(f ${test_graphs})
    get_filename_component(f_name ${f} NAME)
    add_test(mac-${f_name} mac -i ${f} -s)
  endforeach()
  foreach(f ${test_graphs_real})
    get_filename_component(f_name ${f} NAME)
    add_test(mac-${f_name} mac -i ${f} -s --resolve-unmappable --unmap-single-chunks)
  endforeach()
endif()
